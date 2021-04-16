# mrOpen.R : Mark-recapture for open populations.
# Jan Galkowski, extensions of ideas by R. Tanaka.
# Begun 11th February 2020, 
# last changed 19th February 2020.

library(segmented)

tanaka<- function(ratios, covariate, intercept=0, 
                  covariateLabel="number of allocated marks through sample (k-1)",
                  tag="[capture-recapture experiment]", 
                  ftag="TanakaX-", truePop=NULL, show=TRUE,
                  nsd=1, ptCex=1.1, ptCol="blue")
{
  if (0 == intercept)
  {
    # (Would like to make "singular.ok" a variable, but can't, as apparently
    #  it gets frozen out of context when the segmented is called, and 
    #  yields a value error when it is looked up.)
    fitOverall<- lm(ratios ~ covariate + 0, singular.ok=FALSE,
                    model=TRUE, qr=TRUE, x=TRUE, y=TRUE, method="qr")
    slope<- fitOverall$coefficients
    summ<- summary(fitOverall)
    stderrOverall<- sqrt(diag(vcov(fitOverall)))
    solution<- round(1/slope)
    solutionLow<- round(1/(slope+nsd*stderrOverall))
    solutionHigh<- round(1/(slope-nsd*stderrOverall))
  } else
  {
    # (Would like to make "singular.ok" a variable, but can't, as apparently
    #  it gets frozen out of context when the segmented is called, and 
    #  yields a value error when it is looked up.)
    fitOverall<- lm(ratios ~ covariate + 1, 
                    singular.ok=FALSE,
                    model=TRUE, qr=TRUE, x=TRUE, y=TRUE, method="qr")
    slope<- fitOverall$coefficients[2]
    intercept<- fitOverall$coefficients[1]
    summ<- summary(fitOverall)
    stderrOverall<- sqrt(diag(vcov(fitOverall)))
    solution<- round((1-intercept)/slope)
    solutionLow<- round((1-(intercept-stderrOverall[1]))/(slope+nsd*stderrOverall[2]))
    solutionHigh<- round((1-(intercept+stderrOverall[1]))/(slope-nsd*stderrOverall[2]))
  }
  #
  if (show)
  {
    tSNow<- timingStamp()
    pdf(onefile=TRUE, title=sprintf("%s, %s", tag, tSNow),
        file=constructFilenameFrom(root=ftag, suffix=".pdf", stamp=tSNow)$full,
        width=22, height=17, pointsize=21)
    #
    plot(x=covariate, y=ratios, type="n", 
         xlim=c(0,(1.2*solutionHigh)), 
         ylim=c(0,1), xlab=covariateLabel, 
         ylab="proportion of captured in sample k which are marked",
         main=sprintf("%s\n(extended Tanaka method)", tag), cex.sub=0.9, 
         sub=sprintf("(Bounds are %.2f standard errors. Intersections are offset because of rounding to integer.)", nsd))
    abline(a=0, b=slope, lty=6, lwd=2, col="navy")
    abline(h=1, col="red", lwd=1, lty=1)
    abline(v=solution, lwd=2, col="darkred", lty=1)
    if (is.numeric(truePop))
    {
      abline(v=truePop, lty=2, lwd=2, col="darkgreen")
      text(x=truePop, y=0.4, pos=4, offset=2, label=sprintf("true population size %.0f", truePop))
      arrows(x0=truePop+220, y0=0.4, lwd=6, length=0.2, x1=truePop+15, y1=0.4, code=2, col="green")
    }
    print(c(slope, nsd, stderrOverall))
    abline(a=0, b=(slope+nsd*stderrOverall), col="darkgrey", lwd=3)
    abline(a=0, b=(slope-nsd*stderrOverall), col="darkgrey", lwd=3)
    abline(v=solutionLow, col="darkslategrey", lwd=1, lty=3)
    abline(v=solutionHigh, col="darkslategrey", lwd=1, lty=3)
    points(covariate, ratios, pch=21, 
           col=ptCol, bg=ptCol, cex=ptCex)
    text(x=solutionHigh, y=0.7, pos=4, label=sprintf("high bound %.0f", solutionHigh))
    text(x=solutionLow, y=0.3, pos=2, label=sprintf("low bound %.0f", solutionLow))
    text(x=solutionLow, y=0.5, pos=2, label=sprintf("solution %.0f", solution))
    arrows(x0=solutionLow-50, y0=0.5, lwd=6, length=0.2, x1=solution-15, y1=0.5, code=2, col="red")
    dev.off()
    print(summ)
  }
  #
  return(list(fit=fitOverall, slope=slope, summary=summ, stderr=stderrOverall,
              solution=solution, solution.low=solutionLow, 
              solution.high=solutionHigh))
}

reassignSegments<- function(tooSmall, beginEnd)
{
  stopifnot( is.vector(tooSmall) && (0 < length(tooSmall)) )
  stopifnot( is.matrix(beginEnd) )
  stopifnot( 2 == ncol(beginEnd) )
  stopifnot( all(beginEnd[,2] >= beginEnd[,1]) )
  nB<- nrow(beginEnd)
  stopifnot( 2 <= nB )
  # Could be more sophisticated, but here simply assigning elided
  # segments all to the previous one, except in the case where an
  # elided segment is the very first.
  beginEndEdit<- beginEnd
  cat("at start:\n")
  print(beginEndEdit)
  while (1 == tooSmall[1])
  {
    beginEndEdit[2,1]<- beginEndEdit[1,1]
    beginEndEdit<- beginEndEdit[(2:nB),]
    nB<- nB-1
    tooSmall<- tooSmall - 1
    tooSmall<- tooSmall[-1]
    cat(sprintf("beginEnd reduced to %.0f rows.\n", nB))
  }
  cat("after initial:\n")
  print(beginEndEdit)
  if (0 < length(tooSmall))
  {
    while (0 < length(tooSmall))
    {
      k<- tooSmall[1]
      stopifnot( 1 < k )
      beginEndEdit[(k-1),2]<- beginEndEdit[k,2]
      beginEndEdit<- beginEndEdit[-k,]
      tooSmall<- tooSmall[-1]
      nB<- nB-1
      cat(sprintf("beginEnd reduced to %.0f rows.\n", nB))
      print(beginEndEdit)
      cat("\n")
      print(tooSmall)
      cat("\n")
    }
  }
  return(beginEndEdit)
}

popEstimatesPerSegmentFrom<- function(nS, ratios, previouslyAllocated, 
                                      segmentingStages, intercept=0,
                                      tag="",
                                      ftag="",
                                      nsd=1,
                                      ptCex=1.1,
                                      ptCol="blue",
                                      eliminateSmallSegments=FALSE)
{
  stopifnot( all( c("psi", "id.group") %in% names(segmentingStages) ) )
  psi<- segmentingStages$psi
  id.group<- 1+segmentingStages$id.group
  nSegments<- 1+nrow(psi)
  breakIndices<- c(1,((3:(nS-1))[which(as.logical(diff(id.group)))]),(nS-1))
  nB<- length(breakIndices)
  #
  beginEnd<- matrix(c((1+breakIndices[1:(nB-1)]), breakIndices[2:nB]), 
                    nrow=(nB-1), ncol=2, byrow=FALSE)
  nB<- nrow(beginEnd)
  if (eliminateSmallSegments)
  {
    # Reassociate to be sure segments have at least 5 points apiece
    tooSmall<- which(0 < apply(X=beginEnd, MARGIN=1, FUN=function(r) (5 > (1+r[2]-r[1]))))
    if (0 < length(tooSmall))
    {
      cat("Consolidating segments: ")
      print(tooSmall)
      cat(" because some less than 5 points.\n")
      beginEnd<- reassignSegments(tooSmall, beginEnd)
      nB<- nrow(beginEnd)
    }
  }
  #
  print(beginEnd)
  #
  fittings<- vector(mode="list", length=nB)
  for (kB in (1:nB))
  {
    beginning<- beginEnd[kB,1]
    ending<- beginEnd[kB,2]
    slice<- seq(beginning, ending, 1)
    ratiosSlice<- ratios[slice]
    previouslyAllocatedSlice<- previouslyAllocated[slice]
    cL<- sprintf(
        "number of allocated marks from sample stages %.0f through %.0f",
                 beginning, ending)
    fO<- tanaka(ratios=ratiosSlice, 
                covariate=previouslyAllocatedSlice,
                intercept=intercept, 
                covariateLabel=cL,
                tag=sprintf("%s, segment %-.0f", tag, kB), 
                ftag=sprintf("%s-segment-%-.0f", ftag, kB),
                nsd=nsd, ptCex=ptCex, ptCol=ptCol)
    fittings[[kB]]<- fO
  }
  return(list(groups=id.group, breaks=breakIndices, fits=fittings))
}

# Need to borrow entirety of Professor Muggeo's confint.segmented
# because a call to "signif" occasionally faults. This is a result
# of faults in one of the calls to the "splinefun" built-in function, 
# which complains that the number of points being submitted to them
# is too small. This results in a final variable being coerced 
# to character, and "signif" faults because it expects numeric.

source("myConfInt.R")

tanakaX<- function(marksMarked, show=TRUE, nsd=1, intercept=0, 
                   ptCex=1.1, ptCol="blue", extend=TRUE,
                   tag="[capture-recapture experiment]", 
                   ftag="TanakaX-", Kbreaks=10, 
                   method="score", 
                   maxIterationsOuter=400,
                   maxIterationsInner=200,
                   truePop=NULL, conf=0.8, eliminateSmallSegments=FALSE)
{
  stopifnot( is.matrix(marksMarked) )
  stopifnot( 3 == ncol(marksMarked) )
  stopifnot( all( c("markedInSample", "sampleSize", "allocatedDistinctMarks") == 
                  colnames(marksMarked) ) )
  #
  segmenting<- NA
  confidenceSegmenting<- NA
  segmentingStages<- NA
  confidenceSegmentingStages<- NA
  popsBySegment<- NA
  #
  nS<- nrow(marksMarked)
  ratios<- marksMarked[(2:nS), "markedInSample"]/marksMarked[(2:nS),"sampleSize"]
  previouslyAllocated<- marksMarked[(1:(nS-1)),"allocatedDistinctMarks"]
  #
  fO<- tanaka(ratios=ratios, covariate=previouslyAllocated,
              intercept=intercept, show=show,
              covariateLabel="number of allocated marks through sample (k-1)",
              tag=tag, ftag=ftag, nsd=nsd, ptCex=ptCex, ptCol=ptCol)
  fitOverall<- fO$fit
  #
  if (extend)
  {
    cat("extend == TRUE, segmenting ...\n")
    ctrl<- seg.control(n.boot=0, display=TRUE, tol=1e-05, 
                       it.max=maxIterationsOuter, 
                       fix.npsi=FALSE,
                       K=Kbreaks, quant=TRUE, 
                       maxit.glm=maxIterationsInner, 
                       h=1, 
                       size.boot=NULL, jt=FALSE, 
                       nonParam=TRUE, random=TRUE, 
                       seed=randomizeSeed(external=TRUE), 
#                      fn.obj="sum(x$residuals^2)", 
                       fn.obj=NULL, 
                       conv.psi=FALSE,
                       alpha=.02, min.step=.0001,
                       digits=4,
#                      powers=c(1,1), # internal ... do not touch
                       last=TRUE) 
#                      stop.if.error=NULL, # do not touch
#                      gap=FALSE) # do not touch
#
    if (0 == intercept)
    {
      segmenting<- segmented(obj=fitOverall, 
#                            seg.Z, # missing because only univariate covariate
                             psi=NA, 
#                            npsi,  # ignored because psi specified
                             control=ctrl,
                             model=TRUE, 
                             keep.class=TRUE)  
    } else
    {
      candidateBreaks<- round(hdquantile(previouslyAllocated, probs=c(0.30, 0.70), names=FALSE))
      segmenting<- segmented(obj=fitOverall, 
                             seg.Z=~covariate, # missing because only univariate covariate
                             psi=list(covariate=candidateBreaks), 
#                            npsi,  # ignored because psi specified
                             control=ctrl,
                             model=TRUE, 
                             keep.class=TRUE)  
      
    }
#   confidenceSegmenting<- confint(object=segmenting, 
    confidenceSegmenting<- MyConfintSegmented(object=segmenting, 
#                                  parm, # omit
                                   level=conf, 
#                                  method=c("delta", "score", "gradient"),
                                   method=method,
                                   rev.sgn=FALSE, 
                                   var.diff=FALSE, 
                                   digits=4,
                                   is=FALSE
                                )
    cat("Confidence intervals report from segmenting ...\n")
    print(confidenceSegmenting)
    print( summary(object=segmenting, short=FALSE, var.diff=FALSE, 
                   p.df="K", signif.stars=TRUE) )
    if (show)
    {
      cat("Plotting segmentation ...\n")
      tSNow<- timingStamp()
      pdf(onefile=TRUE, title=sprintf("%s confidence by marked, %s", tag, tSNow),
          file=constructFilenameFrom(root=sprintf("%s-confidence-by-marked", ftag), 
                                     suffix=".pdf", stamp=tSNow)$full,
          width=22, height=17, pointsize=21)
      plot(x=segmenting, 
#          term, # omitted
           add=FALSE, res=TRUE, conf.level=conf, interc=FALSE,
#          link=TRUE, # ignored
           res.col="green", rev.sgn=FALSE, const=0, shade=TRUE, rug=TRUE,
           dens.rug=TRUE, dens.col = grey(0.8), transf=I, isV=FALSE, is=FALSE,
           var.diff=FALSE, p.df="K", .vcov=NULL)
      dev.off()
    }
  }
  #
  if (show && extend)
  {
    #
    tSNow<- timingStamp()
    pdf(onefile=TRUE, title=sprintf("%s, %s", 
                               sprintf("%s, extend close-up", tag), tSNow),
        file=constructFilenameFrom(root=sprintf("%s-xcu", ftag), suffix=".pdf", 
                                   stamp=tSNow)$full,
        width=22, height=17, pointsize=21)
    #
    stagesWithData<- 2:nS
    plot(x=stagesWithData, y=ratios, type="n", 
         xlab="number of stages or samples", 
         ylab="proportion of captured in sample k which are marked",
         main=sprintf("%s\n(extended Tanaka method)", tag))
    points(stagesWithData, ratios, pch=21, 
           col=ptCol, bg=ptCol, cex=ptCex)
    dev.off()
    #
    # Fit to stages
    if (0 == intercept)
    {
      fitToStages<- lm(ratios ~ stagesWithData + 0, singular.ok=FALSE,
                       model=TRUE, qr=TRUE, x=TRUE, y=TRUE, method="qr")
    } else
    {
      fitToStages<- lm(ratios ~ stagesWithData + 1, 
                       singular.ok=FALSE,
                       model=TRUE, qr=TRUE, x=TRUE, y=TRUE, method="qr")
    }
    #
    if (0 == intercept)
    {
      segmentingStages<- segmented(obj=fitToStages, 
#                            seg.Z, # missing because only univariate covariate
                             psi=NA, 
#                            npsi,  # ignored because psi specified
                             control=ctrl,
                             model=TRUE, 
                             keep.class=TRUE)  
    } else
    {
      candidateBreaks<- round(hdquantile(stagesWithData, probs=c(0.30, 0.70), names=FALSE))
      segmentingStages<- segmented(obj=fitToStages, 
                             seg.Z=~stagesWithData, 
                             psi=candidateBreaks,
#                            npsi,  # ignored because psi specified
                             control=ctrl,
                             model=TRUE, 
                             keep.class=TRUE)  
    
    }
#   confidenceSegmentingStages<- confint(object=segmentingStages, 
    confidenceSegmentingStages<- MyConfintSegmented(object=segmentingStages, 
#                                        parm, # omit
                                         level=conf, 
#                                        method=c("delta", "score", "gradient"),
                                         method=method,
                                         rev.sgn=FALSE, 
                                         digits=4,
                                         var.diff=FALSE, 
                                         is=FALSE
                                      )
    cat("Confidence intervals report from segmenting of stages...\n")
    print(confidenceSegmentingStages)
    print( summary(object=segmentingStages, short=FALSE, var.diff=FALSE, 
                   p.df="K", signif.stars=TRUE) )
    #
    cat("Plotting segmentation of stages ...\n")
    tSNow<- timingStamp()
    pdf(onefile=TRUE, title=sprintf("%s confidence by stages, %s", tag, tSNow),
        file=constructFilenameFrom(root=sprintf("%s-confidence-by-stages", ftag), 
                                   suffix=".pdf", stamp=tSNow)$full,
        width=22, height=17, pointsize=21)
    plot(x=segmentingStages, 
#        term, # omitted
         add=FALSE, res=TRUE, conf.level=conf, interc=FALSE,
#        link=TRUE, # ignored
         res.col="green", rev.sgn=FALSE, const=0, shade=TRUE, rug=TRUE,
         dens.rug=TRUE, dens.col = grey(0.8), transf=I, isV=FALSE, is=FALSE,
         var.diff=FALSE, p.df="K", .vcov=NULL)
    dev.off()
    #
    popsBySegment<- 
       popEstimatesPerSegmentFrom(nS=nS, 
                                  ratios=ratios, 
                                  intercept=intercept,
                                  previouslyAllocated=previouslyAllocated, 
                                  segmentingStages=segmentingStages, 
                                  tag=tag, 
                                  ftag=ftag,
                                  ptCol=ptCol,
                                  ptCex=ptCex,
                                  eliminateSmallSegments=eliminateSmallSegments
                               )
  }
  #
  #
  return(list(N.hat=fO$solution, 
              N.hat.interval=c(fO$solution.low, fO$solution.high),
              fit.global=fitOverall, extend=extend, segmented=segmenting,
              segmented.confidence=confidenceSegmenting,
              segmented.stages=segmentingStages, 
              segmented.stages.confidence=confidenceSegmentingStages,
              population.by.segments=popsBySegment
           ))
}

