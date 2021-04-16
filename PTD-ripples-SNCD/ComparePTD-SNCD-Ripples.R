# ComparePTD-SNCD-Ripples.R
#
# Compare PTD and SNCD on ripple cases.
# Jan Galkowski, bayesianlogic.1@gmail.com, 24th February 2019.
# Last changed 24th February 2019.
# 

library(Hmisc)
library(random)
library(Matrix)
library(shapes)
library(gtools)
library(corrgram)
source("plottableSVG.R")


stopifnot( file.exists("RippleCases.data") )

load("RippleCases.data")

Cases.df<- as.data.frame(t(Cases[2:17,]), stringsAsFactors=FALSE)


imputeMissingInSeries<- function(SERIES, nonnegative=TRUE)
{
  if (any(is.na(SERIES)))
  {
    cat("-divprocs- : Some missing values ...\n")
    N<- ncol(SERIES)
    I.lm<- 1:nrow(SERIES)
    for (j in 1:N)
    {
      if (any(is.na(SERIES[,j])))
      {
        cat(sprintf("Series %.0f has NAs:\n", j))
        print(SERIES[,j], digits=5)
        cat("Imputing by spline interpolation ...\n")
        missing<- which(is.na(SERIES[,j]))
        if (nonnegative)
        {
          splined<- interpSpline(obj1=I.lm, obj2=sqrt(SERIES[,j]), bSpline=TRUE, ord=4L, na.action=na.omit, sparse=TRUE)
          P<- predict(object=splined, x=missing, deriv=0)
          V<- P$y
          V2<- V*V
        } else
        {
          splined<- interpSpline(obj1=I.lm, obj2=SERIES[,j], bSpline=TRUE, ord=4L, na.action=na.omit, sparse=TRUE)
          P<- predict(object=splined, x=missing, deriv=0)
          V2- P$y
          V2<- V*V
        }
        SERIES[missing,j]<- V2
        cat(sprintf("Repaired series %.0f:\n", j))
        print(SERIES[,j], digits=5)
      }
    }
  }
  return(SERIES)
}


numericToStringForCompression<- function(x, y)
{
  n.x<- length(x)
  n.y<- length(y)
# (This the default number of bins for binning from the sm package, but there are
#  two vectors here, and they need a common number of bins.)
  nb<-  max(round(log(n.x)/log(2)+1), round(log(n.y)/log(2)+1))
  Q<- unique(hdquantile(c(x,y), probs=seq(0,1,length.out=1+nb), names=FALSE))
  alphaBet<- c(letters, LETTERS, sapply(X=0:9, FUN=function(n) sprintf("-%.0f", n)))
  m<- length(Q) - 1
  stopifnot( m <= length(alphaBet) )
  if (1 >= m) 
  {
    warning(sprintf("m in numeric-to-string is %.0f, alphabet length %.0f.", m, length(alphaBet)))
    chx<- sample(alphaBet, n.x, replace=TRUE)
    chy<- sample(alphaBet, n.y, replace=TRUE)
  } else
  {
    stopifnot( (1 < m) && (m <= length(alphaBet)) )
    codes<- c("!", mapply(A=rev(alphaBet[1:m]), K=(1:m), 
                   FUN=function(A,K) Reduce(f=function(a,b) paste0(a,b,collapse=NULL), x=rep(A, (1+K)), init="", right=FALSE, accumulate=FALSE)))
    cx<- 1+unclass(cut(x, Q, labels=FALSE))
    cx[which(is.na(cx))]<- 1
    cy<- 1+unclass(cut(y, Q, labels=FALSE))
    cy[which(is.na(cy))]<- 1
    chx<- codes[cx]
    chy<- codes[cy]
  }
  return(list(x=chx, y=chy))
}


divc.NCD <- function(xGiven, yGiven, trans=function(x) x) 
{
  typCompr<- "xz"
  if (is.numeric(xGiven))
  {
    coding<- numericToStringForCompression(x=xGiven, y=yGiven)
    x<- coding$x
    y<- coding$y
  } else
  {
    stopifnot( is.character(xGiven) )
    stopifnot( is.character(yGiven) )
    x<- xGiven
    y<- yGiven
  }
  #
  xy<- c(x,y)
  yx<- c(y,x)
  zero<- length(memCompress("", type=typCompr))
  cx<- length(memCompress(x, type=typCompr)) - zero
  cy<- length(memCompress(y, type=typCompr)) - zero
  cxy<- length(memCompress(xy, type=typCompr)) - zero
  cyx<- length(memCompress(yx, type=typCompr)) - zero
  #
  # Symmetrized NCD of the above.
  mnxy<- min(cx, cy)
  mxxy<- max(cx, cy)
  ncd<- max(0, min(1, ( (cxy - mnxy) + (cyx - mnxy) ) / (2*mxxy) ) )
  #
  return(trans(ncd))
}

divs<- function(SERIES, period=25, nonnegative=TRUE)
{
  stopifnot( is.data.frame(SERIES) ) 
  N<- ncol(SERIES)
  #
  SERIES<- imputeMissingInSeries(SERIES, nonnegative=nonnegative)
  #
  divergences<- Matrix(0, N, N, dimnames=list(NULL, NULL))
  # Since logits are so common in inference, calculate those, too.
  logit.divergences<- Matrix(-Inf, N, N, dimnames=list(NULL, NULL))
  N1<- N-1
  for (i in (1:N1))
  {
    for (j in ((1+i):N))
    {
      d<- divc.NCD(xGiven=SERIES[,i], yGiven=SERIES[,j], trans=function(x) x)
      divergences[i,j]<- d
      divergences[j,i]<- d
      ld<- logit(d)
      logit.divergences[i,j]<- ld
      logit.divergences[j,i]<- ld
    }
    if (0 == (i%%period))
    {
      cat(sprintf("... did %.0f\n", i))
    }
  }
  stopifnot( !is.null(colnames(SERIES)) )
  colnames(divergences)<- colnames(SERIES)
  rownames(divergences)<- colnames(SERIES)
  colnames(logit.divergences)<- colnames(SERIES)
  rownames(logit.divergences)<- colnames(SERIES)
  #
  # Return Matrix objects, leaving conversion to a matrix, a  distance matrix, or a data
  # from to the consumer of the output. Can't anticipate that here. 
  return(list(divergences=divergences, logit.divergences=logit.divergences))
}

divprocs2<- function(SERIES, period=25, comparison=TRUE)
{
  stopifnot( is.data.frame(SERIES) ) 
  N.lm<- nrow(SERIES)
  I.lm<- 1:N.lm
  N<- ncol(SERIES)
  divergencesRho<- Matrix(0, N, N, dimnames=list(NULL, NULL))
  divergencesRhoRMS<- Matrix(0, N, N, dimnames=list(NULL, NULL))
  divergences.dF.RMS<- Matrix(0, N, N, dimnames=list(NULL, NULL))
  divergences.SNCD<- Matrix(0, N, N, dimnames=list(NULL, NULL))
  N1<- N-1
  for (i in (1:N1))
  {
    for (j in ((1+i):N))
    {
      A<- array(dim=c(N.lm, 2, 2))
      A[,1,1]<- I.lm
      A[,2,1]<- SERIES[,i]
      A[,1,2]<- I.lm
      A[,2,2]<- SERIES[,j]
      fpa<- procGPA(x=A, scale=TRUE, reflect=FALSE, eigen2d=FALSE, tol1=1e-5, tol2=1e-5, tangentcoords="residual",
                    proc.output=FALSE, distances=TRUE, pcaoutput=FALSE, alpha=0, affine=FALSE)
      divergencesRho[i,j]<- mean(fpa$rho)
      divergencesRhoRMS[i,j]<- fpa$rmsrho
      divergences.dF.RMS[i,j]<- fpa$rmsd1
      #
      divergencesRho[j,i]<- divergencesRho[i,j]
      divergencesRhoRMS[j,i]<- divergencesRhoRMS[i,j]
      divergences.dF.RMS[j,i]<- divergences.dF.RMS[i,j]
      #
      if (comparison)
      {
        S<- cbind(SERIES[,i], SERIES[,j])
        colnames(S)<- colnames(SERIES)[c(i,j)]
        S<- as.data.frame(S, stringsAsFactors=FALSE)
        d<- divs(SERIES=S, period=100)$divergences[1,2]
        divergences.SNCD[i,j]<- d
        divergences.SNCD[j,i]<- d
      }
    }
    if (0 == (i%%period))
    {
      cat(sprintf("... did %.0f\n", i))
    }
  }
  stopifnot( !is.null(colnames(SERIES)) )
  colnames(divergencesRho)<- colnames(SERIES)
  rownames(divergencesRho)<- colnames(SERIES)
  colnames(divergencesRhoRMS)<- colnames(SERIES)
  rownames(divergencesRhoRMS)<- colnames(SERIES)
  colnames(divergences.dF.RMS)<- colnames(SERIES)
  rownames(divergences.dF.RMS)<- colnames(SERIES)
  #
  if (comparison)
  {
    colnames(divergences.SNCD)<- colnames(SERIES)
    rownames(divergences.SNCD)<- colnames(SERIES)
  } else
  {
    divergences.SNCD<- NA
  }
  #
  # Return Matrix objects, leaving conversion to a matrix, a  distance matrix, or a data
  # from to the consumer of the output. Can't anticipate that here. 
  return(list(divergences=divergencesRho, divergences.rms=divergencesRhoRMS, dF.rms=divergences.dF.RMS, sncd=divergences.SNCD))
}


properAngle<- function(e1, e2)
{
  if (Re(e1) > 0)
  {
    a<- atan(e2/e1)
  } else
  {
    a<- atan(e2/e1) + pi
  }
  return(a)
}


evectorOrderingOf<- function(cmat)
{
  x.eigen<- eigen(cmat)
  x.eigen.v<- x.eigen$vectors
  e1<- x.eigen.v[, 1]
  e2<- x.eigen.v[, 2]
  angles<- as.vector(mapply(e1=x.eigen.v[,1], e2=x.eigen.v[,2], 
                            FUN=properAngle))
  ord<- order(angles)
  evectorAngles<- 180*angles[ord]/pi
  cmatReturn<- cmat[ord, ord]
  return(list(reordered.correlation.matrix=cmatReturn, ordering=ord, is.posdef=all(0 < x.eigen$values), angles=evectorAngles))
}


shapeCorrgram<- function(sc.cor, tag, kind="Corrgram")
{
  #  D. J. Murdoch & E. D. Chow (1996) A Graphical Display of Large Correlation
  # Matrices, The American Statistician, 50:2, 178-180
  # https://doi.org/10.1080/00031305.1996.10474371
  #
  # M. Friendly, E. Kwan, Effect ordering for data displays, 
  # Computational Statistics & Data Analysis 43 (2003) 509-539
  #
  # Michael Friendly (2002) Corrgrams, The American Statistician, 56:4, 316-324,
  # https://doi.org/10.1198/000313002533
  #
  # Calculate eigenvectors, keeping two principal ones
  m<- nrow(sc.cor)
  #
  eOrdering<- evectorOrderingOf(sc.cor)
  reordered.sc.cor<- eOrdering$reordered.correlation.matrix
  #
  AOE.sc.cor<- corrgram(reordered.sc.cor, type="corr", order=TRUE, lower.panel=panel.cor, upper.panel=panel.shade, 
                        main=sprintf("%s for %s", kind, tag), cex.main=2, cex.labels=14)  
  #
  return(list(aoe.cor=AOE.sc.cor, evector.ordered.cor=reordered.sc.cor, evector.angles=eOrdering$angles))
}

adjustLabels<- function(x) 
{
  parts<- strsplit(x=x, split=".", fixed=TRUE)[[1]]
  p1<- parts[1]
  p2<- parts[2]
  return( sprintf("(%.1f,%.1f)", as.numeric(p1)/10, as.numeric(p2)/10) )
}



if (TRUE)
{
  Comparison<- divprocs2(SERIES=Cases.df, comparison=TRUE)
  labelings<- sapply(X=colnames(Comparison$divergences), FUN=adjustLabels)
  colnames(Comparison$divergences)<- labelings
  rownames(Comparison$divergences)<- labelings
  colnames(Comparison$sncd)<- labelings
  rownames(Comparison$sncd)<- labelings
}


#print(Comparison$divergences, digits=4)


if (FALSE)
{
  xtable.PTDs<- xtable(as.matrix(100*Comparison$divergences), caption="100X Procrustes tangent distances for the 16 ripple cases", digits=4)
  print((xtable.PTDs), type="latex", rotate.rownames=FALSE,rotate.colnames=FALSE)

  xtable.SNCDs<- xtable(as.matrix(Comparison$sncd), caption="SNCDs for the 16 ripple cases", digits=2)
  print((xtable.SNCDs), type="latex", rotate.rownames=FALSE,rotate.colnames=FALSE)
}

ttl<-"corrgram-PTD-16cases-100X"
fnx<- constructFilenameFrom(root=ttl, suffix=".pdf")$full
pdf(onefile=TRUE, title=ttl, file=fnx, width=98, height=round(98/2), pointsize=8)
corrgPTD<- shapeCorrgram(sc.cor=100*as.matrix(Comparison$divergences), 
                         kind="Distances",
                         tag="100X Procrustes tangent distances for the 16 ripple cases")
dev.off()

ttl<-"corrgram-SNCD-16cases"
fnx<- constructFilenameFrom(root=ttl, suffix=".pdf")$full
pdf(onefile=TRUE, title=ttl, file=fnx, width=98, height=round(98/2), pointsize=8)
corrgPTD<- shapeCorrgram(sc.cor=as.matrix(Comparison$sncd), kind="Divergences", tag="SNCDs for the 16 ripple cases")
dev.off()







