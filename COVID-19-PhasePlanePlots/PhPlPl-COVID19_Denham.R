# Phase plane plots of COVID-19 deaths. 
#
# PhPlPl-COVI19_Denham.R
#
# Jan Galkowski, 9th June 2020.
# Last changed 30th October 2020.

library(random)
#library(sfsmisc)
library(matrixcalc)
#library(data.table)
library(tictoc)
library(akima)
library(spectral)
library(pspline)
library(EnvStats)
library(plotrix)

source("plottableSVG.R")

#setwd("c:/Users/Jan/Documents/Westwood Statistical Studios/COVID-19-2020/")

randomizeSeed<- function(external=FALSE)
{
  #set.seed(31415)
  # Futz with the random seed
  if (!external)
  {
    rf<-  (as.numeric(Sys.time())*100000)%%99989
    s<- round(rf)
    set.seed(s)
    cat(sprintf("Seed chosen is: %.0f, internal.\n", s))
  } else
  {
    set.seed(randomNumbers(n=1, min=1, max=10000, col=1, base=10, check=TRUE))
  }
  return( sample.int(2000000, size=sample.int(2000, size=1), replace=TRUE)[1] )
}

scrubFunctions<- function()
{
  # To be used as:
  #
  #  (scrubFunctions())()
  #
  X<- objects(envir=globalenv())
  removed<- c()
  for (x in X)
  {
    if ((x != "scrubFunctions") && is.function(get(x, envir=globalenv())))
    {
      rm(list=c(x), envir=globalenv())
      removed<- c(removed, x)
    }
  }
  return(function (){ rm(list=c("scrubFunctions"), envir=globalenv()) ; removed})
}


pause.with.message<- function (message) 
{
  # Pause with a message between plots
  cat(message)
  cat("\n")
  cat("Paused. Press <Enter> to continue...")
  readline()
  invisible()
}


wonkyRandom<- randomizeSeed(external=TRUE)

significantDifference<- function(s, delta)
{
  stopifnot(0 < delta)
  p<- 1+which(delta<diff(s))[1]
  return(p)
}


check<- function(x,y)
{
  stopifnot( is.vector(x) )
  stopifnot( is.vector(y) )
  n<- length(y)
  stopifnot( length(x) == n )
  stopifnot( all(!is.na(x)) )
  stopifnot( all(!is.na(y)) )
  return(n)
}

splinedDerivativesPreface<- function(x, y, nderiv=2)
{
  n<- check(x,y)
  sorted<- sort.int(x=x, index.return=TRUE, method="shell", decreasing=FALSE)
  x<- x[sorted$ix]
  y<- y[sorted$ix]
  norder<- nderiv+3
  P<- smooth.Pspline(x, y, w=rep(1, n), norder=norder, method=3)
  return(P)
}

splinedDerivatives<- function(x, y, nderiv=2)
{
  P<- splinedDerivativesPreface(x=x, y=y, nderiv=nderiv)
  yHat<- as.vector(predict(object=P, xarg=x, nderiv=0))
  yHatDot<- as.vector(predict(object=P, xarg=x, nderiv=1))
  yHatDotDot<- as.vector(predict(object=P, xarg=x, nderiv=2))
  return(list(x=x, y=yHat, ydot=yHatDot, ydotdot=yHatDotDot))
}

predictionIntervalsFrom<- function(x, alpha=0.9)
{
# predIntNpar(x, k = m, m = 1, lpl.rank = ifelse(pi.type == "upper", 0, 1),
#             n.plus.one.minus.upl.rank = ifelse(pi.type == "lower", 0, 1),
#             lb = -Inf, ub = Inf, pi.type = "two-sided")
#
  stopifnot( is.vector(x) )
  n<- length(x)
  estimateOb<- predIntNpar(x=x, k=round(alpha*1000), m=1000, 
                         lpl.rank=round(n/8), n.plus.one.minus.upl.rank=((n+1)-round(7*n/8)),
                         lb=1, ub=Inf, pi.type="two-sided")
  intervalOb<- estimateOb$interval
  bounds<- intervalOb$limits
  ranks<- intervalOb$limit.ranks
  return(bounds)
}

oneBooted<- function(x, residualsCentered, predictedSeries, nderiv)
{
  stopifnot( is.vector(residualsCentered) )
  stopifnot( is.vector(predictedSeries) )
  n<- length(residualsCentered)
  stopifnot( n == length(predictedSeries) )
  #
  eStar<- sample(residualsCentered, size=n, replace=TRUE)
  yStar<- predictedSeries + eStar
  PStar<- splinedDerivativesPreface(x=x, y=yStar, nderiv=nderiv)
  yHat<- as.vector(predict(object=PStar, xarg=x, nderiv=0))
  return(yHat)
}

predictionIntervalsByBootFrom<- function(x, y, alpha=0.9, debug=TRUE, R=10000)
{
  x<- as.vector(x)
  y<- as.vector(y)
  n<- check(x,y)
  nderiv<- 2
  stopifnot( (0 < alpha) && (alpha < 1) )
# stopifnot(0 == (n%%2))
  #
  # Fit to original data.
  yBase<- splinedDerivatives(x=x, y=y, nderiv=nderiv)
  # 
  # Obtain natural surrogates for derivatives in the data. Noisiness is okay.
  sdp1<- splinedDerivativesPreface(x=x[-n], y=diff(y), nderiv=4)
  yd<- as.vector(predict(sdp1, xarg=x, nderiv=0))
  sdp2<- splinedDerivativesPreface(x=x[c((1-n),(-n))], y=diff(diff(y)), nderiv=4)
  ydd<- as.vector(predict(sdp2, xarg=x, nderiv=0))
  #
  # Residuals.
  cat("Residuals.\n")
  tic("Residuals time")
  #
  r<- y - yBase$y
  rCentered<- r - mean(r, na.rm=TRUE)
  #
  rdot<- yd - yBase$ydot
  rdotCentered<- rdot - mean(rdot, na.rm=TRUE)
  rdotdot<- ydd - yBase$ydotdot
  rdotdotCentered<- rdotdot - mean(rdotdot, na.rm=TRUE)
  #
  toc(log=TRUE)
  #
  # Prediction bootstrapping for responses.
  #
  cat("Booting responses.\n")
  tic("Booting responses time")
  #
  bootedResponses<- matrix(0, nrow=R, ncol=n)
  for (b in (1:R))
  {
    bootedResponses[b,]<- oneBooted(x=x, residualsCentered=rCentered, predictedSeries=yBase$y, nderiv=nderiv)
  }
  #
  toc(log=TRUE)
  #
  # Prediction bootstrapping for first derivatives.
  #
  cat("Booting first derivatives.\n")
  tic("Booting first derivatives time")
  #
  bootedDots<- matrix(0, nrow=R, ncol=n)
  for (b in (1:R))
  {
    bootedDots[b,]<- oneBooted(x=x, residualsCentered=rdotCentered, predictedSeries=yBase$ydot, nderiv=nderiv)
  }
  #
  toc(log=TRUE)
  #
  # Prediction bootstrapping for second derivatives.
  #
  cat("Booting second derivatives.\n")
  tic("Booting second derivatives time")
  #
  bootedDotDots<- matrix(0, nrow=R, ncol=n)
  for (b in (1:R))
  {
    bootedDotDots[b,]<- oneBooted(x=x, residualsCentered=rdotdotCentered, predictedSeries=yBase$ydotdot, nderiv=nderiv)
  }
  #
  toc(log=TRUE)
  #
  # Prediction intervals.
  #
  cat("Calculating prediction intervals.\n")
  tic("Prediction intervals time")
  #
  piR<- t(apply(X=bootedResponses, MARGIN=2, FUN=function(column) predictionIntervalsFrom(x=column, alpha=alpha)))
  piD<- t(apply(X=bootedDots, MARGIN=2, FUN=function(column) predictionIntervalsFrom(x=column, alpha=alpha)))
  piDD<- t(apply(X=bootedDotDots, MARGIN=2, FUN=function(column) predictionIntervalsFrom(x=column, alpha=alpha)))
  #
  toc(log=TRUE)
  #
  return(list(x=yBase$x, y=yBase$y, ydot=yBase$ydot, ydotdot=yBase$ydotdot, 
              lo0=(piR[,1]), up0=(piR[,2]), 
              lo1=(piD[,1]), up1=(piD[,2]),
              lo2=(piDD[,1]), up2=(piDD[,2])))
}


estimateQuantityAndTwoDerivativeWithStandardErrors<- function(S.scaled, S.mean, S.sd, alpha=0.9)
{
  S.scaled<- as.vector(S.scaled)
  tocks<- 1:length(S.scaled)
# ########################################
  fit<- predictionIntervalsByBootFrom(x=tocks, y=S.scaled, alpha=alpha)
  #
  # ########################################
  # Package in matrices.
  smoothPack<- cbind(fit$lo0, fit$y, fit$up0)
  fstPack<- cbind(fit$lo1, fit$ydot, fit$up1)
  secPack<- cbind(fit$lo2, fit$ydotdot, fit$up2)
# ########################################
  smoothPack<- S.mean+S.sd*smoothPack
  ########################################
  return(list(smoothed=smoothPack, fst=fstPack, sec=secPack, model=fit))
}

drawVaryingWidthCurve<- function(x, y, ex, ey, col="blue")
{
  stopifnot( is.matrix(ex) )
  stopifnot( is.matrix(ey) )
  stopifnot( 2 == ncol(ex) )
  stopifnot( 2 == ncol(ey) )
  L<- nrow(ex)
  stopifnot( L == nrow(ey) )
  #
  major<- (ex[,2]-ex[,1])/2
  minor<- (ey[,2]-ey[,1])/2
  plotrix::draw.ellipse(x=x, y=y, a=major, b=minor, col=col, border=col)
}



ppp<- function(S, kth, what, tag="phase plane plots", ordinateLabel="cumulative count of deaths", abscissaLabel="date",
               fstLabel="rate of change in number of deaths", secLabel="acceleration in number of deaths",
               svg=TRUE, one="svg1", two="svg2", three="svg3", nCloudShape=100, internal=FALSE, debug=TRUE)
{
  rejectionQuantile<- 0.90
  pointSizeNA<- 10
  # Vector series S
  stopifnot( is.vector(S) )
  stopifnot( 0 < S[1] )
  #
  if (FALSE)
  {
    cat(sprintf("tag : '%s'\n", tag))
    cat(sprintf("what : '%s'\n", what))
    cat(sprintf("ordinate : '%s'\n", ordinateLabel))
    cat(sprintf("fst : '%s'\n", fstLabel))
    cat(sprintf("sec : '%s'\n", secLabel))
    cat(sprintf("one : '%s'\n", one))
    cat(sprintf("two : '%s'\n", two))
    cat(sprintf("three : '%s'\n", three))
  }
  #
  L<- length(S)
  ticks<- names(S)
  dateLabels<- sapply(X=ticks, FUN=function(s) substr(s, 2, nchar(s)))
  R<- 1:L
  #
  # ESTIMATION.
  #
  Smean<- mean(S, na.rm=TRUE)
  Ssd<- sd(S, na.rm=TRUE)
  Smedian<- median(S, na.rm=TRUE)
  S.scaled<- t(matrix(as.vector(unclass(scale(x=S, center=TRUE, scale=TRUE))), nrow=1, ncol=L))
  #
  dynamics<- estimateQuantityAndTwoDerivativeWithStandardErrors(S.scaled=S.scaled, S.mean=Smean, S.sd=Ssd, alpha=0.667)
  cat("Estimation done.\n")
  ########################################
  smoothed<- dynamics$smoothed[,2]
  fst<- dynamics$fst[,2]
  sec<- dynamics$sec[,2]
  peSmoothed<- dynamics$smoothed[,c(1,3)]
  peFst<- dynamics$fst[,c(1,3)]
  peSec<- dynamics$sec[,c(1,3)]
  ########################################
  #
  # PLOT CUMULATIVE QUANTITY VERSUS DAYS IN.
  #
   pluck5<- R[unique(c(seq(2,length(R), 5), length(R)))]
  #
 if (svg)
  {
    fx<- openSVG(root=one, width=24, height=18, antialias="gray")
  }
  plot(R[3:L], S[3:L], type="p", pch=21, cex=2, col="maroon", bg="darkred", main=sprintf("%s %s", what, tag), cex.main=3, xaxt="n", cex.lab=2, 
       ylab=ordinateLabel, xlab=abscissaLabel, family="Times New Roman", cex.axis=2)
  # Uncertainty. Note centered about S.estimated, not S.
  segments(x0=R[], y0=peSmoothed[,1], x1=R[], y1=peSmoothed[,2], col="lightblue", lty=1, lwd=100)
  #
  text(R[3:L], par("usr")[1], labels=dateLabels[2:L], srt=90, cex=2, col="black", font=2, pos=1, offset=2, xpd=TRUE, family="Times New Roman")
  abline(v=pluck5, lty=6, lwd=1, col="grey")
  # Using output of Student local linear trend model directly rather than smoothed.
  lines(R[3:L], smoothed[3:L], lwd=4, lty=1, col="navy")
  points(R[3:L], S[3:L], pch=21, cex=4, col="maroon", bg="darkred")
  if (svg)
  {
    closeSVG(fx)
  } else
  {
    pause.with.message("basis")
  }
  #
  # PLOT CUMULATIVE QUANTITY VERSUS FIRST DERIVATIVE OF QUANTITY.
  #
  if (svg)
  {
    fx<- openSVG(root=two, width=11, height=7, antialias="gray")
  }
  stopifnot( L == length(smoothed) )
  # Uncertainty.
  plot(fst, smoothed, type="l", lwd=2, lty=1, col="navy", main=tag, cex.main=1.5, 
       ylab=ordinateLabel, xlab=fstLabel, family="Times New Roman")
  #
  drawVaryingWidthCurve(x=fst, y=smoothed, ex=peFst, ey=peSmoothed, col="lightblue")
  #
  lines(fst, smoothed, lwd=2, lty=1, col="navy")
  #
  text(fst, smoothed, labels=dateLabels, 
       srt=0, cex=0.9, col="darkred", font=2, pos=1, offset=0.2, family="Times New Roman")
  points(fst, smoothed, pch=7, col="red", bg="red", cex=1, font=2)
  points(data.table::last(fst), data.table::last(smoothed), pch=8, cex=10, col="maroon", font=2)
  par(family="Times New Roman")
  legend("bottomright", border="black", 
         legend="See text for explanation of what uncertainty envelopes mean.",
         cex=0.7, text.font=2, text.col="black", box.lwd=1, box.lty=1)
  if (svg)
  {
    closeSVG(fx)
  } else
  {
    pause.with.message("velocity vs basis")
  }
  #
  # PLOT SECOND DERIVATIVE OF QUANTITY VERSUS FIRST DERIVATIVE OF QUANTITY.
  #
  if (svg)
  {
    fx<- openSVG(root=three, width=11, height=7, antialias="gray")
  }
  #
  stopifnot( L == length(sec) )
  plot(fst, sec, type="l", lwd=2, lty=1, col="navy", main=tag, cex.main=1.5,
       ylab=secLabel, xlab=fstLabel, family="Times New Roman")
  #
  drawVaryingWidthCurve(x=fst, y=sec, ex=peFst, ey=peSec, col="lightblue")
  #
  lines(fst, sec, lwd=2, lty=1, col="navy")
  #
  text(fst, sec, labels=dateLabels, 
       srt=0, cex=0.9, col="darkred", font=2, pos=1, offset=0.2, family="Times New Roman")
  points(fst, sec, pch=7, col="red", bg="red", cex=1, font=2)
  points(data.table::last(fst), data.table::last(sec), pch=8, cex=10, col="maroon", font=2)
  par(family="Times New Roman")
  legend("topleft", border="black", 
         legend="See text for explanation of what uncertainty envelopes mean.", 
         cex=0.7, text.font=2, text.col="black", box.lwd=1, box.lty=1)
  if (svg)
  {
    closeSVG(fx)
  } else
  {
    pause.with.message("velocity vs acceleration")
  }
}

