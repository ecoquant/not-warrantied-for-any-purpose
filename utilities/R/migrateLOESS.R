#
# migrateLOESS.R : Migrating a time-series using a LOESS smoother.
# Jan Galkowski, bayesianlogic.1@gmail.com, 27 January 2019.
# Last changed 27 January 2019.

library(random)
library(splines)

randomizeSeed<- function(external=FALSE)
{
  #set.seed(31415)
  # Futz with the random seed
  if (!external)
  {
    E<- proc.time()["elapsed"]
    names(E)<- NULL
    rf<- E - trunc(E)
    set.seed(round(10000*rf))
  } else
  {
    set.seed(randomNumbers(n=1, min=1, max=10000, col=1, base=10, check=TRUE))
  }
  return( sample.int(2000000, size=sample.int(2000, size=1), replace=TRUE)[1] )
}

wonkyRandom<- randomizeSeed(external=TRUE)

findLoessSpan<- function(x, y, spanRange=c(0.5, 0.7), nb=15, report=FALSE, cell=0.2, iterations=5, mGiven=NULL, tag="nil")
{
#
# (Adapted from "Fit data points with LOESS + cross validation", 
#  https://rpubs.com/mengxu/loess_cv.)
#
  stopifnot( is.vector(spanRange) )
  stopifnot( 2 == length(spanRange))
  stopifnot( all( 0 < spanRange ) )
  stopifnot( is.vector(x) )
  stopifnot( is.vector(y) )
  N<- length(x)
  stopifnot( N == length(y) )
  stopifnot( 5 < N )
  #
  if (is.null(mGiven))
  {
    m<- round(N/10)
  } else
  {
    m<- mGiven
    stopifnot( is.wholenumber(m) && ((N/1000) < m) && (5 < m) )
  }
  #
  if (report)
  {
    cat(sprintf("Seeking LOESS span for '%s', N = %.0f, %.0f bootstrap replicas on range %.3f to %.3f, %.0f length sequences.\n", tag, N, 
                nb, spanRange[1], spanRange[2], m))
  }
  #
  ctrl<- loess.control(surface="direct", statistics="approximate", trace.hat="approximate",
                       cell=cell, iterations=iterations, iterTrace=FALSE)
  localDf<- data.frame(x=x, y=y, stringsAsFactors=FALSE)
  wonkyRandom<- randomizeSeed(external=FALSE)
  #
  MAEs<- rep(NA, nb)
  spans<- rep(NA, nb)
  for (k in (1:nb))
  {
    p<- sample.int((1+N-m), 1)
    x.s<- x[p:(p+m-1)]
    y.s<- y[p:(p+m-1)]
    spans[k]<- runif(1, min=spanRange[1], max=spanRange[2])
    localDf<- data.frame(x=x.s, y=y.s, stringsAsFactors=FALSE)
    trialL<- loess(y ~ x, data=localDf, degree=2, parametric=FALSE, drop.square=FALSE, normalize=FALSE,
                   family="symmetric", method="loess", control=ctrl, span=spans[k])
    # Use mean absolute error.
    V<- predict(object=trialL, newdata=localDf, se=FALSE)
    MAEs[k]<- mean(abs(V - localDf$y))
    if (report && (0 == (k%%50)))
    {
      nz<- which(0 < MAEs[1:k])
      span<- weighted.mean(x=spans[nz], w=1/MAEs[nz])
      cat(sprintf("'%s', at replica %.0f. Interim span %.3f.\n", tag, k, span))
    }
  }
  nz<- which(0 < MAEs)
  span<- weighted.mean(x=spans[nz], w=1/MAEs[nz])
  if (report)
  {
    cat(sprintf("Found span %.3f for tag '%s'\n", span, tag))
  }
  return(span)
}


migrateLOESS<- function(givenJumble, ticks, tag="nil", cell=0.2, iterations=5, nb=20, 
                        mGiven=4*24*2, report=TRUE, useSpan=NA)
{
  #
  # Migrates a time series onto a common time grid. Assumes values are non-negative, 
  # such as counts. Function can be easily changed to relax that assumption, but
  # doesn't support it at present.
  #
  # (useSpan == 0.523 found for some datasets)
  #
  if (report)
  {
    cat(sprintf("Beginning migration for '%s' ...\n", tag))
  }
  #
  if (is.matrix(givenJumble) )
  {
    sorted<- sort(x=givenJumble[,1], decreasing=FALSE, index.return=TRUE, method="shell")
    jumble<- givenJumble[sorted$ix,]
    # Note how the square root of the ordinate is splined.
    splined<- interpSpline(obj1=jumble[,1], obj2=sqrt(jumble[,2]), bSpline=TRUE, ord=4L, na.action=na.fail, sparse=TRUE)
    P<- predict(object=splined, x=ticks, deriv=0)
    Vsplined<- sapply(X=P$y, FUN=function(v) max(v,1))
    #
    acceptableRange<- c(0.50, 0.55)
    ctrl<- loess.control(surface="direct", statistics="approximate", trace.hat="approximate", cell=cell, 
                         iterations=(3+iterations), iterTrace=FALSE)
    localDf<- data.frame(x=ticks, y=Vsplined, stringsAsFactors=FALSE)
    if (is.na(useSpan))
    {
      span<- findLoessSpan(x=ticks, y=Vsplined, spanRange=acceptableRange, report=report, nb=nb, 
                           tag=tag, cell=cell, iterations=iterations, mGiven=mGiven)
    } else
    {
      stopifnot( is.numeric(useSpan) && (acceptableRange[1] < useSpan) && (useSpan < acceptableRange[2]))
      span<- useSpan
    }
    Vsmoothed<- loess(y ~ x, data=localDf, degree=2, parametric=FALSE, drop.square=FALSE, normalize=FALSE,
                      family="symmetric", method="loess", control=ctrl, span=span)
    V<- predict(object=Vsmoothed, newdata=localDf)
    #
    # The interpolated and smoothed square root of the ordinate is then squared to assure positivity.
    V2<- V*V
    names(V2)<- sapply(X=P$x, FUN=function(t) sprintf("%-.0f", t))
  } else
  {
    cat("Series only has one point.\n")
    browser()
    V2<- NA
  }
  return(V2)
}
