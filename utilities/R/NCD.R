# Symmetrized Normalized Comrpession Divergence.
# Jan Galkowski, 3rd December 2018. 
# Last changed 14th January 2019.

library(Matrix) # Data structure for large divergence matrices
library(random) # Source of the random function
library(Hmisc)  # Source of the hdquantile function
library(gtools) # Source of the logit function

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

wonkyRandom<- randomizeSeed(external=FALSE)

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

evenbins <- function(x, bin.count=10, order=TRUE) 
{
  # Function by Mr Flick, https://stackoverflow.com/users/2372064/mrflick
  # from https://stackoverflow.com/posts/24359943/revisions
  bin.size<- rep(length(x) %/% bin.count, bin.count)
  bin.size<- bin.size + ifelse(1:bin.count <= length(x) %% bin.count, 1, 0)
  bin<- rep(1:bin.count, bin.size)
  if (order) 
  {    
    bin <- bin[rank(x,ties.method="random")]
  }
  B<- unclass(factor(bin, levels=1:bin.count, ordered=order))
  attr(B, names(attributes(B)))<- NULL
  return(B)
}

dissVSTR<- function(VSTR, period=25, logitp=FALSE)
{
  stopifnot( is.vector(VSTR) ) 
  N<- length(VSTR)
  zero<- length(memCompress(""))
  ncdf<- function(cx, cy, cxy, cyx) { mnxy<- min(cx,cy) ; mxxy<- max(cx,cy) ; return( max(0, min(1, (cxy + cyx - 2*mnxy)/(2*mxxy) ))) }
  #
  CV<- sapply(X=VSTR, FUN=function(s) {length(memCompress(s)) - zero})
  if ((N>200) & (0 < period))
  {
    cat(sprintf("Preconditioning of %.0f items completed.\n", N))
  }
  #
  if (logitp)
  {
    dInitial<- -Inf
    trans<- logit
  } else
  {
    dInitial<- 0
    trans<- function(x) x
  }
  #
  divergences<- Matrix(dInitial, N, N, dimnames=list(NULL, NULL))
  #
  N1<- N-1
  for (i in (1:N1))
  {
    sx<- VSTR[i]
    cx<- CV[i]
    for (j in ((1+i):N))
    {
      sy<- VSTR[j]
      cy<- CV[j]
      sxy<- sprintf("%s%s", sx, sy)
      syx<- sprintf("%s%s", sy, sx)
      cxy<- length(memCompress(sxy)) - zero
      cyx<- length(memCompress(syx)) - zero
      d<- trans(ncdf(cx, cy, cxy, cyx))
      if (is.nan(d))
      {
        cat("NANs within VSTR. Inspection:\n")
        browser()
      }
      divergences[i,j]<- d
      divergences[j,i]<- d
    }
    if ((0 < period) && (200 < N) && (0 == (i%%period)))
    {
      cat(sprintf("... did %.0f\n", i))
    }
  }
  colnames(divergences)<- names(VSTR)
  rownames(divergences)<- names(VSTR)
  # Return a Matrix object, leaving conversion to a matrix, a  distance matrix, or a data
  # from to the consumer of the output. Can't anticipate that here.
  return(divergences)
}

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


compression.lengths<- function(xGiven, yGiven, type="xz")
{
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
  xx<- c(x,x)
  yy<<-c(y,y)
  xy<- c(x,y)
  yx<- c(y,x)
  stopifnot( is.character(xx) )
  stopifnot( is.character(yy) )
  stopifnot( is.character(xy) )
  stopifnot( is.character(yx) )
  zero<- length(memCompress("", type=type))
  cx<- length(memCompress(x, type=type)) - zero
  cy<- length(memCompress(y, type=type)) - zero
  cxx<- length(memCompress(xx, type=type)) - zero
  cyy<- length(memCompress(yy, type=type)) - zero
  cxy<- length(memCompress(xy, type=type)) - zero
  cyx<- length(memCompress(yx, type=type)) - zero
  return(list(cx=cx, cy=cy, cxx=cxx, cyy=cyy, cxy=cxy, cyx=cyx, csymmetric=(cxy+cyx)/2))
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


