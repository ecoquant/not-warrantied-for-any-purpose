###############################################################################################
# This section of code up to the next group of #### marks is owned and copyrighted by         #
# Jan Galkowski, who has placed it in the public domain, available for any purpose by anyone. #
###############################################################################################
# Last changed 8th September 2018, bayesianlogic.1@gmail.com

library(arsenal) ## for `compare' function

fastMeansFastMomentsPrecondition<- function(X, filtering=I)
{
# # (From Jan Galkowski, ``Fast means, fast moments'', 1984,
# #  IBM Federal Systems Division, Owego, NY. Released into the public domain, 1994.)
# See https://667-per-cm.net/2018/09/06/fast-means-fast-moments-originally-devised-1984/
  stopifnot( is.matrix(X) )
  stopifnot( is.function(filtering) )
  M<- nrow(X)
  stopifnot( M == ncol(X) )
  if (compare(filtering,I)$equal)
  {
    AX1<- apply(X=X, MARGIN=2, FUN=cumsum)
  } else
  {
    # (For instance, filtering == function(v) max(v,0)
    AX1<- apply(X=X, MARGIN=2, FUN=function(column) cumsum(sapply(X=column, FUN=filtering)))
  }
  AX12<- t(apply(X=AX1, MARGIN=1, FUN=cumsum))
  return(AX12)
}


fastMeansFastMomentsBlock<- function(P, iUL, jUL, iLR, jLR)
{
# # (From Jan Galkowski, ``Fast means, fast moments'', 1984,
# #  IBM Federal Systems Division, Owego, NY. Released into the public domain, 1994.)
# See https://667-per-cm.net/2018/09/06/fast-means-fast-moments-originally-devised-1984/
#
# P is the preconditioned AX12 from above.
#
  stopifnot( is.matrix(P) )
  M<- nrow(P)
  stopifnot( M == ncol(P) )
  stopifnot( (1 <= iUL) && (iUL <= M) )
  stopifnot( (1 <= jUL) && (jUL <= M) )
  stopifnot( (1 <= iLR) && (iLR <= M) )
  stopifnot( (1 <= jLR) && (jLR <= M) )
#
  iUL1<- iUL-1
  jUL1<- jUL-1
  iLR1<- iLR-1
  jLR1<- jLR-1
  #
  if (0 == iUL1)
  {
    W.AL<- 0
    W.A<- 0
    if (0 == jUL1)
    {
      W.L<- 0
    } else
    {
      W.L<- P[iLR,jUL1]
    }
  } else if (0 == jUL1)
  {
    W.AL<- 0
    W.L<- 0
    if (0 == iUL1)
    {
      W.A<- 0
    } else
    {
      W.A<- P[iUL1,jLR]
    }
  } else 
  {
    W.AL<- P[iUL1,jUL1]
    W.A<- P[iUL1,jLR]
    W.L<- P[iLR,jUL1]
  }
  #
  W<- P[iLR,jLR] + W.AL - W.A - W.L
  #
  return(W)
}

# Self-test FMFM ...
cat("Fast means, fast moments self-test ...\n")
Z<- matrix(round(runif(100, min=1, max=100)), 10, 10)
Z.P<- fastMeansFastMomentsPrecondition(Z)
stopifnot( sum(Z[1:4,1:5]) == fastMeansFastMomentsBlock(Z.P, 1, 1, 4, 5) )
stopifnot( sum(Z[8:10, 8:9]) == fastMeansFastMomentsBlock(Z.P, 8, 8, 10, 9) )
stopifnot( sum(Z[4:7, 3:5]) == fastMeansFastMomentsBlock(Z.P, 4, 3, 7, 5) )
rm(list=c("Z", "Z.P"))
cat("... Self-test completed.\n")

###############################################################################################
# End of public domain code.                                                                  #
###############################################################################################

randomizeSeed<- function()
{
  #set.seed(31415)
  # Futz with the random seed
  E<- proc.time()["elapsed"]
  names(E)<- NULL
  rf<- E - trunc(E)
  set.seed(round(10000*rf))
# rm(list=c("E", "rf"))
  return( sample.int(2000000, size=sample.int(2000, size=1), replace=TRUE)[1] )
}

wonkyRandom<- randomizeSeed()

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

palette<- function(col, border = "light gray", ...)
{
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
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

is.df.numeric<- function(df)
{
  stopifnot( is.data.frame(df) )
  m<- ncol(df)
  nope<- FALSE
  for (j in (1:m))
  {
    if (!(typeof(df[,j]) %in% c("integer", "double")))
    {
      nope<- TRUE
      break
    }
  }
  return(!nope)
}


dfToMatrix<- function(df)
{
  stopifnot( is.data.frame(df) )
  stopifnot( is.df.numeric(df) )
  CN<- colnames(df)
  RN<- rownames(df)
  m<- ncol(df)
  n<- nrow(df)
  V<- unlist(df)
  names(V)<- NULL
  A<- Matrix(V, nrow=n, ncol=m)
  rownames(A)<- RN
  colnames(A)<- CN
  return(A)
}

# addalpha()
addAlpha <- function(colors, alpha=1.0) {
  # Transparency
  r <- col2rgb(colors, alpha=TRUE)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

weightToAlpha<- function(w, m)
{
  stopifnot( 1 < m )
  return((m*w - 1)/(m - 1))
}

createPDF<- function(tag1="", tag2="", tag3="", grp="", family=NULL)
{
  E<- round(100*proc.time()["elapsed"])
  if (is.null(family))
  {
    pdf(onefile=FALSE, 
        title="",
        family="Helvetica", paper="special", bg="white", compress=TRUE, width=14, height=11)
  } else
  {
    pdf(onefile=TRUE, 
        title=sprintf("%s %s: %-.0f %s", tag1, grp, E, tag3),
        file=sprintf("%s-%s-%-.0f.pdf", tag2, grp, E), 
        family=family, paper="special", bg="white", compress=TRUE, width=14, height=11)
  }
}



##################################################################################
# From: http://gifi.stat.ucla.edu/apl/
#  Ref: https://rpubs.com/deleeuw/158476

# utilities below

is.wholenumber<-
# This is a standard R predicate
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

first <- function(x) {
  return(x[1])
}

butFirst <- function(x) {
  return(x[-1])
}

last <- function(x) {
  return(x[length(x)])
}

butLast <- function(x) {
  return(x[-length(x)])
}

unit <- function(i, n) {
  as.numeric(i == (1:n))
}

#  representation -- will be overwritten by the C version

aplEncode <- function(rrr, base) {
  b <- c(1, butLast(cumprod(base)))
  r <- rep(0, length(b))
  s <- rrr - 1
  for (j in length(base):1) {
    r[j] <- s %/% b[j]
    s <- s - r[j] * b[j]
  }
  return(1 + r)
}

#  base value -- will be overwritten by the C version

aplDecode <- function(ind, base) {
  b <- c(1, butLast(cumprod(base)))
  return(1 + sum(b * (ind - 1)))
}


#  expand vector

aplEXV <- function(x, y) {
  z <- rep(0, length(y))
  m <- which(y == TRUE)
  if (length(m) != length(x))
    stop("Incorrect vector length")
  z[m] <- x
  return(z)
}

#  expand

aplExpand <- function(x, y, axis = 1) {
  if (is.vector(x))
    return(aplEXV(x, y))
  d <- dim(x)
  m <- which(y == TRUE)
  n <- length (y)
  e <- d
  e[axis] <- n
  if (length(m) != d[axis])
    stop("Incorrect dimension length")
  z <- array(0, e)
  for (i in 1:prod(d)) {
    k <- aplEncode (i, d)
    k[axis] <- m[k[axis]]
    z[aplDecode(k, e)] <- x[i]
  }
  return (z)
}

flipTopBottom<- function(A)
{
  stopifnot( is.matrix(A) )
  F<- rev(seq_len(nrow(A)))
  return( A[F,] )
}

flipLeftRight<- function(A)
{
  stopifnot( is.matrix(A) )
  F<- rev(seq_len(ncol(A)))
  return( A[,F] )
}

##################################################################################

# Two samples from numbered population enumerated N, each sized m.
# Combined lower tail probability of getting d or less duplicates from N.
# Jan Galkowski, jgalkows@akamai.com, 22nd December 2017.

pr.D.at.a.time<- function(N, m, d, repetitions=10000)
{
  #
  # Do the calculation both in closed form and by resampling simulation.
  #
  # The simulation simply performs the sampling twice, and calculates how many in common.
  R<- rep(NA, repetitions)
  for (k in (1:repetitions))
  {
    x1<- sample.int(N, m, replace=FALSE)
    x2<- sample.int(N, m, replace=FALSE)
    inCommon<- length(intersect(x1,x2))
    R[k]<- inCommon
  }
  #
  # Theoretical, for comparison:
  #
  # {N}\choose{d} \frac{{N-d}\choose{m-d}}{{N}\choose{m}} \frac{{N-m}\choose{m-d}}{{N}\choose{m}}
  #
  # This is that the d duplicates have {N}\choose{d} ways of being picked. Then,
  # there are {N-d}\choose{m-d} ways from a population of {N}\choose{m} samples for the first sample.
  # There are {N-m}\choose{m-d} ways from a population of {N}\choose{m} samples for the seconds sample.
  #
  # d varies from zero to the given d, and the probabilities are summed.
  #
  log.p<- function(N,m,d) lchoose(N,d) + (lchoose((N-d),(m-d))-lchoose(N,m)) + (lchoose((N-m),(m-d))-lchoose(N,m))
  p<- sum(sapply(X=0:d, FUN=function(j) exp(log.p(N,m,j))))
  cat(sprintf("Theoretical for comparison: %.7f\n", p))
  return(length(which(R<=d))/length(R))
}


dowFrom<- function(dtob)
{
  aDOW<- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
  d<- dayOfWeek(as.timeDate(dtob))
  return( match(d, aDOW, nomatch=0) )
}

isWeekend<- function(dtob)
{
  dow<- dowFrom(dtob)
  stopifnot( dow > 0 )
  return(as.numeric((dow == 6)||(dow == 7)))
}



splitIntoEquilengthed<- function(X, w)
{
  stopifnot( is.vector(X) )
  stopifnot( is.numeric(w) && (0 < w) && (w <= length(X)) )
  # Split a vector X into parts each having a length of w, excepting
  # conceivably the past part, which could have a length less than w.
  S<- split(X,as.numeric(gl(length(X),w,length(X)))) 
  names(S)<- NULL
  return(S)
}

splitIntoRoughlyEquilengthed<- function(X, w)
{
  stopifnot( is.vector(X) )
  stopifnot( is.numeric(w) && (0 < w) )
  w<- min(w, length(X))
  # Split a vector X into parts each having a length of w, excepting
  # conceivably the past part, which could have a length less than w.
  S<- list()
  L<- length(X)
  j<- 1
  k<- w
  i<- 1
  while(k < L)
  {
    S[[i]]<- X[j:k]
    j<- w + j
    k<- min(L, (w + k))
    i<- 1 + i
  }
  if (k == L)
  {
    S[[i]]<- X[j:k]
  }
  names(S)<- NULL
  return(S)
}
