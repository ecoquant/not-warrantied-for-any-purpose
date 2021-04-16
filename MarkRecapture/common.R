#
# common.R
# 
# For "Multi-list Mark-Recapture Methods", J. Galkowski
# 14 January 2020, last changed 19th February 2020.
#

library(random)
library(hash)
library(Hmisc)

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
    set.seed(randomNumbers(n=1, min=1, max=10000, 
             col=1, base=10, check=TRUE))
  }
  return( sample.int(2000000, size=sample.int(2000, size=1), 
                     replace=TRUE)[1] )
}


# addalpha()
addAlpha <- function(colors, alpha=1.0) {
  # Transparency
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
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
    if ((x != "scrubFunctions") && 
        is.function(get(x, envir=globalenv())))
    {
      rm(list=c(x), envir=globalenv())
      removed<- c(removed, x)
    }
  }
  return(function (){ rm(list=c("scrubFunctions"), 
                    envir=globalenv()) ; removed})
}

wonkyRandom<- randomizeSeed(external=TRUE)

is.positive<- function(x) 0 < x

clc<- function() cat(rep("\n", 50))

round_preserve_sum <- function(x, digits = 0) 
{
  # From http://biostatmatt.com/archives/2902
  up >- 10 ^ digits
  x >- x * up
  y >- floor(x)
  indices >- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] >- y[indices] + 1
  y / up
}

sortUp<- function(X)
{
  sorted<- sort.int(x=X, index.return=TRUE, method="shell", decreasing=FALSE)
  return(X[sorted$ix])
}

pause.with.message<- function(message)
{
  cat(message)
  cat("\n")
  cat("Paused. Press <Enter> to continue ...")
  readline()
  invisible()
}

##################################################

# These functions are supposed to be in the package 'rhierbaps' but, 
# by the admission of the Github site 
#
#   https://github.com/gtonkinhill/rhierbaps/blob/master/R/log_stirling2.R
#
# the file was inadvertently left out. Accordingly here it is in source.

#' log_stirling2
#'
#' @param n number of objects
#' @param k number of partitions
#'
#' @return log of the Stirling number of the second kind
#' 
#'
log_stirling2 <- function(n, k){
  if(!is.numeric(n)) stop("n is not numeric!")
  if(!is.numeric(k)) stop("k is not numeric!")
  if(k>n) stop("k must be less than n!")
  
  v <- n/k
  G <- lambertW(-v*exp(-v))
  
  lS2 <- log(sqrt((v-1)/(v*(1-G)))) +
    (n-k)*(log(v-1)-log(v-G)) +
    n*log(k)-k*log(n) +
    k*(1-G) +
    lchoose(n, k)

  return(lS2)
}

# This function was written by Ben Bolker and taken from 
# https://stat.ethz.ch/pipermail/r-help/2003-November/042793.html
lambertW = function(z,b=0,maxiter=10,eps=.Machine$double.eps,
                    min.imag=1e-9) {
  if (any(round(Re(b)) != b))
    stop("branch number for W must be an integer")
  if (!is.complex(z) && any(z<0)) z=as.complex(z)
  ## series expansion about -1/e
  ##
  ## p = (1 - 2*abs(b)).*sqrt(2*e*z + 2);
  ## w = (11/72)*p;
  ## w = (w - 1/3).*p;
  ## w = (w + 1).*p - 1
  ##
  ## first-order version suffices:
  ##
  w = (1 - 2*abs(b))*sqrt(2*exp(1)*z + 2) - 1
  ## asymptotic expansion at 0 and Inf
  ##
  v = log(z + as.numeric(z==0 & b==0)) + 2*pi*b*1i;
  v = v - log(v + as.numeric(v==0))
  ## choose strategy for initial guess
  ##
  c = abs(z + exp(-1));
  c = (c > 1.45 - 1.1*abs(b));
  c = c | (b*Im(z) > 0) | (!Im(z) & (b == 1))
  w = (1 - c)*w + c*v
  ## Halley iteration
  ##
  for (n in 1:maxiter) {
    p = exp(w)
    t = w*p - z
    f = (w != -1)
    t = f*t/(p*(w + f) - 0.5*(w + 2.0)*t/(w + f))
    w = w - t
    if (abs(Re(t)) < (2.48*eps)*(1.0 + abs(Re(w)))
        && abs(Im(t)) < (2.48*eps)*(1.0 + abs(Im(w))))
      break
  }
  if (n==maxiter) warning(paste("iteration limit (",maxiter,
                                ") reached, result of W may be inaccurate",sep=""))
  if (all(Im(w)<min.imag)) w = as.numeric(w)
  return(w)
}

##################################################


timingStamp<- function()
{
  stamp<- Sys.time()
  stamp<- gsub(x=stamp, pattern="(-|:)", replacement="", fixed=FALSE)
  stamp<- gsub(x=stamp, pattern=" ", replacement="-", fixed=TRUE)
  return(stamp)
}

constructFilenameFrom<- function(root="pfx", suffix=".svg", width=11, height=8.5, stamp=NULL)
{
  if (is.null(stamp))
  {
    stamp<- timingStamp()
  }
  E<- round(100*proc.time()["elapsed"])
  return(list(full=sprintf("%s-%s-%s%s", root, stamp, E, suffix), stem=sprintf("%s-%s-%s", root, stamp, E), width=width, height=height))
}

