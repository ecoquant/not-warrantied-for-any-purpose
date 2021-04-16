#
# PostStratify.R
#
# Post-stratification example. Jan Galkowski.
# 15 March 2019, last change 15 March 2019.
#
library(ggridges)
library(ggplot2)
library(reshape2)
library(rlist)
library(random)
library(gtools)
library(Hmisc)
library(MASS)

source("plottableSVG.R")

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

Section1<- c("North", "South", "East", "West")
S1Means<- c(100, 50, 30, 10)
Section2<- c("Red", "Green", "Blue")
S2Means<- c(2, 4, 8)
Section3<- c("Gather1", "Gather2", "Gather3")
S3Means<- c(10, -30, 50)
Treatments<- c("Control1", "Control2", "Control3", "Control4")
TMeans<- c(90, 110, -50, -90)

# Generate parameters for hyperpriors
k<- 1
for (s1 in Section1)
{
  for (s2 in Section2)
  {
    for (s3 in Section3)
    {
      for (ctrl in Treatments)
      {
        i1<- match(s1, Section1, nomatch=0)
        i2<- match(s2, Section2, nomatch=0)
        i3<- match(s3, Section3, nomatch=0)
        i4<- match(ctrl, Treatments, nomatch=0)
        mu<- S1Means[i1]+S3Means[i3]+TMeans[i4]+400
        scale<- mu/S2Means[i2]
        cat(sprintf("shape %.3f, mean %.0f, scale %.3f\n", S2Means[i2], mu, scale))
        v<- rgamma(1, shape=S2Means[i2], scale=scale)
        stopifnot( 0 < v )
        if (1 == k)
        {
          Collection<- data.frame(A1=s1, A2=s2, A3=s3, treatment=ctrl, observed=sum(rpois(3, lambda=v)), stringsAsFactors=FALSE)
        } else
        {
          Collection<- rbind(Collection,
                             data.frame(A1=s1, A2=s2, A3=s3, treatment=ctrl, observed=sum(rpois(3, lambda=v)), stringsAsFactors=FALSE),
                             stringsAsFactors=FALSE, make.row.names=FALSE)
        }
        k<-1 + k
      }
    }
  }
}

rm(list=c("k", "s1", "s2", "s3", "ctrl", "i1", "i2", "i3", "i4", "v", "mu", "scale"))


Censored<- Collection
L.parent<- nrow(Collection)
keep<- rep(TRUE, L.parent)
for (k in (1:L.parent))
{
  if (0.05>runif(1, min=0, max=1))
  {
    keep[k]<- FALSE
  }
  Censored$observed[k]<- Censored$observed[k] - rpois(1, lambda=15)
}

Censored<- Censored[which(keep),]

rm(list=c("k", "keep"))

L<- nrow(Censored)

cat(sprintf("\nParent Collection %.0f cases, Censored %.0f cases.\n", L.parent, L))




