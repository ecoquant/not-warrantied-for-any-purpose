# ripplingFocus.R
# Rippling studies for PTD vs SNCD.
# Jan Galkowski, bayesianlogic.1@gmail.com, 23rd February 2019.
# Last changed 24th February 2019.

library(random)

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

pause.with.message<- function (message) 
{
  # Pause with a message between plots
  cat(message)
  cat("\n")
  cat("Paused. Press <Enter> to continue...")
  readline()
  invisible()
}


x<- 2*pi*seq(0.1, 100.01, 0.1)
layout(matrix(1:16, nrow=4, ncol=4, byrow=TRUE))

y.upper<- 100/2.5 + 4 + 1
k<- 1
Cases<- matrix(NA, nrow=17, ncol=length(x))
Cases[1,]<- x
RN<- c("x")
ColorsToUse<- sample(c(rainbow(n=8, start=0, end=0.1), rainbow(n=8, start=.33, end=0.4), rainbow(n=16, start=.55, end=0.7)), 16, replace=FALSE)
for (level in c(0, 25,10,2.5))
{
  for (wave in c(5, 10, 16, 20))
  {
    if (0 == level)
    {
      y<- 4*sin(x/wave)^2
    } else
    {
      y<- x/(2*pi*level) + 4*sin(x/wave)^2
    }
    plot(x/(2*pi), y, type="l", lwd=2, col=ColorsToUse[k],
         xlab="x", ylab="", main=sprintf("(L: %.1f, W: %.1f)", level, wave), cex.main=1.2, xlim=c(0.1,100), ylim=c(0.1,y.upper))
    Cases[1+k,]<- y
    RN[1+k]<- sprintf("%.0f.%.0f", (10*level), (10*wave))
    k<- 1 + k
  }
}

pause.with.message("splat")

dev.off()

