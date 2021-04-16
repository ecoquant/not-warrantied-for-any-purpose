# Simple example of Schwarz-Arnason: mlmSchwArna01.R
# Last changed 10th February 2020.

library(bayess)
library(Rlab)

source("common.R")

N<- 91
zeta<- 21
H<- matrix(rbern(n=(N*zeta), prob=0.30), nrow=N, ncol=zeta)

#H.seen<- H[sample.int(N, 50),]
H.seen<- H

resolved<- bayess::gibbscap2(nsimu=5000, z=H.seen)

plot(resolved$p, type="l", col="blue", xlab="iterations", ylab="p", lwd=2,
     main="Probability of capture per subject")
abline(h=0.30, lty=6, lwd=1)
pause.with.message("probability of capture")

plot(resolved$phi[,1], type="l", col="blue", xlab="iterations", ylab="p", lwd=2,
     main="Probability of survival in first region")
abline(h=0.30, lty=6, lwd=1)
pause.with.message("probability of survival in first region")


plot(resolved$psi[1,1,], type="l", col="blue", xlab="iterations", ylab="p", lwd=2,
     main="Probability of migration from first region")
abline(h=0.30, lty=6, lwd=1)
pause.with.message("probability of migration from first region")

