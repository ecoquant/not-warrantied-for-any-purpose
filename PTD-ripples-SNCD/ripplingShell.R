# ripplingShell.R
# Rippling studies for PTD vs SNCD.
# Jan Galkowski, bayesianlogic.1@gmail.com, 23rd February 2019.
# Last changed 24th February 2019.

library(Hmisc)
library(Matrix)
library(shapes)
library(geomorph)
library(xtable)

pause.with.message<- function (message) 
{
  # Pause with a message between plots
  cat(message)
  cat("\n")
  cat("Paused. Press <Enter> to continue...")
  readline()
  invisible()
}


divprocs2<- function(SERIES, period=25)
{
  stopifnot( is.data.frame(SERIES) ) 
  N.lm<- nrow(SERIES)
  I.lm<- 1:N.lm
  N<- ncol(SERIES)
  divergencesRho<- Matrix(0, N, N, dimnames=list(NULL, NULL))
  divergencesRhoRMS<- Matrix(0, N, N, dimnames=list(NULL, NULL))
  divergences.dF.RMS<- Matrix(0, N, N, dimnames=list(NULL, NULL))
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
  # Return Matrix objects, leaving conversion to a matrix, a  distance matrix, or a data
  # from to the consumer of the output. Can't anticipate that here. 
  return(list(divergences=divergencesRho, divergences.rms=divergencesRhoRMS, dF.rms=divergences.dF.RMS))
}

divprocsX<- function(SERIES, period=25)
{
  # SERIES is different than divprocs2, for example, because datasets come in 
  # pairs: Abscissa is given in first column, and then ordinate in second. 
  stopifnot( is.data.frame(SERIES) ) 
  N.lm<- nrow(SERIES)
  N2<- ncol(SERIES)
  stopifnot( 0 == (N2%%2) )
  N<- N2/2
  divergencesRho<- Matrix(0, N, N, dimnames=list(NULL, NULL))
  divergencesRhoRMS<- Matrix(0, N, N, dimnames=list(NULL, NULL))
  divergences.dF.RMS<- Matrix(0, N, N, dimnames=list(NULL, NULL))
  N1<- N-1
  for (i in (1:N1))
  {
    for (j in ((1+i):N))
    {
      A<- array(dim=c(N.lm, 2, 2))
      A[,1,1]<- SERIES[,(2*(i-1)+1)]
      A[,2,1]<- SERIES[,(2*i)]
      A[,1,2]<- SERIES[,(2*(j-1)+1)]
      A[,2,2]<- SERIES[,(2*j)]
      fpa<- procGPA(x=A, scale=TRUE, reflect=FALSE, eigen2d=FALSE, tol1=1e-5, tol2=1e-5, tangentcoords="residual",
                    proc.output=FALSE, distances=TRUE, pcaoutput=FALSE, alpha=0, affine=FALSE)
      divergencesRho[i,j]<- mean(fpa$rho)
      divergencesRhoRMS[i,j]<- fpa$rmsrho
      divergences.dF.RMS[i,j]<- fpa$rmsd1
      #
      divergencesRho[j,i]<- divergencesRho[i,j]
      divergencesRhoRMS[j,i]<- divergencesRhoRMS[i,j]
      divergences.dF.RMS[j,i]<- divergences.dF.RMS[i,j]
    }
    if (0 == (i%%period))
    {
      cat(sprintf("... did %.0f\n", i))
    }
  }
  stopifnot( !is.null(colnames(SERIES)) )
  subnames<- colnames(SERIES[,seq(2,N2,2)])
  colnames(divergencesRho)<- subnames
  rownames(divergencesRho)<- subnames
  colnames(divergencesRhoRMS)<- subnames
  rownames(divergencesRhoRMS)<- subnames
  colnames(divergences.dF.RMS)<- subnames
  rownames(divergences.dF.RMS)<- subnames
  #
  # Return Matrix objects, leaving conversion to a matrix, a  distance matrix, or a data
  # from to the consumer of the output. Can't anticipate that here. 
  return(list(divergences=divergencesRho, divergences.rms=divergencesRhoRMS, dF.rms=divergences.dF.RMS))
}

sh1<-  curve(6-(x-0.3)^2, -5, 5, col="blue", asp=1, type="n", n=1001)
x1<- 2*pi*seq(0,100,length.out=length(sh1$x))
sh1$y<- sh1$y + 4*(7/20)*sin(x1/5)^2
sh2<-  curve(13-(x-0.3)^2, -5, 5, col="blue", asp=1, type="n", n=1001)
x2<- 2*pi*seq(0,100,length.out=length(sh2$x))
sh2$y<- 2+sh2$y + 4*(15/20)*sin(x2/10)^2
sh3<-  curve(20-x^2, -5, 5, col="blue", asp=1, type="n", n=1001)
x3<- 2*pi*seq(0,100,length.out=length(sh3$x))
sh3$y<- 4 + sh3$y + 4*sin(x3/10)^2

plot(sh3$x, sh3$y, type="l", col="steelblue", lwd=2, xlab="", ylab="", main="3 bivalve-like edges", cex.main=1.2, ylim=c(0,30))
lines(sh2$x, sh2$y, lwd=2, col="blue")
lines(sh1$x, sh1$y, lwd=2, col="violet")

p1<- which.max(sh1$y)
p2<- which.max(sh2$y)
p3<- which.max(sh3$y)
text(sh1$x[p1], sh1$y[p1], cex=1.5, col="darkviolet", pos=3, labels="Edge1")
text(sh2$x[p2], sh2$y[p2], cex=1.5, col="darkblue", pos=3, labels="Edge2")
text(sh3$x[p3], sh3$y[p3], cex=1.5, col="navy", pos=3, labels="Edge3")

SERIES<- data.frame(Edge1.x=sh1$x, Edge1.y=sh1$y, Edge2.x=sh2$x, Edge2.y=sh2$y, Edge3.x=sh3$x, Edge3.y=sh3$y, stringsAsFactors=FALSE)

PTDs<- divprocsX(SERIES=SERIES)

print(PTDs$divergences, digits=4)

xtable.PTDs<- xtable(as.matrix(PTDs$divergences), caption="Procrustes tangent distances for the 3 bivalve-like edges", digits=4)
print((xtable.PTDs), type="latex", rotate.rownames=FALSE,rotate.colnames=FALSE)


