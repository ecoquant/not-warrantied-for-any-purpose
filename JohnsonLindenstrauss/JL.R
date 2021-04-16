# JL.R Basic Johnson Lindenstrauss.
# Jan Galkowski, bayesianlogic.1@gmail.com
# 20th November 2018.

library(utils)
library(stringr)
library(rlist)
library(random)
library(FRACTION)
#library(matlab)
library(MBESS)
library(mvtnorm)
library(MASS)
library(TDA)

source("C:/builds/R/plottableSVG.R")
#source("plottableSVG.R")

randomizeSeed<- function(external=TRUE)
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

eye<- function(n) matlab::eye(n)

randomProjectionOfOnto<- function(Y, targetDim)
{
  N<- nrow(Y)
  M.X<- ncol(Y)
  L<- log(N)
  # Randomizing matrix
  Rinit<- matrix(rnorm((M.X*targetDim), 0, 1), nrow=M.X, ncol=targetDim)
  R<- t(apply(X=Rinit, MARGIN=1, FUN=function(r) r/sqrt(sum(r*r))))
  stopifnot( all( dim(R) == c(M.X, targetDim) ) )
  #
  # Deviation from identity
  dev<- sqrt(sum((eye(targetDim) - t(R) %*% R)^2))/prod(dim(R))
  cat(sprintf("Deviation per element is %.5f versus %.5f theoretical.\n", dev, 1/targetDim))
  # Map
  Z<- (Y %*% R)/sqrt(targetDim/L)
  return(list(R=R, Z=Z))
}


generateMultivariateNormal<- function(n=1000)
{
  utri<- matrix(c(0, 0.2, 0.001, 0.3, 0.05, 0.01, 0, 
                  0,   0, 0.01,    0, 0.05,    0, 0.01,
                  0,   0,    0,  0.1,  0.2, 0.01, 0,
                  0,   0,    0,    0,  0.01,   0, 0,
                  0,   0,    0,    0,     0, 0.2, 0, 
                  0,   0,    0,    0,     0,   0, 0.1,
                  0,   0,    0,    0,     0,   0, 0  ), 
                nrow=7, ncol=7, byrow=TRUE)/10
  corm<- diag(7) + utri + t(utri)
  sds<- 10*c(.2, .3, .2, 0.05, .1, .2, .3)
  covm<- cor2cov(cor.mat=corm, sd=sds)
  yields<- rmvnorm(n=n, mean=c(10, -20, 10, 5, 0, -10, 30), sigma=covm,
                   method="svd", pre0.9_9994=FALSE)
  return( yields )
}

#modesFromKDE2D<- function(kde, levelsForSets) { }

doGaussianMultivariateCase<- function(n=200, show=TRUE, tag="7-dimensional multivariate Gaussian", n.kde=500, kclu=170)
{
  X<- generateMultivariateNormal(n=n)
  Y<- scale(X, center=TRUE, scale=TRUE)
  onto<- randomProjectionOfOnto(Y=Y, targetDim=2)
  R<- onto$R
  Z<- as.matrix(onto$Z)
  rownames(Z)<- rownames(onto$Z)
  colnames(Z)<- colnames(onto$Z)
  cat("Projected ...\n")
  # Group and cluster
  h<- bandwidth.nrd(as.vector(Z))/10
  cat(sprintf("Bandwidth for density estimate is: %.5f\n", h))
  Z.kde<- kde2d(Z[,1], Z[,2], n=2*n.kde, h=h)
  cat("Kernel density ...\n")
# Z.clust<- clusterTree(X=Z, k=kclu, h=h, density="kde", dist="euclidean", Nlambda=n.kde, printProgress=TRUE)
  Z.clust<- clusterTree(X=Z, k=kclu, h=h, density="knn", dist="euclidean", Nlambda=n.kde, printProgress=TRUE)
  m<- length(Z.clust[["id"]])
  cat("Clusters ...\n")
  # Show if requested
  if (show)
  {
    h0<- openSVG(root="JLmultivariateGaussian-contour", width=11, height=8, pointsize=12, antialias="subpixel")
    image(Z.kde, col=topo.colors(round(5*m/4)), xlim=c(-5,5), ylim=c(-5,5))
    title(main=sprintf("%s : kernel density", tag))
    closeSVG(h0)
    h1<- openSVG(root="JLmultivariateGaussian-clustering", width=11, height=8, pointsize=12, antialias="subpixel")
    # plot clusters
    colorsToUse<- heat.colors(m)
    plot(Z, pch=21, cex=0.8, col="navy", bg="navy", xlab="Projection 1", ylab="Projection 2", xlim=c(-5,5), ylim=c(-5,5))
    for (i in Z.clust[["id"]])
    {
       w<- Z.clust[["DataPoints"]][[i]]
       cat(sprintf("@ %.0f : ", i)) 
       print(w)
       points(matrix(Z[w,],ncol = 2), col = colorsToUse[i], pch = 22, cex = 2)
    }
    title(main=sprintf("%s: cluster labels", tag))
    closeSVG(h1)
  }
  return( list(randomizer=R, projection=Z, kde=Z.kde, tag=tag, clustering=Z.clust) )
}

clG1<- doGaussianMultivariateCase(n=10000, show=TRUE, tag="7-dimensional multivariate Gaussian")


