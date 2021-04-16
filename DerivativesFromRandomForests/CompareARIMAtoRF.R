#setwd("c:/Users/Jan/Documents/publish/blog/")

# AirPassengers ... ranger vs auto.arima
# and derivatives with random forests.
#
# Jan Galkowski, 27th June 2020.

########################################################################
# Setup and configuration

#memory.limit(size=21000)

expander.interp<- 3
look.back<- 18
NcoresToUse<- 2
#Ntrees<- 3000
#Ntrees<- 10000
Ntrees<- 100000
#Ntrees<- 300000


data(AirPassengers)
AirPassengers

library(assertthat)
library(random)
library(forecast)
#library(zoo)
#library(matrixStats)
#library(scales)
#library(plotrix)
library(ranger)
library(tictoc)
library(spectral)
library(akima)
library(bazar)
library(coda)
library(fda.usc)
#library(quantregForest)

pause.with.message<- function (message) 
{
  # Pause with a message between plots
  cat(message)
  cat("\n")
  cat("Paused. Press <Enter> to continue...")
  readline()
  invisible()
}


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

constructFilenameFrom<- function(root="pfx", suffix=".svg")
{
  stamp<- Sys.time()
  stamp<- gsub(x=stamp, pattern="(-|:)", replacement="", fixed=FALSE)
  stamp<- gsub(x=stamp, pattern=" ", replacement="-", fixed=TRUE)
  E<- round(100*proc.time()["elapsed"])
  return(list(full=sprintf("%s-%s-%s%s", root, stamp, E, suffix), stem=sprintf("%s-%s-%s", root, stamp, E)))
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


wonkyRandom<- randomizeSeed(external=TRUE)
#wonkyRandom<- randomizeSeed(external=FALSE)

mF<- function(series, nLags)
{
  N<- length(series)
  M<- (1+nLags):N
  TM<- t(sapply(X=M, FUN=function(k) series[k:(k-nLags)]))
  return(TM)
}

mW<- function(series, nLags)
{
  N<- length(series)
  n<- N-nLags
  M<- (1+nLags):N
  L<- length(M)
  TM<- matrix(NA, nrow=L, ncol=(nLags-1))
  for (k in (2:nLags))
  {
    TM[,(k-1)]<- tail(n=L, zoo::rollmean(series, k=k, align="left"))
  }
  return(TM)
}

matrixFromLagsAndMedianOf<- function(series, med, nLags, givenTocks)
{
  TM<- cbind(mF(series, nLags), mW(series, nLags))
  nr<- nrow(TM)
  adjustedTocks<- tail(n=nr, givenTocks)
  MTM<- cbind(TM, rep(med,nr))
  lagnames<- sapply(X=1:nLags, FUN=function(ell) sprintf("lag%-.0f", ell))
  windownames<- sapply(X=2:nLags, FUN=function(ell) sprintf("win%-.0f", ell))
  colnames(MTM)<- c("response", lagnames, windownames, "median")
  return(list(MTM=MTM, tocks=adjustedTocks))
}

vote<- function(r) length(which(r)) > length(r)/2

# From McElreath's -rethinking- package which otherwise
# cannot be installed. (Don't know why.)
#
# highest posterior density interval, sensu Box and Tiao
# requires coda library
HPDI <- function( samples , prob=0.89 ) {
    # require(coda)
    coerce.list <- c( "numeric" , "matrix" , "data.frame" , "integer" , "array" )
    if ( inherits(samples, coerce.list) ) {
        # single chain for single variable
        samples <- coda::as.mcmc( samples )
    }
    x <- sapply( prob , function(p) coda::HPDinterval( samples , prob=p ) )
    # now order inside-out in pairs
    n <- length(prob)
    result <- rep(0,n*2)
    for ( i in 1:n ) {
        low_idx <- n+1-i
        up_idx <- n+i
        # lower
        result[low_idx] <- x[1,i]
       # upper
        result[up_idx] <- x[2,i]
        # add names
        names(result)[low_idx] <- concat("|",prob[i])
        names(result)[up_idx] <- concat(prob[i],"|")
    }
    return(result)
}

basisSize<- function(L) round(3*sqrt(L))

bootstrapUncertaintyCloudsFrom<- function(P, tocks, base0, base1, base2, nBoot=100, pHPDI=0.67)
{
  assert_that( is.count(nBoot) )
  assert_that( is.scalar(nBoot) )
  assert_that( is.matrix( P ) )
  assert_that( all( c("data", "argvals") %in% names(base0) ) )
  assert_that( all( c("data", "argvals") %in% names(base1) ) )
  assert_that( all( c("data", "argvals") %in% names(base2) ) )
  assert_that( nrow(P) < ncol(P) )
  #
  basePredictions<- base0$data
  fstPredictions<- base1$data
  secPredictions<- base2$data
  #
  nr<- nrow(P)
  nc<- ncol(P)
  baseRange<- matrix(NA, nrow=nr, ncol=nBoot)
  fstRange<- matrix(NA, nrow=nr, ncol=nBoot)
  secRange<- matrix(NA, nrow=nr, ncol=nBoot)
  #
  tic("bootstrapping uncertainty")
  for (kB in (1:nBoot))
  {
    picks<- sample.int(n=nc, size=nc, replace=TRUE)
#   generated<- matrixStats::rowMedians(P[,picks])
    generated<- rowMeans(P[,picks])
    fdFromGenerated<- suppressWarnings(fdata(mdata=matrix(generated, nrow=nr, ncol=1),argvals=tocks))
    splining<- suppressWarnings(fdata.deriv(fdataobj=fdFromGenerated, nderiv=0, method="bspline",
                                  class.out="fdata", nbasis=basisSize(nr)))
    fst<-      suppressWarnings(fdata.deriv(fdataobj=fdFromGenerated, nderiv=1, method="bspline", 
                                             class.out="fdata", nbasis=basisSize(nr)))
    sec<-      suppressWarnings(fdata.deriv(fdataobj=fdFromGenerated, nderiv=2, method="bspline", 
                                               class.out="fdata", nbasis=basisSize(nr)))
    baseRange[,kB]<- splining$data
    fstRange[,kB]<- fst$data
    secRange[,kB]<- sec$data
  }
  toc(log=TRUE)
  #
  baseHPDI<- matrix(NA, nrow=nr, ncol=2)
  fstHPDI<- matrix(NA, nrow=nr, ncol=2)
  secHPDI<- matrix(NA, nrow=nr, ncol=2)
  for (k in (1:nr))
  {
#   baseHPDI[k,]<- HPDI(baseRange[k,], prob=pHPDI)
#   fstHPDI[k,]<- HPDI(fstRange[k,], prob=pHPDI)
#   secHPDI[k,]<- HPDI(secRange[k,], prob=pHPDI)
    # See, for example: https://datascienceplus.com/prediction-interval-the-wider-sister-of-confidence-interval/
    dfBase<- data.frame(v=baseRange[k,], strings.as.factors=FALSE)
    PI.base<-predict(lm(dfBase$v~1), interval="predict", level=pHPDI, type="response") 
    baseHPDI[k,]<- PI.base[1,c(2,3)]
    dfFst<- data.frame(v=fstRange[k,], strings.as.factors=FALSE)
    PI.fst<-predict(lm(dfFst$v~1), interval="predict", level=pHPDI, type="response") 
    fstHPDI[k,]<- PI.fst[1,c(2,3)]
    dfSec<- data.frame(v=secRange[k,], strings.as.factors=FALSE)
    PI.sec<-predict(lm(dfSec$v~1), interval="predict", level=pHPDI, type="response") 
    secHPDI[k,]<- PI.sec[1,c(2,3)]
  }
  #
  return(list(HPDI.base0=baseHPDI, HPDI.base1=fstHPDI, HPDI.base2=secHPDI))
}

########################################################################
# Script execution:


fit<- auto.arima(AirPassengers)
names(fit)
fit
fit$coef
fit$var.coef

S<- as.vector(unclass(AirPassengers))
tocks.original<- as.vector(unclass(time(AirPassengers)))

S.scaled<- scale(S, center=TRUE, scale=TRUE)

S.mean<- mean(S, na.rm=TRUE)
S.sd<- sd(S, na.rm=TRUE)

N.interp<- expander.interp*length(S.scaled)

tocks.interpolated<- seq(min(tocks.original), max(tocks.original), length.out=N.interp)
S.scaled.interpolated<- aspline(x=tocks.original, S.scaled, n=N.interp, ties=mean, method="improved", degree=4)$y


MTMbuild<- matrixFromLagsAndMedianOf(series=S.scaled.interpolated, med=median(S.scaled), 
                                     nLags=expander.interp*look.back, givenTocks=tocks.interpolated)


MTM.df<- as.data.frame(MTMbuild$MTM, stringsAsFactors=FALSE)
MTM.df$tocks<- MTMbuild$tocks
MTM.df$original<- apply(FUN=any, MARGIN=1, 
                         X=min(diff(MTM.df$tocks))/2 > abs(outer(FUN="-", X=MTM.df$tocks, Y=tocks.original)))
rm(list=c("MTMbuild"))


train<- MTM.df

train

cformula<- Reduce(x=colnames(MTM.df)[3:(ncol(MTM.df)-2)], init="response ~ lag1 ", f=function(a,b) sprintf("%s + %s", a, b))

tic("Importance determination")
rfImporta<- ranger(formula=as.formula(cformula),
                     data=train, num.trees=Ntrees, mtry=ceiling(sqrt(ncol(MTM.df))), 
                     importance = "permutation",
                     scale.permutation.importance=TRUE, local.importance=TRUE,
                     write.forest=TRUE, probability=FALSE,
                     min.node.size=5, max.depth=0,
                     replace=TRUE, sample.fraction=1, 
                     case.weights = NULL, class.weights = NULL,
                     splitrule="extratrees", num.random.splits=10,
#                    respect.unordered.factors="partition", 
                     regularization.factor=1, regularization.usedepth=FALSE,
                     keep.inbag=TRUE, inbag=NULL, holdout=FALSE, quantreg=FALSE,
                     oob.error = TRUE,
                     num.threads=NcoresToUse, save.memory=FALSE, seed=NULL,
                     verbose=TRUE, classification=FALSE
                )
toc(log=TRUE)

sorted<- sort.int(x=rfImporta$variable.importance, decreasing=TRUE, method="shell", index.return=TRUE)
cat("raw variable importance:\n")
print(rfImporta$variable.importance[sorted$ix], digits=4)

cat("normalized variable importance`:\n")
rfImporta$normalized.variable.importance<- rfImporta$variable.importance/sum(rfImporta$variable.importance)
sorted<- sort.int(x=rfImporta$normalized.variable.importance, decreasing=TRUE, method="shell", index.return=TRUE)
rfImporta$normalized.variable.importance<- rfImporta$normalized.variable.importance[sorted$ix]
print(rfImporta$normalized.variable.importance, digits=4)
rm(list=c("sorted"))

which.terms<- names(rfImporta$normalized.variable.importance)[which(0.005<rfImporta$normalized.variable.importance)]

rformula<- Reduce(x=which.terms[-1],
                   init=sprintf("%s ~ %s", "response", which.terms[1]), 
                   f=function(a,b) sprintf("%s + %s", a, b))


tic("Training")
rfTrained<- ranger(formula=as.formula(rformula),
                     data=train, num.trees=Ntrees, mtry=ceiling(sqrt(ncol(MTM.df))), importance="none",
                     write.forest=TRUE, probability=FALSE,
                     min.node.size=5, max.depth=0,
                     replace=TRUE, sample.fraction=1, 
                     case.weights = NULL, class.weights = NULL,
                     splitrule="extratrees", num.random.splits=10,
#                    respect.unordered.factors="partition", 
#                    scale.permutation.importance=FALSE, local.importance=FALSE,
                     regularization.factor=1, regularization.usedepth=FALSE,
                     keep.inbag=TRUE, inbag=NULL, holdout=FALSE, quantreg=FALSE,
                     oob.error = TRUE,
                     num.threads=NcoresToUse, save.memory=FALSE, seed=NULL,
                     verbose=TRUE, classification=FALSE
                )
toc(log=TRUE)

slicer.interpolated<- (1+expander.interp*look.back):length(S.scaled.interpolated)
tocks.sliced<- tocks.interpolated[slicer.interpolated]
S.interpolated<- S.mean+ S.sd*S.scaled.interpolated[slicer.interpolated]
predFromTrained<- S.sd*rfTrained$predictions + S.mean

residualsFromTrained<- S.interpolated - predFromTrained

cat("ranger R-squared: ")
print(rfTrained$r.squared, digits=4)

reweightedSDResiduals<- zoo::rollapply(data=residualsFromTrained, width=look.back,
                                       FUN=function(w) sd(win.tukey(w, a=0.5)))

treeFrom.rfTrained.1<- treeInfo(rfTrained, 1)

########################################################################
# Predictions

tic("prediction overall")
rfPredicted<- predict(object=rfTrained, data=train, predict.all=FALSE, 
                       num.trees=rfTrained$num.trees, type="response", 
                       inbag.counts=rfTrained$inbag.counts,
                       seed=NULL, num.threads=NcoresToUse, verbose=TRUE)
toc(log=TRUE)

predFromPred<- predictions(rfPredicted)*S.sd + S.mean

are_equal(length(predFromPred), length(predFromTrained))
are_equal(length(predFromPred), length(slicer.interpolated))

tic("prediction for each tree")
rfPredictedSingletonTrees<- predict(object=rfTrained, data=train, predict.all=TRUE, 
                                     num.trees=rfTrained$num.trees, type="response", 
                                     seed=NULL, num.threads=NcoresToUse, verbose=TRUE,
                                     inbag.counts=rfTrained$inbag.counts)
toc(log=TRUE)

########################################################################
# Derivatives
tic("derivatives setup")
L<- length(slicer.interpolated)

fdFromForest<- suppressWarnings(fdata(mdata=matrix(predictions(rfPredicted), nrow=length(slicer.interpolated), ncol=1),
                                        argvals=tocks.sliced))
toc(log=TRUE)                     

tic("derivatives calculation")
splining<- suppressWarnings(fdata.deriv(fdataobj=fdFromForest, nderiv=0, method="bspline", 
                                           class.out="fdata", nbasis=basisSize(L)))
fst<-      suppressWarnings(fdata.deriv(fdataobj=fdFromForest, nderiv=1, method="bspline", 
                        class.out="fdata", nbasis=basisSize(L)))
sec<-      suppressWarnings(fdata.deriv(fdataobj=fdFromForest, nderiv=2, method="bspline", 
                                           class.out="fdata", nbasis=basisSize(L)))

# No need to rescale since the rescaled data was fit to fdata.
predFromSplined<- as.vector(splining$data)*S.sd + S.mean
dotFromSplined<- as.vector(fst$data) + S.sd

rm(list=c("L"))
toc(log=TRUE)
########################################################################
uClouds<- bootstrapUncertaintyCloudsFrom(P=rfPredictedSingletonTrees$predictions, 
                                         tocks=tocks.sliced,
                                         base0=splining, base1=fst, base2=sec,
                                         nBoot=200, pHPDI=0.9)
HPDI.splining<- uClouds$HPDI.base0
HPDI.fst<- uClouds$HPDI.base1
HPDI.sec<- uClouds$HPDI.base2
rm(list=c("uClouds"))
########################################################################

tic("checkpointing")
ToSave<- constructFilenameFrom(root=sprintf("ArimaToRF-checkpoint-%-.0f-%-.0f-%-.0f",
                                            Ntrees, look.back, expander.interp),
                               suffix=".RData")

save(train, rfTrained, rformula, cformula, tocks.original, S.mean, S.sd, expander.interp, 
     look.back, NcoresToUse, Ntrees, N.interp, S.scaled, tocks.interpolated, S.scaled.interpolated,
     slicer.interpolated,S.interpolated, predFromTrained, residualsFromTrained, 
     reweightedSDResiduals, treeFrom.rfTrained.1, rfPredicted, predFromPred,
     rfPredictedSingletonTrees, fdFromForest, splining, fst, sec, 
     predFromSplined, dotFromSplined,HPDI.splining, HPDI.fst, HPDI.sec,
     file=ToSave$full, compress=TRUE)

toc(log=TRUE)
plot(predFromTrained[which(train$original)], fit$fitted[19:144], type="p", pch=21, col="navy", bg="navy", 
main="ranger versus auto.arima on AirPassengers\n(ranger results on training)", xlab="auto.arima", ylab="ranger")
abline(b=1, a=0, lty=6, col="blue", lwd=2)

pause.with.message("post-training plot")

cat("rfTrained residuals:\n")
rfTrained$residuals<- S.interpolated - predFromTrained
plot(density(rfTrained$residuals))

pause.with.message("density of residuals")

plot(S.interpolated, rfTrained$residuals, type="p", pch=22, col="blue", bg="blue",
xlab="index of sample", ylab="residuals", main="ranger performance on AirPassengers data set")

pause.with.message("ranger residuals")

plot(tocks.original, as.vector(unclass(fit$residuals)), type="p", pch=22, col="blue", bg="blue",
xlab="index of sample", ylab="residuals", main="auto.arima performance on AirPassengers data set")

pause.with.message("auto.arima residuals")

alphamagenta<- scales::alpha("magenta", 0.4)

plot(as.vector(unclass(time(AirPassengers))), S, type="p", pch=22, col="green", lwd=2, bg="green", cex=1.3,
        xlab="time", ylab="passenger counts",
      main="Comparison of auto.arima and random forest predictions for Air Passengers series data\n(via R package ranger)")
points(as.vector(unclass(time(AirPassengers))), fit$fitted, pch=24, col="red", bg="red")
points(tocks.sliced, predFromTrained, pch=21, col="blue",cex=1, bg="blue")
points(tocks.sliced, predFromPred, pch=25, col=alphamagenta, bg=alphamagenta, cex=1)
lines(tocks.sliced, predFromSplined, lwd=2, col="darkred", lty=1)

legend("topleft", col=c("green", "red", "blue", alphamagenta,"darkred"),  pch=c(22,24,21,25,NA), text.font=2,
        pt.bg=c("green", "red", "blue", alphamagenta,NA), pt.cex=c(1.3,1.3,1.3,1.3,NA), bty="n",
        legend=c("data", "auto.arima prediction", "ranger prediction from training", 
                 "ranger prediction, all trees", "smoothing spline derived from ranger predictions"),
        lty=c(NA, NA, NA, NA, 1), 
        lwd=c(NA, NA, NA, NA, 2), 
        seg.len=1)


pause.with.message("verus data")

alphagrey<- scales::alpha("grey", 0.7)
alphanavy<- scales::alpha("navy", 0.4)
alphagreen<- scales::alpha("green", 0.4)

layout(matrix(1:3, nrow=3, ncol=1))
plot(as.vector(unclass(time(AirPassengers))), S, type="p", pch=22, col="green", bg="green", cex=1,
        xlab="time", ylab="passenger counts",
      main="Random forest predictions and uncertainty for Air Passengers series data")
#
base.top<- HPDI.splining[,2]*S.sd + S.mean
base.bot<- HPDI.splining[,1]*S.sd + S.mean
#
points(tocks.sliced, predFromPred, pch=25, col="navy", bg="navy", cex=0.7)
lines(tocks.sliced, base.top, col=alphagrey, lty=1, lwd=4)
lines(tocks.sliced, base.bot, col=alphagrey, lty=1, lwd=4)
lines(tocks.sliced, predFromSplined, lwd=1, col="red", lty=1)
#
fst.top<- HPDI.fst[,2]
fst.bot<- HPDI.fst[,1]
#
plot(x=c(tocks.sliced, tocks.sliced), y=c(fst.top, fst.bot), type="n", 
      xlab="time", ylab="value of first derivative",
      main="Estimated first  derivatives from Random Forest regression\nusing centered and rescaled data")
#
lines(tocks.sliced, fst.top, col="grey", lty=1, lwd=4)
lines(tocks.sliced, fst.bot, col="grey", lty=1, lwd=4)
lines(tocks.sliced, as.vector(fst$data), lwd=1, col="red")
#
sec.top<- HPDI.sec[,2]
sec.bot<- HPDI.sec[,1]
#
plot(x=c(tocks.sliced, tocks.sliced), y=c(sec.top, sec.bot), type="n", 
      xlab="time", ylab="value of second derivative",
      main="Estimated second derivatives from Random Forest regression\nusing centered and rescaled data")
#
lines(tocks.sliced, sec.top, col="grey", lty=1, lwd=4)
lines(tocks.sliced, sec.bot, col="grey", lty=1, lwd=4)
lines(tocks.sliced, as.vector(sec$data), lwd=1, col="red")
#

pause.with.message("derivatives overall")

dev.off()

xlimit<- c(1955, 1956.2)

layout(matrix(1:3, nrow=3, ncol=1))
plot(as.vector(unclass(time(AirPassengers))), S, type="p", pch=22, col="green", bg="green", cex=1,
        xlab="time", ylab="passenger counts", xlim=xlimit,
      main="Random forest predictions and uncertainty for Air Passengers series data")
#
base.top<- HPDI.splining[,2]*S.sd + S.mean
base.bot<- HPDI.splining[,1]*S.sd + S.mean
#
points(tocks.sliced, predFromPred, pch=25, col="navy", bg="navy", cex=0.7)
lines(tocks.sliced, base.top, col=alphagrey, lty=1, lwd=4)
lines(tocks.sliced, base.bot, col=alphagrey, lty=1, lwd=4)
lines(tocks.sliced, predFromSplined, lwd=1, col="red", lty=1)
#
fst.top<- HPDI.fst[,2]
fst.bot<- HPDI.fst[,1]
#
plot(x=c(tocks.sliced, tocks.sliced), y=c(fst.top, fst.bot), type="n", 
      xlab="time", ylab="value of first derivative", xlim=xlimit,
      main="Estimated first  derivatives from Random Forest regression\nusing centered and rescaled data")
#
lines(tocks.sliced, fst.top, col="grey", lty=1, lwd=4)
lines(tocks.sliced, fst.bot, col="grey", lty=1, lwd=4)
lines(tocks.sliced, as.vector(fst$data), lwd=1, col="red", lty=1)
#
sec.top<- HPDI.sec[,2]
sec.bot<- HPDI.sec[,1]
#
plot(x=c(tocks.sliced, tocks.sliced), y=c(sec.top, sec.bot), type="n", 
      xlab="time", ylab="value of second derivative", xlim=xlimit,
      main="Estimated second derivatives from Random Forest regression\nusing centered and rescaled data")
#
lines(tocks.sliced, sec.top, col="grey", lty=1, lwd=4)
lines(tocks.sliced, sec.bot, col="grey", lty=1, lwd=4)
lines(tocks.sliced, as.vector(sec$data), lwd=1, col="red", lty=1)
#

pause.with.message("close-up")

dev.off()

plot(as.vector(unclass(time(AirPassengers))), S, type="p", pch=22, col="green", bg="green", cex=1,
        xlab="time", ylab="passenger counts",
      main="Random forest predictions and uncertainty for Air Passengers series data")
#
base.top<- HPDI.splining[,2]*S.sd + S.mean
base.bot<- HPDI.splining[,1]*S.sd + S.mean
#
points(tocks.sliced, predFromPred, pch=25, col="navy", bg="navy", cex=0.7)
lines(tocks.sliced, base.top, col=alphagrey, lty=1, lwd=4)
lines(tocks.sliced, base.bot, col=alphagrey, lty=1, lwd=4)
lines(tocks.sliced, predFromSplined, lwd=1, col="red", lty=1)
#

legend("topleft", col=c("green", "red", "navy", alphagrey),  pch=c(22,NA,25,NA), text.font=2,
        pt.bg=c("green", NA, "navy", NA), pt.cex=c(1.3,NA,1.3,NA), bty="n",
        legend=c("data", "splined forest prediction", "ranger prediction from forest", 
                 "prediction error from bootstrapped trees"),
        lty=c(NA, 1, NA, 1), 
        lwd=c(NA, 1, NA, 4), 
        seg.len=1)




