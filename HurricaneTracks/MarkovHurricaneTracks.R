# MarkovHurricaneTracks.R
#
# Generating hurricane tracks with a Markov
# spatial process. Jan Galkowski, implementing the
# code and sketches of Arthur Charpentier at 
# http://freakonometrics.hypotheses.org/17113. Professor
# Charpentier credits Dr Christophe Denuse-Baillon for
# the basic idea.
#
# Begun 11th August 2016. Last changed 11th August 2016.
#
# See also Markus Gesmann's
# http://www.magesblog.com/2014/10/visualising-seasonality-of-atlantic.html
#
# And Gaston Sanchez'
# http://rpubs.com/gaston/hurricanes
#

library(XML)
library(maps)
library(ks)
library(RColorBrewer)

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

# This value gets randomized on each call, so the JVP Cauchy mixture won't be exactly the same
# as shown on the blog or in successive runs.  This can be fixed if this behavior is unattractive.
# As it stands, it can be used to demonstrate performance of the blocks fitter for slightly different
# instances of the Cauchys mixture.
wonkyRandom<- randomizeSeed()

pause.with.message<- function(message)
{
  cat(message)
  cat("\n")
  cat("Paused. Press <Enter> to continue ...")
  readline()
  invisible()
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
    if ((x != "scrubFunctions") && is.function(get(x, envir=globalenv())))
    {
      rm(list=c(x), envir=globalenv())
      removed<- c(removed, x)
    }
  }
  return(function (){ rm(list=c("scrubFunctions"), envir=globalenv()) ; removed})
}

extract.track<-function(year=2012,p=TRUE){
loc <-  paste("http://weather.unisys.com/hurricane/atlantic/",year,"/index.php",sep="")
tabs <-  readHTMLTable(htmlParse(loc)) 
storms <-  unlist(strsplit(as.character(
tabs[[1]]$Name),split=" "))
index <-  storms%in%c("Tropical","Storm",paste("Hurricane-",1:6,sep=""),"Depression","Subtropical","Extratropical","Low",paste("Storm-",1:6,sep=""),
"Xxx")
#
nstorms  <-  storms[!index]
# Repair misspellings:
rewrites<- list(c("SIXTEE", "SIXTEEN"), c("FIFTEE", "FIFTEEN"), c("LAUR", "LAURA"), c("CHANTA", "CHANTAL"), c("ANDR", "ANDREA"),
                c("JOSEPH", "JOSEPHINE"), c("FLOY", "FLOYD"), c("KEIT", "KEITH"), c("CHARLI", "CHARLIE"))
for (z in rewrites)
{
  nstorms<- sapply(X=nstorms, FUN=function(v) gsub(x=v, pattern=sprintf("^%s$", z[1]), replacement=z[2], fixed=FALSE, perl=FALSE))
}
#
TRACK<- NULL
for(i in length(nstorms):1){
loc<- paste("http://weather.unisys.com/hurricane/atlantic/",year,"/",nstorms[i],"/track.dat",sep="")
track<- read.fwf(loc,skip=3,widths = c(4,6,8,12,4,6,20))
names(track)<- c("ADV","LAT","LON","TIME","WIND","PR","STAT")
track$LAT<- as.numeric(as.character(track$LAT))
track$LON<- as.numeric(as.character(track$LON))
track$WIND<- as.numeric(as.character(track$WIND))
track$PR<- as.numeric(as.character(track$PR))
track$year<- year
track$name<- nstorms[i]
TRACK<- rbind(TRACK,track)
if(p==TRUE){  cat(year,i,nstorms[i],nrow(track),"\n")}}
return(TRACK)}

# Load the hurricane tracks. This can take a while.
if (TRUE)
{
  m<-  extract.track(2012)
  print(head(m))
  TOTTRACK<-  NULL
  for(y in 2012:1851)
  {
    TOTTRACK<-  rbind(TOTTRACK,extract.track(y))
  }
# (scrubFunctions())()
# rm(list=c("y", "m"))
# save.image("HurricaneTracks.RData")
# stop()
}

stopifnot( exists("TOTTRACK") )

lineColor<- addAlpha(colors=c("blue"), alpha=0.3)
map("world",xlim=c(-100,-30),ylim=c(10,50),col="lightyellow",fill=TRUE, lwd=2)
title(xlab="longitude", ylab="latitude", main="Hurricane Tracks, 1851-2012")
for(n in unique(TOTTRACK$name)){
lines(TOTTRACK$LON[TOTTRACK$name==n],TOTTRACK$LAT[TOTTRACK$name==n],lwd=.5,col=lineColor)}
axis(side = 1, tck = -.01)
axis(side = 2, las = 1, tck = -.01)

pause.with.message("tracks")

U<- TOTTRACK[,c("LON","LAT")]
U<- U[!is.na(U$LON),]
H<- diag(c(.2,.2))
fat<- kde(U,H,xmin=c(min(U[,1]),min(U[,2])),xmax=c(max(U[,1]),max(U[,2])))
z<- fat$estimate
densityColors<- colorRampPalette(brewer.pal(9, "Blues"))(800)
dc<- addAlpha(densityColors, 0.8)
image(x=seq(-100, -30, length.out=151), y=seq(10, 50, length.out=151), z=z, col=dc, xlim=c(-100,-30),ylim=c(10,50), 
      xlab="longitude", ylab="latitude")
map("world",add=TRUE, lwd=2, xlim=c(-80,-40),ylim=c(10,50),col="black",fill=FALSE)
title(main="Density Plot of Hurricanes, 1851-2012")

pause.with.message("distribution")

gridx<- seq(-100,-30) 
gridy<- seq(-10,50)
pasgrid<- 1

mtransition<- function(i,j)
{
  MOVEMENT<- NA
  sumstop<- 0
  p<- 1
  idx<- which( (TOTTRACK$LON>=gridx[i])&(TOTTRACK$LON<gridx[i+1]) &
               (TOTTRACK$LAT>=gridy[j])&(TOTTRACK$LAT<gridy[j+1]) )
  if (length(idx)>0)
  {
    MOVEMENT<- data.frame(LON=rep(NA,length(idx)),LAT=rep(NA,length(idx)),D=rep(NA,length(idx)))
    for (s in 1:length(idx))
    {
      stops<- TRUE
      if ((is.na(TOTTRACK$name[idx[s]+1])==FALSE)&(is.na(TOTTRACK$name[idx[s]])==FALSE))
      {
        stops<- TOTTRACK$name[idx[s]+1]!=TOTTRACK$name[idx[s]]}
       if (stops==TRUE ){sumstop<- sumstop+1}
       if (stops==FALSE)
       {
         d<- (TOTTRACK$LON[idx[s]+1]-TOTTRACK$LON[idx[s]])^2+(TOTTRACK$LAT[idx[s]+1]-TOTTRACK$LAT[idx[s]])^2
         locx<- NA
         locy<- NA
         if ((d<90)& TOTTRACK$LON[idx[s]+1]<0)
         {
           locx<- floor(TOTTRACK$LON[idx[s]+1]/pasgrid)
           locy<- floor(TOTTRACK$LAT[idx[s]+1]/pasgrid)
         }
         MOVEMENT[s,]<- c(locx,locy,d)
       }
    }
    p<- sumstop/length(idx)
  }
  return(list(listemouv=MOVEMENT,probstop=p))
}

ListTransMat<- list()
ListStopMat <- list()
for (i in 1:(length(gridx)-1))
{
  for (j in 1:(length(gridy)-1))
  {
    nom<- abs(gridx[i])*1000+abs(gridy[j])
    M<- mtransition(i,j)
    ListTransMat[[nom]]<-  M$listemouv
    ListStopMat[[nom]]<-  M$probstop
  }
}

n<- nrow(TOTTRACK)
idx<- which(TOTTRACK$name[1:(n-1)]!=TOTTRACK$name[2:n])
STARTINGVALUE<- TOTTRACK[c(1,idx+1),c("LON","LAT")]

image(x=seq(-100, -30, length.out=151), y=seq(10, 50, length.out=151), z=z, col=dc, xlim=c(-100,-30),ylim=c(10,50), 
      xlab="longitude", ylab="latitude")
map("world",add=TRUE, lwd=2, xlim=c(-100,-30),ylim=c(10,50),col="black",fill=FALSE)
title(main="Generated Starting Points of Simulated Trajectories for Hurricanes, 1851-2012")
abline(v=gridx,col="lightgrey",lwd=.5)
abline(h=gridy,col="lightgrey",lwd=.5)
points(STARTINGVALUE$LON,STARTINGVALUE$LAT,
pch=22, col="green", bg="green", cex=0.5)
pause.with.message("generated trajectories")

sim.mouv<- function(LOC)
{
  lon<- LOC[1]
  lat<- LOC[2]
  nom<- as.numeric(abs(lon)*1000+abs(lat))
  if ((1>nom) || (nom>length(ListTransMat)))
  {
    LOC2<- NA
  } else
  {
    M<- ListTransMat[[nom]]
    if (is.null(M))
    {
      LOC2<- NA
    } else
    {
      M<- M[which(!is.na(M[,"LON"])),]
      stop<- ListStopMat[[nom]]
      if (is.null(stop))
      {
        LOC2<- NA 
      } else
      {
        u<- runif(1)
        if(u<=stop) { LOC2<- NA }
        if(u>stop) { LOC2<- M[sample(1:nrow(M),size=1),c("LON","LAT")] }
      }
    }
  }
  return(LOC2)
}

sim.traj<- function()
{
  i<- sample(1:nrow(STARTINGVALUE),size=1)
  LOC<- round(STARTINGVALUE[i,])
  MLOC<- LOC
  stop<- FALSE
  j<- 0
  while(stop==FALSE)
  {
    LOC2<- sim.mouv(LOC)
    if(any(is.na(LOC2))) { stop<- TRUE }
    if(all(!is.na(LOC2))) { LOC<- LOC2; MLOC<- rbind(MLOC,LOC)}
    j<- 1 + j
    if (0 == (j%%100))
    {
      cat(sprintf("%.0f interior iterations in trajectory simulation ...\n", i))
    }
  }
  return(MLOC)
}

map("world",xlim=c(-100,-30),ylim=c(10,50),col="lightyellow",fill=TRUE, lwd=2)
title(xlab="longitude", ylab="latitude")
lineColor<- addAlpha(colors=c("blue"), alpha=0.2)
for(n in unique(TOTTRACK$name)){
lines(TOTTRACK$LON[TOTTRACK$name==n],TOTTRACK$LAT[TOTTRACK$name==n],lwd=.5,col=lineColor)}

nSim<- 10
#simColors<- brewer.pal(nSim, "Paired")
simColors<- rainbow(1000)[sample.int(1000, nSim, replace=FALSE)]
for (kSim in (1:nSim))
{
   simulation<- sim.traj()
#  cat(sprintf("Simulation %.0f ...\n", kSim))
#  print(t(simulation))
   lines(x=simulation[,1], y=simulation[,2], lwd=1, col=simColors[kSim], lty=1)
   points(x=simulation[,1], y=simulation[,2], pch=21, col=simColors[kSim], bg=simColors[kSim], cex=1.1)
}
title(main=sprintf("Hurricane Trajectories, 1851-2012, with %.0f Simulated Ones Atop", nSim))
axis(side = 1, tck = -.01)
axis(side = 2, las = 1, tck = -.01)

pause.with.message("simulated hurricane trajectories")







