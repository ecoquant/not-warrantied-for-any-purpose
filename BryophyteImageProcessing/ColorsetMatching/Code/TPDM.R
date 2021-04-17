#
# Ternary plots for discriminating mosses (TPDM).
# Jan Galkowski, begun 10th April 2021. Last changed 13th April 2021.
#

library(tictoc)
library(random)
library(scales)
library(sigmoid)
library(vegan)
library(imager)
library(concaveman)
library(Ternary)
library(xtable)


source("plottableSVG.R")

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


xcons<- function(x,y) {stopifnot(is.list(x)) ; append(x, list(y))}

codify<- function(rgbMatrix)
{
  stopifnot( is.matrix(rgbMatrix) )
  stopifnot( 3 == ncol(rgbMatrix) )
  coding<- function(r) sprintf("%-.0f|%-.0f|%-.0f", r[1], r[2], r[3])
  codes<- apply(X=rgbMatrix, MARGIN=1, FUN=coding)
  return(codes)
}

jaccardAnalysis<- function(codesList, labels)
{
  L<- length(codesList)
  jaccardCoefficients<- matrix(1, nrow=L, ncol=L)
  colnames(jaccardCoefficients)<- labels
  rownames(jaccardCoefficients)<- labels
  #
  for (i in (1:(L-1)))
  {
    for (j in ((1+i):L))
    {
      coef<- length(intersect(codesList[[i]], codesList[[j]]))/length(union(codesList[[i]], codesList[[j]]))
      jaccardCoefficients[i,j]<- coef
      jaccardCoefficients[j,i]<- coef
    }
  }
  #
  jaccardDistances<- 1 - jaccardCoefficients
  return(list(coefficients=jaccardCoefficients, distances=jaccardDistances))
}

ternaryFrom<- function(imFileList, nsamples=3000, chatty=TRUE, pointsAlphaSetting=0.15, 
                       contourAlphaSetting=0.5, 
                       caption="", alternativeLegends=NULL, plotType="ccm")
{
  #
  L<- length(imFileList)
  stopifnot( (6 > L) && (L > 0) )
  #
  TernaryPlot(alab = "Redder \u2192", blab = "\u2190 Greener", clab = "Bluer \u2192",
              lab.col = c('red', 'darkgreen', 'blue'), lab.font=2, lab.offset=c(0.13, 0.13, 0.13),
              point = 'right', lab.cex = 0.8, grid.minor.lines = 0,
              grid.lty = 'solid', col = "whitesmoke", grid.col = 'grey', 
              axis.col = "black", ticks.col = "black",
              axis.rotate = FALSE,
              padding = 0.08)
  pointsPch<- c(8, 21,22,24,25)
  pointsCex<- rep(1,5)
  colorsToUse<- c("blue", "darkgreen", "darkorange", "cyan", "red")
  alphaColors<- sapply(X=colorsToUse, FUN=function(x) alpha(x, pointsAlphaSetting))
  contourColors<- sapply(X=colorsToUse, FUN=function(x) alpha(x, contourAlphaSetting))
  #
  # (This hackery is needed because Ternary cannot handle an RGB of exactly triple zero.)
  scaler<- function(v) sapply(X=v, FUN=function(x) max(0.001, round(255*x)))
  #
  # (Expand out the Green channel.)
  logisticOf<- function(v, edge=6) round(255*logistic(((edge*(v-127.5))/127.5), k=1, x0=0))
  #
  tic("Ternary plot")
  #
  kthImage<- 1
  convertedImages<- list()
  for (imFile in imFileList)
  {
    stopifnot( file.exists(imFile) )
    imOrig<- load.image(imFile)
    stopifnot( 3 == spectrum(imOrig) )
    W<- width(imOrig)
    H<- height(imOrig)
    imArray<- array(NA, dim=c(W, H, 3), dimnames=list(width=NULL, height=NULL, channel=NULL))
    imArray[,,1]<- logisticOf(scaler(imOrig[,,1,1]), edge=3)
    imArray[,,2]<- logisticOf(scaler(imOrig[,,1,2]), edge=5) #6
    imArray[,,3]<- logisticOf(scaler(imOrig[,,1,3]), edge=3)
    imArrayDim<- dim(imArray)
    #
    if (chatty)
    {
      cat(sprintf("Image '%s' read and converted to array, dimensions (%.0f, %.0f, %.0f).\n", 
                  imFile, imArrayDim[1], imArrayDim[2], imArrayDim[3]))
    }
    #
    # Add data points
    dataPoints<- cbind(as.vector(imArray[,,1]), as.vector(imArray[,,2]), as.vector(imArray[,,3]))
    colnames(dataPoints)<- c("R", "G", "B")
    dataPointsUnique<- dataPoints[which(!duplicated(dataPoints)),]
    #
    # Analyze and report on RGB population.
    #
    whichToSample<- sample.int(nrow(dataPoints), nsamples, replace=TRUE)
    dataPointsSample<- dataPoints[whichToSample,]
    #
    if ("ccm" == plotType)
    {
      dataPointsUnique.xy<- t(apply(X=dataPointsUnique, MARGIN=1, FUN=TernaryCoords))
      ccm<- concaveman(dataPointsUnique.xy, concavity=2.5, length_threshold =0.05)
      polygon(x=ccm[,1], y=ccm[,2], col=NA, lty=1, lwd=3, border=contourColors[[kthImage]])
    } else if ("contour" == plotType)
    {
      TernaryDensityContour(dataPointsSample, resolution=80L, tolerance=-0.2/25, nlevels=18,
                            edgeCorrection=TRUE, col=contourColors[[kthImage]], lty=1, lwd=3,
                            drawlabels=FALSE)
    } else
    {
      stopifnot( "points" == plotType )
      TernaryPoints(dataPointsSample, pch=pointsPch[kthImage], cex=pointsCex[kthImage], 
#                   bg=alphaColors[kthImage], col=alphaColors[kthImage])
                    bg=NA, col=alphaColors[kthImage])
    }
    #
    convertedImages[[kthImage]]<- list(title=imFile, array=imArray, width=W, height=H, places=whichToSample,
                                       sample=dataPointsSample, sample.unique=dataPointsUnique)
    #
    kthImage<- 1+kthImage
  }
  #
  if (0 < nchar(caption))
  {
    title(main=caption, cex=2)
  }
  #
  if (is.null(alternativeLegends) || length(alternativeLegends) != length(imFileList))
  {
    if (("contour" == plotType) || ("ccm" == plotType))
    {
#     legend(locator(1), text.font=3,
      legend(x=.45, y=.65, text.font=3,
             legend = unlist(imFileList), 
             cex = 0.7, lty=1, lwd=7,
             col=colorsToUse[1:L], 
            )
    } else
    {
      stopifnot( "points" == plotType )
#     legend(locator(1), text.font=3,
      legend(x=.45, y=.65, text.font=3,
             legend = unlist(imFileList), 
             cex = 0.7, bty = 'n', pch = pointsPch[1:L], pt.cex = 1, 
             col=colorsToUse[1:L], pt.bg=NA
            )
    
    }
  } else
  {
    if (("contour" == plotType) || ("ccm" == plotType))
    {
#     legend(locator(1), text.font=3,
      legend(x=.45, y=.65, text.font=3,
             legend = alternativeLegends, 
             cex = 0.7, bty = 'n', lty=1, lwd=7,
             col=colorsToUse[1:L]
            )
    } else
    {
      stopifnot( "points" == plotType )
#     legend(locator(1), text.font=3,
      legend(x=.45, y=.65, text.font=3,
             legend = alternativeLegends, 
             cex = 0.7, bty = 'n', pch = pointsPch[1:L], pt.cex = 1, 
             col=colorsToUse[1:L], pt.bg=NA
            )
    }
  }
  #
  rgbCodes<- vector(mode="list", length=L)
  for (j in (1:L))
  {
    rgbCodes[[j]]<- codify(convertedImages[[j]]$sample.unique)
  }
  jac<- jaccardAnalysis(rgbCodes, labels=imFileList)
  if (chatty)
  {
    cat("Jaccard coefficients:\n")
    print(jac$coefficients, digits=3)
    cat("Jaccard distances:\n")
    print(jac$distances, digits=3)
    cat("\n")
    xtableCoefficients<- xtable(jac$coefficients, digits=rep(3,(1+L)), align=rep("c",(1+L)),
                           caption=sprintf("Jaccard coefficients among the %.0f images", L))
    print((xtableCoefficients), type="latex", include.rownames=FALSE, rotate.rownames=FALSE, rotate.colnames=FALSE)
    cat("\n")
    xtableDistances<- xtable(jac$distances, digits=rep(3,(1+L)), align=rep("c",(1+L)),
                           caption=sprintf("Jaccard distances among the %.0f images", L))
    print((xtableDistances), type="latex", include.rownames=FALSE, rotate.rownames=FALSE, rotate.colnames=FALSE)
    cat("\n")
  }
  #
  toc()
  #
  return(list(conversions=convertedImages, jaccard.coefficients=jac$coefficients, jaccard.distances=jac$distances))
}

# Inserting images in a plot
# library(png)
# library(grid)
# 
# mypng = readPNG('homer.png')
# image(volcano)
# grid.raster(mypng, x=.3, y=.3, width=.25) # print homer in ll conrner
# grid.raster(mypng, x=.9, y=.7, width=.5) # print bigger homer in ur corner

if (FALSE)
{ 
  h40a<- ternaryFrom(imFileList=list("M0770455.JPG", "SphagnumGirgensohnii2.jpg", "Thuidium1.jpg", "SphagnumGirgensohnii0.jpg"),
                     chatty=TRUE)
}

if (FALSE)
{
  h40b<- ternaryFrom(imFileList=list("AtrichumAugustatum2.jpg", "AtrichumAugustatum1.jpg", "AtrichumAugustatumManyCapsules.JPG"),
                     chatty=TRUE)
                     
}

if (FALSE)
{
  h40c<- ternaryFrom(imFileList=list("PlagiomniumCuspidatum1.jpg", "PlagiomniumCuspidatum2.jpg", "PlagiomniumCuspidatumCapsules.jpg"),
                     chatty=TRUE)
                     
}

if (FALSE)
{
  h40d<- ternaryFrom(imFileList=list("SphagnumQSite3InstanceC-1.jpg", "SphagnumQSite3InstanceC-2.jpg", 
                                     "SphagnumQSite3InstanceC-3.jpg", "SphagnumQSite3InstanceC-4.jpg", "SphagnumGirgensohnii0.jpg"),
                     chatty=TRUE)
                     
}

if (FALSE)
{
  h40e<- ternaryFrom(imFileList=list("Random1.jpg", "Random2.jpg", "Random3.jpg", "Random4.jpg", "Random5.jpg"),
                     chatty=TRUE)
                     
}


if (FALSE)
{
  h40g<- ternaryFrom(imFileList=list("Aug1.jpg", "Aug2.jpg", "Aug3.jpg", "Aug4.jpg", "Aug5.jpg"),
                     chatty=TRUE)
                     
}

if (FALSE)
{
  t202104120a<- ternaryFrom(imFileList=list("CapsuleTemplate.jpg", "Slice1MasterAtrichumAugustatum.jpg", "Slice2MasterAtrichumAugustatum.jpg", 
                                            "Slice3MasterAtrichumAugustatum.jpg", "Slice4MasterAtrichumAugustatum.jpg"),
                           caption="Assessing degrees of capsule coverage in photographs of Atrichum augustatum",
                           alternativeLegends=c("capsule template", "25% of test image", "50% of test image", "75% of test image", "100% of test image"),
                           chatty=TRUE, plotType="ccm")
}

if (TRUE)
{
  t202104120b<- ternaryFrom(imFileList=list("CapsuleTemplate.jpg", "MasterTestImageAtrichumAugustatum0-NoSetae.jpg", 
                                            "Slice4MasterAtrichumAugustatum.jpg", "AtrichumAugustatum20201204Capsules.jpg"),
                           caption="Assessing degrees of capsule coverage in photographs of Atrichum augustatum",
                           alternativeLegends=c("capsule template", "no capsules or setae", "test image 1", "test image 2"),
                           chatty=TRUE, plotType="ccm")
}

