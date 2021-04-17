#
# Semiautomatic finding and counting moss gametophytes in a community (SFCMGC).
# Jan Galkowski, begun 14th March 2021. Last changed 19th March 2021.
#

library(tictoc)
library(random)
library(imager)


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

allInstancesOfIn<- function(templateFile, imFile, scaleRange=c(1.2, 0.1), scaleIncrement=0.05, 
                            rotationalIncrement=2, thresholdCor=0.85, per=20, chatty=TRUE, debug=TRUE)
{
  stopifnot( is.vector(scaleRange) )
  stopifnot( 2 == length(scaleRange) )
  stopifnot((0 < scaleIncrement) && (scaleIncrement < 1))
  stopifnot((0 < rotationalIncrement) && (rotationalIncrement < 45))
  stopifnot( file.exists(imFile) )
  stopifnot( file.exists(templateFile) )
  #
  imOrig<- load.image(imFile)
  tmOrig<- load.image(templateFile)
  imFlat<- grayscale(imOrig)
  tmFlat<- grayscale(tmOrig)
  if (chatty)
  {
    cat(sprintf("Images '%s' and '%s' read and grayed.\n", imFile, templateFile))
  }
  #
  tic("Scaling with inner rotations")
  seen<- list()
  nseen<- 0
  for (scaling in seq(scaleRange[1], scaleRange[2], (-scaleIncrement)))
  {
     tmScaled<- imresize(im=tmFlat, scale=scaling, interpolation=6)
     tic("Inner rotations")
     rseen<- 0
     for (rotation in seq(0, (360-rotationalIncrement), rotationalIncrement))
     {
        # Cubic interpolation, Neumann boundary.
        tmRotatedScaled<- imrotate(im=tmScaled, angle=rotation, interpolation=2L, boundary=1L)
        imCor<- correlate(im=imFlat, filter=tmRotatedScaled, dirichlet=TRUE, normalise=TRUE) 
        pxCor<- (imCor > thresholdCor)
        if (any(pxCor))
        {
          seen<- xcons(seen, list(px=pxCor, rot=rotation, sca=scaling))
          rseen<- 1 + rseen
          if ((0 == (rseen%%per)) && chatty)
          {
            cat(sprintf("... added %.0f matches ...\n", rseen))
          }
        }
     }
     nseen<- rseen + nseen
     if ((0 < nseen) && chatty)
     {
       cat(sprintf("Rotated at scale %.2f, overall %.0f, this spin %.0f.\n", scaling, nseen, rseen))
#      browser()
     }
     toc()
  }
  if (chatty)
  {
    cat(sprintf("Completed '%s' on '%s', %.0f seen.\n", templateFile, imFile, nseen))
  }
  toc()
  return(list(seen=seen, nseen=nseen, image=imFile, template=templateFile))
}

if (TRUE)
{ 
  h40a<- allInstancesOfIn(templateFile="ThuidiumRTemplateDrawn40.jpg", imFile="IMG_20210310_162058.jpg", 
                         scaleRange=c(1.2, 0.1), scaleIncrement=0.05, rotationalIncrement=5, 
                         thresholdCor=0.85, chatty=TRUE, debug=TRUE)
  h40b<- allInstancesOfIn(templateFile="ThuidiumRTemplateDrawn40.jpg", imFile="IMG_20210310_162058.jpg", 
                         scaleRange=c(1.2, 0.1), scaleIncrement=0.05, rotationalIncrement=2, 
                         thresholdCor=0.85, chatty=TRUE, debug=TRUE)
  h40c<- allInstancesOfIn(templateFile="ThuidiumRTemplateDrawn40.jpg", imFile="IMG_20210310_162058.jpg", 
                         scaleRange=c(1.2, 0.1), scaleIncrement=0.025, rotationalIncrement=2, 
                         thresholdCor=0.85, chatty=TRUE, debug=TRUE)
  h40d<- allInstancesOfIn(templateFile="ThuidiumRTemplateDrawn40.jpg", imFile="IMG_20210310_162058.jpg", 
                         scaleRange=c(1.2, 0.1), scaleIncrement=0.05, rotationalIncrement=10, 
                         thresholdCor=0.85, chatty=TRUE, debug=TRUE)
}

