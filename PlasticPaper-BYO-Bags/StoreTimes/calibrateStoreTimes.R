# calibrateStoreTimes.R
# Calibrate store times for plastic v paper bags survey, WEAC.
# Jan Galkowski, jan@westwood-statistical-studios.org, 18 February 2019.
# Last changed 18 February 2019.

# Data files are from the Google data of 17 February 2019 for these three stores.

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

getADay<- function(fn)
{
  parts<- strsplit(x=fn, split=".", fixed=TRUE)[[1]]
  if (2 == length(parts))
  {
    fn1<- parts[1]
  } else
  {
    stopifnot( 1 == length(parts) )
    fn1<- fn
  }
  aDay<- read.table(file=sprintf("%s.csv", fn1), 
                    sep=",", quote="", col.names=c("X", "Y"), 
                    row.names=NULL, stringsAsFactors=FALSE, check=TRUE, fill=FALSE)
  attr(aDay, "origin")<- fn1
  return(aDay)
}

# > SatRoche
#                X            Y
# 1   36.172332943 503.78194607
# 2   35.514654162 391.31887456
# 3   65.767878077 504.43962485
# 4   65.767878077 485.36694021
# 5   95.363423212 504.43962485
# 6   97.336459555 472.87104338
# 7  122.985932005 459.05978898
# 8  155.869871043 448.53692849
# 9  184.150058617 440.64478312
# 10 213.087924971 432.75263775
# 11 248.602579132 422.88745604
# 12 279.513481829 412.36459555
# 13 305.820633060 408.41852286
# 14 336.073856975 414.33763189
# 15 362.381008206 429.46424385
# 16 398.553341149 449.85228605
# 17 455.113716295 484.70926143
# 18 488.655334115 494.57444314
# 

rationalize<- function(aDay)
{
  origin<- attr(aDay, "origin")
  low<- mean(aDay[c(1,3,5),2])
  highest<- aDay[2,2]
  highs<- aDay[c(4,(6:nrow(aDay))),2]
  diffs<- low - highs
  diffsAsFracScored<- diffs/(low - highest)
  diffsAsFracOfMax<- diffs/max(diffs)
  diffsNormalized<- diffsAsFracScored/sum(diffsAsFracScored)
  return(list(origin=origin, heights=diffs, heights.fraction.scored=diffsAsFracScored, 
              heights.fraction.max=diffsAsFracOfMax, heights.normalized=diffsNormalized))
}

getAndRat<- function(fn)
{
  aDay<- getADay(fn)
  rtDay<- rationalize(aDay)
  rtDay$original<- aDay
  return(rtDay)
}

dataFiles<- c("FridaysLambert.csv",      # 1
              "FridaysRoche.csv",        # 2
              "FridaysWegmans.csv",      # 3
              "MondaysLamberts.csv",     # 4
              "MondaysRoche.csv",        # 5
              "MondaysWegmans.csv",      # 6
              "SaturdaysLambert.csv",    # 7
              "SaturdaysRoche.csv",      # 8
              "SaturdaysWegmans.csv",    # 9 
              "SundaysLamberts.csv",     # 10
              "SundaysRoche.csv",        # 11
              "SundaysWegmans.csv",      # 12
              "ThursdaysLambert.csv",    # 13
              "ThursdaysWegmans.csv",    # 14
              "ThurssdaysRoche.csv",     # 15
              "TuesdaysLambert.csv",     # 16
              "TuesdaysRoche.csv",       # 17
              "TuesdaysWegmans.csv",     # 18
              "WednesdaysLambert.csv",   # 19
              "WednesdaysRoche.csv",     # 20
              "WednesdaysWegmans.csv"    # 21
            )
             
RocheHours<- sapply(X=seq(700, 2100, 100), FUN=function(x) sprintf("%.0f", x))

WegmansHours<- sapply(X=seq(600, 2300, 100), FUN=function(x) sprintf("%.0f", x))

LambertHours<- sapply(X=seq(700, 1900, 100), FUN=function(x) sprintf("%.0f", x))

RocheFiles<- dataFiles[c(5, 17, 20, 15, 2, 8, 11)]
WegmansFiles<- dataFiles[c(6, 18, 21, 14, 3, 9, 12)]
LambertFiles<- dataFiles[c(4, 16, 19, 13, 1, 7, 10)]

allDays<- function(listOfFiles, 
                   hours, whichStore,
                   DOWs=c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"))
{
  stopifnot( whichStore %in% c("RocheBros", "Lamberts", "Wegmans") )
  k<- 1
  byDays<- matrix(NA, nrow=length(DOWs), ncol=length(hours))
  for (f in listOfFiles)
  {
    aDay<- getAndRat(f)
    m<- length(aDay$heights.fraction.scored)
    if (whichStore == "RocheBros")
    {
      if ("Sunday" == DOWs[k])
      {
        sd<- 2
        ed<- ncol(byDays) - 1
      } else
      {
        stopifnot( m == ncol(byDays) )
        sd<- 1
        ed<- ncol(byDays)
      }
    } else if (whichStore == "Lamberts")
    {
      if ("Sunday" == DOWs[k])
      {
        sd<- 1
        ed<- ncol(byDays) - 1
      } else
      {
        stopifnot( m == ncol(byDays) )
        sd<- 1
        ed<- ncol(byDays)
      }
    } else
    {
      stopifnot( m == ncol(byDays) )
      sd<- 1
      ed<- ncol(byDays)
    }
    stopifnot( m == (1+ed-sd) )
    byDays[k,(sd:ed)]<- aDay$heights.fraction.scored
    k<- 1 + k
  }
  rownames(byDays)<- DOWs
  colnames(byDays)<- hours
  return(byDays)
}


RocheBros<- allDays(RocheFiles, hours=RocheHours, whichStore="RocheBros")
Wegmans<- allDays(WegmansFiles, hours=WegmansHours, whichStore="Wegmans")
Lamberts<- allDays(LambertFiles, hours=LambertHours, whichStore="Lamberts")

rm(list = c("RocheFiles", "WegmansFiles", "LambertFiles", "dataFiles", "RocheHours", "WegmansHours", "LambertHours"))

(scrubFunctions())()







