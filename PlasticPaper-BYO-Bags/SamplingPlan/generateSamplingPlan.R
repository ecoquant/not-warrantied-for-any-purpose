# generateSamplingPlan.R
# Produce a sampling plan based upon the visit tables generated from calibrateStoreTimes.R 
# using digitized versions of Google visitation data, for plastic v paper bags survey, WEAC.
# Jan Galkowski, jan@westwood-statistical-studios.org, 18 February 2019.
# Last changed 18 February 2019.

library(random)
library(FRACTION)
library(splines)
library(Hmisc)
library(pspline)
library(reshape2)

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

wonkyRandom<- randomizeSeed(external=FALSE)

stopifnot( exists("RocheBros") )
stopifnot( exists("Wegmans") )
stopifnot( exists("Lamberts") )

# A visit to Roche Bros in Westwood on Saturday, 16 Feb 2019, between 10:36 a.m. and 11:08 a.m.
# found 72 people. 

#  9:00-10:00 RocheBros Saturday 0.40331491713
# 10:00-11:00 RocheBros Saturday 0.49723756906
# 11:00-12:00 RocheBros Saturday 0.56353591160
# 12:00-13:00 RocheBros Saturday 0.63535911602
# 13:00-14:00 RocheBros Saturday 0.72928176796

VisitsSplined<- interpSpline(obj1=c(9.5, 10.5, 11.5, 12.5, 13.5), 
                             obj2=c(0.40331491713, 0.49723756906, 0.56353591160, 0.63535911602, 0.72928176796), 
                             bSpline=TRUE, ord=4L, na.action=na.fail, sparse=TRUE)
PredictScores<- predict(object=VisitsSplined, x=seq(10.63, 11.13, 0.1), deriv=0)$y
PredictScore<- mean(PredictScores)
PredictVisitsPerUnit<- round(72/PredictScore)
cat(sprintf("Predicted number of people for score of unity is: %.0f\n", PredictVisitsPerUnit))

# First stage sampling frame is hours within days, so sampling from all hour-day combinations

RocheBrosTrafficNormalizer<- sum(RocheBros, na.rm=TRUE)
WegmansTrafficNormalizer<- sum(Wegmans, na.rm=TRUE)
LambertsTrafficNormalizer<- sum(Lamberts, na.rm=TRUE)

# Limit need for people to survey from 9:00 a.m. through 10 p.m.

WegmansAcceptable<- Wegmans[,(4:16)]
RocheBrosAcceptable<- RocheBros[,(3:15)]
LambertsAcceptable<- Lamberts[,(3:13)]

RocheBrosTrafficNormalizer.acceptable<- sum(RocheBrosAcceptable, na.rm=TRUE)
WegmansTrafficNormalizer.acceptable<- sum(WegmansAcceptable, na.rm=TRUE)
LambertsTrafficNormalizer.acceptable<- sum(LambertsAcceptable, na.rm=TRUE)
OverallTrafficNormalizer.acceptable<- RocheBrosTrafficNormalizer.acceptable + WegmansTrafficNormalizer.acceptable +
                                      LambertsTrafficNormalizer.acceptable

RocheBros.melted<- melt(RocheBrosAcceptable, na.rm=TRUE)
RocheBros.melted$value<- RocheBros.melted$value/RocheBrosTrafficNormalizer.acceptable

Wegmans.melted<- melt(WegmansAcceptable, na.rm=TRUE)
Wegmans.melted$value<- Wegmans.melted$value/WegmansTrafficNormalizer.acceptable

Lamberts.melted<- melt(LambertsAcceptable, na.rm=TRUE)
Lamberts.melted$value<- Lamberts.melted$value/LambertsTrafficNormalizer.acceptable

# First stage sample per store
N.stage1<- 12
N.stage1.RocheBros<- round(N.stage1*RocheBrosTrafficNormalizer.acceptable/OverallTrafficNormalizer.acceptable)
N.stage1.Wegmans<- round(N.stage1*WegmansTrafficNormalizer.acceptable/OverallTrafficNormalizer.acceptable)
N.stage1.Lamberts<- round(N.stage1*LambertsTrafficNormalizer.acceptable/OverallTrafficNormalizer.acceptable)

cat(sprintf("RocheBros: %.2f\n", N.stage1.RocheBros))
cat(sprintf("Wegmans: %.2f\n", N.stage1.Wegmans))
cat(sprintf("Lamberts: %.2f\n", N.stage1.Lamberts))

RocheBros.stage1<- RocheBros.melted[sample.int(nrow(RocheBros.melted), size=N.stage1.RocheBros, replace=FALSE, prob=RocheBros.melted$value),]
Wegmans.stage1<- Wegmans.melted[sample.int(nrow(Wegmans.melted), size=N.stage1.Wegmans, replace=FALSE, prob=Wegmans.melted$value),]
Lamberts.stage1<- Lamberts.melted[sample.int(nrow(Lamberts.melted), size=N.stage1.Lamberts, replace=FALSE, prob=Lamberts.melted$value),]

