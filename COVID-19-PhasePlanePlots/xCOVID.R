# Extract and build COVID-17 deaths, case, and some recovery data.
# Uses the COVID-19-JHSSE formats.  Pulled from their Github, 
# https://github.com/CSSEGISandData/COVID-19.git.
#
# Jan Galkowski. 2nd May 2020.
#
# Last changed 30th October 2020.
#


fstNonzeroOf<- function(v)
{
  labeling<- names(v)
  z<- as.vector(as.numeric(v))
  fst<- which(0 < z)[1]
  n<- length(z)
  z<- z[fst:n]
  names(z)<- labeling[fst:n]
  return(z)
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

#source("PhPlPl-COVID19_Gibbs.R")
##source("PhPlPl-COVID19_statespacer.R")
#source("PhPlPl-COVID19_KFAS.R")
#source("PhPlPl-COVID19_nloptr.R")
#source("PhPlPl-COVID19_rf.R")
#source("PhPlPl-COVID19_afo.R")
#source("PhPlPl-COVID19_scp.R")
source("PhPlPl-COVID19_Denham.R")

# Global data has a format like this:


#> load("20200502bGlobal.RData")
#> head(tsc[,1:12])
#  Province.State      Country.Region      Lat     Long X1.22.20 X1.23.20 X1.24.20 X1.25.20 X1.26.20 X1.27.20 X1.28.20 X1.29.20
#1           <NA>         Afghanistan  33.0000  65.0000        0        0        0        0        0        0        0        0
#2           <NA>             Albania  41.1533  20.1683        0        0        0        0        0        0        0        0
#3           <NA>             Algeria  28.0339   1.6596        0        0        0        0        0        0        0        0
#4           <NA>             Andorra  42.5063   1.5218        0        0        0        0        0        0        0        0
#5           <NA>              Angola -11.2027  17.8739        0        0        0        0        0        0        0        0
#6           <NA> Antigua and Barbuda  17.0608 -61.7964        0        0        0        0        0        0        0        0
#> head(tsr[,1:12])
#  Province.State      Country.Region      Lat     Long X1.22.20 X1.23.20 X1.24.20 X1.25.20 X1.26.20 X1.27.20 X1.28.20 X1.29.20
#1           <NA>         Afghanistan  33.0000  65.0000        0        0        0        0        0        0        0        0
#2           <NA>             Albania  41.1533  20.1683        0        0        0        0        0        0        0        0
#3           <NA>             Algeria  28.0339   1.6596        0        0        0        0        0        0        0        0
#4           <NA>             Andorra  42.5063   1.5218        0        0        0        0        0        0        0        0
#5           <NA>              Angola -11.2027  17.8739        0        0        0        0        0        0        0        0
#6           <NA> Antigua and Barbuda  17.0608 -61.7964        0        0        0        0        0        0        0        0
#> head(tsd[,1:12])
#  Province.State      Country.Region      Lat     Long X1.22.20 X1.23.20 X1.24.20 X1.25.20 X1.26.20 X1.27.20 X1.28.20 X1.29.20
#1           <NA>         Afghanistan  33.0000  65.0000        0        0        0        0        0        0        0        0
#2           <NA>             Albania  41.1533  20.1683        0        0        0        0        0        0        0        0
#3           <NA>             Algeria  28.0339   1.6596        0        0        0        0        0        0        0        0
#4           <NA>             Andorra  42.5063   1.5218        0        0        0        0        0        0        0        0
#5           <NA>              Angola -11.2027  17.8739        0        0        0        0        0        0        0        0
#6           <NA> Antigua and Barbuda  17.0608 -61.7964        0        0        0        0        0        0        0        0
#>
#

# United States data has a format like this, and lacks recovery data:

#> head(tsd.us[,1:13])
#       UID iso2 iso3 code3 FIPS  Admin2           Province_State Country_Region          Lat         Long_                 Combined_Key Population X1.22.20
#1       16   AS  ASM    16   60                   American Samoa             US -14.27100000 -170.13200000           American Samoa, US      55641        0
#2      316   GU  GUM   316   66                             Guam             US  13.44430000  144.79370000                     Guam, US     164229        0
#3      580   MP  MNP   580   69         Northern Mariana Islands             US  15.09790000  145.67390000 Northern Mariana Islands, US      55144        0
#4      630   PR  PRI   630   72                      Puerto Rico             US  18.22080000  -66.59010000              Puerto Rico, US    2933408        0
#5      850   VI  VIR   850   78                   Virgin Islands             US  18.33580000  -64.89630000           Virgin Islands, US     107268        0
#6 84001001   US  USA   840 1001 Autauga                  Alabama             US  32.53952745  -86.64408227         Autauga, Alabama, US      55869        0
#> head(tsc.us[,1:13])
#       UID iso2 iso3 code3 FIPS  Admin2           Province_State Country_Region          Lat         Long_                 Combined_Key X1.22.20 X1.23.20
#1       16   AS  ASM    16   60                   American Samoa             US -14.27100000 -170.13200000           American Samoa, US        0        0
#2      316   GU  GUM   316   66                             Guam             US  13.44430000  144.79370000                     Guam, US        0        0
#3      580   MP  MNP   580   69         Northern Mariana Islands             US  15.09790000  145.67390000 Northern Mariana Islands, US        0        0
#4      630   PR  PRI   630   72                      Puerto Rico             US  18.22080000  -66.59010000              Puerto Rico, US        0        0
#5      850   VI  VIR   850   78                   Virgin Islands             US  18.33580000  -64.89630000           Virgin Islands, US        0        0
#6 84001001   US  USA   840 1001 Autauga                  Alabama             US  32.53952745  -86.64408227         Autauga, Alabama, US        0        0
#> 
#

# Also, there are annoying differences: In global, for example, Province.State is spelled like that, but
# in United States data is is spelled Province_State. Same for Country.Region versus Country_Region. 
# Also, the values of the Province_State fields are right-justified in a fixed width field having leading 
# blanks. 
#
# Because of this kind of mess, felt it best to write a script to prepare these from the JSSE data 
# rather than doing it by hand. None of this, of course, can compensate for other typical irregularities,
# such as incomplete counting of deaths (because causes have not formally been determined: the USA doesn't
# automatically do autopsies) and irregular points of reporting where, say, deaths are accumulated and then
# reported all at once.
#

#Countries<- c("Iceland", "Switzerland", "Germany", "Sweden", "China", "US", "Taiwan*", "United Kingdom")
# ("China" excluded because shows singularity in solving smoothing equations.)
Countries<- c("Iceland", "Switzerland", "Germany", "Sweden", "US", "Taiwan*", "United Kingdom")
#Countries<- c("Iceland", "Switzerland", "Sweden", "US", "Taiwan*", "United Kingdom")
#Countries<- c("Germany")

USStates<- c("New York", "Washington", "Michigan", "Massachusetts", "Georgia", "Florida", "Tennessee")
#USStates<- c("New York", "Massachusetts", "Michigan")

rowsOfHavingCountries<- function(tsx, countriesList=Countries, countryColumnName="Country.Region")
{
  rowList<- list()
  for (country in countriesList)
  {
    p<- grep(pattern=country, x=tsx[,countryColumnName], value=FALSE)
    stopifnot( 0 < length(p) )
    rowList[[country]]<- tsx[p,]
  }
  return(rowList)
}

rowsOfHavingStates<- function(tsx, statesList=USStates, stateColumnName="Province_State")
{
  rowList<- list()
  for (state in statesList)
  {
    p<- grep(pattern=state, x=tsx[,stateColumnName], value=FALSE)
    stopifnot( 0 < length(p) )
    rowList[[state]]<- tsx[p,]
  }
  return(rowList)
}

collapseCountsToCumulative<- function(rowList, countsBeginAtColumn=5)
{
  # countsBeginAtColumn=13 for states
  lapply(X=rowList, FUN=function(dataset) fstNonzeroOf(colSums(dataset[, (countsBeginAtColumn:ncol(dataset))])))
}

source("includeCOVID-19Datasets.R")

# Setup the lot.

tsd.list<- rowsOfHavingCountries(tsx=tsd)
tsd.cum<- collapseCountsToCumulative(rowList=tsd.list, countsBeginAtColumn=5)
tsc.list<- rowsOfHavingCountries(tsx=tsc)
tsc.cum<- collapseCountsToCumulative(rowList=tsc.list, countsBeginAtColumn=5)
tsr.list<- rowsOfHavingCountries(tsx=tsr)
tsr.cum<- collapseCountsToCumulative(rowList=tsr.list, countsBeginAtColumn=5)

tsd.us.list<- rowsOfHavingStates(tsx=tsd.us)
tsd.us.cum<- collapseCountsToCumulative(rowList=tsd.us.list, countsBeginAtColumn=13)
tsc.us.list<- rowsOfHavingStates(tsx=tsc.us)
tsc.us.cum<- collapseCountsToCumulative(rowList=tsc.us.list, countsBeginAtColumn=13)

recodePlace<- function(place)
{
  places<- c("United Kingdom", "New York", "Washington", "Michigan", "Massachusetts", "Georgia", "Florida", "Taiwan*",
             "China", "Sweden", "Germany", "Tennessee", place)
  recoding<- c("UK", "NY", "WA", "MI", "MA", "GA", "FL", "TW", "CN", "SE", "DE", "TN", place)
  # The following by construction cannot fail. The "match" function scans left-to-right.
  return( recoding[match(place, places)] )
}

genTag<- function(place, kind) sprintf("Plots of %s COVID-19 %s", recodePlace(place), kind)
  
genOrdinate<- function(kind) sprintf("cumulative count of %s", kind)

genFst<- function(kind) sprintf("rate of change in number of %s", kind)

genSec<- function(kind) sprintf("acceleration in number of %s", kind)

genOne<- function(place, kind) sprintf("cumulative-%s-%s", recodePlace(place), kind)

genTwo<- function(place, kind) sprintf("cumulative-%s-%s-vs-rate", recodePlace(place), kind)

genThree<- function(place, kind) sprintf("%s-%s-rate-vs-acceleration", recodePlace(place), kind)

pp<- function(cum, place, kind, kth, what)
     { print(place) ;
       ppp(S=cum, what=what, ordinateLabel=genOrdinate(kind), kth=kth, 
           fstLabel=genFst(kind), secLabel=genSec(kind), svg=TRUE, one=genOne(place, kind),
           two=genTwo(place, kind), three=genThree(place,kind))
     }

#stop("noplot")

if (TRUE)
{
  invisible(mapply(FUN=function(cumulative, place, kth) 
                         pp(cum=cumulative, place=place, kind="deaths", kth=kth, what=place),
                   cumulative=tsd.cum, place=names(tsd.cum), kth=(1:length(tsd.cum))))
#  invisible(mapply(FUN=function(cumulative, place, kth) 
#                         pp(cum=cumulative, place=place, kind="cases", kth=kth, what=place),
#                   cumulative=tsc.cum, place=names(tsc.cum), kth=(1:length(tsc.cum))))
# invisible(mapply(FUN=function(cumulative, place, kth) 
#                        pp(cum=cumulative, place=place, kind="recovered", kth=kth, what=place),
#                  cumulative=tsr.cum, place=names(tsr.cum), kth=(1:length(tsr.cum))))
}

if (TRUE)
{
  invisible(mapply(FUN=function(cumulative, place, kth) 
                         pp(cum=cumulative, place=place, kind="deaths", kth=kth, what=place),
                   cumulative=tsd.us.cum, place=names(tsd.us.cum), kth=(1:length(tsd.us.cum))))
}

if (FALSE)
{
  invisible(mapply(FUN=function(cumulative, place, kth) 
                         pp(cum=cumulative, place=place, kind="cases", kth=kth, what=place),
                   cumulative=tsc.us.cum, place=names(tsc.us.cum), kth=(1:length(tsc.us.cum))))
}





