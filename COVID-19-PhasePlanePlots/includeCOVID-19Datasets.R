# includeCOVID-19Datasets.R
# Revised for RStudio, 15th September 2020, 
# bayesianlogic.1@gmail.com


rtCOVID<- function(fn)
              read.table(file=fn, sep=",", header=TRUE, quote="\"", row.names=NULL, stringsAsFactors=FALSE, check=TRUE, fill=FALSE)

tsc<- rtCOVID("./replica-covid-19-jhsse/time_series_covid19_confirmed_global.csv")
tsd<- rtCOVID("./replica-covid-19-jhsse/time_series_covid19_deaths_global.csv")
tsr<- rtCOVID("./replica-covid-19-jhsse/time_series_covid19_recovered_global.csv")

tsc.us<- rtCOVID("./replica-covid-19-jhsse/time_series_covid19_confirmed_US.csv")
tsd.us<- rtCOVID("./replica-covid-19-jhsse/time_series_covid19_deaths_US.csv")


#tsc<- rtCOVID("time_series_covid19_confirmed_global.csv")
#tsd<- rtCOVID("time_series_covid19_deaths_global.csv")
#tsr<- rtCOVID("time_series_covid19_recovered_global.csv")

#tsc.us<- rtCOVID("time_series_covid19_confirmed_US.csv")
#tsd.us<- rtCOVID("time_series_covid19_deaths_US.csv")

