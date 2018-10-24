##Script to find the test reference year at each grid cell for a provided netcdf file

library(ncdf4)
library(PCICt)
library(data.table)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

find.closest.month <- function(month,data,mon.ts,mon.indices,yrs,trys) { 
  print(month)
  mon.ix <- mon.indices[,as.numeric(month)]
  mon.vl <- which(mon.ix)  
  
  mon.data <- data[mon.ix]
  long.cdf <- ecdf(mon.data)

  cdf.diffs <- rep(0,length(yrs))

  plot(mon.data,long.cdf(mon.data),main=month.abb[as.numeric(month)],xlab='units',ylab='CDF')

  for (i in seq_along(yrs)) {
    yr.mon.ix <- mon.vl[mon.ts[mon.ix] %chin% paste0(yrs[i],'-',month)]
    yr.mon <- data[yr.mon.ix]
    short.cdf <- ecdf(yr.mon)
    points(yr.mon,short.cdf(yr.mon),col='blue')
  }

  slct.ix <- mon.ts %in% paste0(trys[as.numeric(month)],'-',month)
  slct.mon <- data[slct.ix]
  slct.cdf <- ecdf(slct.mon)
  points(mon.data,long.cdf(mon.data))
  points(slct.mon,slct.cdf(slct.mon),col='red')

}

##----------------------------------------------------------------------------------------------

var.name <- 'insol'
coord.ix <- c(10,10)

read.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/hourly/'

load('/storage/data/climate/downscale/RCM/CRCM5/reconfig/daily/year_data/cwec.tmy.years.from.daily.RData')


data.file <- paste0(read.dir,var.name,'_hour_WC011_ERA-Interim+CRCM5_historical_19800101-20141231.nc')

months <- sprintf('%02d',1:12)
mlen <- length(months)

nc <- nc_open(data.file)
time <- netcdf.calendar(nc)
mon.ts <- format(time,'%Y-%m')
mon.fac <- format(time,'%m')
yrs <- levels(as.factor(format(time,'%Y')))
n.lat <- nc$dim$rlat$len ##Latitude Length
n.lon <- nc$dim$rlon$len ##Longitude Length
n.time <- nc$dim$time$len ##Longitude Length

mon.indices <- matrix(NA,nrow=n.time,ncol=mlen)
for (i in 1:mlen) {
  mon.indices[,i] <- mon.fac %chin% months[i]
}

ymlen <- length(unique(mon.ts))
yr.mons <- unique(mon.ts)
yr.mon.indices <- matrix(NA,nrow=n.time,ncol=ymlen)
for (i in 1:ymlen) {
  yr.mon.indices[,i] <- mon.ts %chin% yr.mons[i]
}

par(mfrow=c(3,4))

series <- ncvar_get(nc,var.name,start=c(coord.ix,1),count=c(1,1,-1))
trys <- tmy.years[coord.ix[1],coord.ix[2],]

yr.match <- unlist(lapply(months,find.closest.month,series,mon.ts,mon.indices,yrs,trys) )



