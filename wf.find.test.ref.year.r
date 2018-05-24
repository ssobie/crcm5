##Script to find the test reference year at each grid cell for a provided netcdf file

Rprof('wx.out')

library(ncdf4)
library(PCICt)
library(doParallel)
registerDoParallel(cores=4)
library(foreach)
library(data.table)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

find.closest.month <- function(month,data,mon.ts,mon.indices,yrs) { 
  print(month)
  mon.ix <- mon.indices[,as.numeric(month)]
  mon.vl <- which(mon.ix)  
  
  mon.data <- data[mon.ix]
  long.cdf <- ecdf(mon.data)

  cdf.diffs <- rep(0,length(yrs))

##  plot(mon.data,long.cdf(mon.data))

  for (i in seq_along(yrs)) {
    ##yr.mon.ix <- grep(paste0(yrs[i],'-',month),mon.ts)
    ##yr.mon.ix <- mon.ts %in% paste0(yrs[i],'-',month)
    ##yr.mon.ix <- mon.ts %chin% paste0(yrs[i],'-',month)
    yr.mon.ix <- mon.vl[mon.ts[mon.ix] %chin% paste0(yrs[i],'-',month)]
    ##print(which(yr.mon.ix)-test.mon.ix)
    yr.mon <- data[yr.mon.ix]
    short.cdf <- ecdf(yr.mon)
    cdf.diffs[i] <- mean(abs(long.cdf(yr.mon) - short.cdf(yr.mon)))
##    points(yr.mon,short.cdf(yr.mon),col='blue')
  }

  rv <- yrs[which.min(cdf.diffs)]

##  slct.ix <- mon.ts %in% paste0(rv,'-',month)
##  slct.mon <- data[slct.ix]
##  slct.cdf <- ecdf(slct.mon)
##  points(mon.data,long.cdf(mon.data))
##  points(slct.mon,slct.cdf(slct.mon),col='red')

  return(rv)
}

find.all.months <- function(months,data,mon.ts,mon.indices,yrs) { 
  rv <- rep(0,length(months))
  for (m in seq_along(months)) {
    month <- months[m]
    ##mon.ix <- mon.fac %in% month
    ##mon.ix <- mon.fac %chin% month
    mon.ix <- mon.indices[,m]
    mon.vl <- which(mon.ix)  

    mon.data <- data[mon.ix]
    long.cdf <- ecdf(mon.data)

    cdf.diffs <- rep(0,length(yrs))

    for (i in seq_along(yrs)) {
      ##yr.mon.ix <- mon.ts %in% paste0(yrs[i],'-',month)
      ##yr.mon.ix <- mon.ts %chin% paste0(yrs[i],'-',month)
      yr.mon.ix <- mon.vl[mon.ts[mon.ix] %chin% paste0(yrs[i],'-',month)]
      yr.mon <- data[yr.mon.ix]
      short.cdf <- ecdf(yr.mon)
      cdf.diffs[i] <- mean(abs(long.cdf(yr.mon) - short.cdf(yr.mon)))
    }
    rv[m] <- yrs[which.min(cdf.diffs)]
  }

  return(rv)
}

##----------------------------------------------------------------------------------------------

read.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/hourly/'

tas.file <- paste0(read.dir,'tas_hour_WC011_ERA-Interim+CRCM5_historical_19800101-20141231.nc')

months <- sprintf('%02d',1:12)
mlen <- length(months)

nc <- nc_open(tas.file)
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



##par(mfrow=c(3,4))

tas <- ncvar_get(nc,'tas',start=c(1,1,1),count=c(1,1,-1))

yr.match <- unlist(lapply(months,find.closest.month,tas,mon.ts,mon.indices,yrs) )

##browser()

input.nc <- nc
varname <- 'tas'

years.array <- array(0,c(n.lon,n.lat,12))

for (j in 1:n.lat) { ##n.lon) {
  print(paste0('Latitude: ',j,' of ',n.lat))
  input.subset <- ncvar_get(input.nc,varname,start=c(1,j,1),count=c(-1,1,-1))
  flag <- is.na(input.subset[,1])
  input.list <- lapply(seq_len(nrow(input.subset)), function(k) input.subset[k,])
  rm(input.subset)

  years.slct <- foreach(
                  input=input.list,                 
                  .export=c('find.all.months','months','mon.ts','mon.indices','yrs')
                           ) %dopar% {
                              objects <- find.all.months(months=months,data=input,mon.ts=mon.ts,mon.indices=mon.indices,yrs=yrs)
                         }
  years.matrix <- matrix(unlist(years.slct),nrow=n.lon,ncol=length(months),byrow=T)
  years.array[,j,] <- years.matrix  
}



Rprof(NULL)