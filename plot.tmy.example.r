##Script to plot the EPW file series

library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

pull.time.series <- function(var.name,lonc,latc,input.file) {

  print(input.file)              
  nc <- nc_open(input.file)
  time.series <- netcdf.calendar(nc)

  rlon <- ncvar_get(nc,'rlon')
  lat <- ncvar_get(nc,'lat')
  
  lon.ix <- 23 ##which.min(abs(lonc-lon))
  lat.ix <- 13 ##which.min(abs(latc-lat))

  data <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))

  day.fac <- as.factor(format(time.series,'%m-%d'))
  day.data <- tapply(data,day.fac,mean)
  mon.fac <- as.factor(format(time.series,'%m'))
  mon.data <- tapply(data,mon.fac,mean)
  nc_close(nc)
  return(list(hour=data,day=day.data,month=mon.data))

}

##----------------------------------------------------
var.name <- 'tas'
lonc <- 235.60
latc <- 48.93

read.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/tmy_files/'

past.file <- paste0(read.dir,'tas_CWEC_TMY_CanESM2+CRCM5_historical_1981-2010.nc')
proj.file <- paste0(read.dir,'tas_CWEC_TMY_CanESM2+CRCM5_historical_2021-2050.nc')
morph.file <- paste0(read.dir,'morphed_tas_CWEC_TMY_ERA-Interim+CanESM2+CRCM5_historical+rcp85_2021-2050.nc')

nc <- nc_open(past.file)
time.series <- netcdf.calendar(nc)
hours <- as.Date(format(time.series,'%Y-%m-%d %h'))
days <- as.Date(levels(as.factor(format(time.series,'%Y-%m-%d'))))

months <- as.Date(paste0(levels(as.factor(format(time.series,'%Y-%m'))),'-01'))
nc_close(nc)

past.series <- pull.time.series(var.name,lonc,latc,past.file)
proj.series <- pull.time.series(var.name,lonc,latc,proj.file)
morph.series <- pull.time.series(var.name,lonc,latc,morph.file)

par(mfrow=c(3,2))

plot(hours,past.series$hour,type='l',lwd=2,ylim=range(c(past.series$hour,morph.series$hour)))
lines(hours,morph.series$hour,type='l',lwd=2,col='red')
plot(hours,morph.series$hour-past.series$hour,type='l',lwd=2,col='green')

plot(days,past.series$day,type='l',lwd=2,ylim=range(c(past.series$day,morph.series$day)))
lines(days,morph.series$day,type='l',lwd=2,col='red')
plot(days,morph.series$day-past.series$day,type='l',lwd=2,col='green')

plot(months,past.series$month,type='l',lwd=2,ylim=range(c(past.series$month,morph.series$month)))
lines(months,morph.series$month,type='l',lwd=2,col='red')
plot(months,morph.series$month-past.series$month,type='l',lwd=2,col='green')


