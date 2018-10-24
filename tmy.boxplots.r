##Script to plot the EPW file series

library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##Replace box edges with 10th and 90th percentiles
new.box <- function(at,data,add,boxwex,axes,col) { 
  bb <- boxplot(data,plot=FALSE)
  ##bb$stats[c(1,5),] <- quantile(data,c(0.1,0.9),na.rm=TRUE)
  bxp(z=bb,at=at,add=add,boxwex=boxwex,axes=axes,col=col)
}


pull.time.series <- function(var.name,lonc,latc,input.file,interval=NULL) {

  print(input.file)              
  nc <- nc_open(input.file)
  time.series <- netcdf.calendar(nc)

  rlon <- ncvar_get(nc,'rlon')
  lat <- ncvar_get(nc,'lat')
  
  lon.ix <- 23 ##which.min(abs(lonc-lon))
  lat.ix <- 13 ##which.min(abs(latc-lat))

  if (!is.null(interval)) {
    bnds <- strsplit(interval,'-')[[1]]
    st <- head(grep(bnds[1],time.series),1)
    en <- tail(grep(bnds[2],time.series),1)
    cnt <- en-st+1
    day.fac <- as.factor(format(time.series[st:en],'%m-%d'))     
    mon.fac <- as.factor(format(time.series[st:en],'%m'))     
    data <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,st),count=c(1,1,cnt))
  } else {  
    day.fac <- as.factor(format(time.series,'%m-%d'))     
    mon.fac <- as.factor(format(time.series,'%m'))     
    data <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  }

  day.vals <- tapply(data,day.fac,function(x){return(x)})
  mon.vals <- tapply(data,mon.fac,function(x){return(x)})
  nc_close(nc)
  return(day.vals)
}

##----------------------------------------------------
var.name <- 'tas'
lonc <- 235.60
latc <- 48.93

read.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/tmy_files/'

past.file <- paste0(read.dir,'tas_CWEC_TMY_CanESM2+CRCM5_historical_1981-2010.nc')
proj.file <- paste0(read.dir,'tas_CWEC_TMY_CanESM2+CRCM5_historical_2021-2050.nc')
morph.file <- paste0(read.dir,'morphed_by_week.nc')

nc <- nc_open(past.file)
time.series <- netcdf.calendar(nc)
hours <- as.Date(format(time.series,'%Y-%m-%d %h'))
days <- as.Date(levels(as.factor(format(time.series,'%Y-%m-%d'))))
 
months <- as.Date(paste0(levels(as.factor(format(time.series,'%Y-%m'))),'-01'))
nc_close(nc)

##past.series <- pull.time.series(var.name,lonc,latc,past.file)
proj.series <- pull.time.series(var.name,lonc,latc,proj.file)
morph.series <- pull.time.series(var.name,lonc,latc,morph.file)

rcm.file <- "/storage/data/climate/downscale/RCM/CRCM5/reconfig/hourly/tas_hour_WC011_CanESM2+CRCM5_historical+rcp85_19800101-20501231.nc"
full.series <- pull.time.series(var.name,lonc,latc,rcm.file,'2021-2050')
past.series <- pull.time.series(var.name,lonc,latc,rcm.file,'1981-2010')

yvals <- round(range(c(unlist(proj.series),unlist(morph.series),unlist(full.series))))
plot(c(),xlim=c(1,12),ylim=c(-20,20),main='Boxplots',
    xlab='',ylab='Tas Data',cex.axis=2,cex.lab=2,cex.main=2.5,axes=F,yaxs='i')
  axis(1,at=1:12,1:12,cex.axis=2)
  axis(2,at=seq(yvals[1],yvals[2],by=2),seq(yvals[1],yvals[2],by=2),cex.axis=2)

for (i in 1:12) {
    new.box(at=i-0.3,past.series[[i+12]],add=T,boxwex=0.25,axes=F,col='green')
    new.box(at=i-0.1,morph.series[[i+12]],add=T,boxwex=0.25,axes=F,col='green')
    new.box(at=i+0.1,proj.series[[i+12]],add=T,boxwex=0.25,axes=F,col='yellow')    
    new.box(at=i+0.3,full.series[[i+12]],add=T,boxwex=0.25,axes=F,col='red')
}

  abline(v=seq(0.5,31,by=1))
##  grid(ny=NULL,nx=NA,col='gray',lwd=1)
  box(which='plot')


browser()

par(mar=c(10,5,5,5))
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


