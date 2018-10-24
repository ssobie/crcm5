##Script to plot the EPW file series

library(ncdf4)
library(PCICt)
library(zoo)

sub.by.time <- function(var.name,lonc,latc,interval,input.file,gcm,read.dir) {

  print(input.file)              
  nc <- nc_open(paste(read.dir,gcm,'/',input.file,sep=''))
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.series <- origin.pcict + time.values*86400
  years <- format(time.series,'%Y')
  yrs <- strsplit(interval,'-')[[1]]

  feb.flag <- grep('-02-29',time.series)

  new.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal='365_day')
                           
  new.time0 <- (as.PCICt(format(time.series[1],'%Y-%m-%d'),cal='365_day') - new.origin)/86400
  if (length(feb.flag)==0) {
     new.values <- seq(as.numeric(new.time0),by=1,length.out=length(time.values))
  } else {
     new.values <- seq(as.numeric(new.time0),by=1,length.out=length(time.values[-feb.flag]))
  }

  new.series <- new.origin + new.values*86400
  print(range(time.series))
  print(range(new.series))

  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  
  lon.ix <- which.min(abs(lonc-lon))
  lat.ix <- which.min(abs(latc-lat))

  data.raw <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  data <- data.raw
  if (length(feb.flag!=0)) {
    data <- data.raw[-feb.flag]
  }

  years <- format(new.series,'%Y')
  st <- head(grep(yrs[1],years),1)
  en <- tail(grep(yrs[2],years),1)

  mon.data <- tapply(data[st:en],as.factor(format(new.series[st:en],'%m')),mean)

  nc_close(nc)
  return(mon.data)
}


##**************************************************************************************

epw.dir <- '/storage/home/ssobie/code/repos/bc-projected-weather/bcweather/tests/data/' 
present.epw.file <- 'CAN_BC_ABBOTSFORD-A_1100031_CWEC.epw'
future.epw.file <- 'abbotsford_test.epw'

epw.present.data <- read.csv(paste0(epw.dir,present.epw.file),skip=8,header=F,as.is=T)
epw.future.data <- read.csv(paste0(epw.dir,future.epw.file),skip=8,header=F,as.is=T)

##Create one year of daily dates
dates <- as.Date(paste('1999',sprintf('%02d',epw.present.data[,2]),sprintf('%02d',epw.present.data[,3]),sep='-'))

epw.daily.mean <- tapply(epw.present.data[,7],as.factor(format(dates,'%m-%d')),mean)
epw.daily.max <- tapply(epw.present.data[,7],as.factor(format(dates,'%m-%d')),max)
epw.daily.min <- tapply(epw.present.data[,7],as.factor(format(dates,'%m-%d')),min)

dy.dates <- as.Date(unique(format(dates,'%Y-%m-%d')))

epw.mon.mean <- tapply(epw.daily.mean,as.factor(format(dy.dates,'%m')),mean)
epw.mon.max <- tapply(epw.daily.max,as.factor(format(dy.dates,'%m')),mean)
epw.mon.min <- tapply(epw.daily.min,as.factor(format(dy.dates,'%m')),mean)

epw.day.anoms <- epw.present.data[,7]*0

for (d in 1:365) {
   ix <- format(dates,'%j') == sprintf('%03d',d)  
   epw.day.anoms[ix] <- epw.present.data[ix,7] - epw.daily.mean[d]
}

ix <- 7

plot.dir <- '/storage/data/projects/rci/data/cas/wx_files/'
plot.file <- paste0(plot.dir,'abbotsford.epw.tas.jan.png')

if (1==0) {
png(plot.file,width=1200,height=400)
par(mar=c(4.5,4.5,4,2))

plot(1:31,epw.daily.mean[1:31],type='l',lwd=4,xlab='Julian Day',ylab='Daily Mean Temperature (\u00B0C)',
     cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Abbotsford Weather File January Temperature')
abline(h=0,col='gray')
lines(1:31,epw.daily.mean[1:31],lwd=4)

box(which='plot')
dev.off()
}

if (1==0) {

lon <- -122.36
lat <- 49.03

gcm.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'

gcm.list <- c('ACCESS1-0',
              'CanESM2',
              'CCSM4',
              'CNRM-CM5',
              'CSIRO-Mk3-6-0',
              'GFDL-ESM2G',
              'HadGEM2-CC',
              'HadGEM2-ES',
              'inmcm4',
              'MIROC5',
              'MPI-ESM-LR',
              'MRI-CGCM3')

deltas <- matrix(0,nrow=12,ncol=12)
alphas <- matrix(0,nrow=12,ncol=12)

for (g in seq_along(gcm.list)) {
    gcm <- gcm.list[g]
    scen.files <- list.files(path=paste0(gcm.dir,gcm),pattern='rcp85')
    var.files <- scen.files[grep("tasmax",scen.files)]
    past.tx.file <- var.files[grep("1951-2000",var.files)]
    proj.tx.file <- var.files[grep("2001-2100",var.files)]

    print('Past TASMAX')
    past.tx.data <- sub.by.time(var.name='tasmax',lonc=lon,latc=lat,
                            interval='1971-2000',
                            input.file=past.tx.file,gcm=gcm,read.dir=gcm.dir)
    print('Proj TASMAX')
    proj.tx.data <- sub.by.time(var.name='tasmax',lonc=lon,latc=lat,
                            interval='2041-2070',
                            input.file=proj.tx.file,gcm=gcm,read.dir=gcm.dir)

    var.files <- scen.files[grep("tasmin",scen.files)]
    past.tn.file <- var.files[grep("1951-2000",var.files)]
    proj.tn.file <- var.files[grep("2001-2100",var.files)]

    print('Past TASMIN')
    past.tn.data <- sub.by.time(var.name='tasmin',lonc=lon,latc=lat,
                            interval='1971-2000',
                            input.file=past.tn.file,gcm=gcm,read.dir=gcm.dir)
    print('Proj TASMIN')
    proj.tn.data <- sub.by.time(var.name='tasmin',lonc=lon,latc=lat,
                            interval='2041-2070',
                            input.file=proj.tn.file,gcm=gcm,read.dir=gcm.dir)

    delta_tasmax <- proj.tx.data - past.tx.data
    delta_tasmin <- proj.tn.data - past.tn.data
    delta_tas <- (proj.tx.data+proj.tn.data)/2 - (past.tx.data+past.tn.data)/2
    deltas[g,] <- delta_tas
    print('Delta')
    print(delta_tas)

    alpha <- (delta_tasmax - delta_tasmin) / (epw.mon.max - epw.mon.min)
    print('Alpha')
    print(alpha)        
    alphas[g,] <- alpha
}

}


ens.deltas <- apply(deltas,2,mean)
ens.alphas <- apply(alphas,2,mean)

##morphed.tas <- epw.present.data[,7]*0
morphed.tas <- matrix(0,nrow=12,ncol=length(epw.present.data[,7]*0))

for (g in 1:12) {
  for (m in 1:12) {
     ix <- format(dates,'%m') == sprintf('%02d',m)  
     morphed.tas[g,ix] <- epw.present.data[ix,7] + deltas[g,m] + alphas[g,m]*epw.day.anoms[ix]
     ##morphed.tas[g,ix] <- deltas[g,m] + alphas[g,m]*epw.day.anoms[ix]
  }
}

daily.morphed.tas <- t(apply(morphed.tas,1,function(x,dates){tapply(x,as.factor(format(dates,'%m-%d')),mean)},dates))

hour.dates <- strftime(paste('1999-',sprintf('%02d',epw.present.data[,2]),'-',sprintf('%02d',epw.present.data[,3]),' ', 
                          sprintf('%02d',epw.present.data[,4]),':00:00', sep=''),format='%Y-%m-%d %H:%M:%S')

mn <- '07'
hr.dates <- hour.dates[grep(paste0('-',mn,'-'),hour.dates)]
hr.epw <- epw.present.data[grep(paste0('-',mn,'-'),dates),7]
hr.morph <- morphed.tas[,grep(paste0('-',mn,'-'),dates)]
hr.anoms <- epw.day.anoms[grep(paste0('-',mn,'-'),dates)]

if (1==0) {
plot.file <- paste0(plot.dir,'abbotsford.jul.hourly.tas.png')

png(plot.file,width=1200,height=400)
par(mar=c(4.5,4.5,4,2))
plot(1:length(hr.morph),hr.morph,type='l',lwd=4,xlab='Day',ylab='Temperature (\u00B0C)',
     cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Abbotsford July Hourly Temperature',ylim=c(10,35),col='white',axes=F)
axis(1,at=seq(0,744,24)[-32],labels=format(dy.dates[grep(paste0('-',mn,'-'),dy.dates)],'%d'),cex.axis=1.5,cex=1.5)
axis(2,at=seq(-15,35,5),labels=seq(-15,35,5),cex.axis=1.5,cex=1.5)
abline(h=seq(-15,35,5),lty=2,col='gray',lwd=2)
abline(h=0,col='gray')
lines(1:length(hr.morph),hr.epw,lwd=4,col='orange')
lines(1:length(hr.morph),hr.morph,lwd=4,col='red')

legend('topleft',leg=c('Future','Present'),col=c('red','orange'),pch=15,cex=1.5)

box(which='plot')
dev.off()
}

if (1==0) {
plot.file <- paste0(plot.dir,'examples.abbotsford.jan.diurnal.tas.png')
png(plot.file,width=1200,height=800)
par(mar=c(4.5,4.5,4,2))
plot(1:24,hr.morph[1:24],type='l',lwd=4,xlab='Day',ylab='Temperature (\u00B0C)',
     cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Abbotsford Jan Hourly Temperature',ylim=c(-15,15),col='white',axes=F)
axis(1,at=1:24,labels=1:24,cex.axis=1.5,cex=1.5)
axis(2,at=seq(-15,35,5),labels=seq(-15,35,5),cex.axis=1.5,cex=1.5)
abline(h=seq(-15,35,5),lty=2,col='gray',lwd=2)
abline(h=0,col='gray')

slen <- length(unique( format(as.Date(hr.dates),'%d')))
for (dy in 1:slen) {
  ix <- sprintf('%02d',dy) == format(as.Date(hr.dates),'%d')
  lines(1:24,hr.epw[ix],lwd=2)
}

box(which='plot')
dev.off()

}


plot.file <- paste0(plot.dir,'examples.abbotsford.jul1.diurnal.tas.png')
png(plot.file,width=1200,height=800)
par(mar=c(4.5,4.5,4,2))
par(mfrow=c(3,1))

plot(1:24,hr.epw[1:24],type='l',lwd=4,xlab='Hour',ylab='Temperature (\u00B0C)',
     cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Abbotsford Jul 1 Hourly Temperature',ylim=c(10,26),col='white',axes=F)
axis(1,at=1:24,labels=1:24,cex.axis=1.5,cex=1.5)
axis(2,at=seq(-40,40,2),labels=seq(-40,40,2),cex.axis=1.5,cex=1.5)
abline(h=seq(-40,40,2),lty=2,col='gray',lwd=2)
lines(1:24,hr.epw[1:24],lwd=4)
box(which='plot')

plot(1:24,alphas[2,as.numeric(mn)]*hr.anoms[1:24],type='l',lwd=4,xlab='Hour',ylab='Temperature (\u00B0C)',
     cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Abbotsford Jul Alphas',ylim=c(-0.8,0.6),col='white',axes=F)
axis(1,at=1:24,labels=1:24,cex.axis=1.5,cex=1.5)
axis(2,at=seq(-1,1,0.25),labels=seq(-1,1,0.25),cex.axis=1.5,cex=1.5)
abline(h=seq(-1,1,0.25),lty=2,col='gray',lwd=2)
for (i in 1:12) {
  lines(alphas[i,as.numeric(mn)]*hr.anoms[1:24],lwd=2)
}
box(which='plot')

plot(1:24,hr.morph[2,1:24],type='l',lwd=4,xlab='Hour',ylab='Temperature (\u00B0C)',
     cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Abbotsford Jul 1 Hourly Morphed',ylim=c(10,26),col='white',axes=F)
axis(1,at=1:24,labels=1:24,cex.axis=1.5,cex=1.5)
axis(2,at=seq(-40,40,2),labels=seq(-40,40,2),cex.axis=1.5,cex=1.5)
abline(h=seq(-40,40,2),lty=2,col='gray',lwd=2)
for (i in 1:12) {
  lines(hr.morph[i,1:24],lwd=2)
}

box(which='plot')
dev.off()





browser()

plot.file <- paste0(plot.dir,'abbotsford.tas.anomalies.smoothed.11.png')


png(plot.file,width=1200,height=400)
par(mar=c(4.5,4.5,4,2))

plot(1:365,apply(daily.morphed.tas,2,mean),type='l',lwd=4,xlab='Julian Day',ylab='Daily Mean Temperature Anomalies (\u00B0C)',
     cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Abbotsford Temperature Anomalies',ylim=c(0,6.2),col='white')
abline(h=0,col='gray')
lines(6:360,rollmean(apply(daily.morphed.tas,2,mean),11),lwd=4,col='orange')
lines(6:360,rollmean(apply(daily.morphed.tas,2,quantile,0.1),11),lwd=4,col='gold')
lines(6:360,rollmean(apply(daily.morphed.tas,2,quantile,0.9),11),lwd=4,col='red')

box(which='plot')
dev.off()

