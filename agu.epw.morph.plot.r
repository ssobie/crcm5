##Script to plot the EPW file series

library(ncdf4)
library(PCICt)
library(zoo)
library(scales)
source('/storage/home/ssobie/code/repos/crcm5/read.write.epw.r',chdir=T)

##**************************************************************************************

epw.dir <- '/storage/data/projects/rci/weather_files/wx_files/'
morph.dir <- '/storage/data/projects/rci/weather_files/wx_files/morphed_files/tas_only/'

present.epw.file <- 'CAN_BC_VANCOUVER-INTL-A_1108395_CWEC.epw'
future.2020s <- 'MORPHED_ROLL21_TAS_VANCOUVER-INTL-A_1108395_2011-2040_CWEC.epw'
future.2050s <- 'MORPHED_ROLL21_TAS_VANCOUVER-INTL-A_1108395_2041-2070_CWEC.epw'
future.2080s <- 'MORPHED_ROLL21_TAS_VANCOUVER-INTL-A_1108395_2071-2100_CWEC.epw'

plot.dir <- '/storage/data/projects/rci/weather_files/'

epw.present <- read.epw.file(epw.dir,present.epw.file)
epw.2020s <- read.epw.file(morph.dir,future.2020s)
epw.2050s <- read.epw.file(morph.dir,future.2050s)
epw.2080s <- read.epw.file(morph.dir,future.2080s)

##Create one year of daily dates


dates <- as.Date(paste('1999',sprintf('%02d',epw.present$data[,2]),sprintf('%02d',epw.present$data[,3]),sep='-'))
hour.dates <- strftime(paste('1999-',sprintf('%02d',epw.present$data[,2]),'-',sprintf('%02d',epw.present$data[,3]),' ', 
                       sprintf('%02d',epw.present$data[,4]),':00:00', sep=''),format='%Y-%m-%d %H:%M:%S')
dy.dates <- as.Date(unique(format(dates,'%Y-%m-%d')))

epw.tas <- epw.present$data[,7]
##   epw.agg.mean <- epw.daily.mean <- make_average_series(epw.tas,dates,method='daily',agg.fxn=mean)
##   epw.day.anoms <- epw.tas*0
##   ##Daily
##   for (d in 1:365) {
##      ix <- format(dates,'%j') == sprintf('%03d',d)
##      epw.day.anoms[ix] <- epw.tas[ix] - epw.daily.mean[d]
##   }
  ##--------------
  mn <- '-06-|-07-|-08-'
  hr.dates <- hour.dates[grep(mn,hour.dates)]
  hr.dy <- dy.dates[grep(mn,dy.dates)]
  hr.epw <- epw.present$data[grep(mn,dates),7]

  hr.2020s <- epw.2020s$data[grep(mn,dates),7]
  hr.2050s <- epw.2050s$data[grep(mn,dates),7]
  hr.2080s <- epw.2080s$data[grep(mn,dates),7]

##  hr.morph <- morphed.tas[,grep(mn,dates)]
##  hr.anoms <- epw.day.anoms[grep(mn,dates)]
  
  mn.max <- tapply(hr.epw,as.factor(format(as.Date(hr.dates),'%m-%d')),max)
  mn.min <- tapply(hr.epw,as.factor(format(as.Date(hr.dates),'%m-%d')),min)
  mn.mean <- tapply(hr.epw,as.factor(format(as.Date(hr.dates),'%m-%d')),mean)

  mn.max.2080s <- tapply(hr.2080s,as.factor(format(as.Date(hr.dates),'%m-%d')),max)
  mn.min.2080s <- tapply(hr.2080s,as.factor(format(as.Date(hr.dates),'%m-%d')),min)
  mn.mean.2080s <- tapply(hr.2080s,as.factor(format(as.Date(hr.dates),'%m-%d')),mean)


  plot.file <- paste0(plot.dir,'agu.examples.abbotsford.both.diurnal.tas.png')
  png(plot.file,width=1000,height=600)
  par(mfrow=c(2,1))
  par(mar=c(4.5,4.5,4,2))
  days <- format(as.Date(hr.dy),'%b-%d')
  plot(1:length(days),mn.mean,type='l',lwd=4,xlab='Day',ylab='Temperature (\u00B0C)',xaxs='i',
       cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Vancouver Summer Temperature',ylim=c(9,33),col='white',axes=F)
  dx <- c(1,10,20,31,40,50,60,71,81,92)
  axis(1,at=dx,labels=days[dx],cex.axis=1.5,cex=1.5)
  axis(2,at=seq(-40,40,4),labels=seq(-40,40,4),cex.axis=1.5,cex=1.5)
  abline(h=seq(-40,40,4),lty=2,col='gray',lwd=1)
  polygon(c(1:length(days),rev(1:length(days))),
          c(mn.min,rev(mn.max)),col=alpha('gray',0.4),border=alpha('gray',0.4))
  polygon(c(1:length(days),rev(1:length(days))),
          c(mn.min.2080s,rev(mn.max.2080s)),col=alpha('red',0.3),border=alpha('red',0.3))
  lines(1:length(days),mn.mean.2080s,lwd=3,col='red')
  lines(1:length(days),mn.mean,lwd=3)
  legend('topleft',leg=c('Current','Morphed'),col=c('black','red'),cex=1.75,pch=15)
  box(which='plot')

  ##------------------------------
  mn <- '-12-|-01-|-02-'
  hr.dates <- hour.dates[grep(mn,hour.dates)]
  hr.dates <- c(hr.dates[1416:2160],hr.dates[1:1415])
  hr.dy <- as.Date(unique(format(as.Date(hr.dates),'%Y-%m-%d')))
  hr.epw <- epw.present$data[grep(mn,dates),7]
  hr.epw <- c(hr.epw[1416:2160],hr.epw[1:1415])
  days <- unique(format(as.Date(hr.dy),'%b-%d'))

  hr.2020s <- epw.2020s$data[grep(mn,dates),7]
  hr.2050s <- epw.2050s$data[grep(mn,dates),7]
  hr.2080s <- epw.2080s$data[grep(mn,dates),7]
  hr.2080s <- c(hr.2080s[1416:2160],hr.2080s[1:1415])

##  hr.morph <- morphed.tas[,grep(mn,dates)]
##  hr.anoms <- epw.day.anoms[grep(mn,dates)]
  
  mn.max <- tapply(hr.epw,as.factor(format(as.Date(hr.dates),'%m-%d')),max)
  mn.min <- tapply(hr.epw,as.factor(format(as.Date(hr.dates),'%m-%d')),min)
  mn.mean <- tapply(hr.epw,as.factor(format(as.Date(hr.dates),'%m-%d')),mean)

  mn.max.2080s <- tapply(hr.2080s,as.factor(format(as.Date(hr.dates),'%m-%d')),max)
  mn.min.2080s <- tapply(hr.2080s,as.factor(format(as.Date(hr.dates),'%m-%d')),min)
  mn.mean.2080s <- tapply(hr.2080s,as.factor(format(as.Date(hr.dates),'%m-%d')),mean)

  plot(1:length(days),mn.mean,type='l',lwd=4,xlab='Day',ylab='Temperature (\u00B0C)',xaxs='i',
       cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Vancouver Winter Temperature',ylim=c(-5,18),col='white',axes=F) ##ylim=c(9,33)
  dx <- c(1,10,20,32,41,51,63,71,82,90)
  axis(1,at=dx,labels=days[dx],cex.axis=1.5,cex=1.5)
  axis(2,at=seq(-40,40,4),labels=seq(-40,40,4),cex.axis=1.5,cex=1.5)
  abline(h=seq(-40,40,4),lty=2,col='gray',lwd=1)
  polygon(c(1:length(days),rev(1:length(days))),
          c(mn.min,rev(mn.max)),col=alpha('gray',0.4),border=alpha('gray',0.4))
  polygon(c(1:length(days),rev(1:length(days))),
          c(mn.min.2080s,rev(mn.max.2080s)),col=alpha('red',0.3),border=alpha('red',0.3))
  lines(1:length(days),mn.mean.2080s,lwd=3,col='red')
  lines(1:length(days),mn.mean,lwd=3)
  box(which='plot')







  dev.off()



browser()


