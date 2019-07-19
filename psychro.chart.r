##Script to plot psychro chart
library(scales)
source('read.write.epw.r')
source('rpckg.psy.chart.r')

hum_ratio_rel_hum <- function(rel.hum, temp.air, hum.ratio.max, alt = 0) {
  # Humidity ratio (eq. 25 and 12, solved for W)
  rel.hum <- rel.hum / 100
  W_s <- sat_hum_ratio(temp.air, alt)
  p <- bar_press(alt)
  p_ws <- sat_w_press(temp.air)

  W <- (-rel.hum) * W_s * (p - p_ws) / (rel.hum * p_ws - p)

  W <- ifelse(W > hum.ratio.max, NA, W)

  return(W)

}

##--------------------------------------------------------------
epw.dir <- '/storage/data/projects/rci/weather_files/wx_files/'
morph.dir <- '/storage/data/projects/rci/weather_files/wx_files/morphed_files/tas_only/'

present.epw.file <- 'CAN_BC_VANCOUVER-INTL-A_1108395_CWEC.epw'
future.2020s <- 'MORPHED_ROLL21_TAS_VANCOUVER-INTL-A_1108395_2011-2040_CWEC.epw'
future.2050s <- 'MORPHED_ROLL21_TAS_VANCOUVER-INTL-A_1108395_2041-2070_CWEC.epw'
future.2080s <- 'MORPHED_ROLL21_TAS_VANCOUVER-INTL-A_1108395_2071-2100_CWEC.epw'

epw.present <- read.epw.file(epw.dir,present.epw.file)
epw.2020s <- read.epw.file(morph.dir,future.2020s)
epw.2050s <- read.epw.file(morph.dir,future.2050s)
epw.2080s <- read.epw.file(morph.dir,future.2080s)

tas.past <-  epw.present$data[,7]
tas.2020s <- epw.2020s$data[,7]
tas.2050s <- epw.2050s$data[,7]
tas.2080s <- epw.2080s$data[,7]

rh.past <-  epw.present$data[,9]
rh.2020s <- epw.2020s$data[,9]
rh.2050s <- epw.2050s$data[,9]
rh.2080s <- epw.2080s$data[,9]

hr.past <- hum_ratio_rel_hum(rh.past,tas.past,hum.ratio.max=0.04)
hr.2020s <- hum_ratio_rel_hum(rh.2020s,tas.2020s,hum.ratio.max=0.04)
hr.2050s <- hum_ratio_rel_hum(rh.2050s,tas.2050s,hum.ratio.max=0.04)
hr.2080s <- hum_ratio_rel_hum(rh.2080s,tas.2080s,hum.ratio.max=0.04)

rh.lines <- hum_ratio_rel_hum(rep(100,200),seq(-15,40,length.out=200),hum.ratio.max=0.04)

plot.dir <- '/storage/data/projects/rci/weather_files/'
plot.file <- paste0(plot.dir,'agu.single.psychro2.png')
png(plot.file,width=1000,height=400)
if (1==0) {
par(mar=c(5.1,2.1,2.1,4.1))
par(mfrow=c(2,1))
plot(c(),xlim=c(-5,35),ylim=c(0,0.02),main='',xaxs='i',yaxs='i',
         xlab='Dry-bulb Temperature (\u00B0C)',ylab='',axes=F,cex.lab=1.75)
axis(1,at=seq(-20,40,5),seq(-20,40,5),cex.axis=1.75)
axis(4,at=seq(0,0.02,0.005),seq(0,0.02,0.005),cex.axis=1.75)
mtext("Humidity Ratio (kgw/kga)", side=4, line=2.5,cex=1.75)
lines(seq(-15,40,length.out=200),rh.lines)
for (i in seq(10,90,10)) {
  lines(seq(-15,40,length.out=200),hum_ratio_rel_hum(rep(i,200),seq(-15,40,length.out=200),hum.ratio.max=0.04))
}

for (j in seq(-5,30,5)) {
  lines(c(j,j),c(0.0,hum_ratio_rel_hum(100,j,hum.ratio.max=0.04)),lwd=1) 
}

##points(tas.2080s,hr.2080s,col='red',pch=16)
##points(tas.2050s,hr.2050s,col='orange',pch=16)
##points(tas.2020s,hr.2020s,col='gold',pch=16)
points(tas.past,hr.past,col='black',pch=16)

text(21.5,0.0195,'Rel. Hum.  100%',cex=1.3)
text(25.6,0.0195,'90%',cex=1.3)
text(27.7,0.0195,'80%',cex=1.3)
text(29.9,0.0195,'70%',cex=1.3)
text(32.7,0.0195,'60%',cex=1.3)
text(34.2,0.018,'50%',cex=1.3)
text(34.2,0.0145,'40%',cex=1.3)
text(34.2,0.011,'30%',cex=1.3)
text(34.2,0.0075,'20%',cex=1.3)
text(34.2,0.0038,'10%',cex=1.3)

##lines(c(20,25),rep(0.012,2),col='green',lwd=4)
##lines(c(20,25),rep(0.006,2),col='green',lwd=4)
##lines(c(20,20),c(0.006,0.012),col='green',lwd=4)
##lines(c(25,25),c(0.006,0.012),col='green',lwd=4) 

box(which='plot')


}
par(mar=c(5.1,2.1,1.1,4.1))
plot(c(),xlim=c(-5,35),ylim=c(0,0.02),main='',xaxs='i',yaxs='i',
         xlab='Dry-bulb Temperature (\u00B0C)',ylab='',axes=F,cex.lab=1.75)
axis(1,at=seq(-20,40,5),seq(-20,40,5),cex.axis=1.75)
axis(4,at=seq(0,0.02,0.005),seq(0,0.02,0.005),cex.axis=1.75)
mtext("Humidity Ratio (kgw/kga)", side=4, line=2.5,cex=1.75)
lines(seq(-15,40,length.out=200),rh.lines)
for (i in seq(10,90,10)) {
  lines(seq(-15,40,length.out=200),hum_ratio_rel_hum(rep(i,200),seq(-15,40,length.out=200),hum.ratio.max=0.04))
}
for (j in seq(-5,30,5)) {
  lines(c(j,j),c(0.0,hum_ratio_rel_hum(100,j,hum.ratio.max=0.04)),lwd=1) 
}

points(tas.2080s,hr.2080s,col='red',pch=16)
##points(tas.2050s,hr.2050s,col='orange',pch=16)
##points(tas.2020s,hr.2020s,col='gold',pch=16)
points(tas.past,hr.past,col=alpha('black',0.2),pch=16,cex=1.5)

text(21.5,0.0195,'Rel. Hum.  100%',cex=1.3)
text(25.6,0.0195,'90%',cex=1.3)
text(27.7,0.0195,'80%',cex=1.3)
text(29.9,0.0195,'70%',cex=1.3)
text(32.7,0.0195,'60%',cex=1.3)
text(34.2,0.018,'50%',cex=1.3)
text(34.2,0.0145,'40%',cex=1.3)
text(34.2,0.011,'30%',cex=1.3)
text(34.2,0.0075,'20%',cex=1.3)
text(34.2,0.0038,'10%',cex=1.3)

##lines(c(20,25),rep(0.012,2),col='green',lwd=4)
##lines(c(20,25),rep(0.006,2),col='green',lwd=4)
##lines(c(20,20),c(0.006,0.012),col='green',lwd=4)
##lines(c(25,25),c(0.006,0.012),col='green',lwd=4)

##legend('topleft',leg=c('Past','2080s','Comfort Zone'),col=c('black','red','green'),cex=1.75,pch=16)
legend('topleft',leg=c('Past','2080s'),col=c('black','red'),cex=1.75,pch=16)

box(which='plot')


dev.off()

