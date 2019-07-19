##Script to make cleaner version of energy model output for AGU poster


hr.cool.step2 <- c(8,10,13,15)
hr.heat.step2 <- c(42,37,31,24)
hr.step2 <- rbind(hr.heat.step2,hr.cool.step2)
colnames(hr.step2) <- c('Past','2020s','2050s','2080s')

hr.cool.step4 <- c(13,14,15,16)
hr.heat.step4 <- c(8,7,6,5)
hr.step4 <- rbind(hr.heat.step4,hr.cool.step4)
colnames(hr.step4) <- c('Past','2020s','2050s','2080s')

png('/storage/data/projects/rci/weather_files/high.low.rise.energy.png',width=1000,height=600)


par(mfrow=c(2,2))
par(mar=c(3.2,5.1,4.1,2.1))
plot(c(),xlim=c(0,5),ylim=c(0,60),main='High-Rise Energy Use (Step 2)',xlab='',ylab='Energy Use Intensity (kWh m-2 y-1)',axes=F,yaxs='i',cex.main=2,cex.axis=1.75,cex.lab=1.75)
abline(h=seq(0,70,10),col='gray',lty=2)
barplot(hr.step2,col=c('orange','blue'),add=TRUE,cex.axis=1.75,cex.names=1.75)
legend('topright',leg=c('Cooling','Heating'),col=c('blue','orange'),pch=15,bg='white',cex=1.75)
box(which='plot')


plot(c(),xlim=c(0,5),ylim=c(0,60),main='High-Rise Energy Use (Step 4)',xlab='',ylab='Energy Use Intensity (kWh m-2 y-1)',axes=F,yaxs='i',cex.main=2,cex.axis=1.75,cex.lab=1.75)
abline(h=seq(0,70,10),col='gray',lty=2)
barplot(hr.step4,col=c('orange','blue'),add=TRUE,cex.axis=1.75,cex.names=1.75)
legend('topright',leg=c('Cooling','Heating'),col=c('blue','orange'),pch=15,bg='white',cex=1.75)
box(which='plot')

lr.heat.step2 <- c(44,38,33,28)
lr.heat.step3 <- c(28,24,21,18)
lr.heat.step4 <- c(15,14,12,7)
lr.heat <- rbind(lr.heat.step2,lr.heat.step3,lr.heat.step4)
colnames(lr.heat) <- c('Past','2020s','2050s','2080s')
rownames(lr.heat) <- c('Step 2','Step 3','Step 4') 


plot(c(),xlim=c(0,5),ylim=c(0,60),main='Low-Rise Heating',xlab='',ylab='Energy Use Intensity (kWh m-2 y-1)',axes=F,yaxs='i',cex.main=2,cex.axis=1.75,cex.lab=1.75)
abline(h=seq(0,70,10),col='gray',lty=2)
barplot(lr.heat,col=c('red','orange','gold'),add=TRUE,beside=T,width=0.3,cex.names=1.75,cex.axis=1.75)
legend('topright',leg=c('Step 2','Step 3','Step 4'),col=c('red','orange','gold'),pch=15,bg='white',cex=1.75)
box(which='plot')

lr.cool.step2 <- c(400,500,1000,1350,1810)
lr.cool.step3 <- c(400,510,1020,1440,1900)
lr.cool.step4 <- c(600,750,1490,1750,2300)
lr.cool <- rbind(lr.cool.step2,lr.cool.step3,lr.cool.step4)
colnames(lr.cool) <- c('Past','Now','2020s','2050s','2080s')
rownames(lr.cool) <- c('Step 2','Step 3','Step 4') 


plot(c(),xlim=c(0,5),ylim=c(0,2500),main='Low-Rise Unmet Cooling',xlab='',ylab='Unmet Cooling Hours',axes=F,yaxs='i',cex.main=2,cex.axis=1.75,cex.lab=1.75)
abline(h=seq(0,2500,500),col='gray',lty=2)
barplot(lr.cool,col=c('lightblue','royalblue','blue'),add=TRUE,beside=T,width=0.245,cex.names=1.75,cex.axis=1.75)
legend('topleft',leg=c('Step 2','Step 3','Step 4'),col=c('lightblue','royalblue','blue'),pch=15,bg='white',cex=1.75)
box(which='plot')

dev.off()