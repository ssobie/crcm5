##Script to copy step code figure
library(scales)

plot.dir <- '/storage/data/projects/rci/weather_files/'
plot.file <- paste0(plot.dir,'agu.step.figure.png')
png(plot.file,width=1000,height=500)

par(mar=c(0,0,0,0))
par(bg='gray94')
plot(c(),xlim=c(0,120),ylim=c(0,110),main='',xlab='',ylab='',axes=F,xaxs='i',yaxs='i')

text(6.5,105,'2017',font=2,col='dimgray',cex=2)
text(96.5,105,'2032',font=2,col='darkgreen',cex=2)
lines(c(12,87),c(105,105),col='darkgray',lwd=3,lty=2)
arrows(87,105,90,105,col='darkgray',lwd=3)
polygon(c(0,35,45,0),c(1,1,19,19),col=alpha('gray',0.75),border=alpha('gray',0.75))
polygon(c(0,50,60,0),c(21,21,39,39),col=alpha('orange',0.75),border=alpha('orange',0.75))
polygon(c(0,65,75,0),c(41,41,59,59),col=alpha('gold',0.75),border=alpha('gold',0.75))
polygon(c(0,80,90,0),c(61,61,79,79),col=alpha('olivedrab',0.75),border=alpha('olivedrab',0.75))
polygon(c(0,100,100,0),c(81,81,99,99),col=alpha('green',0.75),border=alpha('green',0.75))
abline(h=10,col='darkgray')
abline(h=20,col='darkgray')
abline(h=40,col='darkgray')
abline(h=60,col='darkgray')
abline(h=80,col='darkgray')
abline(h=100,col='darkgray')

text(15,5,'BC Building Code',font=2,col='white',cex=2)
text(65,14,'Enhanced Compliance',font=3,cex=2,col='dimgray')
text(6.5,15,'STEP',font=2,col='white',cex=2)
text(15,15,'1',font=2,col='white',cex=4)

text(6.5,30,'STEP',font=2,col='white',cex=2)
text(15,30,'2',font=2,col='white',cex=4)

text(6.5,50,'STEP',font=2,col='white',cex=2)
text(15,50,'3',font=2,col='white',cex=4)

text(6.5,70,'STEP',font=2,col='white',cex=2)
text(15,70,'4',font=2,col='white',cex=4)

text(6.5,90,'STEP',font=2,col='white',cex=2)
text(15,90,'5',font=2,col='white',cex=4)

polygon(c(101,119,119,101),c(0,0,20,20),col=alpha('green',0.15),border=alpha('green',0.15))
polygon(c(101,119,119,101),c(20,20,40,40),col=alpha('green',0.30),border=alpha('green',0.30))
polygon(c(101,119,119,101),c(40,40,60,60),col=alpha('green',0.45),border=alpha('green',0.45))
polygon(c(101,119,119,101),c(60,60,80,80),col=alpha('green',0.60),border=alpha('green',0.60))
polygon(c(101,119,119,101),c(80,80,90,90),col=alpha('green',0.75),border=alpha('green',0.75))
polygon(c(101,119,110,101),c(90,90,100,90),col=alpha('green',0.75),border=alpha('green',0.75))

text(110,5,'Average',cex=2,font=2,col='dimgray')
text(110,15,'Improved',cex=2,font=2,col='dimgray')
text(110,30,'10% Better',cex=2,font=2,col='dimgray')
text(110,50,'20% Better',cex=2,font=2,col='dimgray')
text(110,70,'40% Better',cex=2,font=2,col='dimgray')
text(110,90,'Net Zero',cex=2,font=2,col='dimgray')
box(which='plot')



dev.off()