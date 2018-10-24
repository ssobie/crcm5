##Script to plot comparisons of the RCM Climatologies and 
##morphed "weather file" degree hours.

source('/storage/home/ssobie/code/repos/crcm5/plot.configured.crcm5.r',chdir='T')
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir='T')

rcm.class.breaks <- get.class.breaks

##----------------------------------------------------------------------------
##Degree Day Values
dd <- function(tmean,tbase) {
  g <- tmean - tbase
  days <- sum(g[g>0], na.rm=T)*1 ##For three-hourly RCM values
  return(round(days))
}

gdd<-function(data,fac){tapply(data,fac, dd, tbase=5)}   ##Growing degree days
cdd<-function(data,fac){tapply(data,fac, dd, tbase=18)}  ##Cooling degree days
hdd<-function(data,fac){tapply(-data,fac,dd, tbase=-18)} ##Heating degree days
fdd<-function(data,fac){tapply(-data,fac,dd, tbase=0)} ##Freezing degree days

dd.fxns <- list(gdd=gdd,
                cdd=cdd,
                hdd=hdd,
                fdd=fdd)
                
make.plots <- function(gcm,var.name) {

  dd.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/degree_hours/'
  tmy.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/'

  present.tmy.file <- paste0(tmy.dir,'/tas_CWEC_TMY_',gcm,'+CRCM5_historical_19810101-20101231.nc')
  ptf.nc <- nc_open(present.tmy.file)
  ptf.tas <- ncvar_get(ptf.nc,'tas')
  present.time <- netcdf.calendar(ptf.nc)  
  nc_close(ptf.nc)

  future.tmy.file <-  paste0(tmy.dir,'/tas_CWEC_TMY_',gcm,'+CRCM5_historical_20210101-20501231.nc')
  ftf.nc <- nc_open(future.tmy.file)
  ftf.tas <- ncvar_get(ftf.nc,'tas')
  future.time <- netcdf.calendar(ftf.nc)  
  nc_close(ftf.nc)

  morphed.tmy.file <- paste0(tmy.dir,'/morphed_tas_CWEC_TMY_ERA-Interim+',gcm,'+CRCM5_historical+rcp85_20210101-20501231.nc')
  mtf.nc <- nc_open(morphed.tmy.file)
  mtf.tas <- ncvar_get(mtf.nc,'tas')
  morphed.time <- netcdf.calendar(mtf.nc)  
  nc_close(mtf.nc)

  fx <- dd.fxns[[var.name]]    
  present.fac <- as.factor(format(present.time,'%m'))
  present.tmy.dd <- aperm(apply(ptf.tas,c(1,2),fx,present.fac),c(2,3,1))
  future.fac <- as.factor(format(future.time,'%m'))
  future.tmy.dd <- aperm(apply(ftf.tas,c(1,2),fx,future.fac),c(2,3,1))
  morphed.fac <- as.factor(format(morphed.time,'%m'))
  morphed.tmy.dd <- aperm(apply(mtf.tas,c(1,2),fx,morphed.fac),c(2,3,1))

  ##---------------------
  ##Degree Days from RCM
  dd.file <- paste0(dd.dir,var.name,'_BC_WC011_',gcm,'+CRCM5_historical+rcp85_19800101-20501231.nc')
  dd.nc <- nc_open(dd.file)
  dd.data <- ncvar_get(dd.nc,var.name)*1
  dd.time <- netcdf.calendar(dd.nc)  

  pst <- head(grep('1981',dd.time),1)
  pen <- tail(grep('2010',dd.time),1)
  fst <- head(grep('2021',dd.time),1)
  fen <- tail(grep('2050',dd.time),1)

  dd.past.fac <- as.factor(format(dd.time[pst:pen],'%m'))
  dd.past.clim <- aperm(apply(dd.data[,,pst:pen],c(1,2),function(x,y){tapply(x,y,mean,na.rm=T)},dd.past.fac),c(2,3,1))
  dd.proj.fac <- as.factor(format(dd.time[fst:fen],'%m'))
  dd.proj.clim <- aperm(apply(dd.data[,,fst:fen],c(1,2),function(x,y){tapply(x,y,mean,na.rm=T)},dd.proj.fac),c(2,3,1))
  ##---------------------

  for (mn in 1:12) { ##1:12) {

    shared.range <- range(c(dd.past.clim[,,mn],dd.proj.clim[,,mn],
                            present.tmy.dd[,,mn],future.tmy.dd[,,mn],
                            morphed.tmy.dd[,,mn]))

    anom.range <- range(c(dd.proj.clim[,,mn]-dd.past.clim[,,mn],future.tmy.dd[,,mn]-present.tmy.dd[,,mn])) 
    plot.file <- paste0('/storage/data/projects/rci/data/RCM/CRCM5/plots/bc/degree_hours/',var.name,'.',gcm,'.',month.abb[mn],'.rcm.png')
    png(plot.file,width=1500,height=500)
    par(mfrow=c(1,3))
    par(mar=c(3,3,3,5))
    par(oma=c(1,1,0,0))
    par(mgp=c(2,1,0))
    make.crcm5.anom.plot(dd.past.clim[,,mn],dd.nc,var.name=var.name,type='past',shared.range=shared.range,
                         plot.title=paste0(toupper(var.name),' ',gcm,' ',month.abb[mn],' Past RCM'))
    rv <- make.crcm5.anom.plot(dd.proj.clim[,,mn],dd.nc,var.name=var.name,type='future',shared.range=shared.range,
                         plot.title=paste0(toupper(var.name),' ',gcm,' ',month.abb[mn],' Future RCM'))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)                       

    rv <- make.crcm5.anom.plot(dd.proj.clim[,,mn]-dd.past.clim[,,mn],dd.nc,var.name=var.name,type='anomaly',shared.range=anom.range,
                       plot.title=paste0(toupper(var.name),' ',gcm,' ',month.abb[mn],' RCM Anomalies'))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    dev.off()

    ##shared.range <- range(c(present.tmy.dd[,,mn],future.tmy.dd[,,mn]))
    plot.file <- paste0('/storage/data/projects/rci/data/RCM/CRCM5/plots/bc/degree_hours/',var.name,'.',gcm,'.',month.abb[mn],'.tmy.png')
    png(plot.file,width=1500,height=500)
    par(mfrow=c(1,3))
    par(mar=c(3,3,3,5))
    par(oma=c(1,1,0,0))
    par(mgp=c(2,1,0))
    make.crcm5.anom.plot(present.tmy.dd[,,mn],dd.nc,var.name=var.name,type='past',shared.range=shared.range,
                         plot.title=paste0(toupper(var.name),' ',gcm,' ',month.abb[mn],' Past TMY'))
    rv <- make.crcm5.anom.plot(future.tmy.dd[,,mn],dd.nc,var.name=var.name,type='future',shared.range=shared.range,
                         plot.title=paste0(toupper(var.name),' ',gcm,' ',month.abb[mn],' Future TMY'))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)                       
    rv <- make.crcm5.anom.plot(future.tmy.dd[,,mn]-present.tmy.dd[,,mn],dd.nc,var.name=var.name,type='anomaly',shared.range=anom.range,
                         plot.title=paste0(toupper(var.name),' ',gcm,' ',month.abb[mn],' TMY Anomalies'))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    dev.off()

    shared.range <- range(c(present.tmy.dd[,,mn],morphed.tmy.dd[,,mn]))
    plot.file <- paste0('/storage/data/projects/rci/data/RCM/CRCM5/plots/bc/degree_hours/',var.name,'.',gcm,'.',month.abb[mn],'.morphed.png')
    png(plot.file,width=1500,height=500)
    par(mfrow=c(1,3))
    par(mar=c(3,3,3,5))
    par(oma=c(1,1,0,0))
    par(mgp=c(2,1,0))
    make.crcm5.anom.plot(present.tmy.dd[,,mn],dd.nc,var.name=var.name,type='past',shared.range=shared.range,
                         plot.title=paste0(toupper(var.name),' ',gcm,' ',month.abb[mn],' Past TMY'))
    rv <- make.crcm5.anom.plot(morphed.tmy.dd[,,mn],dd.nc,var.name=var.name,type='future',shared.range=shared.range,
                         plot.title=paste0(toupper(var.name),' ',gcm,' ',month.abb[mn],' Morphed Future TMY'))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)                       
    rv <- make.crcm5.anom.plot(morphed.tmy.dd[,,mn]-present.tmy.dd[,,mn],dd.nc,var.name=var.name,type='anomaly',shared.range=anom.range,
                         plot.title=paste0(toupper(var.name),' ',gcm,' ',month.abb[mn],' Morphed Anomalies TMY'))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    dev.off()

    shared.range <- range(c(present.tmy.dd[,,mn]-dd.past.clim[,,mn],
                            future.tmy.dd[,,mn]-dd.proj.clim[,,mn],
                            morphed.tmy.dd[,,mn]-dd.proj.clim[,,mn]))

    plot.file <- paste0('/storage/data/projects/rci/data/RCM/CRCM5/plots/bc/degree_hours/',var.name,'.',gcm,'.',month.abb[mn],'.differences.png')
    png(plot.file,width=1500,height=500)
    par(mfrow=c(1,3))
    par(mar=c(3,3,3,5))
    par(oma=c(1,1,0,0))
    par(mgp=c(2,1,0))
    make.crcm5.anom.plot((present.tmy.dd[,,mn]-dd.past.clim[,,mn]), dd.nc,var.name=var.name,type='anomaly',shared.range=shared.range,
                       plot.title=paste0(toupper(var.name),' ',gcm,' ',month.abb[mn],' Past TMY - Past Clim'))
    rv <- make.crcm5.anom.plot((future.tmy.dd[,,mn]-dd.proj.clim[,,mn]), dd.nc,var.name=var.name,type='anomaly',shared.range=shared.range,
                       plot.title=paste0(toupper(var.name),' ',gcm,' ',month.abb[mn],' Future TMY - Future Clim'))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)                       
    rv <- make.crcm5.anom.plot((morphed.tmy.dd[,,mn]-dd.proj.clim[,,mn]), dd.nc,var.name=var.name,type='anomaly',shared.range=shared.range,
                       plot.title=paste0(toupper(var.name),' ',gcm,' ',month.abb[mn],' Morphed TMY - Future Clim'))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    dev.off()

  }

  nc_close(dd.nc)                       
}

##***************************************************************

make.plots('MPI','cdd')