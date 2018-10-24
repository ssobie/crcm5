##Script to plot RCM Anomalies by month

source('/storage/home/ssobie/code/repos/crcm5/plot.configured.crcm5.r',chdir='T')
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir='T')

##----------------------------------------------------------------------------
get.data <- function(var.name,gcm) {

  fs.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/fs_stats/'
  past.file <- paste0(fs.dir,'/HOURLY_FS_STAT_',var.name,'_CWEC_TMY_',gcm,'+CRCM5_historical_20210101-20501231.nc')
  past.nc <- nc_open(past.file)
  past.data <- ncvar_get(past.nc,'fs_stat')

  proj.file <- paste0(fs.dir,'/MORPHED_HOURLY_FS_STAT_',var.name,'_CWEC_TMY_',gcm,'+CRCM5_historical_20210101-20501231.nc')
  proj.nc <- nc_open(proj.file)
  proj.data <- ncvar_get(proj.nc,'fs_stat')

  nc_close(past.nc)
  nc_close(proj.nc)
  rv <- list(past=past.data,proj=proj.data)
  return(rv)
}


make.plots <- function(gcm) {

  fs.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/fs_stats/'
  past.file <- paste0(fs.dir,'/HOURLY_FS_STAT_DWPT_CWEC_TMY_',gcm,'+CRCM5_historical_19810101-20101231.nc')
  past.nc <- nc_open(past.file)

  tas.data <- get.data('TAS',gcm)
  tas.anoms <- tas.data$proj - tas.data$past

  dwpt.data <- get.data('DWPT',gcm)
  dwpt.anoms <- dwpt.data$proj - dwpt.data$past

  wspd.data <- get.data('WSPD',gcm)
  wspd.anoms <- wspd.data$proj - wspd.data$past

  insol.data <- get.data('INSOL',gcm)
  insol.anoms <- insol.data$proj - insol.data$past

  for (i in 1:12) {  
    plot.file <- paste0('/storage/data/projects/rci/data/RCM/CRCM5/plots/bc/MORPHED_HOURLY_FS_components.',gcm,'.',month.abb[i],'.anomaly.png')
    png(plot.file,width=1000,height=900)
    par(mfrow=c(2,2))
    par(mar=c(3,3,3,5))
    par(oma=c(1,1,0,0))
    par(mgp=c(2,1,0))
     
    ##shared.max <- max(c(max(tas.anoms[,,i]),max(dwpt.anoms[,,i]),max(wspd.anoms[,,i]),max(insol.anoms[,,i])))
    ##shared.min <- max(c(min(tas.anoms[,,i]),min(dwpt.anoms[,,i]),min(wspd.anoms[,,i]),min(insol.anoms[,,i])))
    shared.range <- range(c(range(tas.anoms[,,i]),range(dwpt.anoms[,,i]),range(wspd.anoms[,,i]),range(insol.anoms[,,i])))
    print(shared.range)

    make.crcm5.anom.plot(tas.anoms[,,i],past.nc,var.name='tas',type='anomaly',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' FS STAT Anomaly TAS'))                         
    make.crcm5.anom.plot(dwpt.anoms[,,i],past.nc,var.name='tasmax',type='anomaly',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' FS STAT Anomaly DWPT'))                         
    make.crcm5.anom.plot(wspd.anoms[,,i],past.nc,var.name='tasmax',type='anomaly',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' FS STAT Anomaly WSPD'))                         
    rv <- make.crcm5.anom.plot(insol.anoms[,,i],past.nc,var.name='tasmax',type='anomaly',shared.range=shared.range,
                               plot.title=paste0(month.abb[i],' FS STAT Anomaly INSOL'))                         
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    dev.off()

  }
  nc_close(past.nc)
}                         

##***************************************************************

make.plots('CanESM2')