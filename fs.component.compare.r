##Script to plot RCM Anomalies by month

source('/storage/home/ssobie/code/repos/crcm5/plot.configured.crcm5.r',chdir='T')
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir='T')

##----------------------------------------------------------------------------
get.data <- function(var.name,gcm) {

  fs.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/fs_stats/'

  past.file <- paste0(fs.dir,'/FS_STAT_',var.name,'_CWEC_TMY_',gcm,'+CRCM5_historical_20210101-20501231.nc')
  past.nc <- nc_open(past.file)
  past.data <- ncvar_get(past.nc,'fs_stat')

  proj.file <- paste0(fs.dir,'/MORPHED_FS_STAT_',var.name,'_CWEC_TMY_',gcm,'+CRCM5_historical_20210101-20501231.nc')
  proj.nc <- nc_open(proj.file)
  proj.data <- ncvar_get(proj.nc,'fs_stat')

  nc_close(past.nc)
  nc_close(proj.nc)

  rv <- list(past=past.data,proj=proj.data)
  return(rv)
}


make.plots <- function(gcm) {

  fs.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/fs_stats/'
  past.file <- paste0(fs.dir,'/FS_STAT_TASMAX_CWEC_TMY_',gcm,'+CRCM5_historical_20210101-20501231.nc')
  past.nc <- nc_open(past.file)

  tasmax.data <- get.data('TASMAX',gcm)
  tasmax.anoms <- tasmax.data$proj - tasmax.data$past

  tasmin.data <- get.data('TASMIN',gcm)
  tasmin.anoms <- tasmin.data$proj - tasmin.data$past

  tas.data <- get.data('TAS',gcm)
  tas.anoms <- tas.data$proj - tas.data$past

  dwpt.max.data <- get.data('DWPT_MAX',gcm)
  dwpt.max.anoms <- dwpt.max.data$proj - dwpt.max.data$past

  dwpt.min.data <- get.data('DWPT_MIN',gcm)
  dwpt.min.anoms <- dwpt.min.data$proj - dwpt.min.data$past

  dwpt.avg.data <- get.data('DWPT_AVG',gcm)
  dwpt.avg.anoms <- dwpt.avg.data$proj - dwpt.avg.data$past

  wspd.max.data <- get.data('WSPD_MAX',gcm)
  wspd.max.anoms <- wspd.max.data$proj - wspd.max.data$past

  wspd.avg.data <- get.data('WSPD_AVG',gcm)
  wspd.avg.anoms <- wspd.avg.data$proj - wspd.avg.data$past

  insol.data <- get.data('INSOL',gcm)
  insol.anoms <- insol.data$proj - insol.data$past

  for (i in 1:12) {  
    plot.file <- paste0('/storage/data/projects/rci/data/RCM/CRCM5/plots/bc/MORPHED_FS_components.',month.abb[i],'.anomaly.png')
    png(plot.file,width=1000,height=1000)
    par(mfrow=c(3,3))
    par(mar=c(3,3,3,5))
    par(oma=c(1,1,0,0))
    par(mgp=c(2,1,0))

    shared.range <- c(-3,3) ##range(c(range(clim.anoms[,,i]),range(tmy.anoms[,,i]))) ##,range(morph.anoms[,,i])))
    print(shared.range)

    make.crcm5.anom.plot(tasmax.anoms[,,i],past.nc,var.name='tasmax',type='anomaly',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' FS STAT Anomaly TASMAX'))                         
    make.crcm5.anom.plot(tasmin.anoms[,,i],past.nc,var.name='tasmax',type='anomaly',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' FS STAT Anomaly TASMIN'))                         
    make.crcm5.anom.plot(tas.anoms[,,i],past.nc,var.name='tasmax',type='anomaly',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' FS STAT Anomaly TAS'))                         
    make.crcm5.anom.plot(dwpt.max.anoms[,,i],past.nc,var.name='tasmax',type='anomaly',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' FS STAT Anomaly DWPT MAX'))                         
    make.crcm5.anom.plot(dwpt.min.anoms[,,i],past.nc,var.name='tasmax',type='anomaly',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' FS STAT Anomaly DWPT MIN'))                         
    make.crcm5.anom.plot(dwpt.avg.anoms[,,i],past.nc,var.name='tasmax',type='anomaly',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' FS STAT Anomaly DWPT AVG'))                         
    make.crcm5.anom.plot(wspd.max.anoms[,,i],past.nc,var.name='tasmax',type='anomaly',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' FS STAT Anomaly WSPD MAX'))                         
    make.crcm5.anom.plot(wspd.avg.anoms[,,i],past.nc,var.name='tasmax',type='anomaly',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' FS STAT Anomaly WSPD AVG'))                         
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

make.plots('MPI')