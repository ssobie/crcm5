##Script to plot RCM Anomalies by month

source('/storage/home/ssobie/code/repos/crcm5/plot.configured.crcm5.r',chdir='T')
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir='T')

##----------------------------------------------------------------------------

make.plots <- function(gcm,var.name) {

  ##Climatology Anomalies
  anom.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/anomalies/'
  clim.anom.file <- paste0(anom.dir,'/',var.name,'_mean_anomaly_BC_WC011_',gcm,'+CRCM5_historical+rcp85_20210101-20501231.nc')
  clim.nc <- nc_open(clim.anom.file)
  clim.anoms <- ncvar_get(clim.nc,var.name)

  ##TMY Anomalies
  anom.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/anomalies/'
  tmy.anom.file <- paste0(anom.dir,'/',var.name,'_anomaly_CWEC_TMY_',gcm,'+CRCM5_historical_20210101-20501231.nc')
  tmy.nc <- nc_open(tmy.anom.file)
  tmy.anoms <- ncvar_get(tmy.nc,var.name)

  ##Morphed Anomalies
  anom.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/anomalies/'
  morph.anom.file <- paste0(anom.dir,'/morphed_',var.name,'_anomaly_CWEC_TMY_ERA-Interim+',gcm,'+CRCM5_historical+rcp85_20210101-20501231.nc')
  morph.nc <- nc_open(morph.anom.file)
  morph.anoms <- ncvar_get(morph.nc,var.name)

  for (i in 1:12) {  
    plot.file <- paste0('/storage/data/projects/rci/data/RCM/CRCM5/plots/bc/',var.name,'.anomalies.',month.abb[i],'.png')
    png(plot.file,width=3000,height=1000)
    par(mfrow=c(1,3))
    par(mar=c(3,3,3,5))
    par(oma=c(1,1,0,0))
    par(mgp=c(2,1,0))

    shared.range <- range(c(range(clim.anoms[,,i]),range(tmy.anoms[,,i]),range(morph.anoms[,,i])))
    print(shared.range)

    make.crcm5.anom.plot(clim.anoms[,,i],clim.nc,var.name=var.name,type='anomaly',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' ',toupper(var.name),' Anomaly TAS'))                         
    make.crcm5.anom.plot(tmy.anoms[,,i],tmy.nc,var.name=var.name,type='anomaly',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' ',toupper(var.name),' Anomaly TMY'))
    rv <- make.crcm5.anom.plot(morph.anoms[,,i],morph.nc,var.name=var.name,type='anomaly',shared.range=shared.range,
                               plot.title=paste0(month.abb[i],' ',toupper(var.name),' Anomaly Morphed TMY'))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)

    dev.off()

  }


  nc_close(clim.nc)
  nc_close(tmy.nc)
  nc_close(morph.nc)
}                         

##***************************************************************

make.plots('CanESM2','wspd')