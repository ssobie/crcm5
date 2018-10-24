##Script to plot quantiles of TMY values in 30 distributions

source('/storage/home/ssobie/code/repos/crcm5/plot.configured.crcm5.r',chdir='T')
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir='T')

##----------------------------------------------------------------------------

make.plots <- function(var.name,type) {

  ##CanESM2 Quantiles
  qt.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/fs_stats/'
  can.file <- paste0(qt.dir,'/MORPHED_TASMAX_QUANTILE_CWEC_TMY_CanESM2+CRCM5_historical_20210101-20501231.nc')
  can.nc <- nc_open(can.file)
  can.qt <- ncvar_get(can.nc,'tasmax')


  mpi.file <- paste0(qt.dir,'/MORPHED_TASMAX_QUANTILE_CWEC_TMY_MPI+CRCM5_historical_20210101-20501231.nc')
  mpi.nc <- nc_open(mpi.file)
  mpi.qt <- ncvar_get(mpi.nc,'tasmax')
  nc_close(mpi.nc)

  for (i in 1:12) {  
    plot.file <- paste0('/storage/data/projects/rci/data/RCM/CRCM5/plots/bc/quantiles/',
                        var.name,'.',type,'.quantiles.',month.abb[i],'.png')
    png(plot.file,width=1000,height=500)
    par(mfrow=c(1,2))
    par(mar=c(3,3,3,5))
    par(oma=c(1,1,0,0))
    par(mgp=c(2,1,0))

    shared.range <- c(0.75,1) 
    print(shared.range)
    rv <- make.crcm5.anom.plot(can.qt[,,i],can.nc,var.name=var.name,type='past',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' CanESM2 ',toupper(var.name),' ',type,' Quantile'))                         
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)
    rv <- make.crcm5.anom.plot(mpi.qt[,,i],can.nc,var.name=var.name,type='past',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' MPI ',toupper(var.name),' ',type,' Quantile'))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)
    dev.off()

    ##-----------------------------

  }

  nc_close(can.nc)
}                         

##***************************************************************

make.plots('tas','max')