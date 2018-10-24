##Script to plot RCM Anomalies by month

source('/storage/home/ssobie/code/repos/crcm5/plot.configured.crcm5.r',chdir='T')
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir='T')

##----------------------------------------------------------------------------

make.plots <- function(gcm,var.name,type) {

  fac <- 1
  if (var.name=='wspd') {          
    fac <- 3.6
  }

  ##Past Climatologies
  clim.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/climatologies/'
  clim.file <- paste0(clim.dir,'/',var.name,'_',type,'_climatology_BC_WC011_',gcm,'+CRCM5_historical+rcp85_19810101-20101231.nc')
  clim.nc <- nc_open(clim.file)
  past.clim <- ncvar_get(clim.nc,var.name)*fac
  nc_close(clim.nc)

  ##TMY Past Climatologies
  tmy.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/climatologies/'
  tmy.clim.file <- paste0(tmy.dir,'/',var.name,'_',type,'_climatology_CWEC_TMY_',gcm,'+CRCM5_historical_19810101-20101231.nc')
  tmy.nc <- nc_open(tmy.clim.file)
  past.tmy <- ncvar_get(tmy.nc,var.name)*fac
  nc_close(tmy.nc)

  past.diff <- past.tmy - past.clim

  ##Future Climatologies
  clim.file <- paste0(clim.dir,'/',var.name,'_',type,'_climatology_BC_WC011_',gcm,'+CRCM5_historical+rcp85_20210101-20501231.nc')
  clim.nc <- nc_open(clim.file)
  proj.clim <- ncvar_get(clim.nc,var.name)*fac

  ##TMY Future Climatologies
  tmy.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/climatologies/'
  tmy.clim.file <- paste0(tmy.dir,'/',var.name,'_',type,'_climatology_CWEC_TMY_',gcm,'+CRCM5_historical_20210101-20501231.nc')
  tmy.nc <- nc_open(tmy.clim.file)
  proj.tmy <- ncvar_get(tmy.nc,var.name)*fac

  proj.diff <- proj.tmy - proj.clim

  ##Morphed Anomalies
  morph.file <- paste0(tmy.dir,'/morphed_',var.name,'_',type,'_climatology_CWEC_TMY_',gcm,'+CRCM5_historical+rcp85_20210101-20501231.nc')
  morph.nc <- nc_open(morph.file)
  morph.clim <- ncvar_get(morph.nc,var.name)*fac

  morph.diff <- morph.clim - proj.clim

  for (i in 1:12) {  
    plot.file <- paste0('/storage/data/projects/rci/data/RCM/CRCM5/plots/bc/tmy_clim_diffs/',
                        var.name,'.',type,'.',gcm,'.tmy.diffs.',month.abb[i],'.past.png')
    png(plot.file,width=1500,height=500)
    par(mfrow=c(1,3))
    par(mar=c(3,3,3,5))
    par(oma=c(1,1,0,0))
    par(mgp=c(2,1,0))

    diff.range <- range(c(range(past.diff[,,i]),range(proj.diff[,,i]),range(morph.diff[,,i])))
    
    shared.range <-range(c(range(past.clim[,,i]),range(past.tmy[,,i]),
                           range(proj.clim[,,i]),range(proj.tmy[,,i]),
                           range(morph.clim[,,i])))
    print(shared.range)
    rv <- make.crcm5.anom.plot(past.clim[,,i],clim.nc,var.name=var.name,type='past',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' ',gcm,' ',toupper(var.name),' ',type,' Past Clim'))                         
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)
    
    rv <- make.crcm5.anom.plot(past.tmy[,,i],tmy.nc,var.name=var.name,type='past',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' ',gcm,' ',toupper(var.name),' ',type,' Past TMY'))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)
    rv <- make.crcm5.anom.plot(past.diff[,,i],clim.nc,var.name=var.name,type='anomaly',shared.range=diff.range,
                               plot.title=paste0(month.abb[i],' ',gcm,' ',toupper(var.name),' ',type,' Past TMY - Past Clim '))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)
    dev.off()
    ##-----------------------------

    plot.file <- paste0('/storage/data/projects/rci/data/RCM/CRCM5/plots/bc/tmy_clim_diffs/',
                        var.name,'.',type,'.',gcm,'.tmy.diffs.',month.abb[i],'.future.png')
    png(plot.file,width=1500,height=500)
    par(mfrow=c(1,3))
    par(mar=c(3,3,3,5))
    par(oma=c(1,1,0,0))
    par(mgp=c(2,1,0))
    ##shared.range <- range(c(range(proj.clim[,,i]),range(proj.tmy[,,i])))
    rv <- make.crcm5.anom.plot(proj.clim[,,i],clim.nc,var.name=var.name,type='future',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' ',gcm,' ',toupper(var.name),' ',type,' Future Clim'))                         
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)
    rv <- make.crcm5.anom.plot(proj.tmy[,,i],tmy.nc,var.name=var.name,type='future',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' ',gcm,' ',toupper(var.name),' ',type,' Future TMY'))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)
    rv <- make.crcm5.anom.plot(proj.diff[,,i],morph.nc,var.name=var.name,type='anomaly',shared.range=diff.range,
                               plot.title=paste0(month.abb[i],' ',gcm,' ',toupper(var.name),' ',type,' Future TMY - Future Clim '))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    dev.off()

    ##-----------------------------

    plot.file <- paste0('/storage/data/projects/rci/data/RCM/CRCM5/plots/bc/tmy_clim_diffs/',
                        var.name,'.',type,'.',gcm,'.tmy.diffs.',month.abb[i],'.morphed.png')
    png(plot.file,width=1500,height=500)
    par(mfrow=c(1,3))
    par(mar=c(3,3,3,5))
    par(oma=c(1,1,0,0))
    par(mgp=c(2,1,0))
    ##shared.range <- range(c(range(proj.clim[,,i]),range(morph.clim[,,i])))
    rv <- make.crcm5.anom.plot(proj.clim[,,i],clim.nc,var.name=var.name,type='future',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' ',gcm,' ',toupper(var.name),' ',type,' Future Clim'))                         
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)
    rv <- make.crcm5.anom.plot(morph.clim[,,i],tmy.nc,var.name=var.name,type='future',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' ',gcm,' ',toupper(var.name),' ',type,' Morphed TMY'))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)
    rv <- make.crcm5.anom.plot(morph.diff[,,i],morph.nc,var.name=var.name,type='anomaly',shared.range=diff.range,
                               plot.title=paste0(month.abb[i],' ',gcm,' ',toupper(var.name),' ',type,' Morphed TMY - Future Clim '))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    dev.off()
##browser()
  }


  nc_close(clim.nc)
  nc_close(tmy.nc)
  nc_close(morph.nc)
}                         

##***************************************************************

make.plots('MPI','dewpoint','mean')