##Script to plot RCM Anomalies by month

source('/storage/home/ssobie/code/repos/crcm5/plot.configured.crcm5.r',chdir='T')
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir='T')

##----------------------------------------------------------------------------

##TAS, Dewpoint, WSPD, RHS, INSOL
rcm.class.breaks <- function(var.name,type,map.range, manual.breaks=""){
  print('Class Breaks')
  print(map.range)
  map.diff <- diff(map.range)
  class.width <- (map.diff/6)

  if (type=='percent') {
    if (class.width <= 5)
      class.rv <- round(class.width/1)*1
    if (class.width < 10 & class.width > 5)
      class.rv <- round(class.width/5)*5
    if (class.width >= 10)
      class.rv <- round(class.width/10)*10
    class.rv <- min(c(class.rv,50))
    class.rv <- max(c(class.rv,0.5))
  } else {
      if (class.width < 1) {
        class.rv <- ceiling(class.width/0.05)*0.05
        class.rv <- max(c(0.05,class.rv))
      }
      if (class.width < 3) {
        class.rv <- ceiling(class.width/0.1)*0.1
        class.rv <- max(c(0.1,class.rv))
      }
      if (class.width <= 5)
        class.rv <- round(class.width/0.25)*0.25
      if (class.width < 10 & class.width > 5)
        class.rv <- round(class.width/2)*2
      if (class.width >= 10)
        class.rv <- round(class.width/5)*5
      if (class.width >= 100)
        class.rv <- round(class.width/100)*100
      class.rv <- min(c(class.rv,500))
      class.rv <- max(c(class.rv,0.05))
  }

  if (manual.breaks[1] != ""){
    class.breaks <- manual.breaks
  } else {
      map.min <- floor(map.range[1]/class.rv) * class.rv
      map.max <- ceiling(map.range[2]/class.rv) * class.rv
      class.breaks <-seq(from=map.min, to=map.max, by=class.rv)
 }

  return(class.breaks)
}


make.plots <- function(gcm,var.name,type,
                       clim.dir,past.file,future.file,anom.dir,anomaly.file) {
  fac <- 1
  if (var.name=='wspd') {                        
    fac <- 3.6 ##Convert windspeed to km/h
  }

  ##Climatology Anomalies  
  past.nc <- nc_open(paste0(clim.dir,past.file))
  past.clim <- ncvar_get(past.nc,var.name)*fac
  proj.nc <- nc_open(paste0(clim.dir,future.file))
  proj.clim <- ncvar_get(proj.nc,var.name)*fac
  anom.nc <- nc_open(paste0(anom.dir,anomaly.file))
  anom.clim <- ncvar_get(anom.nc,var.name)*fac

  for (i in 1:12) { ##1:12) {  
    plot.file <- paste0('/storage/data/projects/rci/data/RCM/CRCM5/plots/bc/rcm/',var.name,'.',type,'.',gcm,'.clims.anom.',month.abb[i],'.png')
    png(plot.file,width=1500,height=500)
    par(mfrow=c(1,3))
    par(mar=c(3,3,3,5))
    par(oma=c(1,1,0,0))
    par(mgp=c(2,1,0))

    bc.mask.shp <- spTransform(readOGR('/storage/data/projects/rci/data/assessments/bc/shapefiles/',
                               'bc', stringsAsFactors=F, verbose=F),CRS("+init=epsg:3005"))
    past.brick <- make.raster.version(past.nc,past.clim[,,i])
    past.mask <- mask(past.brick,bc.mask.shp)
    proj.brick <- make.raster.version(proj.nc,proj.clim[,,i])
    proj.mask <- mask(proj.brick,bc.mask.shp)
    anom.brick <- make.raster.version(anom.nc,anom.clim[,,i])
    anom.mask <- mask(anom.brick,bc.mask.shp)

    shared.range <- range(c(past.mask@data@max,past.mask@data@min,proj.mask@data@max,proj.mask@data@min))
    print(shared.range)
    rv <- make.crcm5.anom.plot(past.clim[,,i],past.nc,var.name=var.name,type='past',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' ',toupper(var.name),' ',toupper(type),' ',gcm,' Past Clim')) 
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)
    rv <- make.crcm5.anom.plot(proj.clim[,,i],proj.nc,var.name=var.name,type='past',shared.range=shared.range,
                         plot.title=paste0(month.abb[i],' ',toupper(var.name),' ',toupper(type),' ',gcm,' Future Clim')) 
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)
    ##c(anom.mask@data@min,anom.mask@data@max)
    rv <- make.crcm5.anom.plot(anom.clim[,,i],anom.nc,var.name=var.name,type='anomaly',##shared.range=c(-2.1,2.8),
                               plot.title=paste0(month.abb[i],' ',toupper(var.name),' ',toupper(type),' ',gcm,' Anomaly'))
    par(xpd=NA)
    legend('bottomleft', col = "black", legend=rv$labels, pch=22, pt.bg = rv$colours,
            pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=rv$units, xjust=0, cex=1.7)
    par(xpd=FALSE)
    dev.off()

  }
  nc_close(past.nc)
  nc_close(proj.nc)
  nc_close(anom.nc)
}                         

##***************************************************************

gcms <- c('CanESM2','MPI')
var.name <- 'dewpoint'
types <- c('mean')

for (gcm in gcms) {
  print(gcm)
  for (type in types) {
    make.plots(gcm,var.name,type,
             clim.dir='/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/climatologies/',
             past.file=paste0(var.name,'_',type,'_climatology_BC_WC011_',gcm,'+CRCM5_historical+rcp85_19810101-20101231.nc'),
             future.file=paste0(var.name,'_',type,'_climatology_BC_WC011_',gcm,'+CRCM5_historical+rcp85_20210101-20501231.nc'),
             anom.dir='/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/anomalies/',
             anomaly.file=paste0(var.name,'_',type,'_anomaly_BC_WC011_',gcm,'+CRCM5_historical+rcp85_20210101-20501231.nc')) 

  }
}