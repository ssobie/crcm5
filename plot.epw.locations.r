##Script to plot CRCM5 RCM data

library(ncdf4)
library(raster)
library(PCICt)
library(graticule)
library(rgdal)
library(RColorBrewer) # to get the sweet color ramps
##library(scales)

source('/storage/home/ssobie/code/repos/assessments/resource.region.map.support.r',chdir='T')

##--------------------------------------------------
##Coordinates

convert.to.nalcc.coords <- function(lon,lat,nalcc.crs) {

  d <- data.frame(lon=lon,lat=lat)
  coordinates(d) <- c('lon','lat')
  proj4string(d) <- CRS("+init=epsg:4326")
  d.albers <- spTransform(d,nalcc.crs)
  rv <- d.albers@coords
  return(rv)
}

##X-Axis Ticks
get.proj.xaxis <- function(lons,nalcc.crs,plot.window.ylim) {

  y <- seq(0,80,0.1)
  xm <- sapply(lons,rep,length(y))
  S <- apply(xm,2,function(x) {SpatialPoints(cbind(x,y), proj4string = CRS("+proj=longlat +datum=WGS84"))})
  S2<- lapply(S,spTransform, nalcc.crs)
  indices <- lapply(S2,function(x){which.min(abs(x@coords[,'y']-plot.window.ylim[1]))})
  xticks <- mapply(FUN=function(s,i){s@coords[,'x'][i]},S2,indices)
  return(xticks)
}


  ##Y-Axis Ticks
get.proj.yaxis <- function(lats,nalcc.crs,plot.window.xlim) {

  x <- seq(-180,-80,0.1)
  ym <- sapply(lats,rep,length(x))
  S <- apply(ym,2,function(y) {SpatialPoints(cbind(x,y), proj4string = CRS("+proj=longlat +datum=WGS84"))})
  S2<- lapply(S,spTransform, nalcc.crs)
  indices <- lapply(S2,function(x){which.min(abs(x@coords[,'x']-plot.window.xlim[1]))})
  yticks <- mapply(FUN=function(s,i){s@coords[,'y'][i]},S2,indices)
  return(yticks)
}

make.epw.prism.plot <- function(tasmax.file,tasmin.file,type,shared.range=NULL,epw.coords,
                                 plot.title,plot.file) {
  nalcc.crs <- CRS("+init=epsg:3005")  
  raw.tasmax <- subset(brick(tasmax.file),13) ##Annual
  raw.tasmin <- subset(brick(tasmin.file),13) ##Annual

  map.tasmax <- projectRaster(raw.tasmax,crs=nalcc.crs)
  map.tasmin <- projectRaster(raw.tasmin,crs=nalcc.crs)

  map.range <- c(-20,20) ##range(as.matrix(map.tasmax),na.rm=T)
  ##tasmin.range <- range(as.matrix(map.tasmin),na.rm=T)

  nalcc.crs <- CRS("+init=epsg:3005") 
  bounds <- extent(map.tasmax)

  xlim.min <- bounds@xmin
  xlim.max <- bounds@xmax
  ylim.min <- bounds@ymin
  ylim.max <- bounds@ymax

  ##Set plot boundaries
  xlim.adj <- (xlim.max - xlim.min) ##* 0.15
  ylim.adj <- (ylim.max - ylim.min) ##* 0.15

  plot.window.xlim <- c((xlim.min + 0.12*xlim.adj), (xlim.max - 0.01*xlim.adj))
  plot.window.ylim <- c((ylim.min + 0.01*ylim.adj), (ylim.max - 0.1*ylim.adj))

##------------------------------------
##Start Plotting

  color.brewer.yellow.orange.red <- c("#FFFFB2", "#FED976", "#FEB24C", "#FD8D3C", "#F03B20")
  color.brewer.red.orange.yellow <- rev(color.brewer.yellow.orange.red)
  color.brewer.yellow.green.blue <- c("#FFFF66", "#99FF33" , "#00CC66", "#009999","#0080FF","#003366")
  color.brewer.light.blue.to.med.blue <- c("#9ECAE1" , "#4292C6","#0080FF","#003366") ##"#F7FBFF", ###"#DEEBF7", ##First entry

  dark.red.to.yellow <- colorRampPalette(colors=color.brewer.red.orange.yellow, bias=1, space = "Lab", interpolate = "linear")
  light.blue.to.dark.blue <- colorRampPalette(colors=color.brewer.light.blue.to.med.blue, bias=1, space = "Lab", interpolate = "linear")

  class.breaks <- get.class.breaks('tasmax',type,map.range)
  ##tasmin.breaks <- get.class.breaks('tasmin',type,tasmin.range)
  map.class.breaks.labels <- get.class.break.labels(class.breaks)
  ##tasmin.class.breaks.labels <- get.class.break.labels(tasmin.breaks)
  colour.ramp <- get.legend.colourbar(var.name='tasmax',map.range=map.range,
                                      my.bp=0,class.breaks=class.breaks,
                                      type)
  ##tasmin.colour.ramp <- get.legend.colourbar(var.name='tasmin',map.range=tasmin.range,
  ##                                    my.bp=0,class.breaks=tasmin.breaks,
  ##                                    type)

  my.label.units <- leg.label.units('tasmax',type)

  projected.coords <- convert.to.nalcc.coords(epw.coords[,2],
                                              epw.coords[,3],
                                              nalcc.crs)
    
  png(plot.file,width=1400,height=800)
  par(mfrow=c(1,2))
  ##par(mar=c(3,3,3,5))
  ##par(oma=c(1,1,0,0))
  ##par(mgp=c(2,1,0))
  plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
     bg='lightgray',axes=FALSE,
       xlab='Longitude (\u00B0E)',ylab='Latitude (\u00B0N)',main='PRISM Maximum Temperature (1981-2010)',
       cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
  image(map.tasmax, col=colour.ramp,breaks=class.breaks,xlim=plot.window.xlim, ylim=plot.window.ylim, add=TRUE)
  ocean.mask.shp <- readOGR('/storage/data/projects/rci/data/assessments/crd/shapefiles/',
                            'west_coast_ocean', stringsAsFactors=F, verbose=F)
  washington.shp <- readOGR('/storage/data/projects/rci/data/assessments/crd/shapefiles/',
                            'washington', stringsAsFactors=F, verbose=F)
  plot(spTransform(ocean.mask.shp,nalcc.crs),col='gray94',add=TRUE)
  plot(spTransform(washington.shp,nalcc.crs),col='gray94',add=TRUE)
  na.shp <- readOGR('/storage/data/projects/rci/data/assessments/crd/shapefiles/',
                    'north_america_state_provincial_boundaries_pnwa_diss_albers',
                    stringsAsFactors=F, verbose=F)
  plot(na.shp,add=TRUE)
  points(projected.coords[,1],projected.coords[,2],pch=18,cex=2.1,col='white')
  points(projected.coords[,1],projected.coords[,2],pch=18,cex=1.65,col='darkgreen')

  lons <- c(-140.0,-135.0,-130.0,-125.0,-120.0,-115.0,-110.0,-100.0)
  lats <- c(45,47.5,50,52.5,55,57.5,60,65)
  grat <- graticule(lons, lats, proj = nalcc.crs)
  plot(grat,add=TRUE,lty=3,col='gray',lwd=2)

  shape.dir <- '/storage/data/gis/basedata/base_layers/'

  bc.shp <- readOGR(shape.dir, 'bc_province', stringsAsFactors=F, verbose=F)
##                  extent=c(plot.window.xlim,plot.window.ylim))
  plot(bc.shp,add=TRUE)

xtks <- get.proj.xaxis(lons,nalcc.crs,plot.window.ylim)
ytks <- get.proj.yaxis(lats,nalcc.crs,plot.window.xlim)

axis(2,at=ytks,label=lats,cex.axis=1.5)
axis(1,at=xtks,label=lons,cex.axis=1.5)

##  par(xpd=NA)
  legend('topright', col = "black", legend=rev(map.class.breaks.labels), pch=22, pt.bg = rev(colour.ramp),
         pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=my.label.units, xjust=0, cex=1.7)

box(which='plot',lwd=3)

##---------------------------------------------------------
##---------------------------------------------------------


  plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
     bg='lightgray',axes=FALSE,
       xlab='Longitude (\u00B0E)',ylab='Latitude (\u00B0N)',main='PRISM Minimum Temperature (1981-2010)',
       cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
  image(map.tasmin, col=colour.ramp,breaks=class.breaks,xlim=plot.window.xlim, ylim=plot.window.ylim, add=TRUE)
  ocean.mask.shp <- readOGR('/storage/data/projects/rci/data/assessments/crd/shapefiles/',
                            'west_coast_ocean', stringsAsFactors=F, verbose=F)
  plot(spTransform(ocean.mask.shp,nalcc.crs),col='gray94',add=TRUE)
  plot(spTransform(washington.shp,nalcc.crs),col='gray94',add=TRUE)
  na.shp <- readOGR('/storage/data/projects/rci/data/assessments/crd/shapefiles/',
                    'north_america_state_provincial_boundaries_pnwa_diss_albers',
                    stringsAsFactors=F, verbose=F)
  plot(spTransform(na.shp,nalcc.crs),add=TRUE)
  points(projected.coords[,1],projected.coords[,2],pch=18,cex=2.1,col='white')
  points(projected.coords[,1],projected.coords[,2],pch=18,cex=1.75,col='darkgreen')
  lons <- c(-170.0,-150.0,-140.0,-130.0,-120.0,-110.0,-100.0,-90)
  lats <- c(35,40,45,50,55,60,65,70,75,80)
  grat <- graticule(lons, lats, proj = nalcc.crs)
  plot(grat,add=TRUE,lty=3,col='gray',lwd=2)

  shape.dir <- '/storage/data/gis/basedata/base_layers/'

  bc.shp <- readOGR(shape.dir, 'bc_province', stringsAsFactors=F, verbose=F)
##                  extent=c(plot.window.xlim,plot.window.ylim))
  plot(bc.shp,add=TRUE)

xtks <- get.proj.xaxis(lons,nalcc.crs,plot.window.ylim)
ytks <- get.proj.yaxis(lats,nalcc.crs,plot.window.xlim)

axis(2,at=ytks,label=lats,cex.axis=1.5)
axis(1,at=xtks,label=lons,cex.axis=1.5)

##  par(xpd=NA)
  legend('topright', col = "black", legend=rev(map.class.breaks.labels), pch=22, pt.bg = rev(colour.ramp),
         pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=my.label.units, xjust=0, cex=1.7)

box(which='plot',lwd=3)


##rv <- list(labels=map.class.breaks.labels,colours=colour.ramp,units=my.label.units)
##return(rv)
dev.off()

}



##**************************************************************************************************
##**************************************************************************************************

epw.file <- '/storage/data/projects/rci/weather_files/epw.sites.csv'
epw.coords <- read.csv(epw.file,header=FALSE,as.is=TRUE)

tasmax.file <- '/storage/data/climate/PRISM/dataportal/tmax_monClim_PRISM_historical_run1_198101-201012.nc'
tasmin.file <- '/storage/data/climate/PRISM/dataportal/tmin_monClim_PRISM_historical_run1_198101-201012.nc'

plot.file <- paste0('/storage/data/projects/rci/weather_files/epw.locations.png')
make.epw.prism.plot(tasmax.file,tasmin.file,type='past',epw.coords=epw.coords,
                    plot.title='TASMAX',
                    plot.file=plot.file)
