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

netcdf.calendar <- function(nc) {

  time.calendar <- ncatt_get(nc,'t','calendar')$value
  time.units <- ncatt_get(nc, 't', 'units')$value
  time.values <- ncvar_get(nc, 't')
  if(grepl('days', time.units)) time.values <- time.values*86400
  if(grepl('hours', time.units)) time.values <- time.values*3600
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.values <- origin.pcict + time.values
  return(time.values)
}

get.coords <- function(nc) {


  lat <- ncvar_get(nc, "lat")
  lon <- ncvar_get(nc, "lon")
  lon <- ((lon + 180) %% 360) - 180
  prj <- "+init=epsg:4326"

  i <- expand.grid(seq(nc$dim$rlon$len), seq(nc$dim$rlat$len))
  names(i) <- c("xi", "yi")
  xc <- as.vector(nc$dim$rlon$vals[i$xi])
  yc <- as.vector(nc$dim$rlat$vals[i$yi])
  ll <- as.vector(apply(i, 1, function(x) {lat[x[1], x[2]]}))
  ln <- as.vector(apply(i, 1, function(x) {lon[x[1], x[2]]}))

  nc.grid <- data.frame(lon=ln, lat=ll, xc=xc, yc=yc, i)
  coordinates(nc.grid) <- c("lon", "lat")
  proj4string(nc.grid) <- CRS(prj)
  spatial.coords <- nc.grid@coords
  return(spatial.coords)
}

make.raster.version <- function(nc,data) {

  spatial.coords <- get.coords(nc)

  ##Create empty raster
  e <- extent(spatial.coords)

  rs <- raster(e,ncol=nc$dim$rlon$len,nrow=nc$dim$rlat$len)
  ##Rasterize data matrix using gridded spatial coordinates
  ts <- rasterize(spatial.coords,rs,data,fun=mean)

  ##Assign WGS84 projection
  ts@crs <- CRS("+init=epsg:4326")

  ##Reproject to North America Lambert Conformal Conic
  nalcc.crs <- CRS("+init=epsg:3005") ##CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'")
  map.raster <- projectRaster(ts,crs=nalcc.crs)
  return(ts) ##map.raster)
}



make.crcm5.anom.plot <- function(data,nc,var.name,type,shared.range=NULL,
                                 plot.title,plot.file) {

  map.raster <- make.raster.version(nc,data)
  map.range <- range(data,na.rm=T)
  if (!is.null(shared.range)) {
     map.range <- shared.range
  }

  nalcc.crs <- CRS("+init=epsg:3005") ##CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'")
  bounds <- extent(map.raster)

  xlim.min <- bounds@xmin
  xlim.max <- bounds@xmax
  ylim.min <- bounds@ymin
  ylim.max <- bounds@ymax

  ##Set plot boundaries
  xlim.adj <- (xlim.max - xlim.min) ##* 0.15
  ylim.adj <- (ylim.max - ylim.min) ##* 0.15

  plot.window.xlim <- c((xlim.min + 0.2*xlim.adj), (xlim.max - 0.05*xlim.adj))
  plot.window.ylim <- c((ylim.min + 0.15*ylim.adj), (ylim.max - 0.235*ylim.adj))

##------------------------------------
##Start Plotting

  color.brewer.yellow.orange.red <- c("#FFFFB2", "#FED976", "#FEB24C", "#FD8D3C", "#F03B20")
  color.brewer.red.orange.yellow <- rev(color.brewer.yellow.orange.red)
  color.brewer.yellow.green.blue <- c("#FFFF66", "#99FF33" , "#00CC66", "#009999","#0080FF","#003366")
  color.brewer.light.blue.to.med.blue <- c("#9ECAE1" , "#4292C6","#0080FF","#003366") ##"#F7FBFF", ###"#DEEBF7", ##First entry

  dark.red.to.yellow <- colorRampPalette(colors=color.brewer.red.orange.yellow, bias=1, space = "Lab", interpolate = "linear")
  light.blue.to.dark.blue <- colorRampPalette(colors=color.brewer.light.blue.to.med.blue, bias=1, space = "Lab", interpolate = "linear")

  ##if (var.name=='wspd')
  ##  class.breaks <- c(floor(map.range[1]),seq(-0.5,0.5,0.05),ceiling(map.range[2]))

  if (var.name=='tas') 
     class.breaks <-  c(-5,seq(-2,2,0.25),10) 

  class.breaks <- get.class.breaks(var.name,type,map.range) ## ##c(-6,seq(0,8,1)) ##c(-25,seq(-10,10,2),25) ##
  map.class.breaks.labels <- get.class.break.labels(class.breaks)
  colour.ramp <- get.legend.colourbar(var.name=var.name,map.range=map.range,
                                      my.bp=0,class.breaks=class.breaks,
                                      type)
    my.label.units <- leg.label.units(var.name,type)
    
##  png(plot.file,width=1000,height=1000)
##  par(mar=c(3,3,3,5))
##  par(oma=c(1,1,0,0))
##  par(mgp=c(2,1,0))
  plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
     bg='lightgray',axes=FALSE,
       xlab='Longitude (\u00B0E)',ylab='Latitude (\u00B0N)',main=plot.title,
       cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
  ##plot(map.raster,add=TRUE)
  image(map.raster, col=colour.ramp,breaks=class.breaks,xlim=plot.window.xlim, ylim=plot.window.ylim, add=TRUE)

  ocean.mask.shp <- readOGR('/storage/data/projects/rci/data/assessments/crd/shapefiles/',
                            'west_coast_ocean', stringsAsFactors=F, verbose=F)
  plot(spTransform(ocean.mask.shp,nalcc.crs),col='gray94',add=TRUE)

  lons <- c(-170.0,-150.0,-140.0,-130.0,-120.0,-110.0,-100.0,-90)
  lats <- c(35,40,45,50,55,60,65,70,75,80)
  grat <- graticule(lons, lats, proj = nalcc.crs)
  plot(grat,add=TRUE,lty=3,col='gray',lwd=2)

  shape.dir <- '/storage/data/gis/basedata/base_layers/'

  can.shp <- readOGR(shape.dir, 'north_america_state_provincial_boundaries', stringsAsFactors=F, verbose=F)
##                  extent=c(plot.window.xlim,plot.window.ylim))
  plot(spTransform(can.shp,nalcc.crs),add=TRUE)


xtks <- get.proj.xaxis(lons,nalcc.crs,plot.window.ylim)
ytks <- get.proj.yaxis(lats,nalcc.crs,plot.window.xlim)

axis(2,at=ytks,label=lats,cex.axis=1.5)
axis(1,at=xtks,label=lons,cex.axis=1.5)

##  par(xpd=NA)
##  legend('bottomleft', col = "black", legend=map.class.breaks.labels, pch=22, pt.bg = colour.ramp,
##         pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=my.label.units, xjust=0, cex=1.7)

box(which='plot',lwd=3)
rv <- list(labels=map.class.breaks.labels,colours=colour.ramp,units=my.label.units)
return(rv)
##dev.off()

}



##**************************************************************************************************
##**************************************************************************************************

##gcm <- 'CanESM2'

##rcm.file <- paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/degree_hours/test.cdd.nc')
##nc <- nc_open(rcm.file)
##data <- ncvar_get(nc,'cdd')

##plot.file <- paste0('/storage/data/projects/rci/data/RCM/CRCM5/plots/cdd.test.png')
##make.crcm5.anom.plot(data,nc,var.name='cdd',type='past',
##                     plot.title='CDD TEST',
##                     plot.file=plot.file)
