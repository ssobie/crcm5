##Script to map the downscaled output

library(graticule)


##Plotting script for Capital Regional District

##------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------

convert.to.ne.coords <- function(lon,lat,ne.crs) {
 
  d <- data.frame(lon=lon,lat=lat)
  coordinates(d) <- c('lon','lat')
  proj4string(d) <- CRS("+init=epsg:4326")
  d.albers <- spTransform(d,CRS(ne.crs))
  rv <- d.albers@coords
  return(rv)
}

reg.ds.maps <- function(box.data,region,region.range,box.range,
                        var.name,type,ds.type,region.shp,shp.buffer,
                        plot.file,plot.title,
                        coords=NULL,set.breaks=NULL,proj,
                        overlays=NULL,
                        leg.loc='topright',width=800,height=800,
                        shared.range=NULL,shared.box=NULL,draft=TRUE) { ##width=1000,height=800) {

  ne.crs <- "+init=epsg:3005"

  box.data <-projectRaster(box.data,crs=CRS("+init=epsg:3005"))
  bounds <- extent(box.data)
  xlim.min <- bounds@xmin
  xlim.max <- bounds@xmax
  ylim.min <- bounds@ymin
  ylim.max <- bounds@ymax

  white.box <- box.data - box.data
  white.box[is.na(white.box)] <- 0
  
  ##Metadata


  ##Set plot boundaries
  xlim.adj <- (xlim.max - xlim.min) * 0.02 ##0.025
  ylim.adj <- (ylim.max - ylim.min) * 0.02 ##0.025 

  plot.window.xlim <- c((xlim.min + xlim.adj*6), (xlim.max + xlim.adj*0.25))
  plot.window.ylim <- c((ylim.min - ylim.adj), (ylim.max + ylim.adj))
  plot.window.ylim <- c((ylim.min + ylim.adj*2), (ylim.max + ylim.adj*0.25))


  e <- extent(c(plot.window.xlim,plot.window.ylim))

  map.range <- region.range
  print(map.range)
  
  if (!is.null(shared.range))
    map.range <- shared.range

  if (!is.null(shared.box))
    box.range <- shared.box
  class.breaks <- get.class.breaks(var.name,type,map.range,manual.breaks='')

  if (!is.null(set.breaks)) {
    old.breaks <- get.class.breaks(var.name,type,map.range,manual.breaks='')
    class.breaks <- set.breaks
    colour.subset <- class.breaks %in% old.breaks
    colour.ramp <- get.legend.colourbar(var.name=var.name,map.range=map.range,
                                        my.bp=0,class.breaks=class.breaks,
                                        type)    
  } else {
    class.breaks <- get.class.breaks(var.name,type,map.range,manual.breaks='')
    colour.ramp <- get.legend.colourbar(var.name=var.name,map.range=map.range,
                                        my.bp=0,class.breaks=class.breaks,
                                        type)    
  }

  map.class.breaks.labels <- get.class.break.labels(class.breaks)
  print(map.class.breaks.labels)

  ##------------------------------------------------------------------------------------------------
  ##Fix the colour.ramp to include greater or less than breaks for the bounding
  ##box if needed
  ##Both
  ##Note: for ranges that are very narrow (i.e. 1 unit or less) or box.ranges and map.ranges that are very close, this
  ##doesn't work as well. The legend repeats intervals as this gets rounded similar values.
  if ((box.range[1] < map.range[1]) & (box.range[2] > map.range[2]) & (class.breaks[1] != 0)) {
    dx <- diff(class.breaks)[1]
    class.breaks <- c(floor(box.range[1]/dx)*dx,class.breaks,ceiling(box.range[2]/dx)*dx)
    colour.ramp <- get.legend.colourbar(var.name=var.name,map.range=box.range,
                                        my.bp=0,class.breaks=class.breaks,
                                        type)        
    map.class.breaks.labels <- get.class.break.labels(class.breaks,lesser.sign=TRUE,greater.sign=TRUE)
  } else {  
    ##Greater than
    if (box.range[2] > map.range[2]) {
      dx <- diff(class.breaks)[1]
      class.breaks <- c(class.breaks,ceiling(box.range[2]/dx)*dx)
      colour.ramp <- get.legend.colourbar(var.name=var.name,map.range=box.range,
                                          my.bp=0,class.breaks=class.breaks,
                                          type)        
      map.class.breaks.labels <- get.class.break.labels(class.breaks,greater.sign=TRUE)
    }
    ##Less than
    if ((box.range[1] < map.range[1]) & (class.breaks[1] !=0) ) {
      dx <- diff(class.breaks)[1]
      class.breaks <- c(floor(box.range[1]/dx)*dx,class.breaks)
      colour.ramp <- get.legend.colourbar(var.name=var.name,map.range=box.range,
                                          my.bp=0,class.breaks=class.breaks,
                                          type)        
      map.class.breaks.labels <- get.class.break.labels(class.breaks,lesser.sign=TRUE)

    }
  }

##  class.breaks <- c(800,1000,1200,1400,1600,1800,2000,3000,3400,3800,4000,4200,4400,6000) ##Metro Van prism precip
##  class.breaks <- c(0,seq(6,12,by=1.5),12.5,13,13.5,14,14.5,15) ##Metro Van prism tasmax
##  class.breaks <- c(-6,seq(-3,4.5,by=1.5),5,5.5,6,6.5,7,10)
##  class.breaks <- c(0,seq(110,180,by=5),1000)
##  class.breaks <- c(0,seq(100,200,by=10),300,400,500,1000)
  if (var.name=='snowdepth' & !grepl('(anomaly|percent)',type)) {
    class.breaks <- c(0,0.5,1,2,3,4,6,8,1000) ##10,1000)
    map.class.breaks.labels <- get.class.break.labels(class.breaks,greater.sign=TRUE,lesser.sign=FALSE)
  } else if (var.name=='snowdepth' & grepl('(anomaly)',type)) {
    class.breaks <- c(-100,-5,-4,-3,-2,-1,-0.5,-0.4,-0.3,-0.2,-0.1,0)
    map.class.breaks.labels <- get.class.break.labels(class.breaks,greater.sign=FALSE,lesser.sign=TRUE)
  } else {
##    map.class.breaks.labels <- get.class.break.labels(class.breaks,greater.sign=TRUE,lesser.sign=TRUE)
  } 

  ##colour.ramp <- rev(make.fixed.colour.ramps('temp'))
  ##class.breaks <- c(14000,6600,3200,1600,1000,700,500,0) ##For Annual Precipitation
  ##class.breaks <- c(200,90,80,70,60,50,40,0) ##For Precip Return periods (tentative)
  ##class.breaks <- c(100,40,36,32,28,24,22,0) ##Tasmax Return Periods
  ##class.breaks <- c(0,-36,-40,-44,-48,-52,-56,-100) ##Tasmin Return Periods
  ##class.breaks <- c(500,240,200,160,120,80,40,0) ##R95 Precipitation
  ##class.breaks <- c(40,18,16,14,12,10,8,0) ##RX1Day Precipitation
  ##class.breaks <- c(360,300,250,200,150,100,50,0) ##Frost Free Days
  ##class.breaks <- c(50,45,40,35,30,20,10,0) ##TX90p
  ##class.breaks <- c(50,45,40,35,30,20,10,0) ##TN10p Precipitation
  
  ##class.breaks <- c(0,25,50,75,100,200,400,600,800,1000) ##Winter Precipitation
  ##class.breaks <- c(3.5,4,4.5,5,5.5,6,6.5,7,7.5,8)
  ##print('Class Labels')
  map.class.breaks.labels <- get.class.break.labels(class.breaks,lesser.sign=FALSE,greater.sign=TRUE)
  colour.ramp <- get.legend.colourbar(var.name=var.name,map.range=box.range,
                                      my.bp=0,class.breaks=class.breaks,
                                      type)        

  map.class.breaks.labels <- rev(map.class.breaks.labels)
  class.breaks <- rev(class.breaks)
  colour.ramp < rev(colour.ramp)


  lons <- c(-130,-128,-126,-124,-122,-120,-118)
  lats <- c(54,55,56,57,58,59,60)
  grat <- graticule(lons, lats, proj = ne.crs)
  labs <- graticule_labels(lons = lons, lats = lats, xline = -129, yline = 54, proj = CRS(ne.crs))

  fj.coords <- convert.to.ne.coords(lon=-120.84773,lat=56.25600,ne.crs=ne.crs) ##Fort St John
  fn.coords <- convert.to.ne.coords(lon=-122.68722,lat=58.80452,ne.crs=ne.crs) ##Fort Nelson
  dc.coords <- convert.to.ne.coords(lon=-120.23144,lat=55.76058,ne.crs=ne.crs) ##Dawson Creek
  cw.coords <- convert.to.ne.coords(lon=-121.62955,lat=55.70771,ne.crs=ne.crs) ##Chetwynd
  pc.coords <- convert.to.ne.coords(lon=-120.13359,lat=55.71034,ne.crs=ne.crs) ##Pouce Coupe
  nr.coords <- convert.to.ne.coords(lon=-122.64943,lat=59.40487,ne.crs=ne.crs) ##Northern Rockies
  tr.coords <- convert.to.ne.coords(lon=-120.99367,lat=55.12597,ne.crs=ne.crs) ##Tumbler Ridge
  ##------------------------------------------------------------------------------------------------

  ##Set up plot image
  png(file=plot.file,width=width,height=height,bg='gray94') ##,width=500,height=1100,bg='white')
  par(mar=c(6,6,7,6))    
  plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
     bg='lightgray',axes=FALSE,
       xlab='Longitude (\u00B0E)',ylab='Latitude (\u00B0N)',main=plot.title,
       cex.axis=2,cex.lab=2,cex.main=2.1)
  ##plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
  ##     bg='lightgray',
  ##     xlab='',ylab='',main=plot.title,
  ##     cex.axis=2,cex.lab=2,cex.main=2.2,axes=FALSE)
  axis(1,at=unclass(labs@coords)[1:7,1],label=lons,cex.axis=2)  
  axis(2,at=unclass(labs@coords)[8:14,2],label=lats,cex.axis=2)  

  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='lightgray')
  
  ##First plot the entire rectangular region with lighter transparency
  if (var.name=='snowdepth') {
    image(white.box, col='white',breaks=c(0,1),ribbon=FALSE,xlim=plot.window.xlim, ylim=plot.window.ylim, add=TRUE)
  }

  ##First plot the entire rectangular region with lighter transparency
  image(box.data, col=colour.ramp,breaks=class.breaks,xlim=plot.window.xlim, ylim=plot.window.ylim, add=TRUE)   
  ##Add the baselayers for the plot box
  ##bc.roads <- spTransform(get.basedata.shape('bc_hwy_geo','basedata/BC_Roads/'),CRS(proj))
  ##city.dir <- '/storage/data/gis/basedata/BC_Cities/' ##'/storage/data/projects/rci/data/forestry/regional_summaries/shapefiles/'
  


  ##Cover the extra area to diminish its effect
  ##rect(region.shp@bbox[1,1],region.shp@bbox[2,1],region.shp@bbox[1,2],region.shp@bbox[2,2],col="#D3D3D378")
  ##rect(plot.window.xlim[1],plot.window.ylim[1],plot.window.xlim[2],plot.window.ylim[2],col="#D3D3D34B")

  ##Add the region raster data
  ##image(my.raster, col=colour.ramp,breaks=class.breaks,ribbon=FALSE,xlim=plot.window.xlim, ylim=plot.window.ylim, add=TRUE)  

  ##-------------------------------------------------------------------------------------------------
  ##Add the region overlays to plot
  
  ##Overlays common to all regions    
  bc.overlay <- 'h_land_WGS84'
  coast.overlay <- 'h_coast_WGS84'
  shape.dir <- '/storage/data/projects/rci/data/assessments/northeast/shapefiles/'
  bc.shp <- readOGR(shape.dir, bc.overlay, stringsAsFactors=F, verbose=F)
  coast.shp <- readOGR(shape.dir, 'west_coast_ocean', stringsAsFactors=F, verbose=F)
  us.shp <- readOGR(shape.dir,'pnw_us_wgs84',stringsAsFactors=F, verbose=F)
  rivers.shp <- readOGR(shape.dir,'h_rivers_WGS84',stringsAsFactors=F, verbose=F)
  nr.shp <- readOGR(shape.dir,'northern_rockies_muni',stringsAsFactors=F, verbose=F)
  tr.shp <- readOGR(shape.dir,'tumbler_ridge',stringsAsFactors=F, verbose=F)
  
  plot(spTransform(coast.shp,CRS(ne.crs)),add=TRUE,col='lightgray')                
  ##plot(spTransform(rivers.shp,CRS(ne.crs)),add=TRUE,col='lightblue')
  plot(spTransform(us.shp,CRS(ne.crs)),add=TRUE,col='gray')
  plot(spTransform(bc.shp,CRS(ne.crs)),add=TRUE)
  plot(spTransform(region.shp,CRS(ne.crs)),add=TRUE,lwd=4)
  plot(spTransform(nr.shp,CRS(ne.crs)),add=TRUE,lwd=2,lty=2)
  plot(spTransform(tr.shp,CRS(ne.crs)),add=TRUE,lwd=2,lty=2)
  plot(grat,add=TRUE,lty=3,col='gray')

  if (draft) {
    text(x = grconvertX(0.5, from = "npc"),  # align to center of plot X axis
         y = grconvertY(0.5, from = "npc"), # align to center of plot Y axis
         labels = "DRAFT", # our watermark
         cex = 10, font = 2, # large, bold font - hard to miss
         col = rgb(1, 1, 1, .4), # translucent (0.2 = 20%) red color
         srt = 45) # srt = angle of text: 45 degree angle to X axis
  }
  
  ##------------------------------------------------------
  ##Plot the grid lines
#  v.cell <- my.raster@grid@cellsize[1]
#  h.cell <- my.raster@grid@cellsize[2]
#  abline(v=spatial.coords[,1]+v.cell/2,lty=2,col='black')
#  abline(h=spatial.coords[,2]+h.cell/2,lty=2,col='black')
  ##------------------------------------------------------ 

  my.label.units <- leg.label.units(var.name,type)
  
  cx <- unclass(fj.coords)[1]
  cy <- unclass(fj.coords)[2]
  points(cx,cy,pch=17,cex=2,col='black')
  shadowtext(cx,cy+15000,'Fort St. John',adj=4,pos=3,cex=1.25,col='black',bg='white',r=0.1)

  cx <- unclass(fn.coords)[1]
  cy <- unclass(fn.coords)[2]
  points(cx,cy,pch=17,cex=2,col='black')
  shadowtext(cx,cy-15000,'Fort Nelson',adj=4,pos=1,cex=1.25,col='black',bg='white',r=0.1)

  cx <- unclass(dc.coords)[1]
  cy <- unclass(dc.coords)[2]
  points(cx,cy,pch=17,cex=2,col='black')
  shadowtext(cx-8000,cy+10000,'Dawson Creek',adj=4,pos=3,cex=1.25,col='black',bg='white',r=0.1)

  cx <- unclass(cw.coords)[1]
  cy <- unclass(cw.coords)[2]
  points(cx,cy,pch=17,cex=2,col='black')
  shadowtext(cx-8000,cy,'Chetwynd',adj=4,pos=2,cex=1.25,col='black',bg='white',r=0.1)

  cx <- unclass(pc.coords)[1]
  cy <- unclass(pc.coords)[2]
  points(cx,cy,pch=17,cex=2,col='black')
  shadowtext(cx-10000,cy-12000,'Pouce Coupe',adj=4,pos=1,cex=1.25,col='black',bg='white',r=0.1)

  cx <- unclass(nr.coords)[1]
  cy <- unclass(nr.coords)[2]
  text(cx,cy,'Northern Rockies Regional Municipality',font=2,adj=4,pos=1,cex=1.25,col='black',bg='white')

  cx <- unclass(tr.coords)[1]
  cy <- unclass(tr.coords)[2]
  text(cx,cy,'Tumbler Ridge',font=2,adj=4,pos=1,cex=1.25,col='black',bg='white')


  par(xpd=NA)
  legend(leg.loc, col = "black", legend=map.class.breaks.labels, pch=22, pt.bg = rev(colour.ramp),
         pt.cex=2.0, y.intersp=0.8, title.adj=0.2, title=my.label.units, xjust=0, cex=1.7)

  box(which='plot',lwd=3)

  dev.off()


}
##-------------------------------------------------------


   

