##Script to plot comparisons of the RCM Climatologies and 

library(ncdf4)
library(raster)
library(rgdal)

source('/storage/home/ssobie/code/repos/crcm5/plot.configured.crcm5.r',chdir='T')
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir='T')

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

find.year.match <- function(dd.mon,dd.tmy,dd.time) {
  ddim <- dim(dd.mon)  
 
  close.year <- matrix(NA,nrow=ddim[1],ncol=ddim[2])
  yrs <- as.numeric(format(dd.time,'%Y'))

  for (i in 1:ddim[1]) {
    for (j in 1:ddim[2]) {
      dix <- which.min(abs(dd.mon[i,j,]-dd.tmy[i,j]))          
      close.year[i,j] <- yrs[dix]
    }
  }    
  return(close.year)  
}

find.dd.rank <- function(dd.mon,dd.tmy) {
  ddim <- dim(dd.mon)  
  dd.rank <- matrix(NA,nrow=ddim[1],ncol=ddim[2])

  for (i in 1:ddim[1]) {
    for (j in 1:ddim[2]) {
      if (sum(dd.mon[i,j,] == 0) > 10) {
        dd.rank[i,j] <- NA
      } else {
        m.cdf <- ecdf(dd.mon[i,j,])                
        dd.rank[i,j] <- round(m.cdf(dd.tmy[i,j])*30)
      }
    }
  }    
  
  return(dd.rank)
}
               
make.scatter.plot <- function(gcm,var.name) {

  shape.dir <- '/storage/data/projects/rci/data/assessments/bc/shapefiles/'
  bc.shp <- readOGR(shape.dir, 'bc', stringsAsFactors=F, verbose=F)
  albers.crs <- CRS("+init=epsg:4326")   
  mask.shp <- spTransform(bc.shp,albers.crs)

  dd.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/degree_hours/'
  tmy.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/'

  present.tmy.file <- paste0(tmy.dir,'/tas_CWEC_TMY_',gcm,'+CRCM5_historical_19810101-20101231.nc')
  ptf.nc <- nc_open(present.tmy.file)
  ptf.tas <- ncvar_get(ptf.nc,'tas')
##  pfc.rs <- make.raster.version(ptf.nc,ptf.tas)
##  mask.shp <- spTransform(bc.shp,albers.crs)
##  pfs.ms <- mask(pfc.rs,mask.shp)
##  browser()

  present.time <- netcdf.calendar(ptf.nc)  

  nc_close(ptf.nc)

  future.tmy.file <-  paste0(tmy.dir,'/tas_CWEC_TMY_',gcm,'+CRCM5_historical_20210101-20501231.nc')
  ftf.nc <- nc_open(future.tmy.file)
  ftf.tas <- ncvar_get(ftf.nc,'tas')
  future.time <- netcdf.calendar(ftf.nc)  
  nc_close(ftf.nc)

  morphed.tmy.file <- paste0(tmy.dir,'/morphed_tas_CWEC_TMY_',gcm,'+CRCM5_historical+rcp85_20210101-20501231.nc')
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

  dd.past <- dd.data[,,pst:pen]
  dd.past.time <- dd.time[pst:pen]
  dd.proj <- dd.data[,,fst:fen]
  dd.proj.time <- dd.time[pst:pen]

  ##---------------------
  plot.file <- paste0('/storage/data/projects/rci/data/RCM/CRCM5/plots/bc/degree_hours/',
                      var.name,'.',gcm,'.all.months.scatter.plot.png')
##  png(plot.file,width=1500,height=500)
  par(mfrow=c(2,3))
  par(mar=c(3,3,3,5))
  par(oma=c(1,1,0,0))
  par(mgp=c(2,1,0))
 
  for (mn in 5:9) {
    mon <- sprintf('%02d',mn)  
    mn.px <- format(dd.past.time,'%m') == mon
    dd.past.mon <- dd.past[,,mn.px]
    dd.past.tmy <- present.tmy.dd[,,mn]

    dd.past.years <- find.year.match(dd.past.mon,dd.past.tmy,dd.past.time[mn.px])    
    dd.past.ranks <- find.dd.rank(dd.past.mon,dd.past.tmy)    

    mn.fx <- format(dd.proj.time,'%m') == mon
    dd.proj.mon <- dd.proj[,,mn.fx]
    dd.proj.tmy <- future.tmy.dd[,,mn]
    dd.morph.tmy <- morphed.tmy.dd[,,mn]

    dd.proj.years <- find.year.match(dd.proj.mon,dd.proj.tmy,dd.proj.time[mn.fx])    
    dd.proj.ranks <- find.dd.rank(dd.proj.mon,dd.proj.tmy)    

    dd.morph.years <- find.year.match(dd.proj.mon,dd.morph.tmy,dd.proj.time[mn.fx])    
    dd.morph.ranks <- find.dd.rank(dd.proj.mon,dd.morph.tmy)    
    
    past.rs <- make.raster.version(ptf.nc,dd.past.ranks)
    ms.past.ranks <- as.matrix(mask(past.rs,mask.shp))
    proj.rs <- make.raster.version(ptf.nc,dd.proj.ranks)
    ms.proj.ranks <- as.matrix(mask(proj.rs,mask.shp))
    morph.rs <- make.raster.version(ptf.nc,dd.morph.ranks)
    ms.morph.ranks <- as.matrix(mask(morph.rs,mask.shp))
    
##    shared.max <- round(max(c(ms.past.ranks,ms.proj.ranks,ms.morph.ranks),na.rm=T)/100)*100
    hp <- hist(ms.past.ranks,breaks=seq(0,30,1),plot=FALSE)
    hf <- hist(ms.proj.ranks,breaks=seq(0,30,1),plot=FALSE)
    hm <- hist(ms.morph.ranks,breaks=seq(0,30,1),plot=FALSE)
    shared.max <- max(c(hp$counts,hf$counts,hm$counts))
    ##browser()
    hp <- hist(ms.past.ranks,breaks=seq(0,30,1),ylim=c(0,shared.max),
         col=alpha('blue',0.5),border=alpha('blue',0.4),
         main=month.abb[mn],xlab='Ranks',ylab='Frequency')
    par(lwd=2)
    hf <- hist(ms.proj.ranks,breaks=seq(0,30,1),add=T,border='red')
    hm <- hist(ms.morph.ranks,breaks=seq(0,30,1),add=T,border='green')
    box(which='plot')     
    rect(22,0.8*par('usr')[4],28,0.99*par('usr')[4],col='gray94',border='lightgray')
    text(25,0.95*par('usr')[4],round(mean(ms.past.ranks,na.rm=T),2),col='blue')
    text(25,0.9*par('usr')[4],round(mean(ms.proj.ranks,na.rm=T),2),col='red')
    text(25,0.85*par('usr')[4],round(mean(ms.morph.ranks,na.rm=T),2),col='black')

##    browser()
  }
##  dev.off()
  nc_close(dd.nc)                       
}

##***************************************************************

make.scatter.plot('MPI','cdd')