library(ncdf4)
library(PCICt)

sub.by.time <- function(read.dir,read.file,interval,write.dir) {

  nc <- nc_open(paste(read.dir,read.file,sep=''))
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.series <- origin.pcict + time.values*3600
  years <- format(time.series,'%Y')
  yrs <- strsplit(interval,'-')[[1]]

  st <- head(grep(substr(yrs[1],1,4),years),1)-1
  en <- tail(grep(substr(yrs[2],1,4),years),1)-1

  sub.file <- 'sub.nc'
  system(paste('ncks --overwrite -d time,',st,',',en,' ',read.dir,read.file,' ',write.dir,sub.file,sep=''))
  print(paste('ncks --overwrite -d time,',st,',',en,' ',read.dir,read.file,' ',write.dir,sub.file,sep=''))
  nc_close(nc)
  rv <- paste0(write.dir,sub.file)
  return(rv)
}

make.rcm.climatologies <- function(read.dir,read.file,interval,write.dir,write.file,fxn) {
    sub.file <- sub.by.time(read.dir,read.file,interval=interval,write.dir)
    work <- paste0('cdo -s -O ymonmean ',sub.file,' ',write.dir,write.file)
    print(work)
    system(work)
    system(paste('rm ',sub.file,sep=''))
}

make.anomalies <- function(read.dir,past.file,proj.file,write.dir,write.file) {
  work <- paste('cdo -s -O sub ',read.dir,proj.file,' ',read.dir,past.file,' ',write.dir,write.file,sep='')
  print(work)
  system(work)
}

make.tmy.climatologies <- function(read.dir,read.file,write.dir,write.file,fxn) {
  work1 <- paste0('cdo -O day',fxn,' ',read.dir,read.file,' ',write.dir,'tmp.nc')
  print(work1)
  system(work1)

  work2 <- paste0('cdo -s -O monmean ',write.dir,'tmp.nc ',write.dir,write.file)                        
  print(work2)
  system(work2)

  system(paste('rm ',write.dir,'tmp.nc',sep=''))

}

##-------------------------------------------------------------
##RCM Climatologies and Anomalies
rcm.clims <- function(gcm,var.name,interval) {

  read.dir <-  '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/'
  write.dir  <- paste0(read.dir,'climatologies/')   

  all.files <- list.files(path=read.dir,pattern=gcm)
  read.file <- all.files[grep(var.name,all.files)]
  print(read.file)
  time.write <- gsub(pattern='[0-9]{8}-[0-9]{8}',replacement=interval,read.file)
  write.file <- gsub(pattern='_day_',replacement='_climatology_',time.write)
  print(write.file)
  make.rcm.climatologies(read.dir,read.file,interval,write.dir,write.file)
}

rcm.anoms <- function(gcm,var.name,past.int,proj.int) {

  read.dir <-  '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/climatologies/'
  write.dir  <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/anomalies/'

  all.files <- list.files(path=read.dir,pattern=gcm)
  var.files <- all.files[grep(paste0(var.name,'_climatology_'),all.files)]
  future.file <- var.files[grep(proj.int,var.files)]
  past.file <- var.files[grep(past.int,var.files)]
  write.file <- gsub(pattern='_climatology_',replacement='_anomaly_',future.file)
  make.anomalies(read.dir,past.file,future.file,write.dir,write.file)
}


##-------------------------------------------------------------
##TMY Climatologies and Anomalies
tmy.clims <- function(gcm,var.name,interval,fxn,morphed=FALSE) {

  read.dir <-  '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/'
  write.dir  <- paste0(read.dir,'climatologies/')   

  read.file <- paste0(var.name,'_CWEC_TMY_',gcm,'+CRCM5_historical_',interval,'.nc')
  if (morphed)
     read.file <- paste0('morphed_',var.name,'_CWEC_TMY_',gcm,'+CRCM5_historical+rcp85_',interval,'.nc')
  print(read.file)
  write.file <- gsub(pattern=paste0(var.name,'_'),replacement=paste0(var.name,'_',fxn,'_climatology_'),read.file)
  print(write.file)
  make.tmy.climatologies(read.dir,read.file,write.dir,write.file,fxn)
}

tmy.anoms <- function(gcm,var.name,past.int,proj.int,fxn,morphed=FALSE) {

  read.dir <-  '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/climatologies/'
  write.dir  <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/anomalies/'
  past.file <- paste0(var.name,'_',fxn,'_climatology_CWEC_TMY_',gcm,'+CRCM5_historical_',past.int,'.nc')
  future.file <- paste0(var.name,'_',fxn,'_climatology_CWEC_TMY_',gcm,'+CRCM5_historical_',proj.int,'.nc')
  if (morphed)
    future.file <- paste0('morphed_',var.name,'_',fxn,'_climatology_CWEC_TMY_',gcm,'+CRCM5_historical+rcp85_',proj.int,'.nc')
  write.file <- gsub(pattern='_climatology_',replacement='_anomaly_',future.file)
  make.anomalies(read.dir,past.file,future.file,write.dir,write.file)

}

##**************************************************************

##--------------------------------------------------------------

##var.names <- c('dewpoint_max','dewpoint_mean','dewpoint_min',
##               'tas_max','tas_mean','tas_min',
##               'huss_mean','insol_mean','psl_mean',
##               'rhs_mean','wspd_mean','wspd_max')


if (1==0) {
 gcms <- c('CanESM2','MPI')
 for (gcm in gcms) {
   var.name <- 'dewpoint_mean'
   rcm.clims(gcm=gcm,var.name=var.name,interval='19810101-20101231')
   rcm.clims(gcm=gcm,var.name=var.name,interval='20210101-20501231')
   rcm.anoms(gcm=gcm,var.name=var.name,past.int='19810101-20101231',proj.int='20210101-20501231')
 }
}

if (1==1) {
gcms <- c('CanESM2','MPI')
var.name <- 'dewpoint'
fxn <- 'mean'
for (gcm in gcms) {
  tmy.clims(gcm=gcm,var.name=var.name,interval='19810101-20101231',fxn)
  tmy.clims(gcm=gcm,var.name=var.name,interval='20210101-20501231',fxn)
  tmy.clims(gcm=gcm,var.name=var.name,interval='20210101-20501231',fxn,morphed=TRUE)
  tmy.anoms(gcm=gcm,var.name=var.name,past.int='19810101-20101231',proj.int='20210101-20501231',fxn)
  tmy.anoms(gcm=gcm,var.name=var.name,past.int='19810101-20101231',proj.int='20210101-20501231',fxn,morphed=TRUE)
}
}
