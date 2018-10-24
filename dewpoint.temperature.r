##Script to calculate dewpoint temperatures

library(ncdf4)
library(PCICt)
library(doParallel)
registerDoParallel(cores=4)
library(foreach)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/crcm5/add.crcm5.metadata.r',chdir=T)


##Dew Point Temperature
dew.point.temp <- function(psl,huss) {
     epsilon  <<- 0.622      ##kgw/kga    ::kg water / kg dry air
     T.zero   <<- 273        ##K          ::Zero in Kelvin
     R.vapor  <<- 461        ##J K-1 kg-1 ::Gas constant for water vapor
     lh.vape  <<- 2.501*10^6 ## J kg-1    ::latent heat of vaporization at 273 K
     sat.vape <<- 0.611      ##kPa        ::Saturation vapour pressure at 273 K

     ##Vapor Pressure
     vape.press <- huss * (psl/10) / epsilon

     ##Dry bulb temp
     temp.dew <- (1/T.zero - ((R.vapor/lh.vape) * log(vape.press/sat.vape)))^-1

     return(temp.dew - 273)
}


##----------------------------------------------------------------------------------------------

make.new.dwpt.file <- function(gcm,freq,rcp,
                               psl.file,huss.file,
                               tmp.base) {
  new.var <- 'dewpoint'
  nc <- nc_open(paste0(tmp.base,psl.file),write=FALSE)
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.series <- origin.pcict + time.values*3600

  yst <-  gsub('-','',format(head(time.series,1),'%Y-%m-%d'))
  yen <-  gsub('-','',format(tail(time.series,1),'%Y-%m-%d'))

  gcm.name <- switch(gcm,
                     ERAI='ERA-Interim',
                     CanESM2='CanESM2',
                     MPI='MPI')

  write.file <- paste0(new.var,'_',freq,'_BC_WC011_',gcm.name,'+CRCM5_historical+',rcp,'_',yst,'-',yen,'.nc')

  ##Attributes to retain
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')

  rlon <- ncvar_get(nc,'rlon')
  rlat <- ncvar_get(nc,'rlat')

  var.atts <- ncatt_get(nc,'psl')

  n.lon <- length(rlon)
  n.lat <- length(rlat)

  ##--------------------------------------------------------------
  ##Create new netcdf file

  x.geog <- ncdim_def('rlon', 'degrees', rlon)
  y.geog <- ncdim_def('rlat', 'degrees', rlat)
  t.geog <- ncdim_def('time', time.units, time.values,
                      unlim=FALSE, calendar=time.calendar)

  lon.geog <- ncvar_def('lon', units='degrees_east', dim=list(x.geog, y.geog))
  lat.geog <- ncvar_def('lat', units='degrees_north',dim=list(x.geog, y.geog))
  var.geog <- ncvar_def(new.var, units=var.atts$units, dim=list(x.geog, y.geog,t.geog),
                        missval=1.e+20)
  proj.geog <- ncvar_def('projection',units='',dim=list(),prec='char')

  hist.nc <- nc_create(paste(tmp.base,write.file,sep=''), list(lon.geog,lat.geog,proj.geog,var.geog),h_minfree=1048570)

  ##Add Attributes
  time.atts <- get.time.atts(time.calendar,time.units)
  time.names <- names(time.atts)
  for (t in 1:length(time.atts))
    ncatt_put(hist.nc,varid='time',attname=time.names[t],attval=time.atts[[t]])

  standard.atts <- get.standard.atts()
  standard.names <- names(standard.atts)
  print('Lon names')
  lon.names <- names(standard.atts$lon)
  for (j in 1:length(standard.atts$lon))
    ncatt_put(hist.nc,varid='lon',attname=lon.names[j],attval=standard.atts$lon[[j]])
  print('Lat names')
  lat.names <- names(standard.atts$lat)
  for (j in 1:length(standard.atts$lat))
    ncatt_put(hist.nc,varid='lat',attname=lat.names[j],attval=standard.atts$lat[[j]])
  ncvar_put(hist.nc,'lon',lon)
  ncvar_put(hist.nc,'lat',lat)

  print('RLon names')
  rlon.names <- names(standard.atts$rlon)
  for (j in 1:length(standard.atts$rlon))
    ncatt_put(hist.nc,varid='rlon',attname=rlon.names[j],attval=standard.atts$rlon[[j]])
  print('RLat names')
  rlat.names <- names(standard.atts$rlat)
  for (j in 1:length(standard.atts$rlat))
    ncatt_put(hist.nc,varid='rlat',attname=rlat.names[j],attval=standard.atts$rlat[[j]])
  print('Proj names')
  proj.names <- names(standard.atts$proj)
  for (j in 1:length(standard.atts$proj))
    ncatt_put(hist.nc,varid='projection',attname=proj.names[j],attval=standard.atts$proj[[j]])

  print('Global Atts')
  global.atts <- get.global.atts(gcm,freq,rcp)
  global.names <- names(global.atts)
  for (g in 1:length(global.atts))
    ncatt_put(hist.nc,varid=0,attname=global.names[g],attval=global.atts[[g]])
    print('Variable Atts')
  var.atts <- get.variable.atts('dewpoint')
  varnames <- names(var.atts)
  for (j in 1:length(var.atts))
    ncatt_put(hist.nc,varid=new.var,attname=varnames[j],attval=var.atts[[j]])

  ##Clear extraneous history
  ncatt_put(hist.nc,varid=0,attname='history',attval='')

  nc_close(hist.nc)
  nc_close(nc)

  return(write.file)

}



##-----------------------------------------------------------------------------------------

calculate.dewpoint.temp <- function(psl.file,huss.file,dewpoint.file,tmp.base) {

  psl.nc <- nc_open(paste0(tmp.base,psl.file))
  huss.nc <- nc_open(paste0(tmp.base,huss.file))
  dwpt.nc <- nc_open(paste0(tmp.base,dewpoint.file),write=TRUE)

  time <- netcdf.calendar(psl.nc)
  yrs <- levels(as.factor(format(time,'%Y')))
  n.lat <- psl.nc$dim$rlat$len ##Latitude Length
  n.lon <- psl.nc$dim$rlon$len ##Longitude Length
  n.time <- psl.nc$dim$time$len ##Longitude Length

  for (j in 1:n.lat) { 
    print(paste0('Latitude: ',j,' of ',n.lat))
    psl.subset <- ncvar_get(psl.nc,'psl',start=c(1,j,1),count=c(-1,1,-1))
    psl.list <- lapply(seq_len(nrow(psl.subset)), function(k) psl.subset[k,])
    rm(psl.subset)

    huss.subset <- ncvar_get(huss.nc,'huss',start=c(1,j,1),count=c(-1,1,-1))
    huss.list <- lapply(seq_len(nrow(huss.subset)), function(k) huss.subset[k,])
    rm(huss.subset)

    dewpoints <- foreach(
                       psl=psl.list,                 
                       huss=huss.list,                 
                       .export=c('dew.point.temp')
                             ) %dopar% {
                                objects <- dew.point.temp(psl=psl,huss=huss)
                                    }
    dwpt.matrix <- matrix(unlist(dewpoints),nrow=n.lon,ncol=n.time,byrow=T)
    ncvar_put(dwpt.nc,varid='dewpoint',vals=dwpt.matrix,
                      start=c(1,j,1),count=c(-1,1,n.time))
  }
  nc_close(psl.nc)
  nc_close(huss.nc)
  nc_close(dwpt.nc)

}




##----------------------------------------------------------------------------------------------

if (1==1) {
  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
}

##gcm <- 'ERAI'
##tmpdir <- '/local_temp/ssobie/crcm5/'

read.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/hourly/'
tmp.dir <- tmpdir
print(tmp.dir)
print(tmpdir)
if (!file.exists(tmpdir)) {
  dir.create(tmpdir,recursive=T)
}

psl.file <- paste0('psl_hour_BC_WC011_',gcm,'+CRCM5_historical+rcp85_19800101-20501231.nc')
huss.file <- paste0('huss_hour_BC_WC011_',gcm,'+CRCM5_historical+rcp85_19800101-20501231.nc')

if (1==1) {
print('Copying PSL file')
file.copy(from=paste0(read.dir,psl.file),to=tmp.dir,overwrite=T)
print('Copying HUSS file')
file.copy(from=paste0(read.dir,huss.file),to=tmp.dir,overwrite=T)

dewpoint.file <- make.new.dwpt.file(gcm,freq='hour','rcp85',
                                    psl.file,huss.file,
                                    tmp.dir)
calculate.dewpoint.temp(psl.file,huss.file,dewpoint.file,tmp.dir)

file.copy(from=paste0(tmp.dir,dewpoint.file),to=read.dir,overwrite=T)


file.remove(paste0(tmp.dir,psl.file))
print('Done with PSL')
file.remove(paste0(tmp.dir,huss.file))
print('Done with HUSS')
file.remove(paste0(tmp.dir,dewpoint.file))
print('Done with Dewpoint')

}


