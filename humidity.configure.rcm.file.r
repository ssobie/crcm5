##Script to configure raw CRCM5 humidity into a useful format.

library(ncdf4)
library(PCICt)
library(udunits2)

##--------------------------------------------------------------
##Metadata
source('/storage/home/ssobie/code/repos/crcm5/add.crcm5.metadata.r',chdir=T)
##--------------------------------------------------------------

correct.rcm.time.series <- function() {
  time.calendar <- "proleptic_gregorian"                        
  time.units <- "hours since 1980-01-01 00:00:00"   
  new.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                         cal=time.calendar)
  new.time <- seq(from=as.Date('1980-01-01'),by='day',to=as.Date('2050-12-31'))
  test.values <- seq(0,by=3,length.out=25933*8)
  test.series <- new.origin + test.values*3600      
  feb.flag <- grep('*-02-29 *',test.series)

  rv <- list(values=test.values[-feb.flag],
             series=test.series[-feb.flag],
             units="hours since 1980-01-01 00:00:00",
             calendar="proleptic_gregorian")
}

##--------------------------------------------------------------
make.new.netcdf.file <- function(gcm,varname,freq,rcp,
                                 existing.file,
                                 tmp.base) {
  ##Vancouver Island Subset
  ##lon.ix <- 115:155
  ##lat.ix <- 85:125
  ##BC Subset
  ##lon.ix <- 100:200
  ##lat.ix <- 79:230

  new.var <- get.new.var.name(varname)                                 

  nc <- nc_open(paste0(tmp.base,existing.file),write=FALSE)

  time.data <- correct.rcm.time.series() 

  yst <-  gsub('-','',format(head(time.data$series,1),'%Y-%m-%d'))
  yen <-  gsub('-','',format(tail(time.data$series,1),'%Y-%m-%d'))

  gcm.name <- switch(gcm,
                     ERAI='ERA-Interim',
                     CanESM2='CanESM2',
                     MPI='MPI')

  ###write.file <- paste0(new.var,'_',freq,'_BC_WC011_',gcm.name,'+CRCM5_historical',rcp,'_',yst,'-',yen,'.nc')
write.file <- paste0(new.var,'_',freq,'_WC011_',gcm.name,'+CRCM5_historical',rcp,'_',yst,'-',yen,'.nc')
 
  ##Attributes to retain
  lon <- aperm(ncvar_get(nc,'lon'),c(2,1))
  lat <- aperm(ncvar_get(nc,'lat'),c(2,1))

  rlon <- ncvar_get(nc,'rlon')
  rlat <- ncvar_get(nc,'rlat')  
     
  var.units <- switch(varname,
                      HU="kg/kg",
                      HR='fraction')
  n.lon <- length(rlon)
  n.lat <- length(rlat)

  ##--------------------------------------------------------------
  ##Create new netcdf file

  x.geog <- ncdim_def('rlon', 'degrees', rlon)###[lon.ix])
  y.geog <- ncdim_def('rlat', 'degrees', rlat)###[lat.ix])

  t.geog <- ncdim_def('time', time.data$units, time.data$values,
                      unlim=TRUE, calendar=time.data$calendar)

  lon.geog <- ncvar_def('lon', units='degrees_east', dim=list(x.geog, y.geog))
  lat.geog <- ncvar_def('lat', units='degrees_north',dim=list(x.geog, y.geog)) 
  var.geog <- ncvar_def(new.var, units=var.units, dim=list(x.geog, y.geog,t.geog),
                        missval=1.e+20)
  proj.geog <- ncvar_def('projection',units='',dim=list(),prec='char') 

  hist.nc <- nc_create(paste(tmp.base,write.file,sep=''), list(lon.geog,lat.geog,proj.geog,var.geog),h_minfree=1048570)

  ##Add Attributes
  time.atts <- get.time.atts(time.data$calendar,time.data$units)
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
  var.atts <- get.variable.atts(new.var)
  varnames <- names(var.atts)
  for (j in 1:length(var.atts))
    ncatt_put(hist.nc,varid=new.var,attname=varnames[j],attval=var.atts[[j]])

  ##Clear extraneous history
  ncatt_put(hist.nc,varid=0,attname='history',attval='')

  ncvar_put(hist.nc,varid='lon',vals=lon)###[lon.ix,lat.ix])
  ncvar_put(hist.nc,varid='lat',vals=lat)###[lon.ix,lat.ix])

  nc_close(hist.nc)  
  nc_close(nc)  

  return(write.file)
}

move.data.to.new.file <- function(gcm,varname,nc,
                                  file,write.file,
                                  tmp.base) {
  hist.nc <- nc_open(paste0(tmp.base,file),write=FALSE)
  new.var <- get.new.var.name(varname)                                 

  time.atts <- ncatt_get(hist.nc,'t')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(hist.nc,'t')
  time.len <- length(time.values)
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.series <- origin.pcict + time.values*3600
  mulc <- 3600
  if (any(grepl('*-02-29 *',time.series))) {
    print('Found Feb 29')
    feb.flag <- grep('*-02-29 *',time.series)
    time.values[feb.flag] <- time.values[feb.flag] - 24
    fixed.series <- origin.pcict + time.values*mulc
    time.series <- fixed.series
  }
  full.series <- correct.rcm.time.series()  
  ix <- which(full.series$series %in% time.series)
  ist <- ix[1]
  icnt <- tail(ix,1)-ist + 1
  ###var.subset <- ncvar_get(hist.nc,varname,start=c(85,115,22,1),count=c(41,41,1,-1))
  ###var.subset <- ncvar_get(hist.nc,varname,start=c(79,100,22,1),count=c(152,101,1,-1))
  var.subset <- ncvar_get(hist.nc,varname,start=c(1,1,22,1),count=c(-1,-1,1,-1))

  ncvar_put(nc,varid=new.var,vals=aperm(var.subset,c(2,1,3)),
            start=c(1,1,ist),count=c(-1,-1,icnt))     
  nc_close(hist.nc)  
}

##**************************************************************************************
##-------------------------------------------------------------------------


if (1==1) {
  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
}

###  gcm <- 'CanESM2'
###  varname <- 'HU'
###  interval <- '1980-2050'
###  freq <- 'hour'

##  write.file <- 'tasmax_day_WC011_ERA-Interim+CRCM5_historical+rcp85_198001-20141231.nc'  

  rcp <- "+rcp85"
  if (gcm=='ERAI') {
   rcp <- ''
  }

  ##tmp.base <- tmpdir
  tmp.base <- paste('/local_temp/ssobie/crcm5/',varname,'/',sep='')
  
  ##data.dir <- paste0('/storage/data/climate/downscale/CMIP5_delivery/CRCM5/MPI/',varname,'/')
  data.dir <- paste0('/storage/data/climate/downscale/RCM/CRCM5/',gcm,'/',varname,'/')
  ###write.dir <- paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/hourly/')
  write.dir <- paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/',gcm,'/')

  if (!file.exists(tmp.base))
       dir.create(tmp.base,recursive=TRUE)

  starting.file <- paste0(gcm,'_WC011_modified_T5_1980-2050.nc')

  print('Copying file to temp')
  print(paste0('/storage/data/climate/downscale/RCM/CRCM5/',gcm,'/',starting.file))
  file.copy(from=paste0('/storage/data/climate/downscale/RCM/CRCM5/',gcm,'/',starting.file),to=tmp.base,overwrite=TRUE)

  ##-------------------------------------------------    
  ##Past File           
  print('Creating new file')                                           
  write.file <- make.new.netcdf.file(gcm,varname,freq,rcp,
                                     starting.file,
                                     tmp.base)
  
  ##Loop over the set of hourly humidity files
  file.list <- list.files(path=data.dir,pattern=varname)

  nc <- nc_open(paste0(tmp.base,write.file),write=TRUE)

  for (file in file.list) {
    file.copy(from=paste0(data.dir,file),to=tmp.base,overwrite=TRUE) 
    print('Copying')
    print(paste0(data.dir,file))

    move.data.to.new.file(gcm,varname,nc,
                          file,write.file,
                          tmp.base)
    print('File removal')
    file.remove(paste0(tmp.base,file))
  }


  nc_close(nc)

  print('Copying file back')     
  file.copy(from=paste0(tmp.base,write.file),to=write.dir,overwrite=TRUE)

  file.remove(paste0(tmp.base,starting.file))
  file.remove(paste0(tmp.base,write.file))
##-------------------------------------------------------------------------


