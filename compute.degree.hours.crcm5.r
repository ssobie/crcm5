##Script to calculate the degree-day indices from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file

##Updated version from compute.climdex.bccaq.r
##This computes all the climdex variables

library(ncdf4)
library(PCICt)

library(doParallel)
registerDoParallel(cores=2) 

##--------------------------------
##Degree Day Values
dd <- function(tmean,tbase) {
  g <- tmean - tbase
  days <- sum(g[g>0], na.rm=T)
  return(round(days))
}

fd <- function(tmin) {
  days <- sum(tmin>0,na.rm=T)
  return((days))
}

s3 <- function(tmax) {
  days <- sum(tmax>30,na.rm=T)
  return((days))
}

gdd<-function(data,fac){tapply(data,fac, dd, tbase=5)}   ##Growing degree days
cdd<-function(data,fac){tapply(data,fac, dd, tbase=18)}  ##Cooling degree days
hdd<-function(data,fac){tapply(-data,fac,dd, tbase=-18)} ##Heating degree days
fdd<-function(data,fac){tapply(-data,fac,dd, tbase=0)} ##Freezing degree days
ffd<-function(data,fac){tapply(data,fac,fd)} ##Frost Free days
s30<-function(data,fac){tapply(data,fac,s3)} ##Days over 30 degC

dd.fxns <- list(gdd=gdd,
                cdd=cdd,   
                hdd=hdd,
                fdd=fdd,
                s30=s30)

create.base.files <- function(degree.name,tas.file,gcm,scenario,
                              tmp.dir) {

  write.clim.name <- gsub('tas_hour',degree.name,tas.file)
  print(write.clim.name)
  nc <- nc_open(paste0(tmp.dir,tas.file),write=FALSE)
  
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.start <- as.Date(strsplit(time.units, ' ')[[1]][3])
  past.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                          cal=time.calendar)
  time.values <- ncvar_get(nc,'time')
  time.series <- format(past.origin + (time.values)*3600,'%Y-%m-%d %H')

  ##years.ix <- grep('*-01-01 00',time.series)
  ##years <- time.values[years.ix]
  ##dates <- as.numeric(years)
  months.ix <- grep('*-*-01 00',time.series)
  dates <- time.values[months.ix]

  ##Attributes to retain
  rlon <- ncvar_get(nc,'rlon')
  rlat <- ncvar_get(nc,'rlat')  
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')  

  rlon.atts <- ncatt_get(nc,'rlon')
  rlat.atts <- ncatt_get(nc,'rlat')
  lon.atts <- ncatt_get(nc,'lon')
  lat.atts <- ncatt_get(nc,'lat')

  global.atts <- ncatt_get(nc,varid=0)
  
  n.lon <- length(rlon)
  n.lat <- length(rlat)

  ##--------------------------------------------------------------
  ##Create new netcdf file
  x.geog <- ncdim_def('rlon', 'degrees', rlon)
  y.geog <- ncdim_def('rlat', 'degrees', rlat)
  t.geog <- ncdim_def('time', time.units, dates,
                      unlim=TRUE, calendar=time.calendar)
  lon.geog <- ncvar_def('lon', units='degrees_east', dim=list(x.geog, y.geog))
  lat.geog <- ncvar_def('lat', units='degrees_north',dim=list(x.geog, y.geog))

  var.geog <- ncvar_def(degree.name, units='degree hours', dim=list(x.geog, y.geog, t.geog),
                        missval=-32768)  
  file.nc <- nc_create(paste0(tmp.dir,write.clim.name), list(lon.geog,lat.geog,var.geog))
  print(paste0(tmp.dir,write.clim.name))
  ##Loop over subsets of the time series
  ##Past file first
  global.names <- names(global.atts)
  for (g in 1:length(global.atts)) 
    ncatt_put(file.nc,varid=0,attname=global.names[g],attval=global.atts[[g]])

  ##Time attributes
  ncatt_put(file.nc,varid='time',attname='units',attval=time.units)
  ncatt_put(file.nc,varid='time',attname='long_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='standard_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='calendar',attval=time.calendar)  

  ncvar_put(file.nc,varid='lon',vals=lon)
  ncvar_put(file.nc,varid='lat',vals=lat)
  
  lon.names <- names(lon.atts)
  for (j in 1:length(lon.atts))  
    ncatt_put(file.nc,varid='lon',attname=lon.names[j],attval=lon.atts[[j]])
  
  lat.names <- names(lat.atts)
  for (j in 1:length(lat.atts))  
    ncatt_put(file.nc,varid='lat',attname=lat.names[j],attval=lat.atts[[j]])

  ##Climdex Attributes
  ncatt_put(file.nc,varid=degree.name,attname='units',attval='degree hours')
  ncatt_put(file.nc,varid=degree.name,attname='_FillValue',attval=-32768)
  ncatt_put(file.nc,varid=degree.name,attname='standard_name',attval=toupper(degree.name))
  ncatt_put(file.nc,varid=degree.name,attname='long_name',attval=toupper(degree.name))

  nc_close(file.nc)
  nc_close(nc)
}


##---------------------------------------------------------------

degree.days.for.model <- function(degree.names,
                                  tas.file,gcm,scenario,tmp.dir) {

  degree.files <- list.files(path=tmp.dir,pattern='dd_',full.name=TRUE)
  print(degree.files)
  clim.ncs <- lapply(degree.files,nc_open,write=TRUE)

  ##--------------------------------------------------------------
  
  tas.nc <- nc_open(paste0(tmp.dir,tas.file),write=FALSE) 
  tas.units <- ncatt_get(tas.nc,'tas')$units
  temp.offset <- 0
  if (tas.units == 'K')
    temp.offset <- 273
  
  ##Attributes to retain
  rlon <- ncvar_get(tas.nc,'rlon')
  rlat <- ncvar_get(tas.nc,'rlat')  
  n.lon <- length(rlon)
  n.lat <- length(rlat)

  ##Combine the dates
  time.atts <- ncatt_get(tas.nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  past.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                          cal=time.calendar)
  tas.values <- ncvar_get(tas.nc,'time')
  tas.dates <- past.origin + tas.values*3600
  ##yearly.fac <- as.factor(format(tas.dates,'%Y'))
  monthly.fac <- as.factor(format(tas.dates,'%Y-%m'))

  ##--------------------------------------------------------------
  ##Compute degree day values and load into newly created netcdf
  for (i in 1:n.lon) {
    print(paste('Lon: ',i,' in ',n.lon,sep=''))
    tas.subset <- ncvar_get(tas.nc,'tas',start=c(i,1,1),count=c(1,-1,-1))-temp.offset

    flag <- is.na(tas.subset[,1])
    temp.list <- list()

    for (j in 1:n.lat) {
      temp.list[[j]] <- tas.subset[j,]
    } 

    for (k in seq_along(degree.names)) {
      degree.name <- degree.names[k]
      clim.nc <- clim.ncs[[k]]
      fx <- dd.fxns[[degree.names[k]]]
      ptm <- proc.time()
      degree.values <- foreach(
                         temp=temp.list,
                         .export=c('monthly.fac','fx')
                         ) %dopar% { 
                              degree.values <- fx(temp,monthly.fac)
                         }
      ncol <- length(degree.values[[1]])         

      degree.matrix <- matrix(unlist(degree.values),nrow=n.lat,ncol=ncol,byrow=TRUE)
      degree.matrix[flag,] <- NA
      print('Parallel version')
      print(proc.time() - ptm)

      ncvar_put(clim.ncs[[k]],varid=degree.names[k],vals=degree.matrix,
                start=c(i,1,1),count=c(1,-1,-1))
    }##Degree day loop    
  }##Longitude Loop  
  lapply(clim.ncs,nc_close)
}##Function end

##**************************************************************************************
##**************************************************************************************
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------

##Degree Days from the CRCM5 Data
run.crcm5.dd <- function(gcm,degree.names) {
  
  scenario <- 'rcp85'

  data.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/hourly/'
  write.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/degree_hours/'

  ##Move data to local storage for better I/O
  tmp.dir <- '/local_temp/ssobie/crcm5/'

  if (!file.exists(tmp.dir))
    dir.create(tmp.dir,recursive=TRUE)    
  
  degree.names <- sort(degree.names)
 
  rcm.files <- list.files(path=data.dir,pattern=gcm)
  tas.file <- rcm.files[grep('tas_hour',rcm.files)]

  file.copy(from=paste0(data.dir,tas.file),to=tmp.dir,overwrite=T)

  print('Create new files')
  first <- lapply(degree.names,create.base.files,
                  tas.file,gcm,scenario,
                  tmp.dir)

  print('Calculate Degree Days')
  second <- degree.days.for.model(degree.names,
                                  tas.file,gcm,scenario,tmp.dir)
  dd.files <- list.files(path=tmp.dir,pattern='dd_')
  file.copy(from=paste0(tmp.dir,dd.files),to=write.dir,overwrite=TRUE)
  file.remove(paste0(tmp.dir,tas.file))
  file.remove(paste0(tmp.dir,dd.files))
}

##**************************************************************************************

degree.names <- c('cdd',
                  'gdd',
                  'hdd',
                  'fdd')                   
run.crcm5.dd('CanESM2',degree.names)

