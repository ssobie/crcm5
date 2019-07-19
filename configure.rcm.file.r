##Script to calculate return periods from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file
##Modified to add the confidence intervals to the resulting netcdf file

library(ncdf4)
library(PCICt)
library(udunits2)

##--------------------------------------------------------------
##Metadata
source('/storage/home/ssobie/code/repos/crcm5/add.crcm5.metadata.r',chdir=T)
##--------------------------------------------------------------

convert.variable <- function(varname,data) {

  rv <- switch(varname,
               uas=data*0.5144,
               vas=data*0.5144,
               tasmax=data-273,
               tasmin=data-273,
               tas=data-273,
               pr=data/8,
               psl=data,
               ps=data,
               snd=data,
               snc=data,
               runoff=data,       
               discharge=data,
               streamflow=data,
               drainage=data,
               gz=data*10,
               insol=data,
               huss=data,       
               rhs=data,        
               irflux=data,
               latflux=data,    
               giac=data,
               giml=data,
               gld=data,
               glf=data,
               gsab=data,
               gsac=data,
               gsml=data,
               gvol=data,
               gwst=data,
               swater=data,
               sice=data,
               spw=data,
               smass=data,
               smelt=data,      
               asp=data/8,
               swlake=data,
               swriver=data)

  return(rv)               
}


correct.rcm.time.series <- function(nc,freq) {

  time.atts <- ncatt_get(nc,'t')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'t')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.series <- origin.pcict + time.values*3600 ##86400      

  if (freq=='day' & grepl('hour',time.units)) { ##Convert hourly to daily
     time.series <- origin.pcict + time.values*3600
     daily.series <- as.factor(format(time.series,'%Y-%m-%d'))
     time.values <- seq(from=0,by=1,length.out=length(levels(daily.series)))
     time.series <- levels(daily.series)
  } 

  day.units <- "days since 1980-01-01 00:00:00"
  hour.units <- "hours since 1980-01-01 00:00:00"   
  time.units <- switch(freq,
                       day=day.units,
                       hour=hour.units)

  new.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                         cal=time.calendar)
  
  fixed.series <- new.origin + time.values*3600 ##86400      

  rv <- list(values=time.values,
             series=fixed.series,
             units=time.units,
             calendar=time.calendar)
}

##--------------------------------------------------------------
make.new.netcdf.file <- function(gcm,varname,freq,rcp,
                                 existing.file,
                                 tmp.base) {

  new.var <- get.new.var.name(varname)                                 

  nc <- nc_open(paste0(tmp.base,existing.file),write=FALSE)

  time.data <- correct.rcm.time.series(nc,freq) 

  yst <-  gsub('-','',format(head(time.data$series,1),'%Y-%m-%d'))
  yen <-  gsub('-','',format(tail(time.data$series,1),'%Y-%m-%d'))

  gcm.name <- switch(gcm,
                     ERAI='ERA-Interim',
                     CanESM2='CanESM2',
                     MPI='MPI')

  write.file <- paste0(new.var,'_',freq,'_WC011_',gcm.name,'+CRCM5_historical',rcp,'_',yst,'-',yen,'.nc')
 
  ##Attributes to retain
  lon <- aperm(ncvar_get(nc,'lon'),c(2,1))
  lat <- aperm(ncvar_get(nc,'lat'),c(2,1))

  rlon <- ncvar_get(nc,'rlon')
  rlat <- ncvar_get(nc,'rlat')  
     
  var.atts <- ncatt_get(nc,varname)
  
  n.lon <- length(rlon)
  n.lat <- length(rlat)

  ##--------------------------------------------------------------
  ##Create new netcdf file

  x.geog <- ncdim_def('rlon', 'degrees', rlon)
  y.geog <- ncdim_def('rlat', 'degrees', rlat)
  t.geog <- ncdim_def('time', time.data$units, time.data$values,
                      unlim=TRUE, calendar=time.data$calendar)

  lon.geog <- ncvar_def('lon', units='degrees_east', dim=list(x.geog, y.geog))
  lat.geog <- ncvar_def('lat', units='degrees_north',dim=list(x.geog, y.geog)) 
  var.geog <- ncvar_def(new.var, units=var.atts$units, dim=list(x.geog, y.geog,t.geog),
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

  nc_close(hist.nc)  
  nc_close(nc)  

  return(write.file)
}

move.data.to.new.file <- function(gcm,varname,
                                  existing.file,write.file,
                                  tmp.base) {

  new.var <- get.new.var.name(varname)                                 

  nc <- nc_open(paste0(tmp.base,write.file),write=TRUE)
  hist.nc <- nc_open(paste0(tmp.base,existing.file),write=FALSE)
  var.units <- ncatt_get(hist.nc,varname,'units')$value

  lon <- aperm(ncvar_get(hist.nc,'lon'),c(2,1))
  lat <- aperm(ncvar_get(hist.nc,'lat'),c(2,1))

  ncvar_put(nc,varid='lon',vals=lon)
  ncvar_put(nc,varid='lat',vals=lat)

  time.atts <- ncatt_get(hist.nc,'t')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(hist.nc,'t') 
  time.len <- length(time.values)
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  
  ##For Temperature (or similar daily quantities)
  if (grepl('(pr|tas|tasm|giac|giml|gld|glf|gsab|gsac|gsml|gvol|swater|sice)',new.var)) {
    var.dates <- origin.pcict + time.values*86400
    yr.fac <- as.factor(format(var.dates,'%Y'))
    yrs <- levels(yr.fac)
##browser()
    ##Split into loop to handle memory issues effectively
    l.st <- seq(1,time.len,by=400)
    l.en <- c(seq(400,time.len,by=400),time.len)
    l.cnt <- l.en-l.st+1
    for (i in seq_along(l.st)) {
      print(paste0(i,' of ',length(l.st)))
      print(l.st[i])
      var.subset <- ncvar_get(hist.nc,varname,start=c(1,1,1,l.st[i]),count=c(-1,-1,1,l.cnt[i]))
      ncvar_put(nc,varid=new.var,vals=aperm(var.subset,c(2,1,3)),
                start=c(1,1,l.st[i]),count=c(-1,-1,l.cnt[i]))     
    }
  } else { 
      ####For hourly quantities to be converted to daily (i.e. precipitation)
      ##if (grepl('(pr|snd|snc|psl|uas|vas|runoff|discharge|drainage|
      ##        streamflow|gz|insol|tas|huss|rhs|latflux|irflux|
      ##        gwst|spw|smass|smelt|asp|ps|swlake|swriver)',new.var)) {
    var.dates <- origin.pcict + time.values*3600
    mn.fac <- as.factor(format(var.dates,'%Y-%m'))
    mns <- levels(mn.fac)       
    mn.len <- length(mns)
    
    day.dates <- levels(as.factor(format(var.dates,'%Y-%m-%d')))
    
    ##Split into time-based loop to aggregate to daily
    for (i in seq_along(mns)) {
      ltm <- proc.time()
      print(paste0(i,' of ',mn.len))
      ix <- grep(mns[i],var.dates)
      l.st <- ix[1]
      l.cnt <- length(ix)
      print(var.dates[l.st])
      var.raw <- ncvar_get(hist.nc,varname,start=c(1,1,1,l.st),count=c(-1,-1,1,l.cnt))
      var.subset <- convert.variable(new.var,var.raw)

      dy.fac <- as.factor(format(var.dates[ix],'%Y-%m-%d'))

      var.agg <- apply(var.subset,c(1,2),function(x,y){tapply(x,y,mean)},dy.fac)
      dy.dates <- levels(dy.fac)
      d.st <- grep(dy.dates[1],day.dates)
      print(d.st)
      d.cnt <- length(dy.dates)
      ncvar_put(nc,varid=new.var,vals=aperm(var.agg,c(3,2,1)),
                start=c(1,1,d.st),count=c(-1,-1,d.cnt))     
      print('Loop time')
      print(proc.time()-ltm)

    }
  }
  
 

  nc_close(hist.nc)
  nc_close(nc)

}

##**************************************************************************************
##-------------------------------------------------------------------------


if (1==1) {
  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
}

##  gcm <- 'MPI'
##  varname <- 'PR'
##  interval <- '1980-2050'
##  freq <- 'hour'


##  write.file <- 'tasmax_day_WC011_ERA-Interim+CRCM5_historical+rcp85_198001-20141231.nc'  

  rcp <- "+rcp85"
  if (gcm=='ERAI') {
   rcp <- ''
  }

  tmp.base <- tmpdir
  ##tmp.base <- paste('/local_temp/ssobie/crcm5/',varname,'/',sep='')
  
  data.dir <- paste0('/storage/data/climate/downscale/RCM/CRCM5/',gcm,'/')
  write.dir <- paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/',gcm,'/')

  if (!file.exists(tmp.base))
       dir.create(tmp.base,recursive=TRUE)

  existing.file <- paste0(gcm,'_WC011_modified_',varname,'_',interval,'.nc')
  print('Copying file to temp')
  file.copy(from=paste0(data.dir,existing.file),to=tmp.base,overwrite=TRUE)

  ##-------------------------------------------------    
  ##Past File           
  print('Creating new file')                                           
  write.file <- make.new.netcdf.file(gcm,varname,freq,rcp,
                                     existing.file,
                                     tmp.base)
  print('made new file')
  print('Moving data')
  move.data.to.new.file(gcm,varname,
                        existing.file,write.file,
                        tmp.base)
  print('Copying file back')     
  file.copy(from=paste0(tmp.base,write.file),to=write.dir,overwrite=TRUE)

  file.remove(paste0(tmp.base,existing.file))
  file.remove(paste0(tmp.base,write.file))
##-------------------------------------------------------------------------


