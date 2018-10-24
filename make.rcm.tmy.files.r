##Script to calculate dewpoint temperatures

library(ncdf4)
library(PCICt)
library(doParallel)
registerDoParallel(cores=4)
library(foreach)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/crcm5/add.crcm5.metadata.r',chdir=T)

##----------------------------------------------------------------------------------------------

find.tmy.months <- function(input,tmys,input.time,tmy.time,gcm) {
   month.series <- rep(0,length(tmy.time))   

   for (mn in 1:12) {                
     mn.ix <- format(tmy.time,'%m') %in% sprintf('%02d',mn)                  
   
     ##-----------------------------------------
     ##Testing
     if (1==0) {
       date.sub <- 119809:207464 ##2929:90584 ##
       input.future <- input[date.sub]
       mon.ix <- format(input.time[date.sub],'%m') %in% sprintf('%02d',mn)                  
       long.cdf <- ecdf(input.future[mon.ix])
       plot(input.future[mon.ix],long.cdf(input.future[mon.ix]),xlim=xlims[[mn]],pch='*')

       for (yr in 2021:2050) {
          year.slct <- which(format(input.time[date.sub],'%Y-%m') %in% paste0(yr,'-',sprintf('%02d',mn))) 
          short.cdf <- ecdf(input.future[year.slct])
          points(input.future[year.slct],short.cdf(input.future[year.slct]),col='blue',pch='*')
       }
       points(input.future[mon.ix],long.cdf(input.future[mon.ix]),col='green',pch='*')
     }
     ##-----------------------------------------

     mn.yr <- tmys[mn]
     yr.slct <- which(format(input.time,'%Y-%m') %in% paste0(mn.yr,'-',sprintf('%02d',mn))) 
     if (any(grepl('-02-29',input.time[yr.slct]))) {
       leap.ix <- grep('02-29',input.time[yr.slct])
       ##print('Leap Flag')
       ##print(input.time[yr.slct][leap.ix])
       if (gcm=='MPI') {
         yr.slct <- yr.slct[-leap.ix] ## For MPI
         ##yr.slct[leap.ix] <- yr.slct[leap.ix]-24 ##For MPI Isol
         ##print(input.time[yr.slct])
         #browser()
       }
       if (gcm == 'CanESM2') {
         yr.slct[leap.ix] <- yr.slct[leap.ix]-24 ##For CanESM2
         ##browser()
       }
     }
     if (sum(mn.ix) != length(yr.slct)) {
       print('Size mismatch')
       print(input.time[yr.slct])
       print('Wrong Model Settings')       
       ##browser()
     }     
     month.series[mn.ix] <- input[yr.slct]
     ##if (mn==2)
     ##  browser()
     ##slct.cdf <- ecdf(input[yr.slct]) 
     ##points(input[yr.slct],slct.cdf(input[yr.slct]),col='red',pch='*')
     ##fs.cdf <- mean(abs(long.cdf(input[yr.slct]) - slct.cdf(input[yr.slct])))
   }   
   ##browser()
   return(month.series)
}

make.tmy.netcdf.file <- function(gcm,var.name,
                                 input.file,interval,
                                 tmp.base) {

  nc <- nc_open(paste0(tmp.base,input.file),write=FALSE)
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.series <- origin.pcict + time.values*3600

  ##Use 1995 as the example year
  time.95 <- format(time.series,'%Y') %in% '1995'
  tmy.time <- time.values[time.95]

  gcm.name <- switch(gcm,
                     ERAI='ERA-Interim',
                     CanESM2='CanESM2',
                     MPI='MPI')
  yrs <- strsplit(interval,'-')[[1]]
  
  write.file <- paste0(var.name,'_CWEC_TMY_',gcm.name,'+CRCM5_historical_',yrs[1],'0101-',yrs[2],'1231.nc')
  
  ##Attributes to retain
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')

  rlon <- ncvar_get(nc,'rlon')
  rlat <- ncvar_get(nc,'rlat')

  var.atts <- ncatt_get(nc,var.name)

  n.lon <- length(rlon)
  n.lat <- length(rlat)

  ##--------------------------------------------------------------
  ##Create new netcdf file

  x.geog <- ncdim_def('rlon', 'degrees', rlon)
  y.geog <- ncdim_def('rlat', 'degrees', rlat)
  t.geog <- ncdim_def('time', time.units, tmy.time,
                      unlim=FALSE, calendar=time.calendar)

  lon.geog <- ncvar_def('lon', units='degrees_east', dim=list(x.geog, y.geog))
  lat.geog <- ncvar_def('lat', units='degrees_north',dim=list(x.geog, y.geog))
  var.geog <- ncvar_def(var.name, units=var.atts$units, dim=list(x.geog, y.geog,t.geog),
                        missval=1.e+20)
  proj.geog <- ncvar_def('projection',units='',dim=list(),prec='char')

  hist.nc <- nc_create(paste(tmp.base,write.file,sep=''), list(lon.geog,lat.geog,proj.geog,var.geog)) ##,h_minfree=1048570)

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
  global.atts <- get.global.atts(gcm,freq='hour',rcp='rcp85')
  global.names <- names(global.atts)
  for (g in 1:length(global.atts))
    ncatt_put(hist.nc,varid=0,attname=global.names[g],attval=global.atts[[g]])
    print('Variable Atts')
  var.atts <- get.variable.atts(var.name)
  varnames <- names(var.atts)
  for (j in 1:length(var.atts))
    ncatt_put(hist.nc,varid=var.name,attname=varnames[j],attval=var.atts[[j]])

  ##Clear extraneous history
  ncatt_put(hist.nc,varid=0,attname='history',attval='')

  nc_close(hist.nc)
  nc_close(nc)

  return(write.file)

}



##-----------------------------------------------------------------------------------------
add.tmy.values <- function(tmy.file,input.file,var.name,tmy.years,tmp.base,gcm) {

  input.nc <- nc_open(paste0(tmp.base,input.file))
  tmy.nc <- nc_open(paste0(tmp.base,tmy.file),write=TRUE)
  
  input.time <- netcdf.calendar(input.nc)
  tmy.time <- netcdf.calendar(tmy.nc)
  yrs <- levels(as.factor(format(time,'%Y')))
  n.lat <- tmy.nc$dim$rlat$len ##Latitude Length
  n.lon <- tmy.nc$dim$rlon$len ##Longitude Length
  n.time <- tmy.nc$dim$time$len ##Time Length

  lon <- ncvar_get(input.nc,'lon')
  lat <- ncvar_get(input.nc,'lat')

  ncvar_put(tmy.nc,varid='lon',vals=lon)
  ncvar_put(tmy.nc,varid='lat',vals=lat)

  for (j in 1:n.lat) { 
    print(paste0('Latitude: ',j,' of ',n.lat))
    input.subset <- ncvar_get(input.nc,var.name,start=c(1,j,1),count=c(-1,1,-1))
    input.list <- lapply(seq_len(nrow(input.subset)), function(k) input.subset[k,])
    tmy.subset <- tmy.years[,j,]
    tmy.list <- lapply(seq_len(nrow(tmy.subset)), function(k) tmy.subset[k,])
    rm(input.subset)

    ##find.tmy.months(input=input.list[[51]],tmys=tmy.list[[51]],input.time=input.time,tmy.time=tmy.time,gcm=gcm)

    tmy.output <- foreach(
                          input=input.list,                 
                          tmys=tmy.list,                 
                          .export=c('find.tmy.months','input.time','tmy.time','gcm')
                             ) %dopar% {
                                objects <- find.tmy.months(input=input,tmys=tmys,input.time=input.time,tmy.time=tmy.time,gcm=gcm)
                             }
    tmy.matrix <- matrix(unlist(tmy.output),nrow=n.lon,ncol=n.time,byrow=T)
    ##if (j==5)
    ##  browser()


    ##browser()
    ncvar_put(tmy.nc,varid=var.name,vals=tmy.matrix,
                      start=c(1,j,1),count=c(-1,1,n.time))
  }
  nc_close(tmy.nc)
  nc_close(input.nc)
}




##----------------------------------------------------------------------------------------------

if (1==1) {
  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
}

gcm <- 'CanESM2'
varname <- 'dewpoint'
##interval <- '1981-2010'
##tmpdir <- '/local_temp/ssobie/crcm5/'

var.name <- varname
print(varname)
print(gcm)
print(interval)

read.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/hourly/'
write.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/'
tmp.dir <- tmpdir
print(tmp.dir)
print(tmpdir)
if (!file.exists(tmpdir)) {
  dir.create(tmpdir,recursive=T)
}

load(paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/year_data/',
            gcm,'.cwec.tmy.years.from.daily.',interval,'.RData'))

gcm.files <- list.files(path=read.dir,pattern=gcm)
input.file <- gcm.files[grep(var.name,gcm.files)]

file.copy(from=paste0(read.dir,input.file),to=tmp.dir,overwrite=T)

tmy.file <- make.tmy.netcdf.file(gcm,var.name,
                                 input.file,interval,
                                 tmp.dir) 
add.tmy.values(tmy.file,input.file,var.name,tmy.years,tmp.dir,gcm)
file.copy(from=paste0(tmp.dir,tmy.file),to=write.dir,overwrite=T)






