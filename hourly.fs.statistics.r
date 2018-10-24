##Script to calculate the FS-Statistic at each grid cell for a TMY year
##Method follows the Sandia approach for finding test meteo months to form the TMY
##Variables: TMax, TMin, TMean, DewPt Max, DewPt Min, DewPt Mean, Wspd Max, Wspd Mean, Solar 
##Weights (%): 5 ,   5 ,   30 ,    2.5   ,    2.5   ,      5    ,     5   ,     5    ,  40

library(ncdf4)
library(PCICt)
library(doParallel)
registerDoParallel(cores=4)
library(foreach)
library(data.table)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/crcm5/add.crcm5.metadata.r',chdir=T)

make.fs.stat.netcdf.file <- function(gcm,var.name,mflag,
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
  tmy.series <- time.series[time.95]
  tmy.months <- tmy.time[grep('*-*-01 00:00:00',tmy.series)]

  gcm.name <- switch(gcm,
                     ERAI='ERA-Interim',
                     CanESM2='CanESM2',
                     MPI='MPI')
  yrs <- strsplit(interval,'-')[[1]]

  if (mflag) {  
     write.file <- paste0('MORPHED_HOURLY_FS_STAT_',toupper(var.name),'_CWEC_TMY_',gcm.name,'+CRCM5_historical_',yrs[1],'0101-',yrs[2],'1231.nc')
  } else {
     write.file <- paste0('HOURLY_FS_STAT_',toupper(var.name),'_CWEC_TMY_',gcm.name,'+CRCM5_historical_',yrs[1],'0101-',yrs[2],'1231.nc')
  }
  
  ##Attributes to retain
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  
  rlon <- ncvar_get(nc,'rlon')
  rlat <- ncvar_get(nc,'rlat')
  
  n.lon <- length(rlon)
  n.lat <- length(rlat)
  
  ##--------------------------------------------------------------
  ##Create new netcdf file
  
  x.geog <- ncdim_def('rlon', 'degrees', rlon)
  y.geog <- ncdim_def('rlat', 'degrees', rlat)
  t.geog <- ncdim_def('time', time.units, tmy.months,
                      unlim=FALSE, calendar=time.calendar)
  
  lon.geog <- ncvar_def('lon', units='degrees_east', dim=list(x.geog, y.geog))
  lat.geog <- ncvar_def('lat', units='degrees_north',dim=list(x.geog, y.geog))
  var.geog <- ncvar_def('fs_stat', units='degC', dim=list(x.geog, y.geog,t.geog),
                        missval=1.e+20)
  proj.geog <- ncvar_def('projection',units='',dim=list(),prec='char')
  
  hist.nc <- nc_create(paste(tmp.base,write.file,sep=''), list(lon.geog,lat.geog,proj.geog,var.geog)) ##,h_minfree=1048570)

  ncvar_put(hist.nc,'lon',lon)
  ncvar_put(hist.nc,'lat',lat)
 
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
  fs.atts <- list(standard_name = "FS_Statistic",
                  long_name = "FS_Statistic",
                  missing_value = 1.e+20,
                  cell_methods = "time: mean",
                  units = "CDF",
                  coordinates='lon lat')
  var.atts <- fs.atts
  varnames <- names(var.atts)
  for (j in 1:length(var.atts))
    ncatt_put(hist.nc,varid='fs_stat',attname=varnames[j],attval=var.atts[[j]])
  
  ##Clear extraneous history
  ncatt_put(hist.nc,varid=0,attname='history',attval='')
  
  nc_close(hist.nc)
  nc_close(nc)
  
  return(write.file)
  
}

##-------------------------------------------------------------

fs.stat <- function(rcm.data,tmy.data,morph.data,mon.ix,tmy.ix,rcm.time,mn) {
  mon.data <- rcm.data[mon.ix]
  long.cdf <- ecdf(mon.data)

  tmy.mon <- tmy.data[tmy.ix]
  short.cdf <- ecdf(tmy.mon)
  morph.mon <- morph.data[tmy.ix]
  morph.cdf <- ecdf(morph.mon)

  fs.cdf <- mean(abs(long.cdf(tmy.mon) - short.cdf(tmy.mon)))
  xlims <- list(c(-40,5),c(-40,10),c(-35,10),c(-25,10),
                c(-20,15), c(-15,20),c(-10,25),c(-10,30),
                c(-15,25),c(-25,20),c(-35,15),c(-40,10))
  ##fs.cdf <- max(tmy.mon) - max(mon.data) ##
  ##fs.cdf <- long.cdf(max(tmy.mon)) ##Finds the quantile of the monthly maximum relative to the 30 years
  if (1==0) {
    plot(sort(mon.data),long.cdf(sort(mon.data)),type='l',xlim=xlims[[mn]],lwd=4)
    for (yr in 2021:2050) {
       year.slct <- which(format(rcm.time,'%Y-%m') %in% paste0(yr,'-',sprintf('%02d',mn)))
       rcm.cdf <- ecdf(rcm.data[year.slct])
       lines(sort(rcm.data[year.slct]),sort(rcm.cdf(rcm.data[year.slct])),col='gray')
       print(paste0(yr,': ',mean(abs(long.cdf(rcm.data[year.slct]) - short.cdf(rcm.data[year.slct])))))
       ##print(cbind(tmy.mon,rcm.data[year.slct]))
       ##if (mn==7 & yr==2043)
       ##  browser()
    }  
    lines(sort(mon.data),long.cdf(sort(mon.data)),lwd=4)
    lines(sort(tmy.mon),sort(short.cdf(tmy.mon)),col='red',lwd=4)
    lines(sort(morph.mon),sort(morph.cdf(morph.mon)),col='green',lwd=4)

    print('FS')
    print(fs.cdf)
  }   
  
  return(fs.cdf)
}

###compute.fs.stat <- function(months,
###                            rcm.tas,rcm.dwpt,rcm.wspd,rcm.insol,
###                            tmy.tas,tmy.dwpt,tmy.wspd,tmy.insol,
###                            morph.tas,morph.dwpt,morph.wspd,morph.insol,
###                            mon.indices,tmy.indices,rcm.time) { 

compute.fs.stat <- function(months,
                            rcm.data,tmy.data,morph.data,
                            mon.indices,tmy.indices,rcm.time) { 

  rv <- rep(0,length(months))
  all.vals <- matrix(0,nrow=12,ncol=4) ## vector(mode='list',length=12)
  
  ###par(mfrow=c(3,4))

  for (m in seq_along(months)) {
    ###print(paste0('Month: ',m))
    month <- months[m]
    mon.ix <- mon.indices[,m]
    tmy.ix <- tmy.indices[,m]
 
    cdf.diffs <- fs.stat(rcm.data,tmy.data,morph.data,mon.ix,tmy.ix,rcm.time,m)
    ##tas.cdf.diffs <- fs.stat(rcm.tas,tmy.tas,morph.tas,mon.ix,tmy.ix,rcm.time,m)
    ##dwpt.cdf.diffs <- fs.stat(rcm.dwpt,tmy.dwpt,morph.dwpt,mon.ix,tmy.ix,rcm.time,m)
    ##wspd.cdf.diffs <- fs.stat(rcm.wspd,tmy.wspd,morph.wspd,mon.ix,tmy.ix)
    ##insol.cdf.diffs <- fs.stat(rcm.insol,tmy.insol,morph.insol,mon.ix,tmy.ix)

    if (1==0) {
      all.vals[m,1] <- tas.cdf.diffs
      all.vals[m,2] <- dwpt.cdf.diffs
      all.vals[m,3] <- wspd.cdf.diffs
      all.vals[m,4] <- insol.cdf.diffs
    }    
     rv[m] <- cdf.diffs*100 
  }
  ###browser()
  return(rv)
}

##----------------------------------------------------------------------------------------------
separate.into.list <- function(nc,var.name,j,yst,cnt) {
    data.subset <- ncvar_get(nc,var.name,start=c(1,j,yst),count=c(-1,1,cnt))
    data.list <- lapply(seq_len(nrow(data.subset)), function(k) data.subset[k,])
    rm(data.subset)
    return(data.list)
}

find.fs.stats <- function(rcm.files,tmy.files,morph.file,fs.file,var.name,interval,tmp.dir) {

  months <- sprintf('%02d',1:12)
  mlen <- length(months)
  fs.nc <- nc_open(paste0(tmp.dir,fs.file),write=TRUE)
  rcm.tas.nc <- nc_open(paste0(tmp.dir,rcm.files$dewpoint))
  tmy.tas.nc <- nc_open(paste0(tmp.dir,tmy.files$dewpoint))
  morph.tas.nc <- nc_open(paste0(tmp.dir,morph.files$dewpoint))

  rcm.ncs <- lapply(rcm.files,function(x){nc_open(paste0(tmp.dir,x))})
  tmy.ncs <- lapply(tmy.files,function(x){nc_open(paste0(tmp.dir,x))})
  morph.ncs <- lapply(morph.files,function(x){nc_open(paste0(tmp.dir,x))})
  
  rcm.time <- netcdf.calendar(rcm.tas.nc)

  bnds <- strsplit(interval,'-')[[1]]
  yst <- head(grep(bnds[1],rcm.time),1)
  yen <- tail(grep(bnds[2],rcm.time),1)
  cnt <- yen-yst+1
  time <- rcm.time[yst:yen]
  
  mon.ts <- format(time,'%Y-%m')
  mon.fac <- format(time,'%m')
  yr.fac <- format(time,'%Y')
  yrs <- levels(as.factor(format(time,'%Y')))

  n.lat <- rcm.tas.nc$dim$rlat$len ##Latitude Length
  n.lon <- rcm.tas.nc$dim$rlon$len ##Longitude Length
  n.time <- length(time)

  fs.stats.array <- array(0,c(n.lon,n.lat,mlen))

  ylen <- length(yrs)
  yr.indices <- matrix(NA,nrow=n.time,ncol=ylen)
  for (i in 1:ylen) {
    yr.indices[,i] <- yr.fac %chin% yrs[i]
  }

  mon.indices <- matrix(NA,nrow=n.time,ncol=mlen)
  for (i in 1:mlen) {
    mon.indices[,i] <- mon.fac %chin% months[i]
  }

  tmy.time <- netcdf.calendar(tmy.tas.nc)
  tmy.months <- format(tmy.time,'%m')
  tmy.indices <- matrix(NA,nrow=length(tmy.time),ncol=mlen)
  for (i in 1:mlen) {
    tmy.indices[,i] <- tmy.months %chin% months[i]
  }

  for (j in 1:n.lat) { ###Change to 35 for example
    print(paste0('Latitude: ',j,' of ',n.lat))

    rcm.tas.list <- separate.into.list(rcm.ncs$tas,'tas',j,yst,cnt)
    rcm.dwpt.list <- separate.into.list(rcm.ncs$dewpoint,'dewpoint',j,yst,cnt)
    rcm.wspd.list <- separate.into.list(rcm.ncs$wspd,'wspd',j,yst,cnt)
    rcm.insol.list <- separate.into.list(rcm.ncs$insol,'insol',j,yst,cnt)

    tmy.tas.list <- separate.into.list(tmy.ncs$tas,'tas',j,1,-1)
    tmy.dwpt.list <- separate.into.list(tmy.ncs$dewpoint,'dewpoint',j,1,-1)
    tmy.wspd.list <- separate.into.list(tmy.ncs$wspd,'wspd',j,1,-1)
    tmy.insol.list <- separate.into.list(tmy.ncs$insol,'insol',j,1,-1)

    morph.tas.list <- separate.into.list(morph.ncs$tas,'tas',j,1,-1)
    morph.dwpt.list <- separate.into.list(morph.ncs$dewpoint,'dewpoint',j,1,-1)
    morph.wspd.list <- separate.into.list(morph.ncs$wspd,'wspd',j,1,-1)
    morph.insol.list <- separate.into.list(morph.ncs$insol,'insol',j,1,-1)

    rcm.data.list <- switch(var.name,
                            tas=rcm.tas.list,dewpoint=rcm.dwpt.list,
                            wspd=rcm.wspd.list,insol=rcm.insol.list)
    tmy.data.list <- switch(var.name,
                            tas=tmy.tas.list,dewpoint=tmy.dwpt.list,
                            wspd=tmy.wspd.list,insol=tmy.insol.list)
    morph.data.list <- switch(var.name,
                              tas=morph.tas.list,dewpoint=morph.dwpt.list,
                              wspd=morph.wspd.list,insol=morph.insol.list)

    if (1==0) {
    compute.fs.stat(months=months,
                    rcm.tas=rcm.tas.list[[51]],rcm.dwpt=rcm.dwpt.list[[51]],
                    rcm.wspd=rcm.wspd.list[[51]],rcm.insol=rcm.insol.list[[51]],
                    tmy.tas=tmy.tas.list[[51]],tmy.dwpt=tmy.dwpt.list[[51]],
                    tmy.wspd=tmy.wspd.list[[51]],tmy.insol=tmy.insol.list[[51]],
                    morph.tas=morph.tas.list[[51]],morph.dwpt=morph.dwpt.list[[51]],
                    morph.insol=morph.insol.list[[51]],
                    mon.indices=mon.indices,tmy.indices=tmy.indices,rcm.time=time)
    browser()
    }
    if (1==0) {
    fs.months <- foreach(
                  rcm.tas=rcm.tas.list,rcm.dwpt=rcm.dwpt.list,                 
                  rcm.wspd=rcm.wspd.list,rcm.insol=rcm.insol.list,
                  tmy.tas=tmy.tas.list,tmy.dwpt=tmy.dwpt.list,              
                  tmy.wspd=tmy.wspd.list,tmy.insol=tmy.insol.list,
                  morph.tas=morph.tas.list,morph.dwpt=morph.dwpt.list,
                  morph.wspd=morph.wspd.list,morph.insol=morph.insol.list,
                  .export=c('compute.fs.stat','months','mon.indices','tmy.indices')
                           ) %dopar% {
                              objects <- compute.fs.stat(months=months,
                                                         rcm.tas=rcm.tas,rcm.dwpt=rcm.dwpt,
                                                         rcm.wspd=rcm.wspd,rcm.insol=rcm.insol,
                                                         tmy.tas=tmy.tas,tmy.dwpt=tmy.dwpt,
                                                         tmy.wspd=tmy.wspd,tmy.insol=tmy.insol,
                                                         morph.tas=morph.tas,morph.dwpt=morph.dwpt,
                                                         morph.wspd=morph.wspd,morph.insol=morph.insol,
                                                         mon.indices=mon.indices,tmy.indices=tmy.indices)
                           }
    }
    fs.months <- foreach(
                  rcm.data=rcm.data.list,
                  tmy.data=tmy.data.list,
                  morph.data=morph.data.list,
                  .export=c('compute.fs.stat','months','mon.indices','tmy.indices')
                           ) %dopar% {
                              objects <- compute.fs.stat(months=months,
                                                         rcm.data=rcm.data,
                                                         tmy.data=tmy.data,
                                                         morph.data=morph.data,
                                                         mon.indices=mon.indices,tmy.indices=tmy.indices)
                           }
                          
    fs.matrix <- matrix(unlist(fs.months),nrow=n.lon,ncol=length(months),byrow=T)

    ncvar_put(fs.nc,varid='fs_stat',vals=fs.matrix,
                      start=c(1,j,1),count=c(-1,1,mlen))

  }
  nc_close(rcm.tas.nc)
  nc_close(fs.nc)
  lapply(rcm.ncs,nc_close)
}

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------


get.file.names <- function(gcm,read.dir,morphed=FALSE,interval=NULL) {
  
  var.list <- c('tas','dewpoint','wspd','insol')
  var.names <- c('tas','dewpoint','wspd','insol')
  file.names <- vector(mode='list',length=length(var.list))
  names(file.names) <- var.names
  rcm.files <- list.files(path=read.dir,pattern=gcm)

  if (morphed) {
    var.list <- paste0('morphed_',var.list)
  } else {
    var.list <- paste0('^',var.list)
  }     

  if (!is.null(interval)) {
    rcm.files <- rcm.files[grep(interval,rcm.files)]
  }

  for (v in seq_along(var.list)) {
    file.names[[var.names[v]]] <- rcm.files[grep(var.list[v],rcm.files)]
  }  
  return(file.names)    
}



##----------------------------------------------------------------------------------------------

if (1==1) {
  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
}

##gcm <- 'MPI'
##interval <- '2021-2050'

var.name <- varname
print(var.name)
print(flag)

tmpdir <- '/local_temp/ssobie/crcm5/'

tmy.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/'
rcm.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/hourly/'
fs.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/fs_stats/'

tmp.dir <- tmpdir

if (!file.exists(tmpdir)) {
  dir.create(tmpdir,recursive=T)
}

bnds <- strsplit(interval,'-')[[1]]
print(bnds)

if (flag=='TRUE')
   mflag <- TRUE
if (flag=='FALSE')
   mflag <- FALSE
                  


rcm.files <- get.file.names(gcm,rcm.dir)
print('RCM Files')
print(rcm.files)
tmy.files <- get.file.names(gcm,tmy.dir,morphed=mflag,interval=paste0(bnds[1],'0101-',bnds[2],'1231'))
print('TMY Files')
print(tmy.files)
morph.files <- get.file.names(gcm,tmy.dir,morphed=mflag,interval=paste0(bnds[1],'0101-',bnds[2],'1231'))
print('MORPH Files')
print(morph.files)

print('Copying RCM files')
for (rcm.file in rcm.files) {
  file.copy(from=paste0(rcm.dir,rcm.file),to=tmp.dir,overwrite=T)
}

print('Copying TMY files')
for (tmy.file in tmy.files) {
  file.copy(from=paste0(tmy.dir,tmy.file),to=tmp.dir,overwrite=T)
}

print('Copying MORPH files')
for (tmy.file in morph.files) {
  file.copy(from=paste0(tmy.dir,morph.files),to=tmp.dir,overwrite=T)
}

input.file <- rcm.files$tas

fs.file <- make.fs.stat.netcdf.file(gcm,var.name,mflag,
                                 input.file,interval,
                                 tmp.dir)
##fs.file <- 'TEST.nc'
find.fs.stats(rcm.files,tmy.files,morph.files,fs.file,var.name,interval,tmp.dir)

file.copy(from=paste0(tmp.dir,fs.file),to=fs.dir,overwrite=T)

print('Removing RCM files')
for (rcm.file in rcm.files) {
  file.remove(paste0(tmp.dir,rcm.file))
}

print('Removing TMY files')
for (tmy.file in tmy.files) {
  file.remove(paste0(tmp.dir,tmy.file))
}

print('Removing MORPH files')
for (morph.file in morph.files) {
  file.remove(paste0(tmp.dir,morph.file))
}