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

make.fs.stat.netcdf.file <- function(gcm,input.file,interval,
                                     tmp.base) {
  var.name <- 'fs_stat'
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
  tmy.months <- tmy.time[grep('*-*-01 ',tmy.series)]

  gcm.name <- switch(gcm,
                     ERAI='ERA-Interim',
                     CanESM2='CanESM2',
                     MPI='MPI')
  yrs <- strsplit(interval,'-')[[1]]
  
  write.file <- paste0('MORPHED_FS_STAT_WSPD_MAX_CWEC_TMY_',gcm.name,'+CRCM5_historical_',yrs[1],'0101-',yrs[2],'1231.nc')
  
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
  var.geog <- ncvar_def(var.name, units='degC', dim=list(x.geog, y.geog,t.geog),
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
    ncatt_put(hist.nc,varid=var.name,attname=varnames[j],attval=var.atts[[j]])
  
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
  xlims <- list(c(-20,10),c(-20,10),c(-15,15),c(-10,20),
                c(-5,20), c(0,25),c(0,25),c(0,30),
                c(0,25),c(-5,20),c(-20,15),c(-25,10))
  ##fs.cdf <- max(tmy.mon) - max(mon.data) ##
  ##fs.cdf <- long.cdf(max(tmy.mon)) ##Finds the quantile of the monthly maximum relative to the 30 years
  if (1==1) {
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

compute.fs.stat <- function(months,
                            rcm.tasmax,rcm.tasmin,rcm.tas,
                            rcm.dwpt.max,rcm.dwpt.min,rcm.dwpt.avg,
                            rcm.wspd.max,rcm.wspd.avg,rcm.insol,
                            tmy.tasmax,tmy.tasmin,tmy.tas,
                            tmy.dwpt.max,tmy.dwpt.min,tmy.dwpt.avg,
                            tmy.wspd.max,tmy.wspd.avg,tmy.insol,
                            morph.tasmax,morph.tasmin,morph.tas,
                            morph.dwpt.max,morph.dwpt.min,morph.dwpt.avg,
                            morph.wspd.max,morph.wspd.avg,morph.insol,
                            mon.indices,tmy.indices,rcm.time) { 
  rv <- rep(0,length(months))
  all.vals <- matrix(0,nrow=12,ncol=10) ## vector(mode='list',length=12)
  
  par(mfrow=c(3,4))

  for (m in seq_along(months)) {
    ###print(paste0('Month: ',m))
    month <- months[m]
    mon.ix <- mon.indices[,m]
    tmy.ix <- tmy.indices[,m]
 
    tas.cdf.diffs <- fs.stat(rcm.tas,tmy.tas,morph.tas,mon.ix,tmy.ix,rcm.time,m)
    ##tasmax.cdf.diffs <- fs.stat(rcm.tasmax,tmy.tasmax,mon.ix,tmy.ix)
    ##tasmin.cdf.diffs <- fs.stat(rcm.tasmin,tmy.tasmin,mon.ix,tmy.ix)

    ##dwpt.max.cdf.diffs <- fs.stat(rcm.dwpt.max,tmy.dwpt.max,mon.ix,tmy.ix)
    ##dwpt.min.cdf.diffs <- fs.stat(rcm.dwpt.min,tmy.dwpt.min,mon.ix,tmy.ix)
    ##dwpt.avg.cdf.diffs <- fs.stat(rcm.dwpt.avg,tmy.dwpt.avg,mon.ix,tmy.ix)

    ##wspd.max.cdf.diffs <- fs.stat(rcm.wspd.max,tmy.wspd.max,mon.ix,tmy.ix)
    ##wspd.avg.cdf.diffs <- fs.stat(rcm.wspd.avg,tmy.wspd.avg,mon.ix,tmy.ix)

    ##insol.cdf.diffs <- fs.stat(rcm.insol,tmy.insol,mon.ix,tmy.ix)

    if (1==0) {
      all.vals[m,1] <- 0.3*tas.cdf.diffs
      all.vals[m,2] <- 0.05*tasmax.cdf.diffs
      all.vals[m,3] <- 0.05*tasmin.cdf.diffs
      all.vals[m,4] <- 0.025*dwpt.max.cdf.diffs
      all.vals[m,5] <- 0.025*dwpt.min.cdf.diffs
      all.vals[m,6] <- 0.05*dwpt.avg.cdf.diffs
      all.vals[m,7] <- 0.05*wspd.max.cdf.diffs
      all.vals[m,8] <- 0.05*wspd.avg.cdf.diffs
      all.vals[m,9] <- 0.4*insol.cdf.diffs


  
    fs.sum <- 0.050 * tasmax.cdf.diffs + 0.050 * tasmin.cdf.diffs + 0.300 * tas.cdf.diffs +
              0.025 * dwpt.max.cdf.diffs + 0.025 * dwpt.min.cdf.diffs + 0.050 * dwpt.avg.cdf.diffs +
              0.050 * wspd.max.cdf.diffs + 0.050 * wspd.avg.cdf.diffs +
              0.400 * insol.cdf.diffs
    }    
     rv[m] <- 0.3*tas.cdf.diffs*100 ##fs.sum ##tasmin.cdf.diffs ##fs.sum ##
  }
  return(rv)
}

##----------------------------------------------------------------------------------------------
separate.into.list <- function(nc,var.name,j,yst,cnt) {
    data.subset <- ncvar_get(nc,var.name,start=c(1,j,yst),count=c(-1,1,cnt))
    data.list <- lapply(seq_len(nrow(data.subset)), function(k) data.subset[k,])
    rm(data.subset)
    return(data.list)
}


find.fs.stats <- function(rcm.files,tmy.files,morph.file,fs.file,interval,tmp.dir) {

  months <- sprintf('%02d',1:12)
  mlen <- length(months)
  fs.nc <- nc_open(paste0(tmp.dir,fs.file),write=TRUE)
  rcm.tas.nc <- nc_open(paste0(tmp.dir,rcm.files$dewpoint.max))
  tmy.tas.nc <- nc_open(paste0(tmp.dir,tmy.files$dewpoint.max))
  morph.tas.nc <- nc_open(paste0(tmp.dir,morph.files$dewpoint.max))

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

  for (j in 1:n.lat) { 
    print(paste0('Latitude: ',j,' of ',n.lat))

    rcm.tas.list <- separate.into.list(rcm.ncs$tas,'tas',j,yst,cnt)
    rcm.tasmax.list <- separate.into.list(rcm.ncs$tasmax,'tas',j,yst,cnt)
    rcm.tasmin.list <- separate.into.list(rcm.ncs$tasmin,'tas',j,yst,cnt)

    rcm.dwpt.max.list <- separate.into.list(rcm.ncs$dewpoint.max,'dewpoint',j,yst,cnt)
    rcm.dwpt.min.list <- separate.into.list(rcm.ncs$dewpoint.min,'dewpoint',j,yst,cnt)
    rcm.dwpt.avg.list <- separate.into.list(rcm.ncs$dewpoint.avg,'dewpoint',j,yst,cnt)

    rcm.wspd.max.list <- separate.into.list(rcm.ncs$wspd.max,'wspd',j,yst,cnt)
    rcm.wspd.avg.list <- separate.into.list(rcm.ncs$wspd.avg,'wspd',j,yst,cnt)
    rcm.insol.list <- separate.into.list(rcm.ncs$insol,'insol',j,yst,cnt)

    tmy.tas.list <- separate.into.list(tmy.ncs$tas,'tas',j,1,-1)
    tmy.tasmax.list <- separate.into.list(tmy.ncs$tasmax,'tas',j,1,-1)
    tmy.tasmin.list <- separate.into.list(tmy.ncs$tasmin,'tas',j,1,-1)

    tmy.dwpt.max.list <- separate.into.list(tmy.ncs$dewpoint.max,'dewpoint',j,1,-1)
    tmy.dwpt.min.list <- separate.into.list(tmy.ncs$dewpoint.min,'dewpoint',j,1,-1)
    tmy.dwpt.avg.list <- separate.into.list(tmy.ncs$dewpoint.avg,'dewpoint',j,1,-1)

    tmy.wspd.max.list <- separate.into.list(tmy.ncs$wspd.max,'wspd',j,1,-1)
    tmy.wspd.avg.list <- separate.into.list(tmy.ncs$wspd.avg,'wspd',j,1,-1)
    tmy.insol.list <- separate.into.list(tmy.ncs$insol,'insol',j,1,-1)

    morph.tas.list <- separate.into.list(morph.ncs$tas,'tas',j,1,-1)
    morph.tasmax.list <- separate.into.list(morph.ncs$tasmax,'tas',j,1,-1)
    morph.tasmin.list <- separate.into.list(morph.ncs$tasmin,'tas',j,1,-1)

    morph.dwpt.max.list <- separate.into.list(morph.ncs$dewpoint.max,'dewpoint',j,1,-1)
    morph.dwpt.min.list <- separate.into.list(morph.ncs$dewpoint.min,'dewpoint',j,1,-1)
    morph.dwpt.avg.list <- separate.into.list(morph.ncs$dewpoint.avg,'dewpoint',j,1,-1)

    morph.wspd.max.list <- separate.into.list(morph.ncs$wspd.max,'wspd',j,1,-1)
    morph.wspd.avg.list <- separate.into.list(morph.ncs$wspd.avg,'wspd',j,1,-1)
    morph.insol.list <- separate.into.list(morph.ncs$insol,'insol',j,1,-1)


    compute.fs.stat(months=months,
                    rcm.tasmax=rcm.tasmax.list[[51]],rcm.tasmin=rcm.tasmin.list[[51]],rcm.tas=rcm.tas.list[[51]],
                    rcm.dwpt.max=rcm.dwpt.max.list[[51]],rcm.dwpt.min=rcm.dwpt.min.list[[51]],rcm.dwpt.avg=rcm.dwpt.avg.list[[51]],
                    rcm.wspd.max=rcm.wspd.max.list[[51]],rcm.wspd.avg=rcm.wspd.avg.list[[51]],rcm.insol=rcm.insol.list[[51]],
                    tmy.tasmax=tmy.tasmax.list[[51]],tmy.tasmin=tmy.tasmin.list[[51]],tmy.tas=tmy.tas.list[[51]],
                    tmy.dwpt.max=tmy.dwpt.max.list[[51]],tmy.dwpt.min=tmy.dwpt.min.list[[51]],tmy.dwpt.avg=tmy.dwpt.avg.list[[51]],
                    tmy.wspd.max=tmy.wspd.max.list[[51]],tmy.wspd.avg=tmy.wspd.avg.list[[51]],tmy.insol=tmy.insol.list[[51]],
                    morph.tasmax=morph.tasmax.list[[51]],morph.tasmin=morph.tasmin.list[[51]],morph.tas=morph.tas.list[[51]],
                    morph.dwpt.max=morph.dwpt.max.list[[51]],morph.dwpt.min=morph.dwpt.min.list[[51]],morph.dwpt.avg=morph.dwpt.avg.list[[51]],
                    morph.wspd.max=morph.wspd.max.list[[51]],morph.wspd.avg=morph.wspd.avg.list[[51]],morph.insol=morph.insol.list[[51]],
                    mon.indices=mon.indices,tmy.indices=tmy.indices,rcm.time=time)
    browser()

    fs.months <- foreach(
                  rcm.tasmax=rcm.tasmax.list,rcm.tasmin=rcm.tasmin.list,rcm.tas=rcm.tas.list,                 
                  rcm.dwpt.max=rcm.dwpt.max.list,rcm.dwpt.min=rcm.dwpt.min.list,rcm.dwpt.avg=rcm.dwpt.avg.list,                 
                  rcm.wspd.max=rcm.wspd.max.list,rcm.wspd.avg=rcm.wspd.avg.list,rcm.insol=rcm.insol.list,
                  tmy.tasmax=tmy.tasmax.list,tmy.tasmin=tmy.tasmin.list,tmy.tas=tmy.tas.list,                 
                  tmy.dwpt.max=tmy.dwpt.max.list,tmy.dwpt.min=tmy.dwpt.min.list,tmy.dwpt.avg=tmy.dwpt.avg.list,                 
                  tmy.wspd.max=tmy.wspd.max.list,tmy.wspd.avg=tmy.wspd.avg.list,tmy.insol=tmy.insol.list,
                  .export=c('compute.fs.stat','months','mon.indices','tmy.indices')
                           ) %dopar% {
                              objects <- compute.fs.stat(months=months,
                                                         rcm.tasmax=rcm.tasmax,rcm.tasmin=rcm.tasmin,rcm.tas=rcm.tas,
                                                         rcm.dwpt.max=rcm.dwpt.max,rcm.dwpt.min=rcm.dwpt.min,rcm.dwpt.avg=rcm.dwpt.avg,
                                                         rcm.wspd.max=rcm.wspd.max,rcm.wspd.avg=rcm.wspd.avg,rcm.insol=rcm.insol,
                                                         tmy.tasmax=tmy.tasmax,tmy.tasmin=tmy.tasmin,tmy.tas=tmy.tas,
                                                         tmy.dwpt.max=tmy.dwpt.max,tmy.dwpt.min=tmy.dwpt.min,tmy.dwpt.avg=tmy.dwpt.avg,
                                                         tmy.wspd.max=tmy.wspd.max,tmy.wspd.avg=tmy.wspd.avg,tmy.insol=tmy.insol,
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
  
  var.list <- c('tas_max','tas_min','tas_mean',
              'dewpoint_max','dewpoint_min','dewpoint_mean',
              'wspd_max','wspd_mean','insol_sum')
  var.names <- c('tasmax','tasmin','tas',
              'dewpoint.max','dewpoint.min','dewpoint.avg',
              'wspd.max','wspd.avg','insol')
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

if (1==0) {
  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
}

gcm <- 'MPI'
interval <- '2021-2050'

tmpdir <- '/local_temp/ssobie/crcm5/'

tmy.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/daily/'
rcm.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/'
fs.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/fs_stats/'

tmp.dir <- tmpdir

if (!file.exists(tmpdir)) {
  dir.create(tmpdir,recursive=T)
}

rcm.files <- get.file.names(gcm,rcm.dir)
tmy.files <- get.file.names(gcm,tmy.dir,morphed=FALSE,interval='20210101-20501231')
morph.files <- get.file.names(gcm,tmy.dir,morphed=TRUE,interval='20210101-20501231')

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

fs.file <- make.fs.stat.netcdf.file(gcm,
                                 input.file,interval,
                                 tmp.dir)

find.fs.stats(rcm.files,tmy.files,morph.files,fs.file,interval,tmp.dir)

file.copy(from=paste0(tmp.dir,fs.file),to=fs.dir,overwrite=T)

print('Removing RCM files')
for (rcm.file in rcm.files) {
  file.remove(paste0(tmp.dir,rcm.file))
}

print('Removing TMY files')
for (tmy.file in tmy.files) {
  file.remove(paste0(tmp.dir,tmy.file))
}






