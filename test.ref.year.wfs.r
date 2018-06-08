##Script to find the test reference year at each grid cell for a provided netcdf file

Rprof('wx.out')

library(ncdf4)
library(PCICt)
library(doParallel)
registerDoParallel(cores=4)
library(foreach)
library(data.table)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

find.all.months <- function(months,data,mon.ts,mon.indices,yrs) { 
  rv <- rep(0,length(months))
  for (m in seq_along(months)) {
    month <- months[m]
    mon.ix <- mon.indices[,m]
    mon.vl <- which(mon.ix)  
    mon.data <- data[mon.ix]
    long.cdf <- ecdf(mon.data)

    cdf.diffs <- rep(0,length(yrs))
    for (i in seq_along(yrs)) {
      yr.mon.ix <- mon.vl[mon.ts[mon.ix] %chin% paste0(yrs[i],'-',month)]
      yr.mon <- data[yr.mon.ix]
      short.cdf <- ecdf(yr.mon)
      cdf.diffs[i] <- mean(abs(long.cdf(yr.mon) - short.cdf(yr.mon)))
    }
    rv[m] <- yrs[which.min(cdf.diffs)]
  }
  return(rv)
}

##----------------------------------------------------------------------------------------------

find.test.ref.years <- function(var.name,var.file) {

  months <- sprintf('%02d',1:12)
  mlen <- length(months)

  input.nc <- nc_open(var.file)
  time <- netcdf.calendar(input.nc)
  mon.ts <- format(time,'%Y-%m')
  mon.fac <- format(time,'%m')
  yrs <- levels(as.factor(format(time,'%Y')))
  n.lat <- input.nc$dim$rlat$len ##Latitude Length
  n.lon <- input.nc$dim$rlon$len ##Longitude Length
  n.time <- input.nc$dim$time$len ##Longitude Length
  years.array <- array(0,c(n.lon,n.lat,mlen))

  mon.indices <- matrix(NA,nrow=n.time,ncol=mlen)
  for (i in 1:mlen) {
    mon.indices[,i] <- mon.fac %chin% months[i]
  }

  ymlen <- length(unique(mon.ts))
  yr.mons <- unique(mon.ts)
  yr.mon.indices <- matrix(NA,nrow=n.time,ncol=ymlen)
  for (i in 1:ymlen) {
    yr.mon.indices[,i] <- mon.ts %chin% yr.mons[i]  
  }

  for (j in 1:n.lat) { 
    print(paste0('Latitude: ',j,' of ',n.lat))
    input.subset <- ncvar_get(input.nc,var.name,start=c(1,j,1),count=c(-1,1,-1))
    input.list <- lapply(seq_len(nrow(input.subset)), function(k) input.subset[k,])
    rm(input.subset)

    years.slct <- foreach(
                  input=input.list,                 
                  .export=c('find.all.months','months','mon.ts','mon.indices','yrs')
                           ) %dopar% {
                              objects <- find.all.months(months=months,data=input,mon.ts=mon.ts,mon.indices=mon.indices,yrs=yrs)
                         }
    years.matrix <- matrix(unlist(years.slct),nrow=n.lon,ncol=length(months),byrow=T)
    years.array[,j,] <- years.matrix  
  }
  nc_close(input.nc)
  return(years.array)
}

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

find.wind.months <- function(wspd,tyrs,iyrs,ryrs,mon.ts,mon.indices,months) {

  rv <- rep(0,length(months))
  for (m in seq_along(months)) {
    month <- months[m]
    mon.ix <- mon.indices[,m]
    mon.vl <- which(mon.ix)  
    mon.wspd <- wspd[mon.ix]
    long.cdf <- ecdf(mon.wspd)
    
    slct.yrs <- unique(c(tyrs[m],iyrs[m],ryrs[m]))
    if (length(slct.yrs) == 1) {
      rv[m] <- slct.yrs
    } else {
      cdf.diffs <- rep(0,length(slct.yrs))
      for (i in seq_along(slct.yrs)) {
        yr.mon.ix <- mon.vl[mon.ts[mon.ix] %chin% paste0(slct.yrs[i],'-',month)]
        yr.mon <- wspd[yr.mon.ix]
        short.cdf <- ecdf(yr.mon)
        cdf.diffs[i] <- mean(abs(long.cdf(yr.mon) - short.cdf(yr.mon)))
      }
      rv[m] <- slct.yrs[which.min(cdf.diffs)]
    }
  }
  return(rv)

}

##Use the selected years from TAS,INSOL,RHS and find the most average wind years

wind.trys <- function(tas.years,insol.years,rhs.years,uas.file,vas.file) {

  months <- sprintf('%02d',1:12)
  mlen <- length(months) 

  uas.nc <- nc_open(uas.file)
  vas.nc <- nc_open(vas.file)
  time <- netcdf.calendar(uas.nc)
  mon.ts <- format(time,'%Y-%m')
  mon.fac <- format(time,'%m')
  yrs <- unique(format(time,'%Y'))
  years <- format(time,'%Y')
  n.lat <- uas.nc$dim$rlat$len ##Latitude Length
  n.lon <- uas.nc$dim$rlon$len ##Longitude Length
  n.time <- uas.nc$dim$time$len ##Longitude Length
  years.array <- array(0,c(n.lon,n.lat,mlen))

  mon.indices <- matrix(NA,nrow=n.time,ncol=mlen)
  for (i in 1:mlen) {
    mon.indices[,i] <- mon.fac %chin% months[i]
  }

  ymlen <- length(unique(mon.ts))
  yr.mons <- unique(mon.ts)
  yr.mon.indices <- matrix(NA,nrow=n.time,ncol=ymlen)
  for (i in 1:ymlen) {
    yr.mon.indices[,i] <- mon.ts %chin% yr.mons[i]  
  }

  for (j in 1:n.lat) { 
    print(paste0('Latitude: ',j,' of ',n.lat))
    uas.subset <- ncvar_get(uas.nc,'uas',start=c(1,j,1),count=c(-1,1,-1))
    vas.subset <- ncvar_get(vas.nc,'vas',start=c(1,j,1),count=c(-1,1,-1))
    wspd.subset <- sqrt(uas.subset^2 + vas.subset^2)
    wspd.list <- lapply(seq_len(nrow(wspd.subset)), function(k) wspd.subset[k,])
    rm(uas.subset)
    rm(vas.subset)

    tas.yr.sub <- tas.years[,j,]
    tas.yr.list <- lapply(seq_len(nrow(tas.yr.sub)), function(k) tas.yr.sub[k,])
    insol.yr.sub <- insol.years[,j,]
    insol.yr.list <- lapply(seq_len(nrow(insol.yr.sub)), function(k) insol.yr.sub[k,])
    rhs.yr.sub <- rhs.years[,j,]
    rhs.yr.list <- lapply(seq_len(nrow(rhs.yr.sub)), function(k) rhs.yr.sub[k,])

    years.slct <- foreach(
                    wspd=wspd.list,                 
                    tyrs=tas.yr.list,                 
                    iyrs=insol.yr.list,                 
                    ryrs=rhs.yr.list,                 
                    .export=c('find.wind.months','mon.ts','mon.indices','months')
                             ) %do% {
                                objects <- find.wind.months(wspd=wspd,tyrs=tyrs,iyrs=iyrs,ryrs=ryrs,mon.ts=mon.ts,mon.indices=mon.indices,months)
                                    }
    years.matrix <- matrix(unlist(years.slct),nrow=n.lon,ncol=length(months),byrow=T)
    years.array[,j,] <- years.matrix  
  }
  nc_close(uas.nc)
  nc_close(vas.nc)
  return(years.array)
}




##----------------------------------------------------------------------------------------------

if (1==1) {
  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
}

##gcm <- 'ERAI'

read.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/ERA-Interim/'
##tmpdir <- '/local_temp/ssobie/crcm5/'
tmp.dir <- tmpdir

if (!file.exists(tmpdir)) {
  dir.create(tmpdir,recursive=T)
}

tas.file <- 'tas_day_WC011_ERA-Interim+CRCM5_historical_19800101-20141231.nc'
insol.file <- 'insol_day_WC011_ERA-Interim+CRCM5_historical_19800101-20141231.nc'
rhs.file <- 'rhs_day_WC011_ERA-Interim+CRCM5_historical_19800101-20141231.nc'

if (1==1) {
print('Copying TAS file')
file.copy(from=paste0(read.dir,tas.file),to=tmp.dir,overwrite=T)
tas.years   <- find.test.ref.years('tas',paste0(tmp.dir,tas.file))
save(tas.years,file='/storage/data/climate/downscale/RCM/CRCM5/reconfig/daily/year_data/tas.years.RData')
file.remove(paste0(tmp.dir,tas.file))
print('Done with TAS')

print('Copying INSOL file')
file.copy(from=paste0(read.dir,insol.file),to=tmp.dir,overwrite=T)
insol.years <- find.test.ref.years('insol',paste0(tmp.dir,insol.file))
save(insol.years,file='/storage/data/climate/downscale/RCM/CRCM5/reconfig/daily/year_data/insol.years.RData')
file.remove(paste0(tmp.dir,insol.file))
print('Done with INSOL')

print('Copying RHS file')
file.copy(from=paste0(read.dir,rhs.file),to=tmp.dir,overwrite=T)
rhs.years   <- find.test.ref.years('rhs',paste0(tmp.dir,rhs.file))
save(rhs.years,file='/storage/data/climate/downscale/RCM/CRCM5/reconfig/daily/year_data/rhs.years.RData')
file.remove(paste0(tmp.dir,rhs.file))
print('Done with RHS')
}



##Find Average Wind Months
print('Copying Wind files')
uas.file <- 'uas_day_WC011_ERA-Interim+CRCM5_historical_19800101-20141231.nc'
vas.file <- 'vas_day_WC011_ERA-Interim+CRCM5_historical_19800101-20141231.nc'
file.copy(from=paste0(read.dir,uas.file),to=tmp.dir,overwrite=T)
file.copy(from=paste0(read.dir,vas.file),to=tmp.dir,overwrite=T)

try.years   <- wind.trys(tas.years,insol.years,rhs.years,
                         paste0(tmp.dir,uas.file),paste0(tmp.dir,vas.file))
save(try.years,file='/storage/data/climate/downscale/RCM/CRCM5/reconfig/daily/year_data/try.years.RData')
file.remove(paste0(tmp.dir,uas.file))
file.remove(paste0(tmp.dir,vas.file))

print('Done with Winds')



Rprof(NULL)