##Script to find the test meteorological year at each grid cell for a provided netcdf file
##Method follows the Sandia approach for finding test meteo months to form the TMY
##Variables: TMax, TMin, TMean, DewPt Max, DewPt Min, DewPt Mean, Wspd Max, Wspd Mean, Solar 
##Weights (%): 5 ,   5 ,   30 ,    2.5   ,    2.5   ,      5    ,     5   ,     5    ,  40

Rprof('wx.out')

library(ncdf4)
library(PCICt)
library(doParallel)
registerDoParallel(cores=4)
library(foreach)
library(data.table)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

fs.stat <- function(data,mon.ix,yrs,mon.vl,mon.ts,month) {
    mon.data <- data[mon.ix]
    long.cdf <- ecdf(mon.data)

    cdf.diffs <- rep(0,length(yrs))
    for (i in seq_along(yrs)) {
      yr.mon.ix <- mon.vl[mon.ts[mon.ix] %chin% paste0(yrs[i],'-',month)]
      yr.mon <- data[yr.mon.ix]
      short.cdf <- ecdf(yr.mon)
      cdf.diffs[i] <- mean(abs(long.cdf(yr.mon) - short.cdf(yr.mon)))
    }
  return(cdf.diffs)
}


find.all.months <- function(months,tas,insol,rhs,mon.ts,mon.indices,yrs) { 
  rv <- matrix(0,nrow=length(months),ncol=3)
  for (m in seq_along(months)) {
    month <- months[m]
    mon.ix <- mon.indices[,m]
    mon.vl <- which(mon.ix)  
    tas.mon.data <- tas[mon.ix]
    insol.mon.data <- insol[mon.ix]
    rhs.mon.data <- rhs[mon.ix]

    tas.cdf.diffs <- fs.stat(tas,mon.ix,yrs,mon.vl,mon.ts,month)
    insol.cdf.diffs <- fs.stat(insol,mon.ix,yrs,mon.vl,mon.ts,month)
    rhs.cdf.diffs <- fs.stat(rhs,mon.ix,yrs,mon.vl,mon.ts,month)

    fs.sum <- tas.cdf.diffs + insol.cdf.diffs + rhs.cdf.diffs
    rv[m,] <- (yrs[order(fs.sum)])[1:3] ##yrs[which.min(fs.sum)]
  }
  return(rv)
}

##----------------------------------------------------------------------------------------------

find.test.ref.years <- function(tas.file,insol.file,rhs.file) {

  months <- sprintf('%02d',1:12)
  mlen <- length(months)

  tas.nc <- nc_open(tas.file)
  insol.nc <- nc_open(insol.file)
  rhs.nc <- nc_open(rhs.file)

  time <- netcdf.calendar(tas.nc)
  mon.ts <- format(time,'%Y-%m')
  mon.fac <- format(time,'%m')
  yrs <- levels(as.factor(format(time,'%Y')))
  n.lat <- tas.nc$dim$rlat$len ##Latitude Length
  n.lon <- tas.nc$dim$rlon$len ##Longitude Length
  n.time <- tas.nc$dim$time$len ##Longitude Length

  year.one.array <- array(0,c(n.lon,n.lat,mlen))
  year.two.array <- array(0,c(n.lon,n.lat,mlen))
  year.three.array <- array(0,c(n.lon,n.lat,mlen))

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

    tas.subset <- ncvar_get(tas.nc,'tas',start=c(1,j,1),count=c(-1,1,-1))
    tas.list <- lapply(seq_len(nrow(tas.subset)), function(k) tas.subset[k,])
    rm(tas.subset)
    insol.subset <- ncvar_get(insol.nc,'insol',start=c(1,j,1),count=c(-1,1,-1))
    insol.list <- lapply(seq_len(nrow(insol.subset)), function(k) insol.subset[k,])
    rm(insol.subset)
    rhs.subset <- ncvar_get(rhs.nc,'rhs',start=c(1,j,1),count=c(-1,1,-1))
    rhs.list <- lapply(seq_len(nrow(rhs.subset)), function(k) rhs.subset[k,])
    rm(rhs.subset)

    years.slct <- foreach(
                  tas=tas.list,                 
                  insol=insol.list,
                  rhs=rhs.list,
                  .export=c('find.all.months','months','mon.ts','mon.indices','yrs')
                           ) %dopar% {
                              objects <- find.all.months(months=months,tas=tas,insol=insol,rhs=rhs,mon.ts=mon.ts,mon.indices=mon.indices,yrs=yrs)
                         }
    year.one.list <- lapply(years.slct,function(x){return(x[,1])}) 
    year.one.matrix <- matrix(unlist(year.one.list),nrow=n.lon,ncol=length(months),byrow=T)
    year.one.array[,j,] <- year.one.matrix  

    year.two.list <- lapply(years.slct,function(x){return(x[,2])}) 
    year.two.matrix <- matrix(unlist(year.two.list),nrow=n.lon,ncol=length(months),byrow=T)
    year.two.array[,j,] <- year.two.matrix  

    year.three.list <- lapply(years.slct,function(x){return(x[,3])}) 
    year.three.matrix <- matrix(unlist(year.three.list),nrow=n.lon,ncol=length(months),byrow=T)
    year.three.array[,j,] <- year.three.matrix  
  }
  nc_close(tas.nc)
  years.array <- list(one=year.one.array,two=year.two.array,three=year.three.array)
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

read.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/daily/'
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
print('Copying INSOL file')
file.copy(from=paste0(read.dir,insol.file),to=tmp.dir,overwrite=T)
print('Copying RHS file')
file.copy(from=paste0(read.dir,rhs.file),to=tmp.dir,overwrite=T)

primary.years   <- find.test.ref.years(paste0(tmp.dir,tas.file),paste0(tmp.dir,insol.file),paste0(tmp.dir,rhs.file))
years.one <- primary.years$one
save(years.one,file='/storage/data/climate/downscale/RCM/CRCM5/reconfig/daily/year_data/years.one.from.daily.RData')
years.two <- primary.years$two
save(years.two,file='/storage/data/climate/downscale/RCM/CRCM5/reconfig/daily/year_data/years.two.from.daily.RData')
years.three <- primary.years$three
save(years.three,file='/storage/data/climate/downscale/RCM/CRCM5/reconfig/daily/year_data/years.three.from.daily.RData')

file.remove(paste0(tmp.dir,tas.file))
print('Done with TAS')
file.remove(paste0(tmp.dir,insol.file))
print('Done with INSOL')
file.remove(paste0(tmp.dir,rhs.file))
print('Done with RHS')

}



##Find Average Wind Months
print('Copying Wind files')
uas.file <- 'uas_day_WC011_ERA-Interim+CRCM5_historical_19800101-20141231.nc'
vas.file <- 'vas_day_WC011_ERA-Interim+CRCM5_historical_19800101-20141231.nc'
file.copy(from=paste0(read.dir,uas.file),to=tmp.dir,overwrite=T)
file.copy(from=paste0(read.dir,vas.file),to=tmp.dir,overwrite=T)

try.years   <- wind.trys(years.one,years.two,years.three,
                         paste0(tmp.dir,uas.file),paste0(tmp.dir,vas.file))
save(try.years,file='/storage/data/climate/downscale/RCM/CRCM5/reconfig/daily/year_data/try.years.from.daily.RData')
file.remove(paste0(tmp.dir,uas.file))
file.remove(paste0(tmp.dir,vas.file))

print('Done with Winds')



Rprof(NULL)