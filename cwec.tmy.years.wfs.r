##Script to find the test meteorological year at each grid cell for a provided netcdf file
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

find.all.months <- function(months,tasmax,tasmin,tas,dwpt.max,dwpt.min,dwpt.avg,
                                   wspd.max,wspd.avg,insol,
                                   mon.ts,mon.indices,yrs) { 
  rv <- matrix(0,nrow=length(months),ncol=5)
  fs.all <- matrix(0,nrow=30,ncol=9)
  for (m in seq_along(months)) {
  m <- 7
    month <- months[m]
    mon.ix <- mon.indices[,m]
    mon.vl <- which(mon.ix)  

    tasmax.cdf.diffs <- fs.stat(tasmax,mon.ix,yrs,mon.vl,mon.ts,month)
    tasmin.cdf.diffs <- fs.stat(tasmin,mon.ix,yrs,mon.vl,mon.ts,month)
    tas.cdf.diffs <- fs.stat(tas,mon.ix,yrs,mon.vl,mon.ts,month)
    dwpt.max.cdf.diffs <- fs.stat(dwpt.max,mon.ix,yrs,mon.vl,mon.ts,month)
    dwpt.min.cdf.diffs <- fs.stat(dwpt.min,mon.ix,yrs,mon.vl,mon.ts,month)
    dwpt.avg.cdf.diffs <- fs.stat(dwpt.avg,mon.ix,yrs,mon.vl,mon.ts,month)
    wspd.max.cdf.diffs <- fs.stat(wspd.max,mon.ix,yrs,mon.vl,mon.ts,month)
    wspd.avg.cdf.diffs <- fs.stat(wspd.avg,mon.ix,yrs,mon.vl,mon.ts,month)
    insol.cdf.diffs <- fs.stat(insol,mon.ix,yrs,mon.vl,mon.ts,month)

    fs.sum <- 0.050 * tasmax.cdf.diffs + 0.050 * tasmin.cdf.diffs + 0.300 * tas.cdf.diffs +
              0.025 * dwpt.max.cdf.diffs + 0.025 * dwpt.min.cdf.diffs + 0.050 * dwpt.avg.cdf.diffs +
              0.050 * wspd.max.cdf.diffs + 0.050 * wspd.avg.cdf.diffs +
              0.400 * insol.cdf.diffs
    fs.all[,1] <- 0.3*tas.cdf.diffs
    fs.all[,2] <- 0.05*tasmax.cdf.diffs
    fs.all[,3] <- 0.05*tasmin.cdf.diffs
    fs.all[,4] <- 0.025*dwpt.max.cdf.diffs
    fs.all[,5] <- 0.025*dwpt.min.cdf.diffs
    fs.all[,6] <- 0.05*dwpt.avg.cdf.diffs
    fs.all[,7] <- 0.05*wspd.max.cdf.diffs
    fs.all[,8] <- 0.05*wspd.avg.cdf.diffs
    fs.all[,9] <- 0.4*insol.cdf.diffs

    a <- barplot(t(fs.all),legend=c('TAS','TX','TN','DX','DN','DA','WX','WA','IN'),cex.axis=1.5,
                 args.legend=list(x='bottomleft'),col=c('red','orange','yellow','green','blue','purple','black','gray','brown'))
    axis(1,at=a,label=2021:2050)
    box(which='plot')


    rv[m,] <- (yrs[order(fs.sum)])[1:5] ##yrs[which.min(fs.sum)]
    browser()
  }
  return(rv)
}

##----------------------------------------------------------------------------------------------

find.test.ref.years <- function(input.files,interval,tmp.dir) {

  months <- sprintf('%02d',1:12)
  mlen <- length(months)

  tasmax.nc <- nc_open(paste0(tmp.dir,input.files$tasmax))
  tasmin.nc <- nc_open(paste0(tmp.dir,input.files$tasmin))
  tas.nc <- nc_open(paste0(tmp.dir,input.files$tas))

  dwpt.max.nc <- nc_open(paste0(tmp.dir,input.files$dewpoint.max))
  dwpt.min.nc <- nc_open(paste0(tmp.dir,input.files$dewpoint.min))
  dwpt.avg.nc <- nc_open(paste0(tmp.dir,input.files$dewpoint.avg))

  wspd.max.nc <- nc_open(paste0(tmp.dir,input.files$wspd.max))
  wspd.avg.nc <- nc_open(paste0(tmp.dir,input.files$wspd.avg))

  insol.nc <- nc_open(paste0(tmp.dir,input.files$insol))

  raw.time <- netcdf.calendar(tas.nc)
  bnds <- strsplit(interval,'-')[[1]]
  yst <- head(grep(bnds[1],raw.time),1)
  yen <- tail(grep(bnds[2],raw.time),1)
  cnt <- yen-yst+1
  time <- raw.time[yst:yen]

  mon.ts <- format(time,'%Y-%m')
  mon.fac <- format(time,'%m')
  yrs <- levels(as.factor(format(time,'%Y')))

  n.lat <- tas.nc$dim$rlat$len ##Latitude Length
  n.lon <- tas.nc$dim$rlon$len ##Longitude Length
  n.time <- length(time)

  year.one.array <- array(0,c(n.lon,n.lat,mlen))
  year.two.array <- array(0,c(n.lon,n.lat,mlen))
  year.three.array <- array(0,c(n.lon,n.lat,mlen))
  year.four.array <- array(0,c(n.lon,n.lat,mlen))
  year.five.array <- array(0,c(n.lon,n.lat,mlen))

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

  for (j in 35:n.lat) { 
    print(paste0('Latitude: ',j,' of ',n.lat))

    tasmax.subset <- ncvar_get(tasmax.nc,'tas',start=c(1,j,yst),count=c(-1,1,cnt))
    tasmax.list <- lapply(seq_len(nrow(tasmax.subset)), function(k) tasmax.subset[k,])
    rm(tasmax.subset)
    tasmin.subset <- ncvar_get(tasmin.nc,'tas',start=c(1,j,yst),count=c(-1,1,cnt))
    tasmin.list <- lapply(seq_len(nrow(tasmin.subset)), function(k) tasmin.subset[k,])
    rm(tasmin.subset)
    tas.subset <- ncvar_get(tas.nc,'tas',start=c(1,j,yst),count=c(-1,1,cnt))
    tas.list <- lapply(seq_len(nrow(tas.subset)), function(k) tas.subset[k,])
    rm(tas.subset)

    dwpt.max.subset <- ncvar_get(dwpt.max.nc,'dewpoint',start=c(1,j,yst),count=c(-1,1,cnt))
    dwpt.max.list <- lapply(seq_len(nrow(dwpt.max.subset)), function(k) dwpt.max.subset[k,])
    rm(dwpt.max.subset)
    dwpt.min.subset <- ncvar_get(dwpt.min.nc,'dewpoint',start=c(1,j,yst),count=c(-1,1,cnt))
    dwpt.min.list <- lapply(seq_len(nrow(dwpt.min.subset)), function(k) dwpt.min.subset[k,])
    rm(dwpt.min.subset)
    dwpt.avg.subset <- ncvar_get(dwpt.avg.nc,'dewpoint',start=c(1,j,yst),count=c(-1,1,cnt))
    dwpt.avg.list <- lapply(seq_len(nrow(dwpt.avg.subset)), function(k) dwpt.avg.subset[k,])
    rm(dwpt.avg.subset)

    wspd.max.subset <- ncvar_get(wspd.max.nc,'wspd',start=c(1,j,yst),count=c(-1,1,cnt))
    wspd.max.list <- lapply(seq_len(nrow(wspd.max.subset)), function(k) wspd.max.subset[k,])
    rm(wspd.max.subset)
    wspd.avg.subset <- ncvar_get(wspd.avg.nc,'wspd',start=c(1,j,yst),count=c(-1,1,cnt))
    wspd.avg.list <- lapply(seq_len(nrow(wspd.avg.subset)), function(k) wspd.avg.subset[k,])
    rm(wspd.avg.subset)

    insol.subset <- ncvar_get(insol.nc,'insol',start=c(1,j,yst),count=c(-1,1,cnt))
    insol.list <- lapply(seq_len(nrow(insol.subset)), function(k) insol.subset[k,])
    rm(insol.subset)

    objects <- find.all.months(months=months,tasmax=tasmax.list[[51]],tasmin=tasmin.list[[51]],tas=tas.list[[51]],
                               dwpt.max=dwpt.max.list[[51]],dwpt.min=dwpt.min.list[[51]],dwpt.avg=dwpt.avg.list[[51]],
                               wspd.max=wspd.max.list[[51]],wspd.avg=wspd.avg.list[[51]],insol=insol.list[[51]],
                               mon.ts=mon.ts,mon.indices=mon.indices,yrs=yrs)


    years.slct <- foreach(
                  tasmax=tasmax.list,                 
                  tasmin=tasmin.list,                 
                  tas=tas.list,                 
                  dwpt.max=dwpt.max.list,                 
                  dwpt.min=dwpt.min.list,                 
                  dwpt.avg=dwpt.avg.list,                 
                  wspd.max=wspd.max.list,
                  wspd.avg=wspd.avg.list,
                  insol=insol.list,
                  .export=c('find.all.months','months','mon.ts','mon.indices','yrs')
                           ) %dopar% {
                              objects <- find.all.months(months=months,tasmax=tasmax,tasmin=tasmin,tas=tas,
                                                                       dwpt.max=dwpt.max,dwpt.min=dwpt.min,dwpt.avg=dwpt.avg,
                                                                       wspd.max=wspd.max,wspd.avg=wspd.avg,insol=insol,
                                                                       mon.ts=mon.ts,mon.indices=mon.indices,yrs=yrs)
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

    year.four.list <- lapply(years.slct,function(x){return(x[,4])}) 
    year.four.matrix <- matrix(unlist(year.four.list),nrow=n.lon,ncol=length(months),byrow=T)
    year.four.array[,j,] <- year.four.matrix  

    year.five.list <- lapply(years.slct,function(x){return(x[,5])}) 
    year.five.matrix <- matrix(unlist(year.five.list),nrow=n.lon,ncol=length(months),byrow=T)
    year.five.array[,j,] <- year.five.matrix  

  }
  nc_close(tas.nc)
  years.array <- list(one=year.one.array,two=year.two.array,three=year.three.array,four=year.four.array,five=year.five.array)
  return(years.array)
}

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

var.spells <- function(yr.mon,quant,comp) {
    thresh <- quantile(yr.mon,quant)
    spells <- yr.mon*0
    if (comp=='greater') {spells[yr.mon >  thresh] <- 1}
    if (comp=='lesser')  {spells[yr.mon <  thresh] <- 1}    
    sqnce <- rle(spells)
    rv <- list(len=max(sqnce$lengths[sqnce$values==1]),
               cnt=sum(sqnce$values))
    return(rv)                 
}

find.tmy.months <- function(tas,insol,yrs1,yrs2,yrs3,yrs4,yrs5,
                   mon.ts,mon.indices,months) {

  rv <- rep(0,length(months))
  for (m in seq_along(months)) {
    month <- months[m]
    mon.ix <- mon.indices[,m]
    mon.vl <- which(mon.ix)  
    mon.tas <- tas[mon.ix]
    mon.insol <- insol[mon.ix]
    
    slct.yrs <- unique(c(yrs1[m],yrs2[m],yrs3[m],yrs4[m],yrs5[m]))
    warm.sqnce.len <- rep(0,length(slct.yrs))
    warm.sqnce.cnt <- rep(0,length(slct.yrs))
    cold.sqnce.len <- rep(0,length(slct.yrs))
    cold.sqnce.cnt <- rep(0,length(slct.yrs))
    dim.sqnce.len <- rep(0,length(slct.yrs))
    dim.sqnce.cnt <- rep(0,length(slct.yrs))

    for (i in seq_along(slct.yrs)) {
      yr.mon.ix <- mon.vl[mon.ts[mon.ix] %chin% paste0(slct.yrs[i],'-',month)]
      tas.yr.mon <- tas[yr.mon.ix]
      insol.yr.mon <- insol[yr.mon.ix]

      warm.spells <- var.spells(tas.yr.mon,0.67,'greater')
      warm.sqnce.len[i] <- warm.spells$len
      warm.sqnce.cnt[i] <- warm.spells$cnt

      cold.spells <- var.spells(tas.yr.mon,0.33,'lesser')
      cold.sqnce.len[i] <- cold.spells$len
      cold.sqnce.cnt[i] <- cold.spells$cnt

      dim.spells <- var.spells(insol.yr.mon,0.33,'lesser')
      dim.sqnce.len[i] <- dim.spells$len
      dim.sqnce.cnt[i] <- dim.spells$cnt
    }
    matrix.mon <- matrix(TRUE,nrow=3,ncol=5)    
    matrix.len <- rbind(warm.sqnce.len,cold.sqnce.len,dim.sqnce.len)    
    matrix.mon[matrix.len==max(matrix.len)] <- FALSE
    matrix.cnt <- rbind(warm.sqnce.cnt,cold.sqnce.cnt,dim.sqnce.cnt)
    matrix.mon[matrix.cnt==max(matrix.cnt)] <- FALSE
    matrix.mon[matrix.cnt==0] <- FALSE
    mon.choose <- which(apply(matrix.mon,2,prod)!=0)
    if (length(mon.choose)==0) {
       rv[m] <- slct.yrs[1] 
    } else {
       rv[m] <- slct.yrs[mon.choose[1]]
    }
  }
  return(rv)
}

##Use the selected years from TAS,INSOL,RHS and find the most average wind years

persistence.tmys <- function(tas.file,insol.file,
                             years.one,years.two,years.three,years.four,years.five) {

  months <- sprintf('%02d',1:12)
  mlen <- length(months) 

  tas.nc <- nc_open(tas.file)
  insol.nc <- nc_open(insol.file)

  time <- netcdf.calendar(tas.nc)
  mon.ts <- format(time,'%Y-%m')
  mon.fac <- format(time,'%m')
  yrs <- unique(format(time,'%Y'))
  years <- format(time,'%Y')
  n.lat <- tas.nc$dim$rlat$len ##Latitude Length
  n.lon <- tas.nc$dim$rlon$len ##Longitude Length
  n.time <- tas.nc$dim$time$len ##Longitude Length
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
    tas.subset <- ncvar_get(tas.nc,'tas',start=c(1,j,1),count=c(-1,1,-1))
    tas.list <- lapply(seq_len(nrow(tas.subset)), function(k) tas.subset[k,])
    rm(tas.subset)
    insol.subset <- ncvar_get(insol.nc,'insol',start=c(1,j,1),count=c(-1,1,-1))
    insol.list <- lapply(seq_len(nrow(insol.subset)), function(k) insol.subset[k,])
    rm(insol.subset)

    yr.one.sub <- years.one[,j,]
    yr.one.list <- lapply(seq_len(nrow(yr.one.sub)), function(k) yr.one.sub[k,])
    yr.two.sub <- years.two[,j,]
    yr.two.list <- lapply(seq_len(nrow(yr.two.sub)), function(k) yr.two.sub[k,])
    yr.three.sub <- years.three[,j,]
    yr.three.list <- lapply(seq_len(nrow(yr.three.sub)), function(k) yr.three.sub[k,])
    yr.four.sub <- years.four[,j,]
    yr.four.list <- lapply(seq_len(nrow(yr.four.sub)), function(k) yr.four.sub[k,])
    yr.five.sub <- years.five[,j,]
    yr.five.list <- lapply(seq_len(nrow(yr.five.sub)), function(k) yr.five.sub[k,])

    years.slct <- foreach(
                    tas=tas.list,                 
                    insol=insol.list,
                    yrs1=yr.one.list,                 
                    yrs2=yr.two.list,                 
                    yrs3=yr.three.list,                 
                    yrs4=yr.four.list,                 
                    yrs5=yr.five.list,                 
                    .export=c('find.tmy.months','mon.ts','mon.indices','months')
                             ) %dopar% {
                                objects <- find.tmy.months(tas=tas,insol=insol,
                                                           yrs1=yrs1,yrs2=yrs2,yrs3=yrs3,yrs4=yrs4,yrs5=yrs5,
                                                           mon.ts=mon.ts,mon.indices=mon.indices,months)
                                    }
    years.matrix <- matrix(unlist(years.slct),nrow=n.lon,ncol=length(months),byrow=T)
    years.array[,j,] <- years.matrix  
  }
  nc_close(tas.nc)
  nc_close(insol.nc)
  return(years.array)
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

read.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/'

tmp.dir <- tmpdir

if (!file.exists(tmpdir)) {
  dir.create(tmpdir,recursive=T)
}

if (1==1) {
gcm.files <- list.files(path=read.dir,pattern=gcm)
tasmax.file <- gcm.files[grep('tas_max',gcm.files)]
tasmin.file <- gcm.files[grep('tas_min',gcm.files)] 
tas.avg.file <- gcm.files[grep('tas_mean',gcm.files)]

dewpoint.max.file <- gcm.files[grep('dewpoint_max',gcm.files)] 
dewpoint.min.file <- gcm.files[grep('dewpoint_min',gcm.files)] 
dewpoint.avg.file <- gcm.files[grep('dewpoint_mean',gcm.files)]

wspd.max.file <- gcm.files[grep('wspd_max',gcm.files)] 
wspd.avg.file <- gcm.files[grep('wspd_mean',gcm.files)] 
insol.file <- gcm.files[grep('insol_sum',gcm.files)] 

input.files <- list(tasmax=tasmax.file,tasmin=tasmin.file,tas=tas.avg.file,
                    dewpoint.max=dewpoint.max.file,dewpoint.min=dewpoint.min.file,dewpoint.avg=dewpoint.avg.file,
                    wspd.max=wspd.max.file,wspd.avg=wspd.avg.file,insol=insol.file)

print('Copying TAS files')
file.copy(from=paste0(read.dir,tasmax.file),to=tmp.dir,overwrite=T)
file.copy(from=paste0(read.dir,tasmin.file),to=tmp.dir,overwrite=T)
file.copy(from=paste0(read.dir,tas.avg.file),to=tmp.dir,overwrite=T)

print('Copying Dewpoint files')
file.copy(from=paste0(read.dir,dewpoint.max.file),to=tmp.dir,overwrite=T)
file.copy(from=paste0(read.dir,dewpoint.min.file),to=tmp.dir,overwrite=T)
file.copy(from=paste0(read.dir,dewpoint.avg.file),to=tmp.dir,overwrite=T)

print('Copying Windspeed files')
file.copy(from=paste0(read.dir,wspd.max.file),to=tmp.dir,overwrite=T)
file.copy(from=paste0(read.dir,wspd.avg.file),to=tmp.dir,overwrite=T)

print('Copying INSOL file')
file.copy(from=paste0(read.dir,insol.file),to=tmp.dir,overwrite=T)

primary.years   <- find.test.ref.years(input.files,interval,tmp.dir)
years.one <- primary.years$one
save(years.one,file=paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/year_data/',
                           gcm,'.years.one.from.daily.',interval,'.RData'))
years.two <- primary.years$two
save(years.two,file=paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/year_data/',
                           gcm,'.years.two.from.daily.',interval,'.RData'))
years.three <- primary.years$three
save(years.three,file=paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/year_data/',
                             gcm,'.years.three.from.daily.',interval,'.RData'))
years.four <- primary.years$four
save(years.four,file=paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/year_data/',
                            gcm,'.years.four.from.daily.',interval,'.RData'))
years.five <- primary.years$five
save(years.five,file=paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/year_data/',
                            gcm,'.years.five.from.daily.',interval,'.RData'))

file.remove(paste0(tmp.dir,tasmax.file))
file.remove(paste0(tmp.dir,tasmin.file))
file.remove(paste0(tmp.dir,tas.avg.file))
print('Done with TAS')
file.remove(paste0(tmp.dir,dewpoint.max.file))
file.remove(paste0(tmp.dir,dewpoint.min.file))
file.remove(paste0(tmp.dir,dewpoint.avg.file))
print('Done with Dewpoint')
file.remove(paste0(tmp.dir,wspd.max.file))
file.remove(paste0(tmp.dir,wspd.avg.file))
print('Done with Windspeed')
file.remove(paste0(tmp.dir,insol.file))
print('Done with INSOL')


}##If block

if (1==1) {

##Find Average Months from TAS and INSOL Persistence
##Step 3
##– For the top five months the persistence of mean dry bulb temperature and daily
##global solar radiation are evaluated.
##– For temperature the number of consecutive warm days (> 67%ile) and
##consecutive cool days (<33 %ile) is calculated.
##– For solar radiation the number of consecutive low radiation days (<33 %ile) is
##calculated.
##– The highest ranking candidate month that meets the following criteria is included 
##in the TMY:
##• The month with the longest run is excluded,
##• The month with the most runs is excluded
##• Any month with zero runs is excluded.

print('Copying TAS,INSOL files')

gcm.files <- list.files(path=read.dir,pattern=gcm)
tas.file <- gcm.files[grep('tas_mean',gcm.files)]
insol.file <- gcm.files[grep('insol_mean',gcm.files)]


file.copy(from=paste0(read.dir,tas.file),to=tmp.dir,overwrite=T)
file.copy(from=paste0(read.dir,insol.file),to=tmp.dir,overwrite=T)

load(paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/year_data/',gcm,'.years.one.from.daily.',interval,'.RData'))
load(paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/year_data/',gcm,'.years.two.from.daily.',interval,'.RData'))
load(paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/year_data/',gcm,'.years.three.from.daily.',interval,'.RData'))
load(paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/year_data/',gcm,'.years.four.from.daily.',interval,'.RData'))
load(paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/year_data/',gcm,'.years.five.from.daily.',interval,'.RData'))

tmy.years   <- persistence.tmys(paste0(tmp.dir,tas.file),paste0(tmp.dir,insol.file),
                                years.one,years.two,years.three,years.four,years.five)
                         
save(tmy.years,file=paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/year_data/',
                           gcm,'.cwec.tmy.years.from.daily.',interval,'.RData'))
file.remove(paste0(tmp.dir,tas.file))
file.remove(paste0(tmp.dir,insol.file))

print('Done with Persistence')

}

