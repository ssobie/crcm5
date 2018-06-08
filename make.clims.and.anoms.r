library(ncdf4)
library(PCICt)

sub.by.time <- function(input.file,interval,read.dir,write.dir) {

  nc <- nc_open(paste(read.dir,input.file,sep=''))
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.series <- origin.pcict + time.values*86400
  years <- format(time.series,'%Y')
  yrs <- strsplit(interval,'-')[[1]]

  st <- head(grep(yrs[1],years),1)-1
  en <- tail(grep(yrs[2],years),1)-1

  sub.file <- 'sub.nc'
  system(paste('ncks --overwrite -d time,',st,',',en,' ',read.dir,input.file,' ',write.dir,sub.file,sep=''))
  print(paste('ncks --overwrite -d time,',st,',',en,' ',read.dir,input.file,' ',write.dir,sub.file,sep=''))
  nc_close(nc)
  rv <- paste0(write.dir,sub.file)
  return(rv)
}

run.climatologies <- function(gcm,var.name,interval) {

  read.dir <-  '/storage/data/climate/downscale/RCM/CRCM5/reconfig/MPI/'
  write.dir  <- paste0(read.dir,'climatologies/')

  ##-----------------------------------------------------------------------------
  all.files <- list.files(path=read.dir,pattern=gcm)
  var.file <- all.files[grep(var.name,all.files)]
  new.int <- paste0(strsplit(interval,'-')[[1]][1],'0101-',strsplit(interval,'-')[[1]][2],'1231')
  time.write <- gsub(pattern='[0-9]{8}-[0-9]{8}',replacement=new.int,var.file)

  avg.file <- gsub(pattern='_hour_',replacement='_climatology_',var.file)
  sub.file <- sub.by.time(var.file,interval=interval,read.dir,write.dir)
  system(paste('cdo -s -O timmean ',sub.file,' ',write.dir,time.write,sep=''))
  system(paste('rm ',sub.file,sep=''))

}

##run.climatologies('ERA-Interim','tas','1981-2010')
run.climatologies('MPI','tas','1981-2010')
run.climatologies('MPI','tas','2021-2050')
