##Script to convert hourly Van Isle files to daily

make.daily.tmy <- function(gcm,avg.list,fxn) {

  ##gcm <- 'MPI'
  ##avg.list <- c('insol','rhs','tas','wspd','dewpoint') ##MEAN
  ##avg.list <- c('tas','wspd','dewpoint') ##MAX
  ##avg.list <- c('tas','dewpoint') ##MIN
  ##avg.list <- 'insol' ##Sum
  ##fxn <- 'sum'

  read.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/'
  write.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/tmy_files/daily/'

  all.files <- list.files(path=read.dir,pattern=gcm) 

  ##Aggregate by averaging
  for (var.name in avg.list) {
      print(var.name)
      tmy.files <- all.files[grep(paste0('^',var.name),all.files)]

      past.file <- tmy.files[grep('2010',tmy.files)]
      output.past <- gsub(paste0(var.name,'_'),paste0(var.name,'_',fxn,'_day_'),past.file)
      work <- paste0('cdo -O day',fxn,' ',read.dir,past.file,' ',write.dir,output.past)
      ##system(work)

      proj.file <- tmy.files[grep('2021',tmy.files)]
      output.proj <- gsub(paste0(var.name,'_'),paste0(var.name,'_',fxn,'_day_'),proj.file)
      work <- paste0('cdo -O day',fxn,' ',read.dir,proj.file,' ',write.dir,output.proj)
      system(work)

      morph.file <-  all.files[grep(paste0('morphed_',var.name),all.files)]
      output.morph <- gsub(paste0(var.name,'_'),paste0(var.name,'_',fxn,'_day_'),morph.file)
      work <- paste0('cdo -O day',fxn,' ',read.dir,morph.file,' ',write.dir,output.morph)
##      system(work)
  }
}


make.daily.rcm <- function(gcm,avg.list,fxn) {

  ##gcm <- 'MPI'
  ##avg.list <- c('insol','rhs','tas','wspd','dewpoint') ##MEAN
  ##avg.list <- c('tas','wspd','dewpoint') ##MAX
  ##avg.list <- c('tas','dewpoint') ##MIN
  ##avg.list <- 'insol' ##Sum
  ##fxn <- 'sum'

  read.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/hourly/'
  write.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/daily/'

  all.files <- list.files(path=read.dir,pattern=gcm) 

  ##Aggregate by averaging
  for (var.name in avg.list) {
      print(var.name)
      rcm.file <- all.files[grep(paste0('^',var.name),all.files)]
      output.rcm <- gsub(paste0(var.name,'_hour_'),paste0(var.name,'_',fxn,'_day_'),rcm.file)
      work <- paste0('cdo -O day',fxn,' ',read.dir,rcm.file,' ',write.dir,output.rcm)
      system(work)
  }
}



##Windspeed
if (1==0) {
all.files <- list.files(path=read.dir,pattern=gcm) 
##Windspeed
uas.file <- all.files[grep('uas',all.files)]
vas.file <- all.files[grep('uas',all.files)]
uas.out <- gsub('hour','day',uas.file)
vas.out <- gsub('hour','day',vas.file)
wspd.hourly <- gsub('uas','wspd',uas.file)

uas2 <- paste0('cdo -O mul ',read.dir,uas.file,' ',read.dir,uas.file,' ',write.dir,'uas2.nc')
print(uas2)
system(uas2)
vas2 <- paste0('cdo -O mul ',read.dir,vas.file,' ',read.dir,vas.file,' ',write.dir,'vas2.nc')
print(vas2)
system(vas2)

wspd2 <- paste0('cdo -O add ',write.dir,'uas2.nc ',write.dir,'vas2.nc ',write.dir,'wspd2.nc')
print(wspd2)
system(wspd2)

##Hourly Windspeed
wspd <- paste0('cdo -O sqrt ',write.dir,'wspd2.nc ',read.dir,wspd.hourly)
print(wspd)
system(wspd)

output.file <- gsub('hour','day',wspd.hourly)
work <- paste0('cdo -O daymean ',read.dir,wspd.hourly,' ',write.dir,output.file)
print(work)
system(work)

file.remove(paste0(write.dir,'wspd2.nc'))
file.remove(paste0(write.dir,'uas2.nc'))
file.remove(paste0(write.dir,'vas2.nc'))

}