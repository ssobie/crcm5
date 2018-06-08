##Script to convert hourly Van Isle files to daily

gcm <- 'MPI'

read.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/hourly/'
write.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/daily/'

avg.list <- c('insol','irflux','psl','rhs','tas','uas','vas')

all.files <- list.files(path=read.dir,pattern=gcm) 

##Aggregate by averaging
for (var.name in avg.list) {
  input.file <- all.files[grep(var.name,all.files)]
  output.file <- gsub('hour','day',input.file)
  work <- paste0('cdo -O daymean ',read.dir,input.file,' ',write.dir,output.file)
  print(work)
  system(work)
}

##Daily total precipitation
input.file <- all.files[grep('pr',all.files)]
output.file <- gsub('hour','day',input.file)
work <- paste0('cdo -O daysum ',read.dir,input.file,' ',write.dir,output.file)
print(work)
system(work)

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

