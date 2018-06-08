##Script to explore the tas anomalies in the Coast Mountains

library(ncdf4)
library(PCICt)

gcm <- 'MPI'
var.name <- 'tas'
rcm.var <- 'TJ'


read.dir <- '/storage/data/climate/downscale/RCM/CRCM5/'

clim.dir <- paste0('/storage/data/climate/downscale/RCM/CRCM5/reconfig/',gcm,'/climatologies/')
anoms.file <- paste0(clim.dir,gcm,'.anoms.nc')

anc <- nc_open(anoms.file)
anoms <- ncvar_get(anc,var.name)
alon <- ncvar_get(anc,'lon')
alat <- ncvar_get(anc,'lat')

anoms.flag <- which(anoms < 0 ,arr.ind=T)

n.lons <- alon[anoms.flag]
n.lats <- alat[anoms.flag]

nc <- nc_open(paste0(read.dir,gcm,'/',gcm,'_WC011_modified_TJ_1980-2050.nc'))
lon <- ncvar_get(nc,'lon')
lat <- ncvar_get(nc,'lat')

  time.atts <- ncatt_get(nc,'t')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'t')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.series <- origin.pcict + time.values*3600


if (1==1) {
##tas.sub <- matrix(0,nrow=nc$dim$t$len,ncol=nrow(anoms.flag))
lon.ixs <- c()
coords <- c()
for (i in 1:nrow(anoms.flag)) {  
    print(i)
  lon.ix <- which((min(abs(lon-n.lons[i])) == abs(lon-n.lons[i])),arr.ind=T)
  lon.ixs <- rbind(lon.ixs,lon.ix)
  print(lon.ix)
  print(n.lons[i])
  print(lon[lon.ix])
  coords <- rbind(coords,c(lon[lon.ix],lat[lon.ix]))
  ##print(n.lats[i])
  ##print(lat[lat.ix])
  ##lat.ix <- which((min(abs(lat-n.lats[i])) == abs(lat-n.lats[i])),arr.ind=T)
  ##tas.sub[,i] <- ncvar_get(nc,'TJ',start=c(lon.ix[1],lon.ix[2],1,1),count=c(1,1,1,-1))
}
}


##  lon.ix <- which((min(abs(lon-n.lons[4])) == abs(lon-n.lons[4])),arr.ind=T)
##  tas.sub <- ncvar_get(nc,'TJ',start=c(lon.ix[1],lon.ix[2],1,1),count=c(1,1,1,-1))




nc_close(nc)
nc_close(anc)