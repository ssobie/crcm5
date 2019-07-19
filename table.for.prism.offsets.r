##Script to adjust EPW file temperature using PRISM

library(ncdf4)

source('/storage/home/ssobie/code/repos/crcm5/read.write.epw.r',chdir=T)


##----------------------------------------------------------

##Return the lon/lat coordinates for the supplied EPW file

get_epw_coordinates <- function(epw.dir,epw.file) {

   epw <- read.epw.file(epw.dir,epw.file)
   epw.header <- epw$header
   epw.first <- strsplit(epw.header[1],',')[[1]]
   lon <- as.numeric(epw.first[8]) ##Fixed location
   lat <- as.numeric(epw.first[7])
   if (lon < -180 | lon > 180) {
      stop('Ill defined longitude coordinate')
   }
   if (lat < 40 | lat > 90) {
      stop('Ill defined latitude coordinate')
   }
       
   rv <- c(lon,lat)
   return(rv)
}

##----------------------------------------------------------

##Create a list of coordinates for all available EPW files

list_of_epw_coordinates <- function(epw.dir,epw.files) {
   epw.coords <- matrix(NA,nrow=length(epw.files),ncol=2)
   for (i in seq_along(epw.files)) {
      epw.coords[i,] <- get_epw_coordinates(epw.dir,epw.files[i])
   }
   return(epw.coords)
}

##----------------------------------------------------------

##Search the available EPW files and find the nearest to the
##supplied coordinates

find_closest_epw_file <- function(coords,
                                  epw.dir='/storage/data/projects/rci/weather_files/wx_files/') 
                         {
  epw.files <- list.files(path=epw.dir,pattern='CWEC.epw')
  epw.coords <- list_of_epw_coordinates(epw.dir,epw.files)
  wx.ix <- which.min(apply((epw.coords - matrix(coords,nrow=nrow(epw.coords),ncol=2,byrow=T))^2,1,sum))
  wx.selected <- epw.files[wx.ix]
  print(wx.selected)
  return(list(dir=epw.dir,file=wx.selected))
}

##----------------------------------------------------------

##Open prism climatologies file

prism_nc <- function(nc,var.name,prism.dir) { 
   prism.filename <- paste0(prism.dir,var.name,'_monClim_PRISM_historical_run1_198101-201012.nc')
   nc <- nc_open(prism.filename)  
   return(nc)
}

##----------------------------------------------------------

##Find PRISM cell containing coordinates

get_prism_indices <- function(nc,coords) {

  lon <- ncvar_get(nc,'lon')
  lon.ix <- which.min(abs(coords[1]-lon))
  lat <- ncvar_get(nc,'lat')
  lat.ix <- which.min(abs(coords[2]-lat))
  rv <- c(lon.ix,lat.ix)
  return(rv)
}

##----------------------------------------------------------

##Calculate the TAS climatologies for the selected PRISM 
##cell

prism_tas <- function(tx.nc,tn.nc,cell,prism.dir) {
 
  tasmax <- ncvar_get(tx.nc,'tmax',start=c(cell,1),count=c(1,1,-1))
  tasmin <- ncvar_get(tn.nc,'tmin',start=c(cell,1),count=c(1,1,-1))
  tas <- (tasmax+tasmin)/2
  return(tas)  
}

##----------------------------------------------------------

##Given coordinates find the nearest weather file and adjust
##the temperature series based on the PRISM climatologies

generate_prism_offset <- function(new.location,lon,lat,epw.dir,prism.dir,tx.nc,tn.nc) {
   coords <- c(lon,lat)
   epw.closest <- find_closest_epw_file(coords,epw.dir)
   epw.name <- strsplit(epw.closest$file,'_')[[1]][3]
   epw.closest.coords <- get_epw_coordinates(epw.closest$dir,
                                             epw.closest$file)
   prism.cell <- get_prism_indices(tx.nc,coords)
   epw.cell <- get_prism_indices(tx.nc,epw.closest.coords)
   prism.loc.tas <- prism_tas(tx.nc,tn.nc,prism.cell,prism.dir)
   prism.loc.seas <- c(prism.loc.tas[1:12],
                       mean(prism.loc.tas[c(1,2,12)]),
                       mean(prism.loc.tas[c(3,4,5)]),
                       mean(prism.loc.tas[c(6,7,8)]),
                       mean(prism.loc.tas[c(9,10,11)]),
                       prism.loc.tas[13])
   prism.epw.tas <- prism_tas(tx.nc,tn.nc,epw.cell,prism.dir)
   prism.epw.seas <- c(prism.epw.tas[1:12],
                       mean(prism.epw.tas[c(1,2,12)]),
                       mean(prism.epw.tas[c(3,4,5)]),
                       mean(prism.epw.tas[c(6,7,8)]),
                       mean(prism.epw.tas[c(9,10,11)]),
                       prism.epw.tas[13])

   prism.diff <- prism.loc.seas - prism.epw.seas
   rv <- list(loc=round(prism.loc.seas,1),
              epw=round(prism.epw.seas,1),
              diff=round(prism.diff,1),
              name=epw.name)
      
   return(rv)
}

##-------------------------------------------------------------
epw.dir <- '/storage/data/projects/rci/weather_files/wx_files/'
prism.dir <- '/storage/data/climate/PRISM/dataportal/'
new.location <- 'Cowichan-Hospital'

tx.filename <- paste0(prism.dir,'tmax_monClim_PRISM_historical_run1_198101-201012.nc')
tx.nc <- nc_open(tx.filename)  
tn.filename <- paste0(prism.dir,'tmin_monClim_PRISM_historical_run1_198101-201012.nc')
tn.nc <- nc_open(tn.filename)  

locations <- list(list(name='Cowichan-Hospital',lon=-123.722775,lat=48.786147),
                  list(name='Royal-Columbian-Hospital',lon=-122.891355,lat=49.226239),
                  list(name='Dogwood',lon=-123.117856,lat=49.218684),
                  list(name='NRGH',lon=-123.97049,lat=49.18498),
                  list(name='Richmond',lon=-123.1469,lat=49.1687),
                  list(name='Lions-Gate',lon=-123.0684,lat=49.3209),
                  list(name='1st-and-Clark',lon=-123.0774,lat=49.2697),
                  list(name='UVic-Residence',lon=-123.30,lat=48.46))

len <- 19
offset.table <- matrix(NA,nrow=length(locations)*4,ncol=len)
loc.seq <- seq(1,length(locations)*4,4)
locs <- 1:length(locations)
for (i in locs) {
  location <- locations[[i]]
  print(location$name)
  j <- loc.seq[i]
  result <- generate_prism_offset(location$name,location$lon,location$lat,epw.dir,prism.dir,tx.nc,tn.nc)
  offset.table[j,1] <- location$name 
  offset.table[j+1,1] <- result$name
  offset.table[j+2,1] <- ''
  offset.table[j,2] <- 'Site'
  offset.table[j+1,2] <- 'CWEC'
  offset.table[j+2,2] <- 'Site-CWEC'
  offset.table[j,3:len] <- result$loc
  offset.table[j+1,3:len] <- result$epw
  offset.table[j+2,3:len] <- result$diff
  offset.table[j+3,1:len] <- ''

}

header <- c('Site','',month.abb,'Winter','Spring','Summer','Fall','Annual')

write.table(rbind(header,offset.table),file='/storage/data/projects/rci/weather_files/wx_files/offsets/epw.offset.table.csv',
            quote=FALSE,row.name=FALSE,col.name=FALSE,sep=',')

nc_close(tx.nc)
nc_close(tn.nc)

