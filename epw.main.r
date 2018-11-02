##Script to plot the EPW file series

library(ncdf4)
library(PCICt)

source('/storage/home/ssobie/code/repos/crcm5/read.write.epw.r',chdir=T)
source('/storage/home/ssobie/code/repos/crcm5/epw.belcher.functions.r',chdir=T)


##------------------------------------------------------------------------------
##Match for EPW fields

get_field_index <- function(var.name) {

   field.names <- c('year', 'month', 'day', 'hour', 'minute',
      'data_source_and_uncertainty_flags', 'dry_bulb_temperature',
      'dew_point_temperature', 'relative_humidity',
      'atmospheric_station_pressure', 'extraterrestrial_horizontal_radiation',
      'extraterrestrial_direct_normal_radition',
      'horizontal_infrared_radiation_intensity', 'global_horizontal_radiation',
      'direct_normal_radiation', 'diffuse_horizontal_radiation',
      'global_horizontal_illuminance', 'direct_normal_illuminance',
      'diffuse_horizontal_illuminance', 'zenith_luminance', 'wind_direction',
      'wind_speed', 'total_sky_cover', 'opaque_sky_cover', 'visibility',
      'ceiling_height', 'present_weather_observation', 'present_weather_codes',
      'precipitable_water', 'aerosol_optical_depth', 'snow_depth',
      'days_since_last_snowfall', 'albedo', 'liquid_precipitation_depth',
      'liquid_precipitation_quantity')
   ix <- grep(var.name,field.names)
}

##------------------------------------------------------------------------------

sub_by_time <- function(var.name,lonc,latc,interval,input.file,gcm,read.dir) {

  print(input.file)              
  nc <- nc_open(paste(read.dir,gcm,'/',input.file,sep=''))
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.series <- origin.pcict + time.values*86400
  years <- format(time.series,'%Y')
  yrs <- strsplit(interval,'-')[[1]]

  feb.flag <- grep('-02-29',time.series)

  new.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal='365_day')
                           
  new.time0 <- (as.PCICt(format(time.series[1],'%Y-%m-%d'),cal='365_day') - new.origin)/86400
  if (length(feb.flag)==0) {
     new.values <- seq(as.numeric(new.time0),by=1,length.out=length(time.values))
  } else {
     new.values <- seq(as.numeric(new.time0),by=1,length.out=length(time.values[-feb.flag]))
  }

  new.series <- new.origin + new.values*86400
  print(range(time.series))
  print(range(new.series))

  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  
  lon.ix <- which.min(abs(lonc-lon))
  lat.ix <- which.min(abs(latc-lat))

  data.raw <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  data <- data.raw
  if (length(feb.flag!=0)) {
    data <- data.raw[-feb.flag]
  }

  years <- format(new.series,'%Y')
  st <- head(grep(yrs[1],years),1)
  en <- tail(grep(yrs[2],years),1)
  if (grepl('HadGEM',gcm) & yrs[2]=='2100') {
    en <- length(years)
  }
  nc_close(nc)
  rv <- list(data=data[st:en],time=new.series[st:en])
  return(rv)
}

##------------------------------------------------------------------------------

##**************************************************************************************

##scenario <- 'rcp85'
##interval <- '2071-2100'

##lon <- -122.36
##lat <- 49.03

##method <- 'seasonal'
##rlen <- ''

if (method!='roll') { 
  rlen <- ''
}

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}
print(scenario)
print(method)
print(rlen)
print(lon)
print(lat)

print(infile)

epw.dir <- '/storage/data/projects/rci/weather_files/wx_files/offsets/'
write.dir <- '/storage/data/projects/rci/weather_files/wx_files/morphed_files/'
##infile ##'CAN_BC_1st_and_Clark_offset_from_VANCOUVER-INTL-A_1108395_CWEC.epw'

present.epw.file <- 'CAN_BC_UVic_residence_offset_from_VICTORIA-UNIVERSITY-CS_1018598_CWEC.epw'
future.epw.file <- paste0('MORPHED_',toupper(method),rlen,'_TAS_CAN_BC_UVic_residence_',interval,'_CWEC.epw')
print(future.epw.file)

full.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')
sub.list <- c('ACCESS1-0','CanESM2','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MRI-CGCM3')

##------------------------------------------------------------------------------


epw.present <- read.epw.file(epw.dir,present.epw.file)

##Create one year of daily dates

epw.morphed.tas <- morph_dry_bulb_temp(epw.present,
                        lon,lat,
                        gcm.list=full.list,
                        gcm.dir='/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/',
                        scenario,interval=interval,
                        method=method,rlen=rlen)


write.epw.file(epw.morphed.tas$data,epw.morphed.tas$header,write.dir,future.epw.file)

##epw.morphed.rhs <- generate_stretched_series(epw.morphed.tas,'relative_humidity','rhs',
##                        lon,lat,
##                        gcm.list=sub.list,
##                        gcm.dir="/storage/data/climate/downscale/CMIP5/building_code/",
##                        scenario,interval=interval,
##                        method=method,rlen=rlen)


##epw.morphed.dwpt <- morph_dew_point_temp(epw.present,lon,lat,gcm.list,
##                        gcm.dir="/storage/data/climate/downscale/CMIP5/building_code/",
##                        scenario,interval,method=method,rlen=rlen)


##epw.morphed.dnr <- generate_stretched_series(epw.present,'direct_normal_radiation','clt',
##                        lon,lat,
##                        gcm.list=sub.list,
##                        gcm.dir="/storage/data/climate/downscale/CMIP5/building_code/",
##                        scenario,interval=interval,
##                        method=method,rlen=rlen)


