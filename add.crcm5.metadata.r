##Script to create spatially complete DQM files that are 5 years in length

ptm <- proc.time()

library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R')

##Get from the DQM file
get.time.atts <- function(cal,units) { 
  time.atts <- list(standard_name = "time",
                    long_name='time',
                    units = units,
                    axis='T',
                    calendar = cal)
  return(time.atts)
}

get.date.bounds <- function(nc) {
  dates <- netcdf.calendar(nc)
  yst <-  gsub('-','',format(head(dates,1),'%Y-%m-%d'))
  yen <-  gsub('-','',format(tail(dates,1),'%Y-%m-%d'))
  rv <- c(yst,yen)
  return(rv)
}

get.new.var.name <- function(var.name) {

  rv <- switch(var.name,
               GWDI='discharge',
               HR='rhs',
               HU='huss',
               PN='psl',
               PR='pr',
               S6='snc',
               SD='snd',
               STFL='streamflow',
               T9='tasmax',
               T5='tasmin',
               TDRA='drainage',
               TJ='tas',
               TRAF='runoff',
               UD='uas',
               VD='vas')                         
  return(rv)
}

get.standard.atts <- function() {
  lon.atts <- list(standard_name="longitude",
                   long_name = "longitude",
                   units = "degrees_east")
                     
  lat.atts <- list(standard_name="latitude",
                   long_name = "latitude",
                   units = "degrees_north")

  rlon.atts <- list(long_name="longitude in oblique transverse grid",
                    standard_name="grid_longitude",
                    units = "degrees",
                    axis = "X")

  rlat.atts <- list(long_name="latitude in oblique transverse grid",
                    standard_name="grid_latitude",
                    units = "degrees",
                    axis = "Y")

  proj.atts <- list(proj = "ob_tran",
                    R = "6370997.0",
                    units = "m" ,
                    lon_0 = "-97.0" ,
                    o_lon_p = "180.0",
                    o_lat_p = "42.5",
                    o_proj = "longlat",
                    lat_0 = "58.04")
  
  rv <- list(lon=lon.atts,
             lat=lat.atts,
             rlon=rlon.atts,
             rlat=rlat.atts,
             proj=proj.atts)
  return(rv)
}

get.variable.atts <- function(var.name) {

  pr.atts <- list(standard_name = "precipitation_flux",
                  long_name = "Precipitation",
                  missing_value = 1.e+20,
                  cell_methods = "time: mean",
                  units = "kg m-2 d-1",
                  coordinates='lon lat')

  tasmax.atts <- list(standard_name = "air_temperature",
                      long_name = "Daily Maximum Near-Surface Air Temperature",
                      units = "degC",
                      missing_value =1.e+20,
                      cell_methods = "time: maximum",
                      coordinates='lon lat')

  tasmin.atts <- list(standard_name = "air_temperature",
                      long_name = "Daily Minimum Near-Surface Air Temperature",
                      units = "degC",
                      missing_value = 1.e+20,
                      cell_methods = "time: minimum",
                      coordinates='lon lat')

  psl.atts <- list(standard_name = "sea_level_pressure",
                      long_name = "Sea Level Pressure",
                      units = "hPa",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')

  uas.atts <- list(standard_name = "eastward_wind",
                      long_name = "Eastward Near-Surface Wind",
                      units = "m s-1",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')

  vas.atts <- list(standard_name = "northward_wind",
                      long_name = "Northward Near-Surface Wind",
                      units = "m s-1",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')

  snd.atts <- list(standard_name = "Snow Depth",
                      long_name = "surface_snow_thickness",
                      units = "m",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')

  snc.atts <- list(standard_name = "Snow Area Fraction",
                      long_name = "surface_snow_area_fraction",
                      units = "%",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')

  runoff.atts <- list(standard_name = "ACCUM. OF TOTAL SURFACE RUNOFF",
                      long_name = "total_surface_runoff",
                      units = "kg m-2 s-1",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')
 
  streamflow.atts <- list(standard_name = "SURF. WATER STREAMFLOW",
                      long_name = "surface_water_streamflow",
                      units = "m^3 s-1",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')
 
  drainage.atts <- list(standard_name = "ACCUM. OF BASE DRAINAGE",
                      long_name = "base_drainage",
                      units = "kg m-2 s-1",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')
  discharge.atts <- list(standard_name = "GROUND WATER DISCHARGE",
                      long_name = "ground_water_discharge",
                      units = "m^3 s-1",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')                      
                     
  var.atts <- switch(var.name,
                     drainage=drainage.atts,
                     discharge=discharge.atts,
                     pr=pr.atts,
                     tasmax=tasmax.atts,
                     tasmin=tasmin.atts,       
                     psl=psl.atts,
                     runoff=runoff.atts,
                     snd=snd.atts,
                     snc=snc.atts,
                     streamflow=streamflow.atts,
                     uas=uas.atts,
                     vas=vas.atts)

                     
}


##Global Attributes
##gcm is the gcm name e.g. CRCM5
##drive.institude is the centre e.g. 
##rcp is the RCP as RCP8.5
get.global.atts <- function(gcm,freq,rcp) {

  drive.centre <- switch(gcm,
                         ERAI='ECMWF',
                         CanESM2='CCCMA',
                         HadGEM2='MOHC',
                         MPI='MPI-M')

  drive.centre.name <- switch(gcm,
                              ERAI='European Centre for Medium‐Range Weather Forecasts',
                              CanESM2='Canadian Centre for Climate Modelling and Analysis, Victoria, BC, Canada',
                              HadGEM2='Met Office Hadley Centre',
                              MPI='Max Planck Institute for Meteorology')

  global.atts <- list(institution="Universite du Quebec a Montreal",
                   contact="http://cnrcwp.ca",
                   Conventions="CF-1.4",
                   institute_id ="UQAM",
                   domain='Western Canada',
                   creation_date="Thu, Apr 19, 2018 7:11 AM",
                   frequency=freq,
                   product="regional climate model output",
                   modeling_realm="atmos",
                   project_id='CNRCWP, NEI',
                   table_id='Table day (10 Jun 2010)',
                   references="Hernández-Díaz L, Laprise R, Sushama L, Martynov A, Winger K, Dugas B (2013) Climate simulation over the CORDEX-Africa domain using the fifth-generation Canadian regional climate model (CRCM5). Clim Dyn 40:1415–1433. doi: 10.1007/s00382-012-1387-z",
                   driving_experiment=paste0('historical',rcp),
                   driving_experiment_id=paste0('historical',rcp),
                   driving_institution = drive.centre.name, ##Full name
                   driving_institute_id = drive.centre, ##Acronym
                   driving_model_id = gcm,
                   driving_realization = "r1i1p1",
                   driving_initialization_method='1',
                   driving_physics_version = '1',
                   title = "Fifth-generation of the Canadian Regional Climate Model")
  return(global.atts)
}

##Filename format is: 'pr_day_BCCAQ2+ANUSPLIN300+ACCESS1-0_historical+rcp45_r1i1p1_19500101-21001231.nc'

