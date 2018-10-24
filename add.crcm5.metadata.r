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
               AD='irflux',
               N3='asp',
               AV='latflux',
               GIAC='giac',
               GIML='giml',
               GLD='gld',
               GLF='glf',
               GSAB='gsab',
               GSAC='gsac',
               GSML='gsml',
               GVOL='gvol',
               GWDI='discharge',
               GWST='gwst',
               GZ='gz',
               HR='rhs',
               HU='huss',
               I1='swater',
               I2='sice',
               I4='spw',
               I5='smass',
               MS='smelt',
               N4='insol',
               P0='ps',
               PN='psl',
               PR='pr',
               S6='snc',
               SD='snd',
               STFL='streamflow',
               SWSL='swlake',
               SWSR='swriver',
               T9='tasmax',
               T5='tasmin',
               TJ='tas',
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


  tas.atts <- list(standard_name = "air_temperature",
                   long_name = "Daily Average Near-Surface Air Temperature",
                   units = "degC",
                   missing_value =1.e+20,
                   cell_methods = "time: mean",
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

  gz.atts <- list(standard_name = "GEOPOTENTIAL HEIGHT 1000mb",
                      long_name = "geopotential_height_1000",
                      units = "m",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')                      

  insol.atts <- list(standard_name = "ACCUMULATION OF SOLAR RADIATION",
                      long_name = "accumulation_solar_radiation",
                      units = "W m-2",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')                      

  huss.atts <- list(standard_name = "SPECIFIC HUMIDITY",
                      long_name = "specific_humidity",
                      units = "kg/kg",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')                      

  rhs.atts <- list(standard_name = "RELATIVE HUMIDITY",
                      long_name = "relative_humidity",
                      units = "fraction",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')                      

  irflux.atts <- list(standard_name = "ACCUMULATION OF FDSI(IR ENERGY FLUX TOWARDS GROUND)",
                      long_name = "ir_energy_flux_towards_ground",
                      units = "W m-2",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')                      

  latflux.atts <- list(standard_name = "ACCUMULATION OF FV(SURFACE LATENT FLUX)",
                      long_name = "surface_latent_flux",
                      units = "W m-2",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')                      

  giac.atts <- list(standard_name = "ACCUMUL. OF GLACIER ICE ACCUMULATION",
                      long_name = "glacier_ice_accumulation",
                      units = "mm weq s-1",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')                      

  giml.atts <- list(standard_name = "ACCUMUL. OF GLACIER ICE MELT",
                      long_name = "glacier_ice_melt",
                      units = "mm weq s-1",
                      missing_value = 1.e+20,
                      cell_methods = "time: mean",
                      coordinates='lon lat')                      

  gld.atts <- list(standard_name = "MEAN GLACIER DEPTH FOR WHOLE GRID BOX",
                   long_name = "glacier_depth",
                   units = "m",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  glf.atts <- list(standard_name = "GLACIER FRACTION WRT WHOLE GRID",
                   long_name = "glacier_fraction",
                   units = "fraction",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  gsab.atts <- list(standard_name = "ACCUMUL. OF SNOW ABLATION ON GLACIER",
                   long_name = "glacier_snow_ablation",
                   units = "mm weq s-1",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  gsac.atts <- list(standard_name = "ACCUMUL. OF SNOW ACCUMUL. ON GLACIER",
                   long_name = "glacier_snow_accumulation",
                   units = "mm weq s-1",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  gsml.atts <- list(standard_name = "ACCUMUL. OF SNOW MELT ON GLACIER",
                   long_name = "glacier_snow_melt",
                   units = "mm weq s-1",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  gvol.atts <- list(standard_name = "GLACIER VOLUME FOR WHOLE GRID BOX",
                   long_name = "glacier_volume",
                   units = "m^3",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  gwst.atts <- list(standard_name = "GROUND WATER STORE",
                   long_name = "ground_water_store",
                   units = "m^3",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  swater.atts <- list(standard_name = "SOIL VOLUMETRIC WATER CONTENTS",
                   long_name = "soil_water_contents",
                   units = "m^3 / m^3",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  sice.atts <- list(standard_name = "SOIL VOLUMETRIC ICE CONTENTS",
                   long_name = "soil_ice_contents",
                   units = "m^3 / m^3",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  spw.atts <- list(standard_name = "WATER IN THE SNOW PACK",
                   long_name = "snow_pack_water",
                   units = "kg m-2",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  smass.atts <- list(standard_name = "SNOW MASS",
                   long_name = "snow_mass",
                   units = "kg m-2",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  smelt.atts <- list(standard_name = "MELTING SNOW FROM SNOWPACK",
                   long_name = "snow_melt",
                   units = "kg m-2 s-1",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  asp.atts <- list(standard_name = "ACCUM. OF SOLID PRECIP. USED BY LAND SURFACE SCHEMES (LAGGS 1 TIME STEP FROM PR)",
                   long_name = "accumulation_of_solid_precip",
                   units = "kg m-2 day-1",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  ps.atts <- list(standard_name = "SURFACE PRESSURE",
                   long_name = "surface_pressure",
                   units = "hPa",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  swlake.atts <- list(standard_name = "SURF. WATER STORE (LAKE)",
                   long_name = "surface_lake_water",
                   units = "m^3",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  swriver.atts <- list(standard_name = "SURF. WATER STORE (RIVER)",
                   long_name = "surface_river_lake",
                   units = "m^3",
                   missing_value = 1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')                      

  dewpoint.atts <- list(standard_name = "dewpoint_temperature",
                   long_name = "Near-Surface Dewpoint Temperature",
                   units = "degC",
                   missing_value =1.e+20,
                   cell_methods = "time: mean",
                   coordinates='lon lat')
                     
  var.atts <- switch(var.name,                    
                     asp=asp.atts,
                     drainage=drainage.atts,
                     discharge=discharge.atts,
                     dewpoint=dewpoint.atts,
                     giac=giac.atts,
                     giml=giml.atts,
                     gld=gld.atts,
                     glf=glf.atts,
                     gsab=gsab.atts,
                     gsac=gsac.atts,
                     gsml=gsml.atts,
                     gvol=gvol.atts,
                     gwst=gwst.atts,
                     gz=gz.atts,
                     huss=huss.atts,
                     insol=insol.atts,
                     irflux=irflux.atts,
                     latflux=latflux.atts,
                     ps=ps.atts,
                     pr=pr.atts,
                     tas=tas.atts,
                     tasmax=tasmax.atts,
                     tasmin=tasmin.atts,       
                     psl=psl.atts,
                     rhs=rhs.atts,
                     runoff=runoff.atts,
                     snd=snd.atts,
                     snc=snc.atts,
                     spw=spw.atts,
                     streamflow=streamflow.atts,
                     swater=swater.atts,                     
                     sice=sice.atts,
                     smass=smass.atts,
                     smelt=smelt.atts,
                     swlake=swlake.atts,
                     swriver=swriver.atts,
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

