##Script to plot the EPW file series

library(zoo)

##------------------------------------------------------------------------------
make_average_series <- function(series,time,method,rlen=5,agg.fxn) {
  method <- 'monthly'
  factor <- switch(method,
                   daily='%m-%d',
                   monthly='%m',
                   roll='%m-%d',
                   seasonal='%m')
  fac <- as.factor(format(time,factor))

  if (method=='daily') {
    agg.daily <- tapply(series,fac,agg.fxn)
  }

  if (method=='monthly') {
     dtime <- as.Date(unique(format(time,'%Y-%m-%d')))
     months <- format(dtime,'%m')
     agg.data <- tapply(series,fac,agg.fxn)   
     agg.daily <- rep(0,365)
     for (mn in 1:12) {
       ix <- months == sprintf('%02d',mn)
       agg.daily[ix] <- agg.data[mn]
     }
  }

  if (method=='seasonal') {
     seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
     fac <- factor(seasons[fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))     
     agg.data <- tapply(series,fac,agg.fxn)   
     daily.seas <- c(rep('DJF',59),rep('MAM',92),rep('JJA',92),rep('SON',91),rep('DJF',31))
     agg.daily <- rep(0,365)
     agg.daily[daily.seas=='DJF'] <- agg.data['DJF']
     agg.daily[daily.seas=='MAM'] <- agg.data['MAM']
     agg.daily[daily.seas=='JJA'] <- agg.data['JJA']
     agg.daily[daily.seas=='SON'] <- agg.data['SON']
  }

  if (method=='roll') {
    agg.daily <- rollmean(agg.data,rlen,fill='extend')
  }    

  return(agg.data)
}

##------------------------------------------------------------------------------
##Separate stretching morphing functions

morph.relative.humidity <- function(epw.series,alpha) {

      if (method=='daily' | method =='roll') {
        for (d in 1:365) {
           ix <- format(dates,'%j') == sprintf('%03d',d)  
           morphed.series[g,ix] <- morphing.fxn(epw.series[ix],alpha[d]) ##epw.series[ix] * alpha[d]
        }
      }

}


##------------------------------------------------------------------------------
##Morph Dry Bulb Temperature using Belcher method

morph_dry_bulb_temp <- function(epw.present,lon,lat,gcm.list,gcm.dir,scenario,interval,
                                method,rlen=NULL,agg.fxn=mean) {

   tas.ix <- get_field_index('dry_bulb_temperature')
   epw.tas <- epw.present$data[,tas.ix]
   dates <- as.Date(paste('1999',sprintf('%02d',epw.present$data[,2]),sprintf('%02d',epw.present$data[,3]),sep='-'))

   epw.agg.mean <- epw.daily.mean <- make_average_series(epw.tas,dates,method='daily',agg.fxn=mean)
   epw.agg.max <- epw.daily.max <-  make_average_series(epw.tas,dates,method='daily',agg.fxn=max)
   epw.agg.min <- epw.daily.min <-  make_average_series(epw.tas,dates,method='daily',agg.fxn=min)

   epw.day.anoms <- epw.tas*0
   ##Daily
   for (d in 1:365) {
      ix <- format(dates,'%j') == sprintf('%03d',d)
      epw.day.anoms[ix] <- epw.tas[ix] - epw.daily.mean[d]
   }
   epw.agg.anoms <- epw.day.anoms

   epw.agg.mean <- make_average_series(epw.tas,dates,method=method,agg.fxn=mean)
   epw.agg.max <-  make_average_series(epw.tas,dates,method=method,agg.fxn=max)
   epw.agg.min <-  make_average_series(epw.tas,dates,method=method,agg.fxn=min)

   if (method == 'monthly') {
     for (d in seq_along(epw.agg.mean)) {
        ix <- format(dates,'%m') == sprintf('%02d',d)
        epw.agg.anoms[ix] <- epw.tas[ix] - epw.agg.mean[d]
     }
   }
   if (method == 'roll') {
     for (d in seq_along(epw.agg.mean)) {
        ix <- format(dates,'%j') == sprintf('%03d',d)
        epw.agg.anoms[ix] <- epw.tas[ix] - epw.agg.mean[d]
     }     
   }
   if (method == 'seasonal') {
     mn.fac <- as.factor(format(dates,'%m'))
     seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
     seas.fac <- factor(seasons[mn.fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))      
     for (d in seq_along(epw.agg.mean)) {
        ix <- seas.fac == (c('DJF', 'MAM', 'JJA', 'SON')[d])
        print(c('DJF', 'MAM', 'JJA', 'SON')[d])
        epw.agg.anoms[ix] <- epw.tas[ix] - epw.agg.mean[d]
     }
   }

   tlen <- length(epw.agg.mean)
   
   deltas <- matrix(0,nrow=length(gcm.list),ncol=tlen)
   alphas <- matrix(0,nrow=length(gcm.list),ncol=tlen)

   morphed.tas <- matrix(0,nrow=length(gcm.list),ncol=length(epw.tas*0))

   for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      scen.files <- list.files(path=paste0(gcm.dir,gcm),pattern=scenario)
      var.files <- scen.files[grep("tasmax",scen.files)]
      past.tx.file <- var.files[grep("1951-2000",var.files)]
      proj.tx.file <- var.files[grep("2001-2100",var.files)]

      past.tx <- sub_by_time(var.name='tasmax',lonc=lon,latc=lat,
                                  interval='1971-2000',
                                  input.file=past.tx.file,gcm=gcm,read.dir=gcm.dir)
      past.tx.agg <- make_average_series(past.tx$data,past.tx$time,
                                         method,rlen,agg.fxn)

      proj.tx <- sub_by_time(var.name='tasmax',lonc=lon,latc=lat,
                                  interval=interval,
                                  input.file=proj.tx.file,gcm=gcm,read.dir=gcm.dir)
      proj.tx.agg <- make_average_series(proj.tx$data,proj.tx$time,
                                         method,rlen,agg.fxn)

      var.files <- scen.files[grep("tasmin",scen.files)]
      past.tn.file <- var.files[grep("1951-2000",var.files)]
      proj.tn.file <- var.files[grep("2001-2100",var.files)]

      past.tn <- sub_by_time(var.name='tasmin',lonc=lon,latc=lat,
                             interval='1971-2000',
                             input.file=past.tn.file,gcm=gcm,read.dir=gcm.dir)
      past.tn.agg <- make_average_series(past.tn$data,past.tn$time,
                                         method,rlen,agg.fxn)

      proj.tn <- sub_by_time(var.name='tasmin',lonc=lon,latc=lat,
                             interval=interval,
                             input.file=proj.tn.file,gcm=gcm,read.dir=gcm.dir)
      proj.tn.agg <- make_average_series(proj.tn$data,proj.tn$time,
                                         method,rlen,agg.fxn)

      ##
      delta_tasmax <- proj.tx.agg - past.tx.agg
      delta_tasmin <- proj.tn.agg - past.tn.agg
      delta_tas <- (proj.tx.agg+proj.tn.agg)/2 - (past.tx.agg+past.tn.agg)/2
      deltas[g,] <- delta_tas
      past_tas <- (past.tx.agg+past.tn.agg)/2
      print('Delta')
      print(delta_tas)

      alpha <- (delta_tasmax - delta_tasmin) / (epw.agg.max - epw.agg.min)
      print('Alpha')
      print(alpha)        
      alphas[g,] <- alpha 

      if (method=='daily' | method =='roll') {
        for (d in 1:tlen) {
           ix <- format(dates,'%j') == sprintf('%03d',d)  
           morphed.tas[g,ix] <- epw.tas[ix] + delta_tas[d] + alpha[d]*epw.agg.anoms[ix]
           ##morphed.tas[g,ix] <- deltas[g,m] + alphas[g,m]*epw.agg.anoms[ix]
        }
      }

      if (method=='monthly') {
        for (m in 1:tlen) {
           ix <- format(dates,'%m') == sprintf('%02d',m)  
           morphed.tas[g,ix] <- epw.tas[ix] + delta_tas[m] + alpha[m]*epw.agg.anoms[ix]
           ##morphed.tas[g,ix] <- deltas[g,m] + alphas[g,m]*epw.agg.anoms[ix]
        }
      }
      if (method=='seasonal') {
        for (s in 1:tlen) {
           ix <- seas.fac == (c('DJF', 'MAM', 'JJA', 'SON')[d])
           morphed.tas[g,ix] <- epw.tas[ix] + delta_tas[s] + alpha[s]*epw.agg.anoms[ix]
        }
      }



   }##GCM loop

   ens.deltas <- apply(deltas,2,mean)
   ens.alphas <- apply(alphas,2,mean)
   ens.morphed.tas <- apply(morphed.tas,2,mean)
   epw.present$data[,tas.ix] <- round(ens.morphed.tas,1)
   return(epw.present)   
}

##------------------------------------------------------------------------------
##Morph Dewpoint Temperature using Belcher method

morph_dew_point_temp <- function(epw.present,lon,lat,gcm.list,gcm.dir,scenario,interval,
                                 method,rlen=NULL,agg.fxn=mean) {

   dwpt.ix <- get_field_index('dew_point_temperature')
   epw.dwpt <- epw.present$data[,dwpt.ix]
   dates <- as.Date(paste('1999',sprintf('%02d',epw.present$data[,2]),sprintf('%02d',epw.present$data[,3]),sep='-'))

   epw.agg.mean <- epw.daily.mean <- make_average_series(epw.dwpt,dates,method='daily',agg.fxn=mean)
   epw.day.anoms <- epw.dwpt*0
   ##Daily
   for (d in 1:365) {
      ix <- format(dates,'%j') == sprintf('%03d',d)
      epw.day.anoms[ix] <- epw.dwpt[ix] - epw.daily.mean[d]
   }
   epw.agg.anoms <- epw.day.anoms

   if (method == 'monthly') {
     epw.agg.mean <- make_average_series(epw.dwpt,dates,method='monthly',agg.fxn=mean)
     for (d in seq_along(epw.agg.mean)) {
        ix <- format(dates,'%m') == sprintf('%02d',d)
        epw.agg.anoms[ix] <- epw.dwpt[ix] - epw.agg.mean[d]
     }
   }
   if (method == 'roll') {
     epw.agg.mean <- make_average_series(epw.dwpt,dates,method='roll',rlen=rlen,agg.fxn=mean)
     for (d in 1:365) {
        ix <- format(dates,'%j') == sprintf('%03d',d)
        epw.agg.anoms[ix] <- epw.dwpt[ix] - epw.agg.mean[d]
     }     
   }

   tlen <- length(epw.agg.mean)
   
   deltas <- matrix(0,nrow=length(gcm.list),ncol=tlen)
   alphas <- matrix(0,nrow=length(gcm.list),ncol=tlen)

   morphed.dwpt <- matrix(0,nrow=length(gcm.list),ncol=length(epw.dwpt*0))

   for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      scen.files <- list.files(path=paste0(gcm.dir,gcm),pattern=scenario)
      var.files <- scen.files[grep("dewpoint",scen.files)]
      dwpt.file <- var.files[grep("19500101-21001231",var.files)]
      past.dwpt <- sub_by_time(var.name='dewpoint',lonc=lon,latc=lat,
                               interval='1971-2000',
                               input.file=dwpt.file,gcm=gcm,read.dir=gcm.dir)
      past.dwpt.agg <- make_average_series(past.dwpt$data,past.dwpt$time,
                                         method,rlen,agg.fxn)
      past.dwpt.sd <- make_average_series(past.dwpt$data,past.dwpt$time,
                                          method,rlen,sd)

      proj.dwpt <- sub_by_time(var.name='dewpoint',lonc=lon,latc=lat,
                                  interval=interval,
                                  input.file=dwpt.file,gcm=gcm,read.dir=gcm.dir)
      proj.dwpt.agg <- make_average_series(proj.dwpt$data,proj.dwpt$time,
                                         method,rlen,agg.fxn)
      proj.dwpt.sd <- make_average_series(proj.dwpt$data,proj.dwpt$time,
                                          method,rlen,sd)

      ##
      delta_dwpt <- proj.dwpt.agg - past.dwpt.agg
      deltas[g,] <- delta_dwpt

      alpha <- proj.dwpt.sd / past.dwpt.sd

      print('Alpha')
      print(alpha)        
      alphas[g,] <- alpha 

      if (method=='daily' | method =='roll') {
        for (d in 1:tlen) {
           ix <- format(dates,'%j') == sprintf('%03d',d)  
           morphed.dwpt[g,ix] <- epw.agg.mean[d] + delta_dwpt[d] + alpha[d]*epw.agg.anoms[ix]
        }
      }

      if (method=='monthly') {
        for (m in 1:tlen) {
           ix <- format(dates,'%m') == sprintf('%02d',m)  
           morphed.dwpt[g,ix] <- epw.agg.mean[m] + delta_dwpt[m] + alpha[m]*epw.agg.anoms[ix]
        }
      }
   }##GCM loop

   ens.deltas <- apply(deltas,2,mean)
   ens.alphas <- apply(alphas,2,mean)
   ens.morphed.dwpt <- apply(morphed.dwpt,2,mean)
   epw.present$data[,dwpt.ix] <- round(ens.morphed.dwpt,1)
   return(epw.present)   
}

##------------------------------------------------------------------------------
generate_horizontal_radiation <- function(epw.present,lon,lat,gcm.list,gcm.dir,scenario,interval,
                                          method,rlen=NULL,agg.fxn=mean) {

   ghr.ix <- get_field_index('global_horizontal_radiation')
   dhr.ix <- get_field_index('diffuse_horizontal_radiation')
   epw.global <- epw.present$data[,ghr.ix]
   epw.diffuse <- epw.present$data[,dhr.ix]
   flag <- epw.global == 0
   diffuse.to.global.ratio <- epw.diffuse / epw.global
   diffuse.to.global.ratio[flag] <- 0

   dates <- as.Date(paste('1999',sprintf('%02d',epw.present$data[,2]),sprintf('%02d',epw.present$data[,3]),sep='-'))
   
   morphed.global <- matrix(0,nrow=length(gcm.list),ncol=length(epw.global*0))
   morphed.diffuse <- matrix(0,nrow=length(gcm.list),ncol=length(epw.diffuse*0))

   gcm.var <- 'rsds'

   for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      
      scen.files <- list.files(path=paste0(gcm.dir,gcm),pattern=scenario)
      var.files <- scen.files[grep(gcm.var,scen.files)]
      gcm.file <- var.files[grep("19500101-21001231",var.files)]
      past.gcm <- sub_by_time(var.name=gcm.var,lonc=lon,latc=lat,
                               interval='1971-2000',
                               input.file=gcm.file,gcm=gcm,read.dir=gcm.dir)
      past.gcm.agg <- make_average_series(past.gcm$data,past.gcm$time,
                                         method,rlen,agg.fxn)

      proj.gcm <- sub_by_time(var.name=gcm.var,lonc=lon,latc=lat,
                                  interval=interval,
                                  input.file=gcm.file,gcm=gcm,read.dir=gcm.dir)
      proj.gcm.agg <- make_average_series(proj.gcm$data,proj.gcm$time,
                                          method,rlen,agg.fxn)

      ##
      alpha <- proj.gcm.agg / past.gcm.agg

      print('Alpha')
      print(alpha)        

      if (method=='daily' | method =='roll') {
        for (d in 1:365) {
           ix <- format(dates,'%j') == sprintf('%03d',d)  
           morphed.global[g,ix] <- epw.global[ix] * alpha[d]
           morphed.diffuse[g,ix] <- (epw.global[ix] * alpha[d]) * diffuse.to.global.ratio[ix]           
        }
      }

      if (method=='monthly') {
        for (m in 1:12) {
           ix <- format(dates,'%m') == sprintf('%02d',m)  
           morphed.global[g,ix] <- epw.global[ix] * alpha[d]
           morphed.diffuse[g,ix] <- (epw.global[ix] * alpha[d]) * diffuse.to.global.ratio[ix]           
        }
      }

      if (method=='seasonal') {
        monthly.fac <- as.factor(format(dates,'%m'))
        seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
        seasonal.fac <- factor(seasons[monthly.fac], levels=c('DJF', 'MAM', 'JJA', 'SON')) 
        for (s in 1:4) {
           ix <- seasonal.fac == unique(seasons)[s]  
           morphed.global[g,ix] <- epw.global[ix] * alpha[d]
           morphed.diffuse[g,ix] <- (epw.global[ix] * alpha[d]) * diffuse.to.global.ratio[ix]           
        }
      }
   }##GCM loop

   ens.morphed.global <- apply(morphed.global,2,mean)
   ens.morphed.diffuse <- apply(morphed.diffuse,2,mean)
   epw.present$data[,ghr.ix] <- round(ens.morphed.global,0)
   epw.present$data[,dhr.ix] <- round(ens.morphed.diffuse,0)
   return(epw.present)   
}

##------------------------------------------------------------------------------


##------------------------------------------------------------------------------
##Morph by stretching using Belcher method

generate_stretched_series <- function(epw.present,epw.var,gcm.var,lon,lat,gcm.list,gcm.dir,scenario,interval,
                                      method,rlen=NULL,agg.fxn=mean) {

   epw.ix <- get_field_index(epw.var)
   morph.fxn <- get.morph.function(epw.var)
   epw.series <- epw.present$data[,epw.ix]
   dates <- as.Date(paste('1999',sprintf('%02d',epw.present$data[,2]),sprintf('%02d',epw.present$data[,3]),sep='-'))
   
   morphed.series <- matrix(0,nrow=length(gcm.list),ncol=length(epw.series*0))

   for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      
      scen.files <- list.files(path=paste0(gcm.dir,gcm),pattern=scenario)
      var.files <- scen.files[grep(gcm.var,scen.files)]
      gcm.file <- var.files[grep("19500101-21001231",var.files)]
      past.gcm <- sub_by_time(var.name=gcm.var,lonc=lon,latc=lat,
                               interval='1971-2000',
                               input.file=gcm.file,gcm=gcm,read.dir=gcm.dir)
      past.gcm.agg <- make_average_series(past.gcm$data,past.gcm$time,
                                         method,rlen,agg.fxn)

      proj.gcm <- sub_by_time(var.name=gcm.var,lonc=lon,latc=lat,
                                  interval=interval,
                                  input.file=gcm.file,gcm=gcm,read.dir=gcm.dir)
      proj.gcm.agg <- make_average_series(proj.gcm$data,proj.gcm$time,
                                          method,rlen,agg.fxn)

      ##
      alpha <- proj.gcm.agg / past.gcm.agg

      if (method=='daily' | method =='roll') {
        for (d in 1:365) {
           ix <- format(dates,'%j') == sprintf('%03d',d)  
           morphed.series[g,ix] <- morphing.fxn(epw.series[ix],alpha[d]) ##epw.series[ix] * alpha[d]
        }
      }

      if (method=='monthly') {
        for (m in 1:12) {
           ix <- format(dates,'%m') == sprintf('%02d',m)  
           morphed.series[g,ix] <- epw.series[ix] * alpha[m]
        }
      }

      if (method=='seasonal') {
        monthly.fac <- as.factor(format(dates,'%m'))
        seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
        seasonal.fac <- factor(seasons[monthly.fac], levels=c('DJF', 'MAM', 'JJA', 'SON')) 
        browser()
        for (s in 1:4) {
           ix <- seasonal.fac == unique(seasons)[s]  
           morphed.series[g,ix] <- epw.series[ix] * alpha[m]
        }
      }
   }##GCM loop

##   ens.deltas <- apply(deltas,2,mean)
##   ens.alphas <- apply(alphas,2,mean)
   ens.morphed.series <- apply(morphed.series,2,mean)
   ens.morphed.series[ens.morphed.series > 100] <- 100
   epw.present$data[,epw.ix] <- round(ens.morphed.series,0)
   return(epw.present)   
}

##------------------------------------------------------------------------------

