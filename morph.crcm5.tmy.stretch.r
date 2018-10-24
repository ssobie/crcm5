##Script to morph the CRCM5 TMY weather years following the same
 ##"shift and stretch" procedure outlined in the Belcher paper.

library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##-----------------------------------------------

get.rcm.alpha <- function(var.name,rcm.file,past.int,proj.int,std=FALSE) {
   nc <- nc_open(rcm.file)
   time <- netcdf.calendar(nc)
   monthly.fac <- as.factor(format(time,'%m'))
   daily.fac <- as.factor(format(time,'%m-%d'))
   yrs <- levels(as.factor(format(time,'%Y')))
   weekly.fac <- as.factor(rep(rep(1:73,each=5),length(yrs)))
   past.bnds <- strsplit(past.int,'-')[[1]]
   proj.bnds <- strsplit(proj.int,'-')[[1]]

   past.st <- head(grep(substr(past.bnds[1],1,4),time),1)
   past.en <- tail(grep(substr(past.bnds[2],1,4),time),1)
   past.cnt <- past.en-past.st+1
   past.data <- ncvar_get(nc,var.name,start=c(1,1,past.st),count=c(-1,-1,past.cnt))      
   past.mean <- apply(past.data,c(1,2),function(x,y){tapply(x,y,mean)},monthly.fac[past.st:past.en])
   ##past.mean <- apply(past.data,c(1,2),function(x,y){tapply(x,y,mean)},daily.fac[past.st:past.en])
   ##past.mean <- apply(past.data,c(1,2),function(x,y){tapply(x,y,mean)},weekly.fac[past.st:past.en])
   past.sd <- apply(past.data,c(1,2),function(x,y){tapply(x,y,sd)},weekly.fac[past.st:past.en])

   proj.st <- head(grep(substr(proj.bnds[1],1,4),time),1)
   proj.en <- tail(grep(substr(proj.bnds[2],1,4),time),1)
   proj.cnt <- proj.en-proj.st+1
   proj.data <- ncvar_get(nc,var.name,start=c(1,1,proj.st),count=c(-1,-1,proj.cnt))
   proj.mean <- apply(proj.data,c(1,2),function(x,y){tapply(x,y,mean)},monthly.fac[proj.st:proj.en])
   ##proj.mean <- apply(proj.data,c(1,2),function(x,y){tapply(x,y,mean)},daily.fac[proj.st:proj.en])
   ##proj.mean <- apply(proj.data,c(1,2),function(x,y){tapply(x,y,mean)},weekly.fac[proj.st:proj.en])
   proj.sd <- apply(proj.data,c(1,2),function(x,y){tapply(x,y,sd)},weekly.fac[proj.st:proj.en])
   delta <- aperm(proj.mean - past.mean,c(2,3,1))
   alpha <- aperm(proj.mean / past.mean,c(2,3,1))
   if (std) {
     alpha <- aperm(proj.sd/past.sd,c(2,3,1))
   }
   nc_close(nc)

   rv <- list(delta=delta,alpha=alpha)
   return(rv)                
}


##Morph by Stretch

morph.variable <- function(var.name,gcm,past.int,proj.int,
                           tmy.dir,rcm.dir) {

   tmy.file <- paste0(tmy.dir,var.name,"_CWEC_TMY_",gcm,"+CRCM5_historical_",past.int,".nc")
   new.tmy.file <- paste0(tmy.dir,"morphed_",var.name,"_CWEC_TMY_",gcm,"+CRCM5_historical+rcp85_",proj.int,".nc")

   file.copy(from=tmy.file,to=new.tmy.file,overwrite=T)
   tmy.nc <- nc_open(tmy.file)           
   new.nc <- nc_open(new.tmy.file,write=T)

   xl <- tmy.nc$dim$rlon$len
   yl <- tmy.nc$dim$rlat$len

   tmy.data <- ncvar_get(tmy.nc,var.name)
   tmy.time <- netcdf.calendar(tmy.nc)
   
   tmy.daily.fac <- as.factor(format(tmy.time,'%Y-%m-%d'))
   tmy.daily.time <- as.Date(levels(tmy.daily.fac))
   tmy.monthly.fac <- as.factor(format(tmy.daily.time,'%m'))
   
   tmy.weekly.fac <- as.factor(rep(1:73,each=40)) ##Divide the year into 73 intervals of 5 days each

   tmy.daily.data <- aperm(apply(tmy.data,c(1,2),function(x,y){tapply(x,y,mean)},tmy.daily.fac),c(2,3,1))

   tmy.weekly.data <- aperm(apply(tmy.data,c(1,2),function(x,y){tapply(x,y,mean)},tmy.weekly.fac),c(2,3,1))

   tmy.monthly.data <- aperm(apply(tmy.daily.data,c(1,2),function(x,y){tapply(x,y,mean)},tmy.monthly.fac),c(2,3,1))

   ##-------------------

   rcm.files <- list.files(path=rcm.dir,pattern=gcm,full.name=TRUE)
   if (var.name=='pr') {
     var.file <- rcm.files[grep(paste0(var.name,'_sum'),rcm.files)] 
   } else {
     var.file <- rcm.files[grep(paste0(var.name,'_mean'),rcm.files)] 
   }

   var.change <- get.rcm.alpha(var.name,var.file,past.int,proj.int)
   var.alpha <- var.change$alpha 
   var.delta <- var.change$delta 
   var.sd <-  get.rcm.alpha(var.name,var.file,past.int,proj.int,std=TRUE)$alpha

   tmy.monthly.fac <- as.factor(format(tmy.time,'%m'))
   tmy.daily.fac <- as.factor(format(tmy.time,'%j'))
   tmy.morphed.data <- tmy.data*0

   ##Monthly Adjustment
   for (mn in 1:12) {
      ix <- tmy.monthly.fac == sprintf('%02d',mn)
      tmp.delta <- array(var.delta[,,mn],c(xl,yl,sum(ix)))
      tmp.month <- array(tmy.monthly.data[,,mn],c(xl,yl,sum(ix)))

      if (var.name=='dewpoint') {
        tmp.alpha <- array(var.sd[,,mn],c(xl,yl,sum(ix)))
        tmy.morphed.data[,,ix] <- tmp.month + tmp.delta + tmp.alpha * (tmy.data[,,ix] - tmp.month)             
      } else {
        tmp.alpha <- array(var.alpha[,,mn],c(xl,yl,sum(ix)))
        tmy.morphed.data[,,ix] <- tmy.data[,,ix] * tmp.alpha
      }
   }


##   ##Weekly Adjustment
##   for (wk in 1:73) {
##      ix <- tmy.weekly.fac == wk
##      tmp.delta <- array(tas.delta[,,wk],c(xl,yl,sum(ix)))
##      tmp.alpha <- array(tas.alpha[,,wk],c(xl,yl,sum(ix)))
##      tmp.week <- array(tmy.weekly.tas[,,wk],c(xl,yl,sum(ix)))
##      tmy.morphed.tas[,,ix] <- tmy.tas[,,ix] + tmp.delta + tmp.alpha * (tmy.tas[,,ix] - tmp.week)             
##   }

   ##Daily Adjustment
##   for (dy in 1:365) {
##      ix <- tmy.daily.fac == paste0(sprintf('%03d',dy))
##      tmp.delta <- array(tas.delta[,,dy],c(xl,yl,sum(ix)))
##      tmp.alpha <- array(tas.alpha[,,dy],c(xl,yl,sum(ix)))
##      tmp.day <- array(tmy.daily.tas[,,dy],c(xl,yl,sum(ix)))
##      tmy.morphed.tas[,,ix] <- tmy.tas[,,ix] + tmp.delta + tmp.alpha * (tmy.tas[,,ix] - tmp.day)             
##   }

   ##Daily Adjustment with sigma alpha
##   for (dy in 1:365) {
##      ix <- tmy.daily.fac == paste0(sprintf('%03d',dy))
##      tmp.delta <- array(tas.delta[,,dy],c(xl,yl,sum(ix)))
##      tmp.alpha <- array(tas.sd[,,dy],c(xl,yl,sum(ix)))
##      tmp.day <- array(tmy.daily.tas[,,dy],c(xl,yl,sum(ix)))
##      tmy.morphed.tas[,,ix] <- tmp.day + tmp.delta + tmp.alpha * (tmy.tas[,,ix] - tmp.day)             
##   }
   
   ncvar_put(new.nc,var.name,tmy.morphed.data)
         
   nc_close(tmy.nc)
   nc_close(new.nc)
browser()
}

##Morph Dewpoint


##Morph by stretching
##(Pr,PSL,Wspd,Insol)       


##-----------------------------------------------
##Generate Future TMY File

var.names <- c('dewpoint')##,'rhs','wspd')
gcms <- c('CanESM2')
past.int <- '19810101-20101231'
proj.int <- '20210101-20501231'

base.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/bc/'
tmy.dir <- paste0(base.dir,'tmy_files/')
rcm.dir <- paste0(base.dir,'daily/')
for (gcm in gcms) {
  for (var.name in var.names) {
    morph.variable(var.name,gcm,past.int,proj.int,tmy.dir,rcm.dir)
  }
}
   