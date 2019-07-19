##Make list of EPW coordinates

epw.dir <- '/storage/data/projects/rci/weather_files/wx_files/'

epw.files <- list.files(path=epw.dir,pattern='CWEC.epw')

coords <- matrix(0,nrow=length(epw.files),ncol=2)
sites <- rep('A',length(epw.files))

for (i in seq_along(epw.files)) {
    present.epw.file <- epw.files[i]
    print(present.epw.file)
    pef.split <- strsplit(present.epw.file,'_')[[1]]
    print(pef.split[3])
    epw.present <- read.epw.file(epw.dir,present.epw.file)
    lat <- as.numeric(strsplit(epw.present$header[1],',')[[1]][7])
    lon <- as.numeric(strsplit(epw.present$header[1],',')[[1]][8])
    print(lon)
    print(lat)
    coords[i,1] <- lon
    coords[i,2] <- lat
    sites[i] <- pef.split[3]
}
    
output <- cbind(sites,coords)
write.table(output,file='/storage/data/projects/rci/weather_files/epw.sites.csv',sep=',',quote=F,row.name=F,col.name=F)