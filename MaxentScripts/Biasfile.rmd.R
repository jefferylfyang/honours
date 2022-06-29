library(dismo) # interface with MaxEnt
library(raster) # spatial data manipulation
library(MASS) # for 2D kernel density function
library(magrittr) # for piping functionality, i.e., %>%
library(maptools) # reading shapefiles
library(rgeos)
library(rJava)

locations = data.frame()
record.list <- list.files("/Users/u6638201/Desktop/Maxent/pass_occurence_data/", pattern = ".csv", full.names = T)  

for (i in 1:7) {
  temp <- fread(record.list[i], select = c("decimalLongitude", "decimalLatitude"))
  locations <- rbind(temp, locations)
}

locations <- na.omit(locations)
colnames(locations) <- c("lon", "lat")
climdat <- brick("clim/wc2.1_2.5m_bio_1.tif")
occur.ras <- rasterize(locations, climdat, 1)
plot(occur.ras)

extent_aus <- c(110, 158, -46, -9)
occur.ras<- crop(occur.ras, extent_aus)
plot(occur.ras)

presences <- which(values(occur.ras) == 1)
pres.locs <- coordinates(occur.ras)[presences, ]

dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(nrow(occur.ras), ncol(occur.ras)))
dens.ras <- raster(dens)
plot(dens.ras)

writeRaster(dens.ras, "biasfile.tif")
