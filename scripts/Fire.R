#fire
library(sf)
library(rgdal)
setwd("~/Desktop/Fire")
FireNSW <- st_read("NPWSFireHistory_18022022.shp")
fire <- readOGR(dsn = '.')
plot(fire)
fire

fire$Year <- substring(fire$Label, 1, 4)
levels(as.factor(fire$Label))
plot(fire[which(fire$Year==2020),])
fire2020 <- fire[which]
fire["Label"== "1970"])



setwd("~/Desktop/Fire/burnedarea")
?list()
stack(c(list.files("~/Desktop/Fire/burnedarea")))
firelist
firestack <- stack(firelist)
fire1 <- "burnedmosaic20200101.tif"
fire2 <- "burnedmosaic20200201.tif"
firelist1 <- c(fire1, fire2)
stack(firelist1)
brick(firelist1)
stack(fire1)
raster(fire1)
raster("~/Desktop/Fire/burnedarea/burnedmosaic20200201.tif")
brick("~/Desktop/Fire/burnedarea/burnedmosaic20200201.tif")
firelist1 <- c("~/Desktop/Fire/burnedarea/burnedmosaic20200201.tif", "~/Desktop/Fire/burnedarea/burnedmosaic20200101.tif")
firelist1
stack(firelist1)
plot(firestack)
dim(firestack)
attributes(firestack)

library(raster)
library(sp)
library(rgdal)
fire <- raster("cropped_bushavg.tif")
extent(fire)
extent_aus <- extent(c(110, 158, -46, -9))
writeRaster(fire, "bushfire_averaged.asc")
