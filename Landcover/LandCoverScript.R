library(ncdf4)
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(maptools)

setwd("~/Desktop/LandCoverData")
ld_management <- nc_open("LUH2_GCB2019_management.nc4")
ld_states <- nc_open("LUH2_GCB2019_states.nc4")
attributes(ld_states$var)
landat <- "LUH2_GCB2019_states.nc4"
variables <- attributes(ld_states$var)





# Save the print(nc) dump to a text file
{
  sink('LUH2_GCB2019_states.txt')
  print(ld_states)
  sink()
}
lon <- ncvar_get(ld_states, "lon")
lon[lon > 180] <- lon[lon > 180] - 360
lat <- ncvar_get(ld_states, "lat", verbose = F)
t <- ncvar_get(ld_states, "time")


fillvalue <- ncatt_get(ld_states, "primf", "_FillValue")
fillvalue
head(ndvi.array)
dim(ndvi.array)
nc_close(ld_states) 
class(ndvi.array)
head(ld_states)


ndvi.array[ndvi.array == fillvalue$value] <- NA
ndvi.slice <- ndvi.array[60 , ,] 
image(ndvi.slice)

r <- raster(t(secdf.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
r <- flip(r, direction='y')
plot(r)



data(wrld_simpl)
plot(wrld_simpl, add = TRUE)

# plot a date
ndvi.slice <- ndvi.array[1111, ,] 
r <- raster(t(ndvi.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
r <- flip(r, direction='y')
plot(r)
data(wrld_simpl)
plot(wrld_simpl, add = TRUE)

#change variable

secdf.array <- ncvar_get(ld_states, "secdf")
fillvalue <- ncatt_get(ld_states, "secdf", "_FillValue")
secdf.array[secdf.array == fillvalue$value] <- NA
secdf.slice <- secdf.array[1150, ,]
r <- raster(t(secdf.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
r <- flip(r, direction='y')
data(wrld_simpl)
extent_aus <- c(110, 158, -46, -9)
zoomed <- crop(r, extent_aus)
png("~/Desktop/LandCoverData/Plots/plot.png")
plot(zoomed)
plot(wrld_simpl, add = TRUE)
dev.off()

pastr.array <- ncvar_get(ld_states, "pastr")

pastr.slice <- pastr.array[, ,]
r <- raster(t(pastr.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
r <- flip(r, direction='y')
plot(r)
data(wrld_simpl)
plot(wrld_simpl, add = TRUE)
dim(pastr.slice)
?plot

# zoom and save to png
pt <- cbind(133,-27)
plot(r); points(pt)
s <- 25
e <- extent(pt[1]-s, pt[1]+s, pt[2]-s, pt[2]+s)
zoom(r, e)
plot(wrld_simpl, add = TRUE)
dev.copy(png, 'plot1')
dev.off()

#write function to make a plot of a certain year and of a certain variable
vsz <- function(var, year) {
  year1 <- year - 849
  var.array <- ncvar_get(ld_states, var)
  fillvalue <- ncatt_get(ld_states, var, "_FillValue")
  var.array[var.array == fillvalue$value] <- NA
  var.slice <- var.array[year1, ,]
  r <- raster(t(var.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  r <- flip(r, direction='y')
  extent_aus <- c(110, 158, -46, -9)
  zoomed <- crop(r, extent_aus)
  png(paste0("~/Desktop/LandCoverData/Plots/",var,year),
      width = 480, height = 480,
      units = "px", pointsize = 12, bg = "white", res = NA)
      #restoreConsole = TRUE)
  plot(zoomed)
  data(wrld_simpl)
  plot(wrld_simpl, add = TRUE)
  dev.off()
}


primf.array <- ncvar_get(ld_states, "primf")
vsz("pastr",1980)


vsz("range", 2015)
vsz("range", 1980)



ycon(2015)

# location, variable name, value, year, birdID,  

secdf.array[3, 600, 10]


#extracting variables at a point
?extract 
rasValue=extract(r,p)
p <- matrix(c(48, -55, 5, 6, 7, 3), nrow = 3, ncol = 2, byrow = T, dimnames = NULL)

plot(r)
rasValue
pointCord <- read.csv("Point_Data.csv")
rasValue=extract(r, pointCord)
coordinates(pointCord)= ~ LON + LAT
combinePointValue=cbind(pointCord,rasValue)
write.table(rasValue,file="combinedPointValue.csv", append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)
extract(pastrbrick, p, layer = 500, nl = 2)
r <- brick(t(pastr.array), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
rc <- read.csv("RandomCord.csv")
extract(pastrbrick, rc, layer = 1131)
rvalues <- extract(pastrbrick, rc)
write.table(rvalues,file="combinedPointValue.csv", append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)


CordTestTable <- cbind(rbind(c(1,2,3),c(2,2,4), c(3,1,7))) #test coordinate table 
secdf.array[CordTestTable] #extract set of data values using coordinates in table from array 

RandomPointsInArray <- randomPoints(r,100000)
write.csv(RandomPointsInArray,"ThousRandomPointsInArray.csv")
RanCordTable <- read.csv("ThousRandPoints.csv")
RanCordTable <- RanCordTable %>% 
  mutate_if(is.numeric, abs)
RanCordTable <- data.matrix(RanCordTable, rownames.force = NA)
write.csv(secdf.array[RanCordTable], "RandCordTable.csv")
