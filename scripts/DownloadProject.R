library(ncdf4)
library(terra) # package for spatraster manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(maptools)
library(magrittr) #piping
library(dplyr)
library(data.table)
library(climates)
library(ntbox)
library(stringr)

setwd("~/Honours")

#Read in tmin, tmax, rain, tavg, of said year 
varlist <- c('tavg', 'tmax', 'tmin', 'rain')
extent_aus <- ext(112, 154, -44,-9)
months <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", 10:12)

nc.open <- function(month, variable, year) { #Function to get raster with inputs that can be iterated
  a = nc_open(paste0("https://dapds00.nci.org.au/thredds/dodsC/gh70/ANUClimate/v2-0/stable/month/", variable,"/", year, "/ANUClimate_v2-0_", variable, "_monthly_", year, month,".nc"))
  b = ncvar_get(a, variable) %>% t() %>% rast(crs = "WGS84") %>% set.ext(extent_aus) %>% as.data.frame()
  return(b)
}

nc.open.xy <- function(month, variable, year) { #Function to get raster with inputs that can be iterated
  a = nc_open(paste0("https://dapds00.nci.org.au/thredds/dodsC/gh70/ANUClimate/v2-0/stable/month/", variable,"/", year, "/ANUClimate_v2-0_", variable, "_monthly_", year, month,".nc"))
  b = ncvar_get(a, variable) %>% t() %>% rast(crs = "WGS84") %>% set.ext(extent_aus) %>% as.data.frame(xy = T)
  return(b)
}

bioclim.xy <- nc.open.xy("01", "tmax", 1960) 
bioclim.xy <- cbind(bioclim.xy$x, bioclim.xy$y) 
colnames(bioclim.xy) <- c("lon", "lat")
yearlist <- c(1980:2018)
yearlist <- as.character(yearlist)
  


for (i in "1980") {
  bioclimdf.list = list()
  for (ii in varlist) {
    temp.df = data.frame(matrix(nrow=7050282, ncol=0))
    temp.df = lapply(months, nc.open, variable = ii, year = i) %>%
      bind_cols()
    bioclimdf.list[[ii]] = temp.df
  }
  bioclim1 <- bioclim2(bioclimdf.list$tmin, bioclimdf.list$tmax, bioclimdf.list$rain, bioclimdf.list$tavg, files.as.inputs = F)
  bioclim1.dataframe <- cbind(bioclim.xy, bioclim1)
  bioclim1.rast <- rast(bioclim1.dataframe, type = "xyz")
  bioclim1.rast[is.na(bioclim1.rast)] <- -9999
  names.list <- names(bioclim1.rast)
}


for (iii in names.list) {
  writeRaster(subset(bioclim1.rast, iii), paste0("~/Honours/data/proj_layers/temp/", iii,".asc"), NAflag = -9999, overwrite = T)
}


maxent_call(
  "~/Desktop/maxent-2",
  features = c("h","p"),
  environmentallayers = "~/Honours/data/clim_data/ANUBIOCLIM_ascs",
  samplesfile = "~/MaxEnt/BirdData/combined_1960_1979_30km_thin1.csv",
  outputdirectory = "outputs/proj_layers/temp",
  projectionlayers = "data/proj_layers/temp",
  extrapolate = T,
  askoverwrite = T,
  randomtestpoints = 25,
  betamultiplier = 0.5,
  biasfile = "~/MaxEnt/biasfile.asc",
  doclamp = T,
  randomseed = T)

dry_matter <- EPFAST$dry

rast("/Users/u6638201/Honours/outputs/proj_layers/temp/ANUBIOCLIM_ascs_bmult_0.5_hp/Malurus_temp.asc") %>%
writeRaster("outputs.tif")
#make dataframes for that year
#make bioclim variables for that year and turn them into rasters

#Rast cbind(temp.df)r calculate the change between that year and the baseline year 

#workflow that produces a projection of a certain year based on the base model 
#produce folder of 1980 - 2015 projeted rasters. 
#Script that calculates that change over years


BaseModBio.list <- list.files("data/clim_data/ANUBIOCLIM_ascs", pattern = "*.asc$" )

#calculate changes in layers over time 
output.file.list <- list.files("outputs/proj_layers/years", full.names = T)
output.rasters <- list() 
bird.data <- read.csv("data/bird_data/Honours Thesis Metadata v2.csv")
bird.data.trimmed <- cbind(bird.data$Longitude, bird.data$Latitude) %>% na.omit()

bird.data$Year <- str_sub(bird.data$Date, start = -4) #get years from dates format

for (i in output.file.list) {
  output.rasters = c(output.rasters, rast(i))
}

bird_values <- matrix(, nrow = 224, ncol = 0)

for (i in 1:41) {
  bird_values = extract(output.rasters[[i]], bird.data.xy) %>%
     cbind(bird_values) 
}

colnames(bird_values) <- c(1980:2020)
bird_values <- cbind(bird.data.ID, bird_values) %>% na.omit
write.csv(bird_values, "outputs/sampleMaxentValues1980_2020.csv")

for (i in 1:244){
  (bird.data[i, 21] - 20:bird.data[i, 21] - 1) #vector of 20 years before specimen collection excluding collection year
  
}


#plot value variablity against lat, lon, bioclim variables 
bioclim.inputs <- list.files("data/proj_layers/temp", full.names = T, pattern = ".asc$")
bioclim.in.rast <- list() 

for (i in bioclim.inputs) {
  bioclim.in.rast = c(bioclim.in.rast, rast(i))
}


bird.data <- read.csv("data/bird_data/Honours Thesis Metadata v2.csv")
bird.data.ID <- cbind(bird.data$EXP.ID, bird.data$Longitude, bird.data$Latitude) %>% na.omit()
colnames(bird.data.ID) <- c("ID", "lon", "lat")
bird.data.xy <- cbind(bird.data$Longitude, bird.data$Latitude) %>% na.omit()
bird_clim_values <- matrix(, nrow = 224, ncol = 0)

for (i in 1:19) {
  bird_clim_values = extract(bioclim.in.rast[[i]], bird.data.xy) %>%
    cbind(bird_clim_values) 
}

bird_var <- read.csv("bird_values.csv")
bird_var <- data.frame(bird_var$.)

bird_var_plotdata <- cbind(bird.data.ID, bird_clim_values, bird_var)

colnames(bird_var_plotdata[,1]) <- "ID"
c(colnames(bird_clim_values), "lon", 'lat')

for (i in c(colnames(bird_clim_values), "lon", 'lat')){
  png(file = paste0("outputs/1995plot/", i, ".png"))
  plot(cbind(bird_var_plotdata[i], bird_var_plotdata$bird_var..), xlab = i, ylab = "sd")
  dev.off
}

#maxent landcover layers 

varlist <- c("pastr", "primf", "primn", "range", "secdf", "secdn","secma", "secmb", "urban")

for (i in varlist) {
  r = rast(paste0("data/land_data/1960-1979Means/", i,".asc"))
  r[is.na(r)] <- -9999
  writeRaster(r, paste0("data/land_data/land_layers/",i,".asc"), overwrite = T, NAflag = -9999)
}

r = rast("data/land_data/1960-1979Means/primf.asc")
         
maxent_call(
  "~/Desktop/maxent-2",
  features = c("h","p"),
  environmentallayers = "data/land_data/land_layers",
  samplesfile = "data/bird_data/combined_1960_1979_30km_thin1.csv",
  outputdirectory = "outputs/proj_layers/temp",
  projectionlayers = "data/land_data/1980-2015_ascs/2000",
  extrapolate = T,
  askoverwrite = T,
  randomtestpoints = 25,
  betamultiplier = 0.5,
  biasfile = "outputs/biasfile.asc",
  doclamp = T,
  randomseed = T,
  replicates = 20,
  replicatetype = "bootstrap")

#### make landcover layers for years of 1980 - 2020
library(ncdf4)
library(terra) # package for spatraster manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(maptools)
library(magrittr) #piping

#Setup environment
setwd("~/Honours")
ld_statesA <- nc_open("data/land_data/LUH2_GCB2019_states.nc4")
ld_statesB <- nc_open("data/land_data/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-MAGPIE-ssp585-2-1-f_gn_2015-2100.nc")


#Set extents
lon <- ncvar_get(ld_statesA, "lon")
lon[lon > 180] <- lon[lon > 180] - 360
lat <- ncvar_get(ld_statesA, "lat", verbose = F)
t <- ncvar_get(ld_statesA, "time")
extent_aus <- ext(110, 158, -46, -9)
world = ext(min(lon), max(lon), min(lat), max(lat))
temp <- rast("/Users/u6638201/MaxEnt/biasfile.asc")

#function to retrieve rasterlayers from array
get.raster <- function(x, y, z) {
  if (x < 2016) {
    a = x - 849
    b = y[a,,] %>% t() %>% rast(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0") %>% flip(direction='v') %>% set.ext(world) %>% crop(extent_aus)  %>% resample(temp, method = "rms") 
    b[is.na(b)] <- -9999
    writeRaster(b, paste0("data/land_data/1980-2015_ascs/", x,"/", z,".asc"), overwrite = T, NAflag = -9999)
  }
  else {
    a = x - 2014
    b = y[,,a] %>% t() %>% rast(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0") %>% flip(direction='v') %>% set.ext(world) %>% crop(extent_aus)  %>% resample(temp, method = "rms") 
    b[is.na(b)] <- -9999
    writeRaster(b, paste0("data/land_data/1980-2015_ascs/", x,"/", z,".asc"), overwrite = T, NAflag = -9999)
  }
  
}


#variable setup
varlist <- c("pastr", "primf", "primn", "range", "secdf", "secdn", "secma", "secmb", "urban")
yearlistA <- c(1980:2015)
yearlistB <- c(2016:2020) 

#for loop to run get rasters from range of years in yearlist
for (i in "secmb") {
  array = ncvar_get(ld_statesA, i)
  lapply(yearlistA, get.raster, y = array, z = i)
}

for (i in varlist) {
  array = ncvar_get(ld_statesB, i)
  lapply(yearlistB, get.raster, y = array, z = i)
}


###########Landcover means

#function to retrieve rasterlayers from array
land.mean <- function(x, y, z) {
  a = x - 849
  b = y[a,,] %>% t() %>% rast(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0") %>% flip(direction='v') %>% set.ext(world) %>% resample(temp, method = "rms") 
  b[is.na(b)] <- -9999
  NAflag(b) <- -9999
  return(b)
}

#variable setup
varlist <- c("pastr", "primf", "primn", "range", "secdf", "secdn","secma", "secmb", "urban")
mean.yearlist <- c(1960:1979)


#for loop to run get rasters from range of years in yearlist
for (i in "pastr") {
  array = ncvar_get(ld_statesA, i)
  rast.list = lapply(mean.yearlist, land.mean, y = array, z = i)
  b = rast(rast.list) %>% mean() 
  b[is.na(b)] <- -9999
  writeRaster(b, paste0("data/land_data/land_layers/", i, ".asc"), overwrite = T, NAflag = -9999)
}
