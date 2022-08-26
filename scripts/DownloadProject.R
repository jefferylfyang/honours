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



