library(ncdf4)
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(factoextra)

setwd("~/Desktop/LandCover")
ld_states <- nc_open("LUH2_GCB2019_states.nc4")

lon <- ncvar_get(ld_states, "lon")
lon[lon > 180] <- lon[lon > 180] - 360
lat <- ncvar_get(ld_states, "lat", verbose = F)
t <- ncvar_get(ld_states, "time")

#Make list of variables relevant to land cover
varlist <- c("pastr", "primf", "primn", "range", "secdf", "secdn", "secma", "secmb", "urban")

#Function to sample a selected year and variable
SampleRaster <- function(var, year) {
  year1 <- year - 849
  var.array <- ncvar_get(ld_states, var)
  fillvalue <- ncatt_get(ld_states, var, "_FillValue")
  var.array[var.array == fillvalue$value] <- NA
  var.slice <- var.array[year1, ,]
  extent_aus <- c(110, 158, -46, -9)
  ras <- raster(t(var.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  ras <- flip(ras, direction='y')
  ran <- sampleRandom(ras, 1000, na.rm = T, extent_aus, F, F, T, T, F)
  write.csv(ran, paste0("~/Desktop/LandCoverData/PCAfiles/RandomSample", var, ".csv"))
}

read.

SampleRaster("primf", 2010)

#loop to extract data for each variable
for (i in varlist){
  SampleRaster(i, 2015)
}

printPCA <- function(year) {
 
  for (i in varlist){
    SampleRaster(i, year)
  }
  setwd("~/Desktop/LandCoverData/PCAfiles/")
  temp = list.files(pattern="*.csv")
  for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))
  d1 <- RandomSamplepastr.csv[[4]]
  d2 <- RandomSampleprimf.csv[[4]]
  d3 <- RandomSampleprimn.csv[[4]]
  d4 <- RandomSamplerange.csv[[4]]
  d5 <- RandomSamplesecdf.csv[[4]]
  d6 <- RandomSamplesecdn.csv[[4]]
  d7 <- RandomSamplesecma.csv[[4]]
  d8 <- RandomSamplesecmb.csv[[4]]
  d9 <- RandomSampleurban.csv[[4]]
  setwd("~/Desktop/LandCoverData/")
  PCAData <- data.frame(d1, d2, d3, d4, d5, d6, d7, d8, d9)
  colnames(PCAData) <- c( "pastr", "primf", "primn", "range", "secdf", "secdn","secma", "secmb", "urban")
  landcover.pca <- prcomp(PCAData, scale = TRUE)
  png(paste0("~/Desktop/LandCoverData/PCAPlots/PCA",year,".png"))
  print({fviz_pca_var(landcover.pca,
                      col.var = "contrib", # Color by contributions to the PC
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      repel = TRUE     # Avoid text overlapping
  )})
  dev.off()
}

PCA <- function(year) {
  setwd("~/Desktop/LandCoverData/PCAfiles/")
  temp = list.files(pattern="*.csv")
  for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))
  d1 <- RandomSamplepastr.csv[[4]]
  d2 <- RandomSampleprimf.csv[[4]]
  d3 <- RandomSampleprimn.csv[[4]]
  d4 <- RandomSamplerange.csv[[4]]
  d5 <- RandomSamplesecdf.csv[[4]]
  d6 <- RandomSamplesecdn.csv[[4]]
  d7 <- RandomSamplesecma.csv[[4]]
  d8 <- RandomSamplesecmb.csv[[4]]
  d9 <- RandomSampleurban.csv[[4]]
  setwd("~/Desktop/LandCoverData/")
  PCAData <- data.frame(d1, d2, d3, d4, d5, d6, d7, d8, d9)
  colnames(PCAData) <- c( "pastr", "primf", "primn", "range", "secdf", "secdn","secma", "secmb", "urban")
  landcover.pca <- prcomp(PCAData, scale = TRUE)


}

PCA(1980)

yearlist <- c(1980:2018)

for (i in yearlist){
  PCA(i)
  
}




