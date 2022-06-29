#script to print all data visualisation fun with loops
variables <- attributes(ld_states$var)
years <- 1980:2018

splice("pastr")

for (i in variables) {
  var.array <- ncvar_get(ld_states, i)
  fillvalue <- ncatt_get(ld_states, i, "_FillValue")
  var.array[var.array == fillvalue$value] <- NA
  for (i in years) {
    var.slice <- var.array[i, ,]
    r <- raster(t(var.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
    r <- flip(r, direction='y')
    extent_aus <- c(110, 158, -46, -9)
    zoomed <- crop(r, extent_aus)
    png(paste0("~/Desktop/LandCoverData/Plots/",var,i),
        width = 480, height = 480,
        units = "px", pointsize = 12, bg = "white", res = NA,
        restoreConsole = TRUE)
    plot(zoomed)
    data(wrld_simpl)
    plot(wrld_simpl, add = TRUE)
    dev.off()
    
  }
}




for (i in variables) {
  i <- ncvar_get(ld_states, i)
  fillvalue <- ncatt_get(ld_states, i, "_FillValue")
  i[var.array[var.array == fillvalue$value] <- NA]
}


###


splice <- function(var) {
  var.array <- ncvar_get(ld_states, var)
  fillvalue <- ncatt_get(ld_states, var, "_FillValue")
  var.array[var.array == fillvalue$value] <- NA
  return(var.array)
  
}

pastr <- splice("pastr")

for (i in years) {
  var.slice <- var.array[i, ,]
  r <- raster(t(var.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  r <- flip(r, direction='y')
  extent_aus <- c(110, 158, -46, -9)
  zoomed <- crop(r, extent_aus)
  png(paste0("~/Desktop/LandCoverData/Plots/",var,i),
      width = 480, height = 480,
      units = "px", pointsize = 12, bg = "white", res = NA,
      restoreConsole = TRUE)
  plot(zoomed)
  data(wrld_simpl)
  plot(wrld_simpl, add = TRUE)
  dev.off()

#apply function to sum up NA in coordinates

#to sum zeros at a given lat and lon across all years combined
sum.zero <- function(x) {
  sum(is.zero(x))
  
}

#to test for zero 
is.zero <- function(x) {
  x == 0 


}


viszero <- function(var) {
  var.array <- ncvar_get(ld_states, var)
  fillvalue <- ncatt_get(ld_states, var, "_FillValue")
  var.array[var.array == fillvalue$value] <- NA
  zerosum <- apply(var.array, c(2,3), sum.zero)
  zerosum <- raster(t(zerosum), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  zerosum <-flip(zerosum, direction='y')
  extent_aus <- c(110, 158, -46, -9)
  zoomed <- crop(zerosum, extent_aus)
  png(paste0("~/Desktop/LandCoverData/ZeroSumPlots/",var),
      width = 480, height = 480,
      units = "px", pointsize = 12, bg = "white", res = NA)
  plot(zoomed)
  data(wrld_simpl)
  plot(wrld_simpl, add = TRUE)
  dev.off()
}

#Make list of variables relevant to land cover
varlist <- list(attributes(ld_states$var))

#Run above function
apply(varlist, 1, sum.zero)

#Plot above
extent_aus <- c(110, 158, -46, -9)
write.csv(RandomSample1000Secdf, file = "RandomSample1000Secdf.csv")
RandomSample1000Secdf <- sampleRandom(r, 1000, na.rm = T, extent_aus, F, F, T, T, F)
sumzero <- apply(secdf.array, c(2,3), sum.zero)

#make arrays
arraylist <- list()
varlist <- c("pastr", "primf", "primn", "range", "secdf", "secdn", "secma", "secmb", "urban")
extent_aus <- c(110, 158, -46, -9)

AddArray <- function(var){
  var.array <- ncvar_get(ld_states, var)
  fillvalue <- ncatt_get(ld_states, var, "_FillValue")
  var.array[var.array == fillvalue$value] <- NA
  return(var.array)
}

#add them to a list
for(i in varlist) {
    variable <- AddArray(i)  ## define a variable to be added
    arraylist[[paste0(i,".array")]] <- variable
  }



#Function to sample a selected year and variable
SampleRaster <- function(vararray, year) {
year1 <- year - 849
var.slice <- vararray[year1, ,]
ras <- raster(t(var.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
ras <- flip(ras, direction='y')
ran <- sampleRandom(ras, 1000, na.rm = T, extent_aus, F, F, T, T, F)
write.csv(ran, paste0("~/Desktop/LandCoverData/PCAfiles/", vararray, year, ".csv"))
}



#loop to extract data for each year
yearlist <- c(1980:2018)
for (i in arraylist) {
  lapply(yearlist, SampleRaster, var.array = i)
 }


#read list of .csv files
setwd("~/Desktop/LandCoverData/PCAfiles/")
temp = list.files(pattern="*.csv")
 for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

#extract column of layer values
d1 <- RandomSamplepastr.csv[[4]]
d2 <- RandomSampleprimf.csv[[4]]
d3 <- RandomSampleprimn.csv[[4]]
d4 <- RandomSamplerange.csv[[4]]
d5 <- RandomSamplesecdf.csv[[4]]
d6 <- RandomSamplesecdn.csv[[4]]
d7 <- RandomSamplesecma.csv[[4]]
d8 <- RandomSamplesecmb.csv[[4]]
d9 <- RandomSampleurban.csv[[4]]

PCAData2000 <- data.frame(d1, d2, d3, d4, d5, d6, d7, d8, d9)
colnames(PCAData2000) <- c( "pastr", "primf", "primn", "range", "secdf", "secdn","secma", "secmb", "urban")
write.csv(PCAData2000, "PCAData2010.csv")

## ACTUAL PCA PART
install.packages("factoextra")
library(factoextra)
landcover.pca <- prcomp(PCAData2000, scale = TRUE)
fviz_eig(landcover.pca)
fviz_pca_ind(landcover.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
# Eigenvalues
eig.val <- get_eigenvalue(landcover.pca)
eig.val

# Results for Variables
landcover.pca <- get_pca_var(landcover.pca)
landcover.pca$coord          # Coordinates
landcover.pca$contrib        # Contributions to the PCs
landcover.pca$cos2           # Quality of representation 
# Results for individuals
landcover.ind <- get_pca_ind(landcover.pca)
landcover.ind$coord          # Coordinates
landcover.ind$contrib        # Contributions to the PCs
landcover.ind$cos2           # Quality of representation 

fviz_pca_var(landcover.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

ggbiplot(landcover.pca, obs.scale = 1, var.scale = 1,
         ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
