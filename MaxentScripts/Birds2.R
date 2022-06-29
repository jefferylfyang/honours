```{r}
library(dismo) # interface with MaxEnt
library(raster) # spatial data manipulation
library(MASS) # for 2D kernel density function
library(magrittr) # for piping functionality, i.e., %>%
library(maptools) # reading shapefiles
library(rgeos)
library(rJava)
library(data.table)
library(dplyr)
library(sp)

```


download and install maxent java
```{r}
utils::download.file(url = "https://raw.githubusercontent.com/mrmaxent/Maxent/master/ArchivedReleases/3.3.3k/maxent.jar", destfile = paste0(system.file("java", package = "dismo"), "/maxent.jar"), mode = "wb") 
```

Read in bird data 
```{r}
birds1 <- fread("lambertiANDassimilis_post1970.csv", select = c("decimalLongitude", "decimalLatitude")) %>% 
  mutate(species = "Malurus", .before = "decimalLongitude") %>% 
  na.omit()

colnames(birds1) <- c("species", "longitude", "latitude")
write.csv(birds1, "birdsdata.csv", row.names = F)

birds2 <- read.csv("assimilis_trimmed.csv") %>% mutate(species = "Malurus assimilis", .before = "lon")
birds <- rbind(birds1, birds2)
write.csv(birds, "birdsdata.csv")
birds <- read.csv("birdsdata.csv")
na.omit(birds)
```
Read in data from preserved specimens

```{r}
birdsPres <- fread("lambertiANDassimilis_19702000_preserved.csv", select = c("decimalLongitude", "decimalLatitude")) %>%
  mutate(species = "Malurus", .before = "decimalLongitude") %>%
  na.omit()
colnames(birdsPres) <- c("species", "longitude", "latitude")
write.csv(birdsPres, "birdsPres.csv", row.names = F)
```


Read in preserved data from 1980-2020

```{r}
lambertiPresDecades <- fread("lambertiPres1980-2020.csv", select = c("decimalLongitude", "decimalLatitude")) %>%
  mutate(species = "Malurus lamberti", .before = "decimalLongitude")
colnames(lambertiPresDecades) <- c("species", "longitude", "latitude")
lambertiPresDecades <- na.omit(lambertiPresDecades)

```

Thinning
``` {r}
library(spThin)

thin(birdsPres, "latitude", "longitude", "species", 1, 2)
subsetbirds1 <- birds1[-c(500:nrow(birds1)), ]
thin(subsetbirds1, "latitude", "longitude", "species", 1, 2, out.dir = "/Users/u6638201/Desktop/Maxent/ThinnedBirdData/")

thin(birdsPres, "latitude", "longitude", "species", 30, 1, out.dir = "/Users/u6638201/Desktop/Maxent/ThinnedBirdData/")

thin(lambertiPresDecades, "latitude", "longitude", "species", 30, 1, out.dir = "/Users/u6638201/Desktop/Maxent/ThinnedBirdData/")
```


Thin and trim bird data

```{r}

filterByProximity <- function(xy, dist, mapUnits = F) {
    #xy can be either a SpatialPoints or SPDF object, or a matrix
    #dist is in km if mapUnits=F, in mapUnits otherwise
    if (!mapUnits) {
        d <- spDists(xy,longlat=T)
    }
    if (mapUnits) {
        d <- spDists(xy,longlat=F)
    }
    diag(d) <- NA
    close <- (d <= dist)
    diag(close) <- NA
    closePts <- which(close,arr.ind=T)
    discard <- matrix(nrow=2,ncol=2)
    if (nrow(closePts) > 0) {
            while (nrow(closePts) > 0) {
                if ((!paste(closePts[1,1],closePts[1,2],sep='_') %in% paste(discard[,1],discard[,2],sep='_')) & (!paste(closePts[1,2],closePts[1,1],sep='_') %in% paste(discard[,1],discard[,2],sep='_'))) {
                discard <- rbind(discard, closePts[1,])
                closePts <- closePts[-union(which(closePts[,1] == closePts[1,1]), which(closePts[,2] == closePts[1,1])),]
                }
            }
        discard <- discard[complete.cases(discard),]
        return(xy[-discard[,1],])
    }
    if (nrow(closePts) == 0) {
        return(xy)
    }
}

birds$species <- NULL
birds <- data.matrix(birds) %>% na.omit
filterByProximity(birds, 2)

```

```{r}
#Octavio's function 
gridbiostats <- function(spoints, extent, reskm=NA, resdg, minSampSize=0){
  if(length(spoints$species)==0){ stop("Error. There must be a field called 'species' in the object spoints")}
  resdg <- ifelse(is.na(reskm), resdg, reskm*0.008333333)
  r <- raster(ext=extent(extent), resolution=resdg)
  newcolumn <- length(spoints@data)
  r[]<- 1:ncell(r)
  #spoints@data$
  cells <- extract(r, spoints) # vector with the cell identities
  r[] <- NA # 0
  Nspecimens <- SpRichness <- Coverage <- Chao <- r
  # tab <- table(cellFromXY(r, spatial))
  # r[as.numeric(names(tab))] <- tab # number of individuals per cell
 
  for(i in 1:length(levels(as.factor(cells)))){
    cell <- as.numeric(levels(as.factor(cells))[i])
    #sp <- subset(spoints, cells==cell)
    sp <- spoints[which(cells==cell),]
   
    # Number of specimens per cell
    Nspecimens[cell] <- nrow(sp)
   
    # Number of species per cell
    sp@data$species <- droplevels(sp@data$species) # This means that the initial file must have a column called "species"
    SpRichness[cell] <- length(levels(sp@data$species))
   
    # Sample Coverage per cell
    # We could produce a table with the species abundance distribution per all the sites (cells),
    # but for the moment let's stick to do the calculations by a single site.
    # We need a vector of number of individuals per species, ¿ranked by abundance?
    tab <- as.data.frame(table(sp@data$species))
    tab$Freq <- tab$Freq-minSampSize
    tab$Freq[which(tab$Freq<1)] <- 1
    #tab <- tab[order(tab ,decreasing=T)]
    Coverage[cell] <- DataInfo(tab$Freq, datatype="abundance")$SC
   
    # Estimated Species Richness per cell
    Chao[cell] <- ChaoSpecies(tab$Freq, datatype='abundance', conf=0.95)$Estimator
   
   
  }
  STACK <- stack(Nspecimens, SpRichness, Coverage, Chao)
  names(STACK) <- c("Nspecimens", "SpRichness", "SampleCoverage", "ChaoEstimatedRichness")
  return(STACK)
}
```

Read in climate data
```{r}
clim_list <- list.files("/Users/u6638201/Desktop/Maxent/clim", pattern = ".tif$", 
                        full.names = T)  
clim_raster <- lapply(clim_list, raster)
clim_raster <- lapply(clim_raster, crop, y = extent_aus)
clim <- raster::stack(clim_list)
for (i in 1:19) {
  writeRaster(clim_raster[[i]], paste0("clim",i,".tif"))
}
```

```{r}

model1 <- maxent(clim, locations, args = "biasfile=dens.ras")

```
```{r}

```



