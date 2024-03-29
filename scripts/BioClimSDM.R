# Species distribution modeling for BioClim

rm(list = ls())
#Note to self running the accum_curve_occ.R requires the objects from this file 
# Load packages
library("sp")
library("raster")
library("maptools")
library("rgdal")
library("dismo")
library("magrittr")

# Download bioclim data and store in bioclim.data variable
bioclim.data <- getData(name = "worldclim",
                        var = "bio",
                        res = 2.5,
                        path = "data/")

#Add landcover to bioclim

bil.files <- list.files(path = "data/wc2-5", 
                          pattern = "*.bil$", 
                          full.names = TRUE)


extent_aus <- c(110, 158, -46, -9)
env.rasters <- lapply(bil.files, raster)
env.rasters <- lapply(env.rasters,crop, y = extent_aus)
clim <- raster::stack(env.rasters)

  
  
# Read in bird observations
obs.data <- read.csv(file = "data/MalurusPresThin30.csv")

# Drop any rows with NAs
obs.data <- obs.data[!is.na(obs.data$latitude), ]

# Only pull out those columns of interest and in the order we want them
obs.data <- obs.data[, c("longitude", "latitude")]

# Determine geographic extent of our data
max.lat = ceiling(max(obs.data$latitude))
min.lat = floor(min(obs.data$latitude))
max.lon = ceiling(max(obs.data$longitude))
min.lon = floor(min(obs.data$longitude))
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))

# Crop the bioclim data to geographic extent of birds
clim <- crop(x = clim, y = geographic.extent)

# Create pseudo-absence, or background, points
# Use the bioclim data files for sampling resolution
bil.files <- list.files(path = "data/wc2-5", 
                        pattern = "*.bil$", 
                        full.names = TRUE)

# We only need one file, so use the first one in the list of .bil files
mask <- raster(bil.files[1])

# Randomly sample points (same number as our observed points)
background <- randomPoints(mask = mask,     # Provides resolution of sampling points
                           n = nrow(obs.data),      # Number of random points
                           ext = geographic.extent, # Spatially restricts sampling
                           extf = 1.25)             # Expands sampling a little bit

# Arbitrarily assign group 2 as the testing data group
testing.group <- 1

# Create vector of group memberships
group.presence <- kfold(x = obs.data, k = 5) # kfold is in dismo package

# Separate observations into training and testing groups
presence.train <- obs.data[group.presence != testing.group, ]
presence.test <- obs.data[group.presence == testing.group, ]

# Repeat the process for pseudo-absence points
group.background <- kfold(x = background, k = 5)
background.train <- background[group.background != testing.group, ]
background.test <- background[group.background == testing.group, ]

# Build a model using training data
bc.model <- bioclim(x = clim, p = presence.train)

# Predict presence from model
predict.presence <- dismo::predict(object = bc.model, 
                                   x = clim, 
                                   ext = geographic.extent)

# Use testing data for model evaluation
bc.eval <- evaluate(p = presence.test,   # The presence testing data
                    a = background.test, # The absence testing data
                    model = bc.model,    # The model we are evaluating
                    x = clim)    # Climatic variables for use by model

# Determine minimum threshold for "presence"
bc.threshold <- threshold(x = bc.eval, stat = "spec_sens")

# Load map data for plotting
data(wrld_simpl)

# Plot base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")

# Only plot areas where probability of occurrence is greater than the threshold
plot(predict.presence > bc.threshold, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "olivedrab"))

# And add those observations
points(x = obs.data$longitude, 
       y = obs.data$latitude, 
       col = "black",
       pch = "+", 
       cex = 0.6)

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()


#SDM output values for acum_curve_occ
SDMOutPuts <- rasterToPoints(predict.presence) %>% 
  data.frame
  
EnvPoints <- rasterToPoints(clim) %>% data.frame 
EnvPoints <- na.omit(EnvPoints)
SDMOutPuts <- cbind(SDMOutPuts, EnvPoints[3:30])


#Extract points to form occ.pnts

coordinates(obs.data)= ~ longitude + latitude

rasValue=extract(clim, obs.data)

ModelValue <- extract(predict.presence, obs.data)

rasValue.df=data.frame(rasValue)
obs.data.df=data.frame(obs.data)

occ.pntsbirds=cbind(obs.data.df, ModelValue, rasValue.df) %>% na.omit


