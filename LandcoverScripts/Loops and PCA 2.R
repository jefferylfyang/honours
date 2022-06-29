var <- ncvar_get(ld_states, "urban")
fillvalue <- ncatt_get(ld_states, "urban", "_FillValue")
var[var == fillvalue$value] <- NA
yearlist <- c(1980:2018)
extent_aus <- c(110, 158, -46, -9)

SampleRaster <- function(year) {
  year1 <- year - 849
  var.slice <- var[year1, ,]
  ras <- raster(t(var.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  ras <- flip(ras, direction='y')
  ran <- sampleRandom(ras, 1000, na.rm = T, extent_aus, F, F, T, T, F)
}

YearSampleList  <- list()
RanSampleList  <- list()

for(i in yearlist) {
  variable <- SampleRaster(i)  ## define a variable to be added
  YearSampleList[[paste0("urban",i,"array")]] <- variable
}

for(i in YearSampleList){
  variable <- `$`(i, 'layer')
  variable <- list(variable)
  RanSampleList <- append(RanSampleList, variable)
}

RanSampleList <- data.frame(RanSampleList)
colnames(RanSampleList) <- yearlist
write.csv(RanSampleList, "~/Desktop/LandCoverData/PCAfiles/urban.csv")



## ACTUAL PCA PART
setwd("~/Desktop/LandCoverData/PCAfiles/")
temp = list.files(pattern="*.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

my_data <- lapply(temp, read.csv)
PCADataList <- list()

names(PCADataList) <- names(varlist)
pca.variance <- list()
for (i in temp) {
  pca <- prcomp(i, scale = T)
  eigs <- pca$sdev^2
  variance.explained <- eigs[1]/sum(eigs)
  pca.list <- append(pca.list, variance.explained)
}
 
pca <- prcomp(PCA1980, scale = T)
eigs <- pca$sdev^2
variance.explained <- eigs[1]/sum(eigs)
pca.variance <- append(pca.variance, variance.explained)

pcavartable <- cbind(yearlist, pca.variance)

library(dplyr)
library(xts)

#extrct column and bind
yearlist2 <- list("X1980", 'X1981', 'X1982', 'X1983', 'X1984', 'X1985', 'X1986', 'X1987', 'X1988', 'X1989', 'X1990', 'X1991', 'X1992',
                  'X1993', 'X1994', 'X1995', 'X1996', 'X1997', 'X1998', 'X1999', 'X2000', 'X2001', 'X2002', 'X2003', 'X2004', 'X2005', 'X2006', 'X2007', 
                  'X2008', 'X2009', 'X2010', 'X2011', 'X2012', 'X2013', 'X2014', 'X2015', 'X2016', 'X2017', 'X2018')

for (i in yearlist2) {
  dat <- tibble::lst(pastr.csv, primf.csv, primn.csv, range.csv, secdf.csv, secdn.csv, secma.csv, secmb.csv, urban.csv) %>%
    purrr::map(as.data.frame) %>%
    purrr::imap(function(dat, ds) select(dat, {{ ds }} := i)) %>%
    bind_cols()
    pca <- prcomp(dat, scale = T)
    eigs <- pca$sdev^2
    variance.explained <- eigs[1]/sum(eigs)
    pca.variance <- append(pca.variance, variance.explained)
}
dat <- tibble::lst(pastr.csv, primf.csv, primn.csv, range.csv, secdf.csv, secdn.csv, secmb.csv, urban.csv) %>%
  purrr::map(as.data.frame) %>%
  purrr::imap(function(dat, ds) select(dat, {{ ds }} := i)) %>%
  bind_cols()


ggbiplot(pca, obs.scale = 1, var.scale = 1,
         ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

plot(primf.csv)
#Make delta of all variables 5 years in intervals from 1980 - 2020
#atlas of living australia 
#species distribution modeling 
#use climate (generic environmental data) accounting for bias 