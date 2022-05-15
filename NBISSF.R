############################################################
############################################################
#### Helena Wehner                                      ####
#### MB2 Project - EAGLE Master University of Wuerzburg ####
#### First Spring Migration of Northern Bald Ibises     ####
#### Exploratory Movement statistics,                   ####
#### Probability Map &                                  ####
#### Step Selection Analysis                            ####
############################################################
############################################################


## Set working directory ##
###########################

setwd('F:/MB2/') # set your own working directory

## Downdload and install needed R-Packages ##
#############################################

if(!require('rgdal')) {
  install.packages('rgdal');
  library(rgdal)
}

if(!require('sp')) {
  install.packages('sp');
  library(sp)
}

if(!require('raster')) {
  install.packages('raster');
  library(raster)
}

if(!require('dplyr')) {
  install.packages('dplyr');
  library(dplyr)
}

if(!require('XML')) {
  install.packages('XML');
  library(XML)
}

if(!require('lubridate')) {
  install.packages('lubridate');
  library(lubridate)
}

if(!require('ggmap')) {
  install.packages('ggmap');
  library(ggmap)
}

if(!require('ggplot2')) {
  install.packages('ggplot2');
  library(ggplot2)
}

if(!require('amt')) {
  install.packages('amt');
  library(amt)
}

if(!require('mapview')) {
  install.packages('mapview');
  library(mapview)
}

if(!require('tmaptools')) {
  install.packages('tmaptools');
  library(tmaptools)
}

if(!require('sf')) {
  install.packages('sf');
  library(sf)
}

if(!require('rgeos')) {
  install.packages('rgeos');
  library(rgeos)
}

if(!require('ggmap')) {
  install.packages('ggmap');
  library(gfmap)
}

if(!require('move')) {
  install.packages('move');
  library(move)
}

if(!require('tidyr')) {
  install.packages('tidyr');
  library(tidyr)
}

if(!require('dismo')) {
  install.packages('dismo');
  library(dismo)
}

if(!require('mgcv')) {
  install.packages('mgcv');
  library(mgcv)
}

if(!require('ellipse')) {
  install.packages('ellipse');
  library(ellipse)
}

if(!require('RStoolbox')) {
  install.packages('RStoolbox');
  library(RStoolbox)
}

if(!require('ggsn')) {
  install.packages('ggsn');
  library(ggsn)
}

source('varImpBiomod.R')

#------------------------------------------------------------------------------------------------------------
## Load all needed Data ##
##########################

## Load Study Area
area <- readOGR('studyArea.2.shp')
str(area)
mapview(area)

## Load Site of breeding area
site <- readOGR('site.shp')

## Load Movement Data as SpatialPointsDataFrame
rupertSP <- readOGR('RupertSP.shp')
cupiSP <- readOGR('CupiSP.shp')
cupiSP2 <- readOGR('cupiAlps.shp')

birds <- rbind(rupertSP, cupiSP)
birds <- intersect(birds, area)
birdsDF <- as.data.frame(birds)

mapview(birds)

#-----------------------------------------------------------------------------------------------------------
### Data Exploration ###
########################
birds <- birds[order(birds$UTC_dateti),]

birdsDf <- as.data.frame(birds)
birdsDf <- birdsDf[!duplicated(birdsDf[,c("UTC_dateti")]),]

# create Move Object
move <- move::move(x = as.numeric(birdsDf$Longitude), y = as.numeric(birdsDf$Latitude),
                   time = as.POSIXct(birdsDf$UTC_dateti,format="%Y-%m-%d %H:%M:%S",tz="UTC"),
                   data = birdsDf, proj = crs(birds))
move <- spTransform(move, center = T)

# Net Squared Displacement (NSD)
# By calculating NSD we are having a first look on the Displacement of the Northern Bald Ibis Movements Data.
# Does it look like they are staying around at one place? Does it seem like they leave their area for a longer?

nsd <- sp::spDistsN1(pts = move, pt = move[1,])

# Net squared Displacement Plot
birdsDf$UTC_dateti <- as.POSIXct(birdsDf$UTC_dateti,format="%Y-%m-%d %H:%M:%S",tz="UTC")
ggplot(birdsDf, aes(x = UTC_dateti, y = nsd))+geom_point()+scale_x_datetime()+
  xlab('Timespan (Apr. 2021 - Jun. 2021)')+
  ylab('Net squared Displacement')+
  ggtitle('Distance to Wintering Site over Time')

# Point density during Migrtion
ggplot(birdsDF, aes(x = coords.x1, y = coords.x2)) + 
  geom_point(size = 1, alpha = 0.1) + 
  coord_equal() + 
  xlab('Longitude') + 
  ylab('Latitude')


# The Northern Bald Ibis is a migratory bird in the NSD plot, his migration and stay in the breeding area during summer
# is mapped.
# In the following analysis we are having a look at his migration behaviour. Which environmental variables do mostly
# influence his migration route?
# Two young Northern Bald Ibises (Rupert, Cupi) which obviusly flew during different times and took different routes,
# especially while crossing the european Alps, have been chosen as example. It has been their first spring migration.

# Two analyses will be carried out: Probability map based on general additive model and step selection functions.

#----------------------------------------------------------------------------------------------------------------------
### Add Environmental Parameters ###
####################################
# potential variables causing higher point density during Migration

# Elevation
# Slope
# Aspect
# Distance to Breeding Site

### Elevation/SRTM ###
srtm1 <- getData('SRTM', lon = 13, lat = 48)
srtm2 <- getData('SRTM', lon = 9, lat = 49)
srtm3 <- getData('SRTM', lon = 9, lat = 44)
srtm4 <- getData('SRTM', lon = 13, lat = 44)

# merge all SRTM files
srtm <- raster::merge(srtm1, srtm2, srtm3, srtm4)
mapview(srtm)

# calculate Slope and Aspect
slope <- raster::terrain(srtm, opt = 'slope', unit = 'degrees')
aspect <- raster::terrain(srtm, opt = 'aspect', unit = 'degrees')

mapview(slope)
mapview(aspect)


### Distance to Breeding Site ###

# shortest Distance (siteD) between Birds and breeding site in Ueberlingen/final destination
# raster solution
siteD <- distanceFromPoints(object = srtm, xy = site)
str(siteD)
mapview(siteD)

# stack all environmental Raster Data
env <- stack(srtm, slope, aspect, siteD)
env2 <- aggregate(env, 10)
names(env) <- c('elevation','slope','aspect', 'siteD')

#-------------------------------------------------------------------------------------------------------------------
### Probability Map ###
#######################


# selecting 10000 random background samples
set.seed(2)
background <- randomPoints(env, 2000, birds)

# select only one presence record in each cell of environmental layer
presence <- gridSample(birds, env, n=1)

# now we combine the presence and background points, adding a column "species" that contains
# the information about presence (1) and background (2)
fulldata <- SpatialPointsDataFrame(rbind(presence, background),
                                   data = data.frame("species" = rep(c(1,0),c(nrow(presence),nrow(background)))),
                                   match.ID = F,
                                   proj4string = CRS(projection(env)))

# add information about environmental conditions at point locations
fulldata@data <- cbind(fulldata@data, raster::extract(env, fulldata))

# test for collinearity between the environmental variables and decide which to take for the model 
# Visual inspection of collinearity
cm <- cor(getValues(env2), use = "complete.obs") # pearson coefficient (default)
cm
colors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white",
            "#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C")   
plotcorr(cm, col = ifelse(abs(cm)>0.8, 'red','blue'), type = 'lower')


# split data set into training and test/validation data
set.seed(2)
fold <- kfold(fulldata, k=5)
traindata <- fulldata[fold != 1,]
testdata <- fulldata[fold == 1, ]

varname <- names(env)

### the Actual Model
### General Additive Model
gammodel <- gam(species ~ s(elevation) +
                  s(slope) + s(aspect)+ s(siteD),
                family = "binomial", data = traindata)

# variable importance
gamimp <- varImpBiomod(model = gammodel, varnames = c("elevation",
                                                      "slope",
                                                      "aspect","siteD"),
                       data = traindata, n = 10)

# Plot Variable Importance of Probability Map

gamimp3 <- 100*gamimp

gamimp2 <- as.data.frame(gamimp3)
rownames(gamimp2) <- c("Elevation","Slope","Aspect","SiteD")
gamimp2 <- tibble::rownames_to_column(gamimp2, "Name")
colnames(gamimp2) <- c("Variable","Percent")

ggplot(gamimp2, aes(x=Variable,y=Percent, fill = Variable))+geom_col()+
  labs(title="Variable Importance", subtitle = "2021 Apr.- Jun.")+
  scale_fill_manual(values = c("green","lightblue","orange","purple"))+
  geom_text(aes(label=round(Percent, digits = 1), vjust=2.5))+
  ylim(0,100)

### Prediction Map
gammap <- raster::predict(env2, gammodel, type="response")
writeRaster(gammap, filename = 'GammapMB2.tiff')
mapview(gammap)


#--------------------------------------------------------------------------------------------------------------------
### Step Selection Functions - Rupert ###
#########################################

### Preparation for Step Selection Function ###

# Clean points from NA times and select relevant cols
rupertSP <- intersect(rupertSP, area)
mapview(rupertSP)
r <- as.data.frame(rupertSP)
p <- data.frame(r$UTC_dateti, r$Latitude, r$Longitude)
names(p) <- c('time','y','x')

p[2:3] <- lapply(p[2:3],as.numeric) 

p$time <- as.POSIXct(p$time,format="%Y-%m-%d %H:%M:%S",tz="UTC")

pm <- move(x=p$x, y=p$y, 
           time=as.POSIXct(p$time, format="%Y-%m-%d %H:%M:%S", tz="UTC"), 
           data=p, proj=CRS("+proj=longlat +ellps=WGS84"), sensor="GPS")


tl <- timeLag(x = pm, units="mins")
summary(tl)

# Filter date range
p <- p %>%
  dplyr::filter(time > as.POSIXct("2021-04-01 00:00:00", tz='UTC'))
p <- p %>%
  dplyr::filter(time < as.POSIXct("2021-06-30 23:59:00", tz='UTC'))


# Convert Movement Data
p_track <- amt::make_track(tbl = p, .x = x, .y = y, .t = time,
                           crs = sp::CRS('+proj=longlat +datum=WGS84 +no_defs '))

# get a summary of the tracks
sum <- summarize_sampling_rate(p_track)
sum

# Format and generate Random Steps
ssf_dat <- p_track %>%
  track_resample(rate = minutes(15), tolerance = seconds(2))%>%
  # Filter a minimum number of steps per track
  amt::filter_min_n_burst(min_n = 3)  %>%
  # Create steps from the bursts/fixes
  steps_by_burst()%>%
  # calculate random steps to compare with the actual steps -> defaults 10 control step
  random_steps(n_control = 15)

# Extract Covariates
ssf_preds <- extract_covariates(ssf_dat, env, where = 'both')%>%
  mutate(elevation_start = scale(elevation_start), 
         elevation_end = scale(elevation_end),
         slope_start = scale(slope_start),
         aspect_end = scale(aspect_end),
         siteD_start = scale(siteD_start), 
         siteD_end = scale(siteD_end),
         cos_ta_ = cos(ta_), 
         log_sl_ = log(sl_)
  ) %>% 
  filter(!is.na(ta_))

### Fit Step Selection Function ###
model <- ssf_preds %>%
  fit_issf(case_ ~ elevation_end + slope_end + aspect_end + siteD_end + sl_ + cos_ta_ + log_sl_ +
             elevation_start:(sl_ + log_sl_ + cos_ta_) + strata(step_id_), model = T)


### Results of Step Selection Function for one Bird - Rupert ### 
summary(model)

SSF_Coeffs <- summary(model)$coefficients
print(SSF_Coeffs)

### Results Rupert ###
coef    exp(coef)    se(coef)          z     Pr(>|z|)
elevation_end           -1.665335e+00 1.891272e-01  0.72678646 -2.2913680 2.194214e-02
slope_end               -2.252336e-01 7.983297e-01  0.05053276 -4.4571806 8.304462e-06
aspect_end               4.787915e-02 1.049044e+00  0.05011509  0.9553840 3.393836e-01
#siteD_end               -1.163486e+02 2.954333e-51 28.59179888 -4.0692986 4.715488e-05
sl_                     -1.099452e+01 1.679348e-05  8.98575551 -1.2235499 2.211221e-01
cos_ta_                 -5.929742e-01 5.526811e-01  0.06130541 -9.6724604 3.947711e-22
log_sl_                  2.768857e-02 1.028075e+00  0.01814300  1.5261294 1.269776e-01
sl_:elevation_start      1.877815e+01 1.429700e+08  7.42820169  2.5279531 1.147297e-02
log_sl_:elevation_start -8.020490e-03 9.920116e-01  0.01724087 -0.4652021 6.417867e-01
cos_ta_:elevation_start -8.720833e-03 9.913171e-01  0.06161632 -0.1415345 8.874477e-01


#-------------------------------------------------------------------------------------------------------------------
### Step Selection Function - Cupi ###
######################################

# Clean points from NA times and select relevant cols
cupiSP <- intersect(cupiSP, area)
r2 <- as.data.frame(cupiSP)
p2 <- data.frame(r2$UTC_dateti, r2$Latitude, r2$Longitude)
names(p2) <- c('time','y','x')

p2[2:3] <- lapply(p2[2:3],as.numeric) 

p2$time <- as.POSIXct(p2$time,format="%Y-%m-%d %H:%M:%S",tz="UTC")

pm2 <- move(x=p2$x, y=p2$y, 
            time=as.POSIXct(p2$time, format="%Y-%m-%d %H:%M:%S", tz="UTC"), 
            data=p2, proj=CRS("+proj=longlat +ellps=WGS84"), sensor="GPS")


tl2 <- timeLag(x = pm2, units="mins")
summary(tl2)

# Filter date range
p2 <- p2 %>%
  dplyr::filter(time > as.POSIXct("2021-04-01 00:00:00", tz='UTC'))
p2 <- p2 %>%
  dplyr::filter(time < as.POSIXct("2021-06-30 23:59:00", tz='UTC'))

# Convert Movement Data
p_track2 <- amt::make_track(tbl = p2, .x = x, .y = y, .t = time,
                            crs = sp::CRS('+proj=longlat +datum=WGS84 +no_defs '))

# get a summary of the tracks
sum2 <- summarize_sampling_rate(p_track2)
sum2

# Format and generate Random Steps
ssf_dat2 <- p_track2 %>%
  track_resample(rate = minutes(15), tolerance = seconds(2))%>%
  # Filter a minimum number of steps per track
  amt::filter_min_n_burst(min_n = 3)  %>%
  # Create steps from the bursts/fixes
  steps_by_burst()%>%
  # calculate random steps to compare with the actual steps -> 15 control step
  random_steps(n_control = 15)%>%
  extract_covariates(env, where = "both") %>%
  mutate(elevation_start = scale(elevation_start), 
         elevation_end = scale(elevation_end),
         slope_start = scale(slope_start),
         slope_end = scale(slope_end),
         aspect_start = scale(aspect_start), 
         aspect_end = scale(aspect_end), 
         siteD_start = scale(siteD_start),
         siteD_end = scale(siteD_end),
         cos_ta_ = cos(ta_), 
         log_sl_ = log(sl_)
  ) %>% 
  filter(!is.na(ta_))

### Fit Step Selection Function ###

model2 <- ssf_dat2 %>%
  fit_issf(case_ ~ elevation_end + slope_end + aspect_end + siteD_end + sl_ + log_sl_ + cos_ta_ +
             elevation_start:(sl_ + log_sl_ + cos_ta_)+
             strata(step_id_), model = T)


### Results of Step Selection Function for one Bird - Rupert ### 
summary(model2)

SSF_Coeffs2 <- summary(model2)$coefficients
print(SSF_Coeffs2)

### Results - Cupi ###
coef    exp(coef)   se(coef)          z     Pr(>|z|)
elevation_end            -1.010144016 3.641665e-01 0.33204082 -3.0422284 2.348336e-03
slope_end                -0.436937862 6.460116e-01 0.07294645 -5.9898440 2.100424e-09
aspect_end               -0.049842532 9.513792e-01 0.03507499 -1.4210278 1.553087e-01
#siteD_end               -40.901437824 1.724770e-18 6.61027236 -6.1875571 6.110371e-10
sl_                       1.026466471 2.791186e+00 2.71327695  0.3783125 7.051985e-01
log_sl_                   0.012494715 1.012573e+00 0.01148456  1.0879576 2.766138e-01
cos_ta_                  -0.027830358 9.725533e-01 0.04475755 -0.6218025 5.340717e-01
sl_:elevation_start       7.773857638 2.377626e+03 2.49443097  3.1164854 1.830208e-03
log_sl_:elevation_start   0.005674556 1.005691e+00 0.01104582  0.5137286 6.074418e-01
cos_ta_:elevation_start  -0.016176775 9.839534e-01 0.04404451 -0.3672825 7.134083e-01

# When having a look on the overall movement map, it is clearly visible, Cupi is taking a different route while crossing 
# the Alps. Therefore this segment will be analysed additionally.

### Cupi - Alps ###
###################
# Clean points from NA times and select relevant cols
cupiSP2 <- readOGR('cupiAlps.shp')
cupiSP2 <- intersect(cupiSP2, area)
r3 <- as.data.frame(cupiSP2)
p3 <- data.frame(r3$UTC_dateti, r3$Latitude, r3$Longitude)
names(p3) <- c('time','y','x')

p3[2:3] <- lapply(p3[2:3],as.numeric) 

p3$time <- as.POSIXct(p3$time,format="%Y-%m-%d %H:%M:%S",tz="UTC")

pm3 <- move(x=p3$x, y=p3$y, 
            time=as.POSIXct(p3$time, format="%Y-%m-%d %H:%M:%S", tz="UTC"), 
            data=p3, proj=CRS("+proj=longlat +ellps=WGS84"), sensor="GPS")


tl3 <- timeLag(x = pm3, units="mins")
summary(tl3)

# Filter date range
p3 <- p3 %>%
  dplyr::filter(time > as.POSIXct("2021-04-01 00:00:00", tz='UTC'))
p3 <- p3 %>%
  dplyr::filter(time < as.POSIXct("2021-06-30 23:59:00", tz='UTC'))

# Convert Movement Data
p_track3 <- amt::make_track(tbl = p3, .x = x, .y = y, .t = time,
                            crs = sp::CRS('+proj=longlat +datum=WGS84 +no_defs '))

# get a summary of the tracks
sum3 <- summarize_sampling_rate(p_track3)
sum3

# Format and generate Random Steps
ssf_dat3 <- p_track3 %>%
  track_resample(rate = minutes(15), tolerance = seconds(2))%>%
  # Filter a minimum number of steps per track
  amt::filter_min_n_burst(min_n = 3)  %>%
  # Create steps from the bursts/fixes
  steps_by_burst()%>%
  # calculate random steps to compare with the actual steps -> 15 control step
  random_steps(n_control = 15)%>%
  extract_covariates(env, where = "both") %>%
  mutate(elevation_start = scale(elevation_start), 
         elevation_end = scale(elevation_end),
         slope_start = scale(slope_start),
         slope_end = scale(slope_end),
         aspect_start = scale(aspect_start), 
         aspect_end = scale(aspect_end), 
         siteD_start = scale(siteD_start),
         siteD_end = scale(siteD_end),
         cos_ta_ = cos(ta_), 
         log_sl_ = log(sl_)
  ) %>% 
  filter(!is.na(ta_))

### Fit Step Selection Function ###

model3 <- ssf_dat3 %>%
  fit_issf(case_ ~ elevation_end + slope_end + aspect_end + siteD_end + sl_ + log_sl_ + cos_ta_ +
             elevation_start:(sl_ + log_sl_ + cos_ta_)+
             strata(step_id_), model = T)


### Results of Step Selection Function for one Bird - Rupert ### 
summary(model3)

SSF_Coeffs3 <- summary(model3)$coefficients
print(SSF_Coeffs3)

# Results Cupi - Alps

coef    exp(coef)    se(coef)          z     Pr(>|z|)
elevation_end            -0.97168534 3.784447e-01  0.75752412 -1.2827121 0.1995929648
slope_end                -2.10163749 1.222561e-01  0.59611923 -3.5255321 0.0004226329
aspect_end                0.06640159 1.068656e+00  0.16210949  0.4096095 0.6820924179
siteD_end                10.97422583 5.835065e+04  3.88357957  2.8258017 0.0047162448
sl_                      -2.20174172 1.106103e-01 14.40764378 -0.1528176 0.8785421085
log_sl_                   0.06992334 1.072426e+00  0.07117133  0.9824650 0.3258708031
cos_ta_                  -0.04909028 9.520952e-01  0.22581343 -0.2173931 0.8279020256
sl_:elevation_start     -36.05265145 2.200556e-16 15.79365403 -2.2827302 0.0224462641
log_sl_:elevation_start   0.07255029 1.075247e+00  0.07414923  0.9784361 0.3278586583
cos_ta_:elevation_start  -0.44068933 6.435926e-01  0.22578262 -1.9518302 0.0509583650


#---------------------------------------------------------------------------------------------------------------------
### Discussion ###
##################

# Distance to breeding site seems to have a high influence on the decision of the Bird where to fly to
# actually reach this site, even when they have to cross the Alpine mountan ridge.
# Elevation seems having an additional slightly lower influence, as can bee seen in the probability map prior to
# the step selection function.
# Step selection functions will need more clearance in the interpretation of results.
# All in all having a look at a probability map, where birds might fly based on static variables is giving a first
# insight, what to expect from more detailed analyses like step selection functions.
# Nevertheless small-scale variables influencing small-scale movement and flight decisions in small mountanious valley
# should be added in a temporally small resolution.

### Outlook ###
###############

# Potential additional variables:

# wind speed (needs to be fitted with resolution as low as possible)
# wind direction (needs to be fitted with resolution as low as possible)
# sun radiation from ground/surface temperature above ground (needs to be fitted with resolution as low as possible)

### Analysis based on:
# https://conservancy.umn.edu/bitstream/handle/11299/218272/AppB_SSF_examples.html?sequence=26#Continuous_Movement_Predictors
# Signer et al. 2019: Animal movement tools (amt): R package for managing tracking data and conducting habitat selection analyses; https://doi.org/10.1002/ece3.4823

### Movement Data: Waldrappteam Conservation & Research; contact: hwehner@waldrapp.eu