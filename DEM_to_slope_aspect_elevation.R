# clear workspace (DO NOT RUN IF YOU WANT THINGS IN YOUR WORKSPACE)
rm(ls=list())
# load libraries
library(raster)
library(rgdal)
library(ggplot2)
library(dplyr)

source("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/rSFSW2_Files/WeatherDB.R")

# this code needs a raster DEM that emconpasses your points
# I got mine off of https://apps.nationalmap.gov/downloader/

# Load DEM raster into R
DEM <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/Climate_Database/DEM_Rasters/Full_UGRB.tif")

# View raster structure
DEM

# plot DEM
plot(DEM)

# read in latitude and longitude points
pnts<-read.csv("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/Climate_Database/Final_UGRB_Current/1_Input/SWRuns_InputMain_UGRB_Points_v12.csv")

# make sure the dataframe includes just the latitude and longitude (TEM cleaning up her dataset)
pnts2<-pnts[,c("X_WGS84","Y_WGS84")]

# make sure there are no NAs
pnts2<-na.omit(pnts2)

# IF you longitude is not negative this will add the negative sign
#  *do not run if your longitude is already negative (e.g. -105.3645)
#pnts2$Long<-pnts2$Long*-1

# convert your lat/long coordinates into a spatial points dataframe
spnts<-SpatialPoints(coords = cbind(pnts2$X_WGS84,pnts2$Y_WGS84),
              proj4string=CRS(as.character("+proj=longlat +datum=NAD83 +no_defs")))

# plot your points on top of the DEM
plot(DEM)
plot(spnts,add=T)

# extract elevation data from DEM (if you get NAs it might mean your points are
#    off the DEM)
elev<-extract(DEM,spnts)

# calculate slope using the DEM and plot
slope <- terrain(DEM, 'slope', unit='degrees', neighbors=8)
plot(slope)

# calculate aspect using the DEM and plot
aspect <- terrain(DEM, 'aspect', unit='degrees', neighbors=8)
plot(aspect)

# create a raster stack and check to make sure all 3 layers are there
stk.pnts<-stack(slope,aspect,DEM)
stk.pnts

# extract the slope, aspect and DEM from the raster stack at your points
stuff<-extract(stk.pnts,spnts)
stuff2<-cbind(stuff,pnts2)

# rename columns
names(stuff2)<-c("slope","aspect","elevation","X_WGS84","Y_WGS84")
#write.csv(stuff2,"C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/Climate_Database/Test/1_Input/slope_aspect_elev.csv")

# merge dataframe with slope aspect and elev with sheet used in climate database
full_df <- merge(stuff2, pnts, by= c("X_WGS84", "Y_WGS84"))

# clean dataframe to match the original format
full_df2 <- full_df[1:9]
names(full_df2)<-c("X_WGS84", "Y_WGS84", "Slope", "Aspect", "ELEV_m", "Label", "site_id", "Include_YN", "WeatherFolder")
final_df <- full_df2[,c(6,7,8,9,1,2,5,3,4)]
final_df <- final_df[order(final_df$site_id),]
final_df["dailyweather_source"] <- "DayMet_NorthAmerica"
final_df['Include_YN_DailyWeather'] <- 1

# overwrite original csv
write.csv(final_df,"C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/Climate_Database/UGRB_Current/1_Input/SWRuns_InputMain_UGRB_Points_v12.csv", row.names = FALSE)
