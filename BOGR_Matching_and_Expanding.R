##################################
# Purpose: example code for using rMultivariateMatching::multivarmatch to 
#     identify sites with similar climate to a set of coordinates
#
#
#
# Rachel R. Renne
# April 27, 2022
#
##################################


# Install and/or load package
#install.packages(c("devtools","rmarkdown"))
#devtools::install_github("DrylandEcology/rMultivariateMatching", build_vignettes = TRUE)
library(rMultivariateMatching)

# Load raster package
library(raster)

# Read in coordinates
## Coordinates of BOGR patches
#bogr <- read.csv("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/Field_Methods/GPS_Points/BOGR_Points_clean.csv")
bogr <- read.csv("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/Field_Methods/GPS_Points/Final_BOGR_Points.csv")

# Read in raster
# This is just an extent of the area you are interested in (should cover all my STEPWAT points)
ugrb <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/Field_Methods/More_Searching/GIS_Data/Soil_Shapefiles/Rasters/soil_size_tif.tif")
# Read in Rachel's raster
# This is the extent of the sagebrush region Rachel made (FROM WHAT SOURCE: Renne et al. in reivew)
art <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/Field_Methods/More_Searching/Rachels_Files/CoreARTR_combined_DayMet_cropped_trimmed.tif")
# Reproject ugrb raster
projugrb <- projectRaster(ugrb, art, method = 'ngb')


# Read in Wyoming-specific data
data(targetcells)
raster_template <- targetcells[[1]]


# Convert coordinates to SpatialPoints
sp_pts <- SpatialPoints(bogr[,c("Longitude","Latitude")], proj4string = crs(projugrb))

# Extract raster to points, get cellnumbers
sp_pts1 <- extract(raster_template, sp_pts, method = "simple", cellnumbers = TRUE)

# Extract data from raster stack "targetcells"
## This is from Thornton et al. 2018 and Hengl et al. 2017 (use ?targetcells)
## Climate data is from 1980-2010
bioclim <- makeInputdata(targetcells)
rownames(bioclim) <- bioclim$cellnumbers

# Restrict data to matching variables of interest
matchingvars <- bioclim[,c("cellnumbers","x","y","bioclim_01","bioclim_10",
                           "bioclim_11","bioclim_12","bioclim_18","bioclim_19", 
                           "sand", "clay")]
## Site Rachel's Multivariate Matching paper (in review, same as before)
# bioclim_01 = Mean Annual Temperature
# bioclim_10 = Mean Temperature of Warmest Quarter
# bioclim_11 = Mean Temperature of Coldest Quarter
# bioclim_12 = Total Annual Precipitation
# bioclim_18 = Precip of warmest quarter
# bioclim_19 = Precipitation of Coldest Quarter

# Figure out which cells correspond to your sites
bogrvars <- matchingvars[rownames(matchingvars) %in% sp_pts1[,1],]

# Determine what range of values should be
fullrange <- c()
for(i in c(4:11)){
  ranges <- range(bogrvars[i])[2] - range(bogrvars[i])[1]
  fullrange <- append(fullrange, ranges)
}

# Create list of range of value for each variable in bogrvars:
#fullrange2 <- c(0.49,0.72,0.492,76.2,4.29,28.07,10.86,6.26)

# Create criteria_list (for both the full range of values and 10% of the full range)
criteria_list100 <- fullrange
criteria_list10 <- fullrange * 0.1
criteria_list50 <- fullrange * 0.5
criteria_list5 <- fullrange * 0.05
criteria_list25 <- fullrange * 0.25




# Use multivarmatch to find similar sites

# Run matching
potbogrsites <- multivarmatch(matchingvars = na.omit(matchingvars), subsetcells = bogrvars[,c("x","y")],
                              matchingvars_id = "cellnumbers", 
                              subsetcells_id = "cellnumbers",
                              criteria = criteria_list5, n_neighbors = 1,
                              raster_template = art,
                              subset_in_target = TRUE, saveraster = FALSE,
                              plotraster = TRUE, addpoints = FALSE)
points(sp_pts1)                    

# Rasterize output (takes a moment)
tosearchraster <- rasterize(x = potbogrsites[,c("x","y")], y = art, field = potbogrsites$matching_quality)


# Look at map (I have it zoomed in to SW Wyoming)
par(mar = c(1,1,1,1))
image(tosearchraster, col = rev(c("#d7191c", "#fdae61", "#abd9e9","#2c7bb6")),
      breaks = c(0, 0.5, 1, 1.5, 5), xaxt = "n", yaxt = "n", bty = "n",
      ylim = c(41.029271, 43.248840), xlim = c(-111.018772, -107.607517))
points(sp_pts, pch = 16)

# Write raster
writeRaster(tosearchraster, filename = "potBOGRsites_5%.tif")

# Let's say you just want the cells that are within 1 unit of matching
# distance from your BOGR sites (note the lat/long here will be centroids of cells)
best_matches <- potbogrsites[potbogrsites$matching_quality <= 5,]

# See how many
dim(best_matches)

# Calculate the area that is matched this well
areas <- area(projugrb)
ptsx <- SpatialPoints(best_matches[,c("x","y")], proj4string = crs(projugrb))
sum(extract(areas, ptsx))
# The output will be in square km






####### Use other Wyoming-specific dataset
# Load targetcells data for Target Cells
data(targetcells)
# Create data frame of potential matching variables for Target Cells
allvars <- makeInputdata(targetcells)





###################################################################################################
#### Interpolation Station #####


#### TO DO ###
# merge stepwatpoints with stepwat output
# Goal: a data frame containing simulation output and column with name 'cellnumbers' (basically should be stepwatvars df but instead of bioclim it has stepwat output)

stepwatpoints_full<- read.csv("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/Field_Methods/GPS_Points/STEPWAT_Points/STEPWAT_Points_Final_with_Cellnumbers.csv")
stepwatpoints<-stepwatpoints_full[-c(11,12,17,18,36),]
colnames(stepwatpoints)[2] <- "Site"
colnames(stepwatpoints)[24] <- "cellnumbers"
sitematch <- stepwatpoints[,c(2,24)]
soil_data <- stepwatpoints[,c(2,8:23)]
clean.points <- stepwatpoints[,-c(1,8:23)]

# Calculate average soil data for each site
#10(0-10)+10(10-20)+10(20-30)+10(30-40)+20(40-60)+20(60-80)+20(80-100)+50(100-150) # based on Hengl et al. 2017 from ?targetcells
soil.averages<-data.frame()
for(i in c(1:length(soil_data$Site))){
  new_clay <- ((10*soil_data$clay_0_10[i])+(10*soil_data$clay_10_20[i])+(10*soil_data$clay_20_30[i])+
                 (10*soil_data$clay_30_40[i])+(20*soil_data$clay_40_60[i])+(20*soil_data$clay_60_80[i])+
                 (20*soil_data$clay_80_100[i])+(50*soil_data$clay_100_150[i]))/150
  new_sand <- ((10*soil_data$sand_0_10[i])+(10*soil_data$sand_10_20[i])+(10*soil_data$sand_20_30[i])+
                 (10*soil_data$sand_30_40[i])+(20*soil_data$sand_40_60[i])+(20*soil_data$sand_60_80[i])+
                 (20*soil_data$sand_80_100[i])+(50*soil_data$sand_100_150[i]))/150
  Site<-soil_data$Site[i]
  soils<-cbind(Site,new_sand, new_clay)
  soil.averages<-rbind(soil.averages,soils)
}

site.type.soil <- merge(soil.averages, stepwatpoints[,c(2,4)], by = "Site")
site.type.soil$Functional[site.type.soil$Functional=="past"]<- "unn"
sand <- aggregate(site.type.soil$new_sand, by=list(site.type.soil$Functional), FUN=mean)
clay <- aggregate(site.type.soil$new_clay, by=list(site.type.soil$Functional), FUN=mean)

# Read in dfs
full_stats <- read.csv("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/Stepwat_Stats_FULL.csv")
current_stats <- full_stats[which(full_stats$GCM == "Current"),c(2,5,7)]
GCM_med<-read.csv("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/FULL_GCM_Median_Biomass.csv")
#colnames(GCM_med)[2]<- "site"
#colnames(current_stats)[1]<-"site"

# Merge
current.merge <- merge(current_stats, clean.points, by = "Site")
GCM.med.merge <- merge(GCM_med, clean.points, by = "Site")
GCM.med.merge <- GCM.med.merge[,-2]


# First do current data
# Convert coordinates to SpatialPoints
sp_pts2 <- SpatialPoints(current.merge[,c("longitude","latitude")], proj4string = crs(targetcells[[1]]))

# Extract raster to points, get cellnumbers
sp_pts2.2 <- extract(targetcells[[1]], sp_pts2, method = "simple", cellnumbers = TRUE)


## Establish matchingvars
bioclim <- makeInputdata(targetcells)
rownames(bioclim) <- bioclim$cellnumbers
matchingvars <- bioclim[,c("cellnumbers","x","y","bioclim_01","bioclim_10",
                           "bioclim_11","bioclim_12","bioclim_18","bioclim_19", 
                           "sand", "clay")]

# Figure out which cells correspond to your sites
stepwatvars <- matchingvars[rownames(matchingvars) %in% sp_pts2.2[,1],]

# Replace soil data in stepwatvars with data from SSURGO
soil_replacement <- merge(soil.averages, sitematch, by = "Site")
stepwatvars.updated <- merge(stepwatvars, soil_replacement, by = "cellnumbers")
stepwatvars.updated["sand"] <- stepwatvars.updated["new_sand"]
stepwatvars.updated["clay"] <- stepwatvars.updated["new_clay"]
stepwatvars.updated <- stepwatvars.updated[,-c(13,14)]


# Add biomass data to climatic variables
stepwat.clim.bmass1<-merge(current.merge, stepwatvars.updated, by = "cellnumbers")
stepwat.clim.bmass1 <- stepwat.clim.bmass1[,-c(5,6,9,10,11,20)]
colnames(stepwat.clim.bmass1)[2]<- "Site"
stepwat.clim.bmass<- stepwat.clim.bmass1[,c(1,2,5,6,3,4,7:14)]

# Determine what range of environmental values should be
fullrange2 <- c()
for(i in c(7:14)){
  ranges2 <- range(stepwat.clim.bmass[i])[2] - range(stepwat.clim.bmass[i])[1]
  fullrange2 <- append(fullrange2, ranges2)
}

criteria_list5.2 <- fullrange2 * 0.05

# Create subsetcells for multivarmatch
subsetcells <- unique(stepwat.clim.bmass[,c("longitude","latitude","cellnumbers")])
rownames(subsetcells) <- subsetcells$cellnumbers

# Run matching on environmental variables
stepwatsites <- multivarmatch(matchingvars = na.omit(matchingvars), subsetcells = subsetcells,
                              matchingvars_id = "cellnumbers", 
                              subsetcells_id = "cellnumbers",
                              criteria = criteria_list5.2, n_neighbors = 1,
                              raster_template = targetcells[[1]],
                              subset_in_target = TRUE, saveraster = TRUE,
                              plotraster = TRUE, addpoints = FALSE)


points(sp_pts2)

#Output rasters for all rgroups
for(i in unique(stepwat.clim.bmass$Rgroup)){
output_results<-stepwat.clim.bmass[which(stepwat.clim.bmass$Rgroup==i),c(1,6)]
rownames(output_results) <- output_results$cellnumbers

# Interpolate simulation output to rasters
interpolatePoints(matches = stepwatsites, output_results = output_results,
                  exclude_poor_matches = FALSE,
                  subset_cell_names = "subset_cell",
                  quality_name = "matching_quality",
                  matching_distance = 50, raster_template = targetcells[[1]],
                  plotraster = TRUE, overwrite = TRUE, filepath = paste0("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/Interpolation_Rasters/Current/",i))

}


################
# Now same thing for future data

# Convert coordinates to SpatialPoints
sp_pts2 <- SpatialPoints(GCM.med.merge[,c("longitude","latitude")], proj4string = crs(targetcells[[1]]))

# Extract raster to points, get cellnumbers
sp_pts2.2 <- extract(targetcells[[1]], sp_pts2, method = "simple", cellnumbers = TRUE)

# Figure out which cells correspond to your sites
stepwatvars <- matchingvars[rownames(matchingvars) %in% sp_pts2.2[,1],]

# Replace soil data in stepwatvars with data from SSURGO
soil_replacement <- merge(soil.averages, sitematch, by = "Site")
stepwatvars.updated <- merge(stepwatvars, soil_replacement, by = "cellnumbers")
stepwatvars.updated["sand"] <- stepwatvars.updated["new_sand"]
stepwatvars.updated["clay"] <- stepwatvars.updated["new_clay"]
stepwatvars.updated <- stepwatvars.updated[,-c(13,14)]

# Add biomass data to climatic variables
stepwat.clim.bmass1<-merge(GCM.med.merge, stepwatvars.updated, by = "cellnumbers")
stepwat.clim.bmass1 <- stepwat.clim.bmass1[,-c(7,8,9,12,13,14,23)]
stepwat.clim.bmass<- stepwat.clim.bmass1[,c(1,2,7,8,3,4,5,6,9:16)]
colnames(stepwat.clim.bmass)[2]<-"Site"

# Determine what range of environmental values should be
fullrange2 <- c()
for(i in c(9:16)){
  ranges2 <- range(stepwat.clim.bmass[i])[2] - range(stepwat.clim.bmass[i])[1]
  fullrange2 <- append(fullrange2, ranges2)
}

criteria_list5.2 <- fullrange2 * 0.05

# Create subsetcells for multivarmatch
subsetcells <- unique(stepwat.clim.bmass[,c("longitude","latitude","cellnumbers")])
rownames(subsetcells) <- subsetcells$cellnumbers

# Run matching on environmental variables
stepwatsites <- multivarmatch(matchingvars = na.omit(matchingvars), subsetcells = subsetcells,
                              matchingvars_id = "cellnumbers", 
                              subsetcells_id = "cellnumbers",
                              criteria = criteria_list5.2, n_neighbors = 1,
                              raster_template = targetcells[[1]],
                              subset_in_target = TRUE, saveraster = TRUE,
                              plotraster = TRUE, addpoints = FALSE)
points(sp_pts2)


#Output rasters for all rgroups
for(p in c("Mid", "Late")){
  for(r in c("RCP45", "RCP85")){
    for(i in unique(stepwat.clim.bmass$Rgroup)){
  output_results<-stepwat.clim.bmass[which(stepwat.clim.bmass$Period==p & stepwat.clim.bmass$RCP==r & stepwat.clim.bmass$Rgroup==i),c(1,8)]
  rownames(output_results) <- output_results$cellnumbers
  
  # Interpolate simulation output to rasters
  interpolatePoints(matches = stepwatsites, output_results = output_results,
                    exclude_poor_matches = FALSE,
                    subset_cell_names = "subset_cell",
                    quality_name = "matching_quality",
                    matching_distance = 50, raster_template = targetcells[[1]],
                    plotraster = TRUE, overwrite = TRUE, filepath = paste0("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/Interpolation_Rasters/",p,"/",r,"/",i))
  
}}}
### Code that made it possible to add cellnumbers to original csv
# cnrast<-rasterFromXYZ(stepwatvars[c(2,3,1)], crs = crs(targetcells[[1]]))
# 
# writeRaster(cnrast, filename = "cell_numbers.tif")
# plot(site_soils_avg[2:3])


############################################################################
## Make plots of rasters
library(ggplot2)
library(raster)
library(rasterVis)
library(rgdal)
library(grid)
library(scales)
library(viridis)
library(terra)
library(ggthemes)
library(colorspace)



ugrb_extent<- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/UGRB.tif")
xlimits<-c(xmin(ugrb_extent),xmax(ugrb_extent))
ylimits<-c(ymin(ugrb_extent),ymax(ugrb_extent))
matchrast <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Matchingquality.tif")
matchrast[matchrast <= 10] <- 1
matchrast[matchrast > 10] <- 0



## Sagebrush Current
sage.current <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/Interpolation_Rasters/Current/sagebrush/M.Bmass.tif")
sage.current.clip <- sage.current * matchrast

# Plotting in ggplot
sage.current.clip.points<-rasterToPoints(sage.current.clip, spatial = TRUE)
sage.current.clip.df <- data.frame(sage.current.clip.points)
sage.current.clip.df2<-sage.current.clip.df[sage.current.clip.df$layer>0,]

sage.current.map<-ggplot(sage.current.clip.df2)+
  geom_tile(aes(x=x, y=y, fill = layer))+
  coord_sf(xlim = xlimits, ylim = ylimits)+
  #scale_fill_gradient(low = "blue", high = "red")+
  labs(fill = "")+
  scale_fill_viridis(na.value = "white")+
  labs(title = expression("Big Sagebrush Biomass (g/m"^2*")"))+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(plot.title = element_text(size=14))+
  theme(legend.position = c(0.9,0.8))+
  theme(legend.key.size = unit(0.53,'cm'))



## C3 Current
c3.current <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/Interpolation_Rasters/Current/p.cool.grass/M.Bmass.tif")
c3.current.clip <- c3.current * matchrast

# Plotting in ggplot
c3.current.clip.points<-rasterToPoints(c3.current.clip, spatial = TRUE)
c3.current.clip.df <- data.frame(c3.current.clip.points)
c3.current.clip.df2<-c3.current.clip.df[c3.current.clip.df$layer>0,]

c3.current.map<-ggplot(c3.current.clip.df2)+
  geom_tile(aes(x=x, y=y, fill = layer))+
  coord_sf(xlim = xlimits, ylim = ylimits)+
  #scale_fill_gradient(low = "blue", high = "red")+
  labs(fill = "")+
  scale_fill_viridis(na.value = "white")+
  labs(title = expression("C3 Grass Biomass (g/m"^2*")"))+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(plot.title = element_text(size=14))+
  theme(legend.position = c(0.9,0.8))+
  theme(legend.key.size = unit(0.53,'cm'))



## C4 Current
c4.current <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/Interpolation_Rasters/Current/p.warm.grass/M.Bmass.tif")
c4.current.clip <- c4.current * matchrast

# Plotting in ggplot
c4.current.clip.points<-rasterToPoints(c4.current.clip, spatial = TRUE)
c4.current.clip.df <- data.frame(c4.current.clip.points)
c4.current.clip.df2<-c4.current.clip.df[c4.current.clip.df$layer>0,]

c4.current.map<-ggplot(c4.current.clip.df2)+
  geom_tile(aes(x=x, y=y, fill = layer))+
  coord_sf(xlim = xlimits, ylim = ylimits)+
  #scale_fill_gradient(low = "blue", high = "red")+
  labs(fill = "")+
  scale_fill_viridis(na.value = "white")+
  labs(title = expression("C4 Grass Biomass (g/m"^2*")"))+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(plot.title = element_text(size=14))+
  theme(legend.position = c(0.9,0.8))+
  theme(legend.key.size = unit(0.53,'cm'))


## Annual Forb Current
af.current <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/Interpolation_Rasters/Current/a.forb/M.Bmass.tif")
af.current.clip <- af.current * matchrast

# Plotting in ggplot
af.current.clip.points<-rasterToPoints(af.current.clip, spatial = TRUE)
af.current.clip.df <- data.frame(af.current.clip.points)
af.current.clip.df2<-af.current.clip.df[af.current.clip.df$layer>0,]

af.current.map<-ggplot(af.current.clip.df2)+
  geom_tile(aes(x=x, y=y, fill = layer))+
  coord_sf(xlim = xlimits, ylim = ylimits)+
  #scale_fill_gradient(low = "blue", high = "red")+
  labs(fill = "")+
  scale_fill_viridis(na.value = "white")+
  labs(title = expression("Annual Forb Biomass (g/m"^2*")"))+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(plot.title = element_text(size=14))+
  theme(legend.position = c(0.9,0.8))+
  theme(legend.key.size = unit(0.53,'cm'))



## Perennial Forb Current
pf.current <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/Interpolation_Rasters/Current/p.forb/M.Bmass.tif")
pf.current.clip <- pf.current * matchrast

# Plotting in ggplot
pf.current.clip.points<-rasterToPoints(pf.current.clip, spatial = TRUE)
pf.current.clip.df <- data.frame(pf.current.clip.points)
pf.current.clip.df2<-pf.current.clip.df[pf.current.clip.df$layer>0,]

pf.current.map<-ggplot(pf.current.clip.df2)+
  geom_tile(aes(x=x, y=y, fill = layer))+
  coord_sf(xlim = xlimits, ylim = ylimits)+
  #scale_fill_gradient(low = "blue", high = "red")+
  labs(fill = "")+
  scale_fill_viridis(na.value = "white")+
  labs(title = expression("Perennial Forb Biomass (g/m"^2*")"))+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(plot.title = element_text(size=14))+
  theme(legend.position = c(0.9,0.8))+
  theme(legend.key.size = unit(0.53,'cm'))


## Shrub Current
shrub.current <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/Interpolation_Rasters/Current/shrub/M.Bmass.tif")
shrub.current.clip <- shrub.current * matchrast

# Plotting in ggplot
shrub.current.clip.points<-rasterToPoints(shrub.current.clip, spatial = TRUE)
shrub.current.clip.df <- data.frame(shrub.current.clip.points)
shrub.current.clip.df2<-shrub.current.clip.df[pf.current.clip.df$layer>0,]

shrub.current.map<-ggplot(shrub.current.clip.df2)+
  geom_tile(aes(x=x, y=y, fill = layer))+
  coord_sf(xlim = xlimits, ylim = ylimits)+
  #scale_fill_gradient(low = "blue", high = "red")+
  labs(fill = "")+
  scale_fill_viridis(na.value = "white")+
  labs(title = expression("Shrub Biomass (g/m"^2*")"))+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(plot.title = element_text(size=14))+
  theme(legend.position = c(0.9,0.8))+
  theme(legend.key.size = unit(0.53,'cm'))+
  theme(plot.margin = margin(0,0,0,0, "cm"))


## Sagebrush RCP45 Mid
sage.45m <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/Interpolation_Rasters/Mid/RCP45/sagebrush/M.M.Bmass.tif")

# Plotting in ggplot
sage.45m.points<-rasterToPoints(sage.45m, spatial = TRUE)
sage.45m.df <- data.frame(sage.45m.points)

sage.45m.map<-ggplot(sage.45m.df)+
  geom_tile(aes(x=x, y=y, fill = M.M.Bmass))+
  coord_sf(xlim = xlimits, ylim = ylimits)+
  #scale_fill_gradient(low = "blue", high = "red")+
  labs(fill = "")+
  scale_fill_continuous_diverging(palette = "Purple-Green")+
  labs(title = expression(Delta ~ "Big Sagebrush Biomass (g/m"^2*")"))+
  xlab("")+a
  ylab("")+
  theme_bw()+
  theme(plot.title = element_text(size=14))+
  theme(legend.position = c(0.89,0.79))+
  theme(legend.key.size = unit(0.53,'cm'))


## C3 RCP45 Mid
c3.45m <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/Interpolation_Rasters/Mid/RCP45/p.cool.grass/M.M.Bmass.tif")

# Plotting in ggplot
c3.45m.points<-rasterToPoints(c3.45m, spatial = TRUE)
c3.45m.df <- data.frame(c3.45m.points)

c3.45m.map<-ggplot(c3.45m.df)+
  geom_tile(aes(x=x, y=y, fill = M.M.Bmass))+
  coord_sf(xlim = xlimits, ylim = ylimits)+
  labs(fill = "")+
  scale_fill_continuous_diverging(palette = "Purple-Green")+
  labs(title = expression(Delta ~ "C3 Grass Biomass (g/m"^2*")"))+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(plot.title = element_text(size=14))+
  theme(legend.position = c(0.89,0.79))+
  theme(legend.key.size = unit(0.53,'cm'))


## C4 RCP45 Mid
c4.45m <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/Interpolation_Rasters/Mid/RCP45/p.warm.grass/M.M.Bmass.tif")

# Plotting in ggplot
c4.45m.points<-rasterToPoints(c4.45m, spatial = TRUE)
c4.45m.df <- data.frame(c4.45m.points)

c4.45m.map<-ggplot(c4.45m.df)+
  geom_tile(aes(x=x, y=y, fill = M.M.Bmass))+
  coord_sf(xlim = xlimits, ylim = ylimits)+
  labs(fill = "")+
  scale_fill_continuous_diverging(palette = "Purple-Green")+
  labs(title = expression(Delta ~ "C4 Grass Biomass (g/m"^2*")"))+
  xlab("")+
  ylab("")+
  theme(plot.title = element_text(size=14))+
  theme_bw()+
  theme(legend.position = c(0.89,0.79))+
  theme(legend.key.size = unit(0.53,'cm'))


## Perennial Forb RCP45 Mid
pf.45m <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/Interpolation_Rasters/Mid/RCP45/p.forb/M.M.Bmass.tif")

# Plotting in ggplot
pf.45m.points<-rasterToPoints(pf.45m, spatial = TRUE)
pf.45m.df <- data.frame(pf.45m.points)

pf.45m.map<-ggplot(pf.45m.df)+
  geom_tile(aes(x=x, y=y, fill = M.M.Bmass))+
  coord_sf(xlim = xlimits, ylim = ylimits)+
  labs(fill = "")+
  scale_fill_continuous_diverging(palette = "Purple-Green")+
  labs(title = expression(Delta ~ "Perennial Forb Biomass (g/m"^2*")"))+
  xlab("")+
  ylab("")+
  theme(plot.title = element_text(size=14))+
  theme_bw()+
  theme(legend.position = c(0.89,0.79))+
  theme(legend.key.size = unit(0.53,'cm'))


## Annual Forb RCP45 Mid
af.45m <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/Interpolation_Rasters/Mid/RCP45/a.forb/M.M.Bmass.tif")

# Plotting in ggplot
af.45m.points<-rasterToPoints(af.45m, spatial = TRUE)
af.45m.df <- data.frame(af.45m.points)

af.45m.map<-ggplot(af.45m.df)+
  geom_tile(aes(x=x, y=y, fill = M.M.Bmass))+
  coord_sf(xlim = xlimits, ylim = ylimits)+
  labs(fill = "")+
  scale_fill_continuous_diverging(palette = "Purple-Green")+
  labs(title = expression(Delta ~ "Annual Forb Biomass (g/m"^2*")"))+
  xlab("")+
  ylab("")+
  theme(plot.title = element_text(size=14))+
  theme_bw()+
  theme(legend.position = c(0.89,0.79))+
  theme(legend.key.size = unit(0.53,'cm'))


## Shrub RCP45 Mid
shrub.45m <- raster("C:/Users/ander/OneDrive - Yale University/Documents/Lauenroth_Lab/STEPWAT/FINAL_OUTPUTS/Interpolation_Rasters/Mid/RCP45/shrub/M.M.Bmass.tif")

# Plotting in ggplot
shrub.45m.points<-rasterToPoints(shrub.45m, spatial = TRUE)
shrub.45m.df <- data.frame(shrub.45m.points)

shrub.45m.map <- ggplot(shrub.45m.df)+
  geom_tile(aes(x=x, y=y, fill = M.M.Bmass))+
  coord_sf(xlim = xlimits, ylim = ylimits)+
  labs(fill = "")+
  scale_fill_continuous_diverging(palette = "Purple-Green")+
  labs(title = expression(Delta ~ "Shrub Biomass (g/m"^2*")"))+
  xlab("")+
  ylab("")+
  theme(plot.title = element_text(size=14))+
  theme_bw()+
  theme(legend.position = c(0.89,0.79))+
  theme(legend.key.size = unit(0.53,'cm'))



# Combine functional types we want
# Combine all plots into one
mainfig<-ggarrange(sage.current.map, c3.current.map, c4.current.map, sage.45m.map, c3.45m.map, c4.45m.map, nrow=2, ncol=3,common.legend = F, labels = c("(A)", "(B)", "(C)","(D)","(E)","(F)"), align = "hv")
#ggexport(all.current.maps,filename = "current.maps.tiff", width = 1024, height = 768)

png(filename = "Interpolation_Main.png", width = 1000, height = 746)
print(mainfig)
dev.off()

# Combine less important functional types
# Combine all plots into one
supfig<-ggarrange(af.current.map, pf.current.map, shrub.current.map, af.45m.map, pf.45m.map, shrub.45m.map, nrow=2, ncol=3,common.legend = F,labels = c("(A)", "(B)", "(C)","(D)","(E)","(F)"), align = "hv")
#ggexport(all.45l.maps,filename = "45l.maps.tiff", width = 1024, height = 768)

png(filename = "Interpolation_Supplemental.png", width = 1000, height = 746)
print(supfig)
dev.off()