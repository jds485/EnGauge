#Script to read USGS streamflow and water quality data from gauges
# Recommended resource - Package Readme: 
#  (located somewhere on your computer after installing the package) dataRetrieval.html
# Function help at: https://github.com/USGS-R/dataRetrieval
#             Blog: https://owi.usgs.gov/R/training-curriculum/usgs-packages/
#        Slideshow: https://owi.usgs.gov/R/dataRetrieval.html#1

#Streamflow Gauge Retrieval Method 1:
#The function whatNWISsites() will return gauge numbers
# with coordinates for a specified state (county, or various other criteria)
# Search Criteria for whatNWISsites(): Table 1 of https://www.waterqualitydata.us/webservices_documentation/
#                                       References to "domain service" indicate to look at Table 2.
#                  characteristicName: https://www.waterqualitydata.us/public_srsnames/
#                   Likely the most common filtering method. 
#                  Example: whatNWISsites(statecode="MD", characteristicName="Phosphorus")
#  Note: their function must use something like grep when characteristicName is specified. 
#        Searching for Nitrogen returned all instances of Nitrogen on that list.
#         If you want a specific code only, specify the parameterCd (parameter code)
#         instead of the characteristic name
# That function will not return full site information. 
# For full site information, use readNWISsite() with a vector of gauge numbers returned from whatNWISsites().
# NOTE: This method seems more reliable, and returned for my ROI more gauges than Method 2

#Streamflow Gauge Retrieval Method 2:
#You can find USGS gauges of any type in your region of interest on this website:
# https://cida.usgs.gov/enddat/dataDiscovery.jsp
# Record the coordinates of the bounding box for your region of interest! They could be important for the next step. 
# Select the data you want and copy the resulting table of values into a txt file
#  txt is important to retain leading zeros on gauge numbers
# If you inhereted a csv or txt file without leading zeros, a script below can add them to NWIS gauges.
#  Other gauge datasets seem to not require leading zeros.

#Coordinates of gauges are not reported in that download :( 
# The USGS function readNWISsite() can look up the coordinates and all site info, given the gauge numbers.
# You can also find coordinates for your gauges on this site if you provide a bounding box of coordinates:
# https://waterdata.usgs.gov/nwis/inventory?search_criteria=lat_long_bounding_box&submitted_form=introduction

#Comparing altitude of gauge vs. DEM:
# Obtained coordinates may include altitude of the gauge, which can be important vs. DEM elevation 
#  (e.g. if there's a cliff at the gauge vs. DEM mean elevation of the pixel)
#  A plot of DEM vs. Gauge reported elevation is made below to visually detect discrepancies
#   (e.g. USGS data reported in m instead of in ft)
#  Ensure that the units and the vertical datum are the same for your DEM and all gauge altitudes
#   DEMs in the US tend to be NAVD88 in m, whereas gauges tend to be referenced to NGVD29 in ft
#   Differences after unit conversion tend to be minor in the US, 
#    except in the West: https://www.ngs.noaa.gov/TOOLS/Vertcon/vertcon.html
#  Altitude datum codes: https://help.waterdata.usgs.gov/code/alt_datum_cd_query?fmt=html
#   collection method codes: https://help.waterdata.usgs.gov/code/alt_meth_cd_query?fmt=html

#Water Quality Retrieval Method 1:
# You can use the whatWQPsites(), using the same query criteria as the whatNWISsites() function.
# The whatWQPsites() function returns all site information, same as the waterqualitydata website
# described in Method 2.

#Water Quality Retrieval Method 2:
# You can find water quality gauges on this website:
# https://www.waterqualitydata.us/portal/
# Download "site data only" to receive a csv file with gauge/site information with coordinates.
# The MonitoringLocationIdentifier field is used to download data for each gauge/site in the script below.

#Data Quality Codes:
# You can find data quality codes for USGS datasets here:
# https://help.waterdata.usgs.gov/codes-and-parameters/codes#discharge_cd
# You should always look at these quality codes, and process your data accordingly.
# This script colors streamflow time series by error code, but doesn't process further than that.
# codes for streamflow: https://help.waterdata.usgs.gov/codes-and-parameters/daily-value-qualification-code-dv_rmk_cd
#                  AND: https://help.waterdata.usgs.gov/codes-and-parameters/instantaneous-and-daily-value-status-codes
#                  Yes, there are 2 separate reference schemes for streamflow data.

#Author Credits:
# Portions of the function for the streamflow download were provided by:
# Caitline Barber (Caitline.Barber@tufts.edu) and Jonathan Lamontagne (Jonathan.Lamontagne@tufts.edu)
# Modified by Jared Smith (js4yd@virginia.edu) in June, 2019, and started git tracking.
# See git commit history for all later contributions

#Fixme: make each method a separate function or script
#Fixme: Also has a component for downloading weather station data 
#       stored on the USGS database (which includes NOAA ACIS data)

#Set directory names----
#Region of interest shapefile
dir_ROI = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\BES-Watersheds-Land-Cover-Analysis"  
#Color functions - from JDS github repo: Geothermal_ESDA
dir_ColFuns = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\USGSGauges"
#EnGauge repository
dir_EnGauge = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\EnGauge\\EnGauge"
#USGS streamflow gauges
dir_sfgauges = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\USGSGauges"
#DEM - specify as a vector of directories if there are multiple tiles to be mosaicked together.
# The directory order has to match the file name order for f_DEM below.
# Fixme: can DEMs be downloaded from a server instead of downloading manually before using this script?
dir_DEM = c("C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\DEM\\USGS_NED_1_n40w077_ArcGrid\\grdn40w077_1",
            'C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\DEM\\USGS_NED_1_n40w078_ArcGrid\\grdn40w078_1')
#Output for processed DEM
dir_DEM_out = 'C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\DEM\\'
#Water quality gauges
dir_wq = 'C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\USGSGauges'

#Set input filenames----
#Region of interest shapefile name
f_ROI = "Watershed_GF"
#Streamflow gauges and site coordinates filenames
f_StreamGaugeData = "BES_USGS_GaugeStations.csv"
#Optional site data. If not provided, will use function readNWISsite() to collect site coordinates
#f_StreamGaugeSites = "USGS_GaugeSites.txt"
#DEM - all separate DEM tiles may be added to this as a vector (e.g. c("w001001.adf", "w001002.adf") )
f_DEM = c("w001001.adf", "w001001.adf")
#Water quality gauges
f_WQgauges = "BES_WaterQualityGaugeStations.csv"

#Set output filenames----
#DEM - geotiff format is the default. You can change in the script below
f_DEM_mosiac = "DEM_mosaic"
#NWIS streamflow gauges in ROI
f_NWIS_ROI_out = 'NWIS_ROI'
#NWIS streamflow gauges in bounding box
f_NWIS_bb_out = 'NWIS_bb'
#Streamflow gauge data processing name appendage
# e.g. 0159384_p
f_sf_processKey = '_p'

#Set project coordinate system----
#This is the coordinate system that all data will be plotted and written in
# It is not the coordinate system of your data (although it could be)
# EPSG codes from: https://spatialreference.org/ref/?page=2
pCRS = '+init=epsg:26918'

#Load libraries and functions----
#USGS function library - note that a more recent version is available through Github
library(dataRetrieval)
#R libraries
library(stringi)
library(stringr)
library(foreach)
library(doParallel)
library(rgdal)
library(GISTools)
library(raster)
library(rlist)
#Color functions for plots (R script from Jared Smith's Geothermal_ESDA Github Repo)
setwd(dir_ColFuns)
source('ColorFunctions.R')
#Functions from repository
setwd(dir_EnGauge)
source('missingDates.R')
source('addZerosToGaugeNames.R')
source('processDEM.R')

#Streamflow----
setwd(dir_sfgauges)
# Method 1: Get gauges using whatNWISsites()----
AllStations_fn = whatNWISsites(statecode = "MD", parameterCd = '00060')
AllStations_fn = readNWISsite(AllStations_fn$site_no)
#Convert to spatial dataframe
# NOTE: Your data may be all one coordinate system, and therefore not need to split into 2 datasets
#       before joining into 1 dataset.
# NOTE: your coordinate system may be different (epsg code)
# Some of the data are NAD27 projection and others are NAD83 projection. Split the dataset to handle each
#Fixme: function for splitting coordinate systems and returning one same-coordinate system file
#       Would require being able to look up epsg codes.
GaugesLocs_NAD27 = AllStations_fn[AllStations_fn$coord_datum_cd == 'NAD27', ]
coordinates(GaugesLocs_NAD27) = c('dec_long_va', 'dec_lat_va')
proj4string(GaugesLocs_NAD27) = CRS('+init=epsg:4267')
GaugesLocs_NAD83 = AllStations_fn[AllStations_fn$coord_datum_cd == 'NAD83', ]
coordinates(GaugesLocs_NAD83) = c('dec_long_va', 'dec_lat_va')
proj4string(GaugesLocs_NAD83) = CRS('+init=epsg:4269')
#Transform to NAD83 UTM Zone 18N
GaugeLocs_NAD27 = spTransform(GaugesLocs_NAD27, CRS(pCRS))
GaugeLocs_NAD83 = spTransform(GaugesLocs_NAD83, CRS(pCRS))
#Join to one dataset again
GaugeLocs_fn = rbind(GaugeLocs_NAD27, GaugeLocs_NAD83)
#Remove separate datasets
rm(GaugeLocs_NAD27, GaugeLocs_NAD83, GaugesLocs_NAD27, GaugesLocs_NAD83)

# Method 2: Read USGS station data from .csv file----
AllStations <- read.csv(f_StreamGaugeData, stringsAsFactors = FALSE)

#Add a leading 0 to the NWIS gauges to look up their values on the server
# NOTE: This step may not be necessary for your dataset.
# NOTE: suppressing warnings for NAs introduced by coercion, which is intended.
#       Users should check that other warnings are not also being suppressed.
AllStations = addZeros(AllStations)

#  Select data for only NWIS gauges----
NWISstations = AllStations[AllStations$Source == 'NWIS',]

#  Add Gauge locations to that dataset---- 
#   Also contains altitudes of the gauges, which should be crosss-checked with DEM data
#   NOTE: you may have to change the commands to match your file.
if(exists("f_StreamGaugeSites")){
  GaugesLocs = read.table(f_StreamGaugeSites, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}else{
  GaugesLocs = readNWISsite(NWISstations$GaugeNum)
}

# Make spatial dataframe
GaugesLocs_NAD27 = GaugesLocs[GaugesLocs$coord_datum_cd == 'NAD27', ]
coordinates(GaugesLocs_NAD27) = c('dec_long_va', 'dec_lat_va')
proj4string(GaugesLocs_NAD27) = CRS('+init=epsg:4267')
GaugesLocs_NAD83 = GaugesLocs[GaugesLocs$coord_datum_cd == 'NAD83', ]
coordinates(GaugesLocs_NAD83) = c('dec_long_va', 'dec_lat_va')
proj4string(GaugesLocs_NAD83) = CRS('+init=epsg:4269')
#Transform to NAD83 UTM Zone 18N
GaugeLocs_NAD27 = spTransform(GaugesLocs_NAD27, CRS(pCRS))
GaugeLocs_NAD83 = spTransform(GaugesLocs_NAD83, CRS(pCRS))
#Join to one dataset again
GaugeLocs = rbind(GaugeLocs_NAD27, GaugeLocs_NAD83)
#Remove separate datasets
rm(GaugeLocs_NAD27, GaugeLocs_NAD83, GaugesLocs, GaugesLocs_NAD27, GaugesLocs_NAD83)

# Add DEM elevation to the gauge datasets----
#  Gather all of the DEM files together and mosaic into one file
DEM = processDEM(dir_DEMs = dir_DEM, f_DEMs = f_DEM, pCRS = pCRS)
#Add the DEM elevation to the gauge dataset
GaugeLocs$ElevDEM = extract(x = DEM, y = GaugeLocs)
GaugeLocs_fn$ElevDEM = extract(x = DEM, y = GaugeLocs_fn)

#  Compare the DEM elevation to the listed elevation----
png('CompareGaugeElev.png', res = 300, units = 'in', width = 5, height = 5)
plot(GaugeLocs$alt_va[which((is.na(GaugeLocs_fn$alt_va) == FALSE) & (is.na(GaugeLocs_fn$ElevDEM) == FALSE))], GaugeLocs$ElevDEM[which((is.na(GaugeLocs_fn$alt_va) == FALSE) & (is.na(GaugeLocs_fn$ElevDEM) == FALSE))]/.3048,
     xlab = 'USGS Reported Elevation (ft)', ylab = 'DEM Elevation (ft)', main = 'Gauge Elevations')
lines(c(-100,1100), c(-100,1100), col = 'red')
dev.off()

png('CompareGaugeElev_fn.png', res = 300, units = 'in', width = 5, height = 5)
plot(GaugeLocs_fn$alt_va[which((is.na(GaugeLocs_fn$alt_va) == FALSE) & (is.na(GaugeLocs_fn$ElevDEM) == FALSE))], GaugeLocs_fn$ElevDEM[which((is.na(GaugeLocs_fn$alt_va) == FALSE)  & (is.na(GaugeLocs_fn$ElevDEM) == FALSE))]/.3048,
     xlab = 'USGS Reported Elevation (ft)', ylab = 'DEM Elevation (ft)', main = 'Gauge Elevations')
lines(c(-100,1100), c(-100,1100), col = 'red')
dev.off()


#Fixme: Correct the elevations for those gauges that are very different than the DEM
# and have collection codes that suggest lower data quality
# Some of the gauges were likely reported in m instead of in ft in the USGS database
#Identify spatially those gauges that are more than X feet different
plot(GaugeLocs)
plot(GaugeLocs[which(abs(GaugeLocs$ElevDEM/.3048 - GaugeLocs$alt_va) >= 50),], col = 'red', add =T)
plot(GaugeLocs[which(abs(GaugeLocs$ElevDEM/.3048 - GaugeLocs$alt_va) >= 100),], col = 'blue', add =T)
plot(GaugeLocs[which(abs(GaugeLocs$ElevDEM/.3048 - GaugeLocs$alt_va) >= 200),], col = 'green', add =T)
plot(GaugeLocs_fn)
plot(GaugeLocs_fn[which(abs(GaugeLocs_fn$ElevDEM/.3048 - GaugeLocs_fn$alt_va) >= 50),], col = 'red', add =T)
plot(GaugeLocs_fn[which(abs(GaugeLocs_fn$ElevDEM/.3048 - GaugeLocs_fn$alt_va) >= 100),], col = 'blue', add =T)
plot(GaugeLocs_fn[which(abs(GaugeLocs_fn$ElevDEM/.3048 - GaugeLocs_fn$alt_va) >= 200),], col = 'green', add =T)


#  Method 2 only: Add coordinates, elevation, and other data to the NWISstations data----
for (i in 1:nrow(NWISstations)){
  #NOTE: if your data are not both character, errors stating the following will appear:
  # Error in data.frame(..., check.names = FALSE) : 
  #  arguments imply differing number of rows: 1, 0
  Ind = which(GaugeLocs$site_no == NWISstations$GaugeNum[i])
  if (length(Ind) > 1){
    print(paste('More than one gauge number matches the uniqueNum gauge ', i, '. Using only first match.'))
  }
  cmb = cbind(NWISstations[i,], GaugeLocs@data[Ind,], GaugeLocs@coords[Ind,][1], GaugeLocs@coords[Ind,][2])
  colnames(cmb) = c(colnames(NWISstations), colnames(GaugeLocs@data), colnames(GaugeLocs@coords))
  if (i == 1){ 
    NewData = cmb
  }else{
    NewData = rbind(NewData, cmb)
  }
}
NWISstations = NewData
rm(i, Ind, NewData, cmb)

# Clip to ROI----
ROI = readOGR(dsn = dir_ROI, layer = f_ROI, stringsAsFactors = FALSE)
ROI = spTransform(ROI, CRS(pCRS))

#Make NWIS stations a spatial dataframe
coordinates(NWISstations) = c('dec_long_va', 'dec_lat_va')
proj4string(NWISstations) = CRS(pCRS)
NWIS_ROI = NWISstations[ROI,]

NWIS_ROI_fn = GaugeLocs_fn[ROI,]

# Plot locations of gauges----
png('StremflowGauges.png', res = 300, units = 'in', width = 6, height = 6)
# All NWIS streamflow gauges in bounding box
plot(NWISstations, pch = 16, col = 'white')
# ROI
plot(ROI, add = TRUE)
# All NWIS streamflow gauges in bounding box, in color
plot(NWISstations, pch = 16, add = TRUE)
# NWIS streamflow gauges in ROI
plot(NWIS_ROI, pch = 16, col = 'red', add = TRUE)
# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 370000, yb = 4346000, len = 700, col = 'black', lab = 'N')
legend('topleft', title = 'Streamflow Stations', legend = c('In ROI', 'Not in ROI'), pch = 16, col = c('red', 'black'))
dev.off()

png('StremflowGauges_fn.png', res = 300, units = 'in', width = 6, height = 6)
# All NWIS streamflow gauges in bounding box
plot(NWISstations, pch = 16, col = 'white')
# ROI
plot(ROI, add = TRUE)
# All NWIS streamflow gauges in bounding box, in color
plot(GaugeLocs_fn, pch = 16, add = TRUE)
# NWIS streamflow gauges in ROI
plot(NWIS_ROI_fn, pch = 16, col = 'red', add = TRUE)
# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 370000, yb = 4346000, len = 700, col = 'black', lab = 'N')
legend('topleft', title = 'Streamflow Stations', legend = c('In ROI', 'Not in ROI'), pch = 16, col = c('red', 'black'))
dev.off()

# Plot reported vs. DEM elevation of gauges within the ROI----
png('CompareGaugeElev.png', res = 300, units = 'in', width = 5, height = 5)
plot(NWIS_ROI$alt_va, NWIS_ROI$ElevDEM/.3048,
     xlab = 'USGS Reported Elevation (ft)', ylab = 'DEM Elevation (ft)', main = 'Gauge Elevations')
lines(c(-100,1100), c(-100,1100), col = 'red')
dev.off()

png('CompareGaugeElev_fn.png', res = 300, units = 'in', width = 5, height = 5)
plot(NWIS_ROI_fn$alt_va, NWIS_ROI_fn$ElevDEM/.3048,
     xlab = 'USGS Reported Elevation (ft)', ylab = 'DEM Elevation (ft)', main = 'Gauge Elevations')
lines(c(-100,1100), c(-100,1100), col = 'red')
dev.off()

#One of these gauge elevations is a lot lower than DEM. Likely that the gauge was reported in m in USGS database
#identify(NWIS_ROI$alt_va, NWIS_ROI$ElevDEM/.3048)

# Statistics to download for each streamflow gauge, if available----
# All codes defined here: https://help.waterdata.usgs.gov/code/stat_cd_nm_query?stat_nm_cd=%25&fmt=html
#Fixme: lookup codes from the website
Stat.minFlow = "00002"
Stat.avgFlow = "00003"
Stat.maxFlow = "00001"
Stat.varFlow = "00010"
Stat.skewFlow = "00013"
#Xth percentile
Stat.P1 = "01010"
Stat.P5 = "01050"
Stat.P25 = "01250"
Stat.P50 = "01500"
Stat.P75 = "01750"
Stat.P95 = "01950"
Stat.P99 = "01990"

#Collect all of the Stat variables into a vector
#Fixme: is there a way to collect all variables that begin Stat. and collect them into a new vector?
Stats = c(Stat.minFlow, Stat.P1, Stat.P5, Stat.P25, Stat.P50, Stat.P75, Stat.P95, Stat.P99, Stat.maxFlow, Stat.avgFlow, Stat.varFlow, Stat.skewFlow)

# Parameters to download for each gauge, if available----
# All codes defined here: https://help.waterdata.usgs.gov/code/parameter_cd_query?fmt=rdb&inline=true&group_cd=%
#Fixme: lookup codes from the website
Par.cfsFlow = "00060"
Par.Nflow = "00600"
Par.Lat = "91110"
Par.Long = "91111"

#Fixme: is there a way to collect all variables that begin Par. and collect them into a new vector?
Pars = c(Par.cfsFlow, Par.Nflow, Par.Long, Par.Lat)

# Download the within-ROI stream gauge data in parallel----
#  Use only the unique gauge numbers in the dataset 
#  (repeats occur when multiple variables are available for a gauge)
uniqueNums = unique(NWIS_ROI_fn$site_no)
cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
a = foreach(i = uniqueNums, .packages = 'dataRetrieval') %dopar% {
  # Read all of the Pars data for the provided station number
  #Fixme: unable to test the use of Stats instead of only specifying statistic codes one by one.
  # Need a site that has more than one of these statistic codes reported to test.
  stationData <- readNWISdv(siteNumbers = i, parameterCd = Pars, statCd = Stats)
  write.table(stationData, 
              paste0(getwd(), '/StreamStat_', i, ".txt"), 
              sep = "\t")
}
stopCluster(cl)
rm(cl)
#Check that the run was successful
if (any(!is.null(unlist(a)))){
  print('DOWNLOAD UNSUCCESSFUL')
}else{
  print('Stream gauge data download complete!')
  rm(a)
}

#  Gather the records for each gauge into a list of dataframes----
StreamStationList = list()
#Also record their error codes for use in plotting later
ErrCodes = vector('character')
#Find all of the streamflow station file indices in directory
Ind_f_StreamStat = list.files()[grep(pattern = 'StreamStat', x = list.files(), ignore.case = FALSE, fixed = TRUE)]
for (i in 1:length(Ind_f_StreamStat)){
  #Read file
  f = read.table(Ind_f_StreamStat[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE, colClasses = c('character', 'character', 'character', 'Date', 'numeric', 'character'))
  #Gather the error codes for dates
  ErrCodes = unique(c(ErrCodes, f$X_00060_00003_cd))
  #Add to list
  StreamStationList = c(StreamStationList, list(f))
}
rm(f, Ind_f_StreamStat, i)

#Name the list by the gauge number
for (i in 1:length(StreamStationList)){
  names(StreamStationList)[i] = StreamStationList[[i]]$site_no[1]
}
rm(i)

# Plot the time series for each gauge, and the eCDF, colored by error code----
for (i in 1:length(StreamStationList)){
  #Assign colors to the error codes
  colCodes = rainbow(length(ErrCodes))
  
  #Find the indices of the codes that are in this dataset
  codes = which(ErrCodes %in% StreamStationList[[i]]$X_00060_00003_cd)
  
  png(paste0('StreamflowTimeseries_', StreamStationList[[i]]$site_no[1],'.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(x = StreamStationList[[i]]$Date, y = StreamStationList[[i]]$X_00060_00003, type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', ylab = 'Daily Mean Streamflow (cfs)', main = paste0('Station #', StreamStationList[[i]]$site_no[1]),
       ylim = c(0, 7000), xlim = c(as.Date("1950-01-01"), as.Date("2020-01-01")))
  #Add colors
  for (c in 1:length(colCodes)){
    par(new = TRUE)
    plot(x = StreamStationList[[i]]$Date[StreamStationList[[i]]$X_00060_00003_cd == ErrCodes[c]], y = StreamStationList[[i]]$X_00060_00003[StreamStationList[[i]]$X_00060_00003_cd == ErrCodes[c]], pch = 16, cex = 0.3,
         col = colCodes[c],
         ylim = c(0, 7000), xlim = c(as.Date("1950-01-01"), as.Date("2020-01-01")), axes = FALSE, xlab = '', ylab = '')
  }
  legend('topleft', legend = ErrCodes[codes], col = colCodes[codes], pch = 16)
  dev.off()
  
  #Fixme: get prettier streamflow exceedance graphs from a package
  png(paste0('StreamflowExceedance_', StreamStationList[[i]]$site_no[1],'.png'), res = 300, units = 'in', width = 6, height = 6)
  qqnorm(StreamStationList[[i]]$X_00060_00003, pch = 1, 
         ylab = 'Daily Mean Streamflow (cfs)', main = paste0('Non-Exceedance Probability for Daily Mean Streamflow \n Station #', StreamStationList[[i]]$site_no[1]))
  dev.off()
}
rm(i, c, colCodes, codes)

# Identify missing dates and fill them into the timeseries----
Fills = FillMissingDates(Dataset = NWIS_ROI_fn, StationList = StreamStationList, Var = 'X_00060_00003')
NWIS_ROI_fn = Fills$Dataset
StreamStationList = Fills$StationList
rm(Fills)

#Fixme: Missing data fill in with numerical value estimates using prediction in ungauged basins methods for large gaps
#Fixme: check for high and low flow outliers in each record, and compare spatially to other gauges on those dates
# Can include both FFA and daily flow outlier analysis

# Make a map of gauge locations colored by their record lengths, corrected for the total amount of missing data----
NWIS_ROI_fn$RecordLength = NWIS_ROI_fn$RecordLengthMinusGaps = NA
for (i in 1:length(StreamStationList)){
  NWIS_ROI_fn$RecordLength[which(as.numeric(NWIS_ROI_fn$site_no) == as.numeric(StreamStationList[[i]]$site_no[1]))] = as.numeric((max(StreamStationList[[i]]$Date) - min(StreamStationList[[i]]$Date)))
}
rm(i)
NWIS_ROI_fn$RecordLengthMinusGaps = NWIS_ROI_fn$RecordLength - NWIS_ROI_fn$MissingData
#in years
NWIS_ROI_fn$RecordLength = NWIS_ROI_fn$RecordLength/365.25
NWIS_ROI_fn$RecordLengthMinusGaps = NWIS_ROI_fn$RecordLengthMinusGaps/365.25

#Color by decades
scaleRange = c(0,70)
scaleBy = 10
Pal = rev(rainbow((scaleRange[2] - scaleRange[1])/scaleBy))
png('StremflowGauges_RecordLengths.png', res = 300, units = 'in', width = 6, height = 6)
plot(ROI)
# Gauges colored by their record lengths
plot(NWIS_ROI_fn, pch = 16, col = colFun(NWIS_ROI_fn$RecordLengthMinusGaps), add = TRUE)
# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 370000, yb = 4346000, len = 700, col = 'black', lab = 'N')
legend('right', title = 'Streamflow Station \n Record Lengths \n (years)', legend = seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy), pch = 16, col = colFun(seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy)), bty = 'n')
dev.off()

# Write streamflow datasets to files----
setwd(dir_DEM_out)
writeRaster(x = DEM, filename = f_DEM_mosiac, format = "GTiff")
setwd(dir_sfgauges)
writeOGR(obj = NWIS_ROI_fn, dsn = getwd(), driver = 'ESRI Shapefile', layer = f_NWIS_ROI_out)
writeOGR(obj = NWISstations, dsn = getwd(), driver = 'ESRI Shapefile', layer = f_NWIS_bb_out)
#Stream gauges had missing dates added to the files. Write new files.
for (i in names(StreamStationList)){
  write.table(StreamStationList[[i]], 
              paste0(getwd(), '/StreamStat_', i, f_sf_processKey, ".txt"), 
              sep = "\t")
}
rm(i)

#Water Quality----
setwd(dir_wq)
# Method 1: Read water quality station data using the USGS function----
#  Find sites that have any N and P water quality data in a state within the ROI
#Phosphorus
phosSites <- whatWQPsites(statecode="MD", characteristicName="Phosphorus")
#Nitrogen
NitroSites <- whatWQPsites(statecode="MD", characteristicName="Nitrogen")

#Make spatial data
#Nitrogen
GaugesLocs_NAD27 = NitroSites[NitroSites$HorizontalCoordinateReferenceSystemDatumName == 'NAD27', ]
coordinates(GaugesLocs_NAD27) = c('LongitudeMeasure', 'LatitudeMeasure')
proj4string(GaugesLocs_NAD27) = CRS('+init=epsg:4267')
GaugesLocs_NAD83 = NitroSites[NitroSites$HorizontalCoordinateReferenceSystemDatumName == 'NAD83', ]
coordinates(GaugesLocs_NAD83) = c('LongitudeMeasure', 'LatitudeMeasure')
proj4string(GaugesLocs_NAD83) = CRS('+init=epsg:4269')
GaugesLocs_WGS84 = NitroSites[NitroSites$HorizontalCoordinateReferenceSystemDatumName == 'WGS84', ]
coordinates(GaugesLocs_WGS84) = c('LongitudeMeasure', 'LatitudeMeasure')
proj4string(GaugesLocs_WGS84) = CRS('+init=epsg:4326')
#Assuming that all the unknown coordinate systems are NAD83
GaugesLocs_U = NitroSites[NitroSites$HorizontalCoordinateReferenceSystemDatumName == 'UNKWN', ]
coordinates(GaugesLocs_U) = c('LongitudeMeasure', 'LatitudeMeasure')
proj4string(GaugesLocs_U) = CRS('+init=epsg:4269')
#Transform to NAD83 UTM Zone 18N
GaugeLocs_NAD27 = spTransform(GaugesLocs_NAD27, CRS(pCRS))
GaugeLocs_NAD83 = spTransform(GaugesLocs_NAD83, CRS(pCRS))
GaugeLocs_WGS84 = spTransform(GaugesLocs_WGS84, CRS(pCRS))
GaugeLocs_U = spTransform(GaugesLocs_U, CRS(pCRS))
#Join to one dataset again
GaugeLocs_WQN = rbind(GaugeLocs_NAD27, GaugeLocs_NAD83, GaugeLocs_WGS84, GaugeLocs_U)
#Remove separate datasets
rm(GaugeLocs_NAD27, GaugeLocs_NAD83, GaugesLocs_NAD27, GaugesLocs_NAD83, GaugeLocs_WGS84, GaugeLocs_U, GaugesLocs_U, GaugesLocs_WGS84)

#Phosphorus
GaugesLocs_NAD27 = phosSites[phosSites$HorizontalCoordinateReferenceSystemDatumName == 'NAD27', ]
coordinates(GaugesLocs_NAD27) = c('LongitudeMeasure', 'LatitudeMeasure')
proj4string(GaugesLocs_NAD27) = CRS('+init=epsg:4267')
GaugesLocs_NAD83 = phosSites[phosSites$HorizontalCoordinateReferenceSystemDatumName == 'NAD83', ]
coordinates(GaugesLocs_NAD83) = c('LongitudeMeasure', 'LatitudeMeasure')
proj4string(GaugesLocs_NAD83) = CRS('+init=epsg:4269')
GaugesLocs_WGS84 = phosSites[phosSites$HorizontalCoordinateReferenceSystemDatumName == 'WGS84', ]
coordinates(GaugesLocs_WGS84) = c('LongitudeMeasure', 'LatitudeMeasure')
proj4string(GaugesLocs_WGS84) = CRS('+init=epsg:4326')
#Assuming that all the unknown coordinate systems are NAD83
GaugesLocs_U = phosSites[phosSites$HorizontalCoordinateReferenceSystemDatumName == 'UNKWN', ]
coordinates(GaugesLocs_U) = c('LongitudeMeasure', 'LatitudeMeasure')
proj4string(GaugesLocs_U) = CRS('+init=epsg:4269')
#Transform to NAD83 UTM Zone 18N
GaugeLocs_NAD27 = spTransform(GaugesLocs_NAD27, CRS(pCRS))
GaugeLocs_NAD83 = spTransform(GaugesLocs_NAD83, CRS(pCRS))
GaugeLocs_WGS84 = spTransform(GaugesLocs_WGS84, CRS(pCRS))
GaugeLocs_U = spTransform(GaugesLocs_U, CRS(pCRS))
#Join to one dataset again
GaugeLocs_WQP = rbind(GaugeLocs_NAD27, GaugeLocs_NAD83, GaugeLocs_WGS84, GaugeLocs_U)
#Remove separate datasets
rm(GaugeLocs_NAD27, GaugeLocs_NAD83, GaugesLocs_NAD27, GaugesLocs_NAD83, GaugeLocs_WGS84, GaugeLocs_U, GaugesLocs_U, GaugesLocs_WGS84)

#Clip to ROI
WQstations_ROI_N = GaugeLocs_WQN[ROI,]
WQstations_ROI_P = GaugeLocs_WQP[ROI,]

# Method 2: Read water quality station data from file----
WQstations = read.csv(f_WQgauges, stringsAsFactors = FALSE)
#Convert to spatial data
coordinates(WQstations) = c('LongitudeMeasure', 'LatitudeMeasure')
#Split dataset according to the coordinate reference systems used
GaugesLocs_NAD27 = WQstations[WQstations$HorizontalCoordinateReferenceSystemDatumName == 'NAD27', ]
proj4string(GaugesLocs_NAD27) = CRS('+init=epsg:4267')
GaugesLocs_NAD83 = WQstations[WQstations$HorizontalCoordinateReferenceSystemDatumName == 'NAD83', ]
proj4string(GaugesLocs_NAD83) = CRS('+init=epsg:4269')
GaugesLocs_WGS84 = WQstations[WQstations$HorizontalCoordinateReferenceSystemDatumName == 'WGS84', ]
proj4string(GaugesLocs_WGS84) = CRS('+init=epsg:4326')
#Transform to NAD83 UTM Zone 18N
GaugeLocs_NAD27 = spTransform(GaugesLocs_NAD27, CRS(pCRS))
GaugeLocs_NAD83 = spTransform(GaugesLocs_NAD83, CRS(pCRS))
GaugeLocs_WGS84 = spTransform(GaugesLocs_WGS84, CRS(pCRS))
#Join to one dataset again
WQGaugeLocs = rbind(GaugeLocs_NAD27, GaugeLocs_NAD83, GaugeLocs_WGS84)
#Remove separate datasets
rm(GaugeLocs_NAD27, GaugeLocs_NAD83, GaugeLocs_WGS84, GaugesLocs_NAD27, GaugesLocs_NAD83, GaugesLocs_WGS84)

#Clip to ROI
WQstations_ROI = WQGaugeLocs[ROI,]

#Select only those sites in the ROI that have nitrogen data
WQstations_ROI_N = WQstations_ROI[WQstations_ROI$MonitoringLocationIdentifier %in% NitroSites$MonitoringLocationIdentifier,]
#Select only those sites in the ROI that have Phosphorus data
WQstations_ROI_P = WQstations_ROI[WQstations_ROI$MonitoringLocationIdentifier %in% phosSites$MonitoringLocationIdentifier,]

# Plot TN and TP sampling locations on a map----
png('TNTPsites.png', res = 300, units = 'in', width = 6, height = 6)
plot(ROI)
plot(WQstations_ROI_N, pch = 16, add = TRUE, col = 'red')
plot(WQstations_ROI_P, pch = 16, add = TRUE, col = 'blue')
plot(WQstations_ROI_P[WQstations_ROI_N,], pch = 16, add = TRUE, col = 'purple')
plot(WQstations_ROI_N[WQstations_ROI_P,], pch = 16, add = TRUE, col = 'purple')
# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 365000, yb = 4346000, len = 700, col = 'black', lab = 'N')
legend('topright', title = 'Water Quality Sites', legend = c('T Nitrogen Only', 'T Phosphorus Only', 'Both'), col = c('red', 'blue', 'purple'), pch = 16)
dev.off()

# Download data for those sites in parallel----
#   Use only the unique gauge numbers in the dataset 
#   (repeats occur when multiple variables are available for a gauge)
uniqueWQNums_N = unique(WQstations_ROI_N$MonitoringLocationIdentifier)
uniqueWQNums_P = unique(WQstations_ROI_P$MonitoringLocationIdentifier)

#Make directories for storing Nitrogen and Phosphorous data
wd_N = paste0(getwd(), '/Nitrogen')
wd_P = paste0(getwd(), '/Phosphorus')
#Fixme: this should check if the directory already exists before making it
dir.create(path = wd_N)
dir.create(path = wd_P)

# NOTE: These downloads occasionally fail when run in parallel and return internal server errors.
# If that happens to you, try running in serial and see if you still get the errors (i.e. change dopar to do).
# I'm not sure what to do if you still get them. Running in serial has worked for me.
# NOTE: the readWQPqw function used here imports all water quality portal data.
#  A similar function called readNWISqw only downloads USGS data.
#  readWQPqw should be more generic, although this repo has not tested if those downloads contain all of the USGS data.
cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
n = foreach(i = uniqueWQNums_N, .packages = 'dataRetrieval') %dopar% {
  setwd(wd_N)
  # Read all of the data for the provided station number
  stationData <- readWQPqw(siteNumbers = i, parameterCd = "")
  write.table(stationData, 
              paste0(getwd(), '/Nitrogen_', i, ".txt"), 
              sep = "\t")
}
p = foreach(i = uniqueWQNums_P, .packages = 'dataRetrieval') %dopar% {
  setwd(wd_P)
  # Read all of the data for the provided station number
  stationData <- readWQPqw(siteNumbers = i, parameterCd = "")
  write.table(stationData, 
              paste0(getwd(), '/Phosphorus_', i, ".txt"), 
              sep = "\t")
}
stopCluster(cl)
rm(cl)
#Check that the run was successful
if (any(!is.null(unlist(n)))){
  print('NITROGEN DOWNLOAD UNSUCCESSFUL')
}else{
  print('Nitrogen gauge data download complete!')
  rm(n)
}
if (any(!is.null(unlist(p)))){
  print('Phosphorus DOWNLOAD UNSUCCESSFUL')
}else{
  print('Phosphorus gauge data download complete!')
  rm(p)
}

# Process Nitrogen Data----
setwd(wd_N)
#Gather the records for each gauge into a list of dataframes
NitroStationList = list()
#Find all of the Nitrogen station file indices in directory
Ind_f_NitroStat = list.files()[grep(pattern = 'Nitrogen_', x = list.files(), ignore.case = FALSE, fixed = TRUE)]
for (i in 1:length(Ind_f_NitroStat)){
  #Read file
  f = read.table(Ind_f_NitroStat[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  #Add to list
  NitroStationList = c(NitroStationList, list(f))
}
rm(f, Ind_f_NitroStat, i)

#Cycle through all of the unique combinations of ResultSampleFractionText and CharacteristicName and return separate text files for each variable
for (i in 1:length(NitroStationList)){
  #Find all of the unique CharacteristicName, ResultSampleFractionText combinations
  us = unique(NitroStationList[[i]][,c("CharacteristicName", "ResultSampleFractionText")])
  #Write separate files for each combination
  for (u in 1:nrow(us)){
    #Check for NA values in us
    if(is.na(us[u,1])){
      u_Ind1 = which(is.na(NitroStationList[[i]][, c("CharacteristicName", "ResultSampleFractionText")[1]]))
    }else{
      u_Ind1 = which(NitroStationList[[i]][, c("CharacteristicName", "ResultSampleFractionText")[1]] == us[u,1])
    }
    if(is.na(us[u,2])){
      u_Ind2 = which(is.na(NitroStationList[[i]][, c("CharacteristicName", "ResultSampleFractionText")[2]]))
    }else{
      u_Ind2 = which(NitroStationList[[i]][, c("CharacteristicName", "ResultSampleFractionText")[2]] == us[u,2])
    }
    #Gather only the indices for this unique combination
    u_Ind = u_Ind1[u_Ind1 %in% u_Ind2]
    
    #Check that ActivityMediaName = Water for stream-only nitrogen data
    w_Ind = which(NitroStationList[[i]][u_Ind,"ActivityMediaName"] != "Water")
    if (length(w_Ind) > 0){
      print(paste('Some of the samples are taken in media ', unique(NitroStationList[[i]]$ActivityMediaName), ' for station ', NitroStationList[[i]]$MonitoringLocationIdentifier[1]))
    }
    #Note: subbing ~ in for / because files will not save with / in the name.
    write.table(NitroStationList[[i]][u_Ind,], 
                paste0(getwd(), '/Nitrogen_', NitroStationList[[i]]$MonitoringLocationIdentifier[1], "_cn", gsub(pattern = "/", x = us[u,1], replacement = "~", fixed = TRUE), '_rt', us[u,2], ".txt"), 
                sep = "\t")
  }
}
rm(us, u, i, u_Ind, u_Ind1, u_Ind2)

#Read in the CharacteristicName for nitrogen measurements only
cn_Nitro = grep(x = list.files(), pattern = '_cnNitrogen', ignore.case = TRUE)
# For now, taking only ResultSampleFractionText = total nitrogen
cn_Nitro1 = grep(x = list.files()[cn_Nitro], pattern = 'Total', ignore.case = TRUE)
#Total nitrogen files only
f_TN = list.files()[cn_Nitro][cn_Nitro1]
rm(cn_Nitro, cn_Nitro1)
#Make a list of all of the total nitrogen gauge datasets
TN = list()
for (i in 1:length(f_TN)){
  f = read.table(f_TN[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  #Place the measurements in chronological order
  f = f[order(f$ActivityStartDate),]
  
  #Check that the units are all the same for each measurement
  us = unique(f$ResultMeasure.MeasureUnitCode)
  if (length(us) > 1){
    print(paste('Warning: more than one measurement unit for station ', f$MonitoringLocationIdentifier[1], '. Fix manually.'))
  }
  
  #Check for detection limits
  dl = unique(f$DetectionQuantitationLimitMeasure.MeasureValue)
  #remove NAs
  dl = dl[!is.na(dl)]
  if(length(dl) != 0){
    print(paste('Warning: detection limits for station', f$MonitoringLocationIdentifier[1], '. Fix manually.'))
  }
  
  #add data to list
  TN = c(TN, 
         list(f))
}
rm(i, f, dl, us)

#  Plot TN timeseries----
for (i in 1:length(TN)){
  png(paste0('TN_Timeseries_', TN[[i]]$MonitoringLocationIdentifier[1],'.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = TN[[i]]$ResultMeasureValue, x = as.Date(TN[[i]]$ActivityStartDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = paste0(TN[[i]]$ResultSampleFractionText[1], " ", TN[[i]]$CharacteristicName[1], " (", TN[[i]]$ResultMeasure.MeasureUnitCode[1], ")"), 
       main = paste0('Station #', TN[[i]]$MonitoringLocationIdentifier[1]),
       ylim = c(0, 10), xlim = c(as.Date("1980-01-01"), as.Date("2010-01-01")))
  dev.off()
}
rm(i)

#Fixme: handle missing data, as with streamflow 
#        this may be unnecessary because regular sampling is not completed for WQ data

#Fixme: Search for flow-normalized outliers

# Process Phosphorus Data----
setwd(wd_P)
#Gather the records for each gauge into a list of dataframes
PhosStationList = list()
#Find all of the phosphorus station file indices in directory
Ind_f_PhosStat = list.files()[grep(pattern = 'Phosphorus_', x = list.files(), ignore.case = FALSE, fixed = TRUE)]
for (i in 1:length(Ind_f_PhosStat)){
  #Read file
  f = read.table(Ind_f_PhosStat[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  #Add to list
  PhosStationList = c(PhosStationList, list(f))
}
rm(f, Ind_f_PhosStat, i)

#Cycle through all of the unique combinations of ResultSampleFractionText and CharacteristicName and return separate text files for each variable
for (i in 1:length(PhosStationList)){
  #Find all of the unique CharacteristicName, ResultSampleFractionText combinations
  us = unique(PhosStationList[[i]][,c("CharacteristicName", "ResultSampleFractionText")])
  #Write separate files for each combination
  for (u in 1:nrow(us)){
    #Check for NA values in us
    if(is.na(us[u,1])){
      u_Ind1 = which(is.na(PhosStationList[[i]][, c("CharacteristicName", "ResultSampleFractionText")[1]]))
    }else{
      u_Ind1 = which(PhosStationList[[i]][, c("CharacteristicName", "ResultSampleFractionText")[1]] == us[u,1])
    }
    if(is.na(us[u,2])){
      u_Ind2 = which(is.na(PhosStationList[[i]][, c("CharacteristicName", "ResultSampleFractionText")[2]]))
    }else{
      u_Ind2 = which(PhosStationList[[i]][, c("CharacteristicName", "ResultSampleFractionText")[2]] == us[u,2])
    }
    #Gather only the indices for this unique combination
    u_Ind = u_Ind1[u_Ind1 %in% u_Ind2]
    
    #Check that ActivityMediaName = Water for stream-only Phosphorus data
    w_Ind = which(PhosStationList[[i]][u_Ind,"ActivityMediaName"] != "Water")
    if (length(w_Ind) > 0){
      print(paste('Some of the samples are taken in media ', unique(PhosStationList[[i]]$ActivityMediaName), ' for station ', PhosStationList[[i]]$MonitoringLocationIdentifier[1]))
    }
    #Note: subbing ~ in for / because files will not save with / in the name.
    write.table(PhosStationList[[i]][u_Ind,], 
                paste0(getwd(), '/Phosphorus_', PhosStationList[[i]]$MonitoringLocationIdentifier[1], "_cn", gsub(pattern = "/", x = us[u,1], replacement = "~", fixed = TRUE), '_rt', us[u,2], ".txt"), 
                sep = "\t")
  }
}
rm(us, u, i, u_Ind, u_Ind1, u_Ind2, w_Ind)

#Read in the CharacteristicName for Phosphorus measurements only
cn_Phos = grep(x = list.files(), pattern = '_cnPhosphorus', ignore.case = TRUE)
# For now, taking only ResultSampleFractionText = total Phosphorous
cn_Phos1 = grep(x = list.files()[cn_Phos], pattern = 'Total', ignore.case = TRUE)
#Total nitrogen files only
f_TP = list.files()[cn_Phos][cn_Phos1]
rm(cn_Phos, cn_Phos1)
#Make a list of all of the total nitrogen gauge datasets
TP = list()
for (i in 1:length(f_TP)){
  f = read.table(f_TP[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  #Place the measurements in chronological order
  f = f[order(f$ActivityStartDate),]
  
  #Check that the units are all the same for each measurement
  us = unique(f$ResultMeasure.MeasureUnitCode)
  if (length(us) > 1){
    print(paste('Warning: more than one measurement unit for station ', f$MonitoringLocationIdentifier[1], '. Fix manually.'))
  }
  
  #Check for detection limits
  dl = unique(f$DetectionQuantitationLimitMeasure.MeasureValue)
  #remove NAs
  dl = dl[!is.na(dl)]
  if(length(dl) != 0){
    print(paste('Warning: detection limits for station ', f$MonitoringLocationIdentifier[1], '. Fix manually.'))
  }
  
  #add data to list
  TP = c(TP, 
         list(f))
}
rm(i, f, dl, us)

#  Plot TP timeseries----
for (i in 1:length(TP)){
  png(paste0('TP_Timeseries_', TP[[i]]$MonitoringLocationIdentifier[1],'.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = TP[[i]]$ResultMeasureValue, x = as.Date(TP[[i]]$ActivityStartDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = paste0(TP[[i]]$ResultSampleFractionText[1], " ", TP[[i]]$CharacteristicName[1], " (", TP[[i]]$ResultMeasure.MeasureUnitCode[1], ")"), 
       main = paste0('Station #', TP[[i]]$MonitoringLocationIdentifier[1]),
       ylim = c(0, 0.5), xlim = c(as.Date("1980-01-01"), as.Date("2010-01-01")))
  
  #Check for and add detection limits
  dl = unique(TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue)
  #remove NAs
  dl = dl[!is.na(dl)]
  if(length(dl) != 0){
    #Add detection limits
    #Lower detection limit
    par(new = TRUE)
    plot(y = TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue[grep(pattern = "Low", x = TP[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)], x = as.Date(TP[[i]]$ActivityStartDate[grep(pattern = "Low", x = TP[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)]), pch = 16, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = c(0, 0.5), xlim = c(as.Date("1980-01-01"), as.Date("2010-01-01")), axes = FALSE,
         col = 'blue')
    #Upper detection limit
    par(new = TRUE)
    plot(y = TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue[grep(pattern = "Up", x = TP[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)], x = as.Date(TP[[i]]$ActivityStartDate[grep(pattern = "Up", x = TP[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)]), pch = 16, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = c(0, 0.5), xlim = c(as.Date("1980-01-01"), as.Date("2010-01-01")), axes = FALSE,
         col = 'red')
    
    legend('topright', title = 'Detection Limits', legend = c('Upper', 'Lower'), col = c('red', 'blue'), pch = 16)
  }
  dev.off()
}
rm(i)

#Fixme: handle missing data, as with streamflow 
#        this may be unnecessary because regular sampling is not completed for WQ data

#Fixme: Search for flow-normalized outliers


#Fixme: handle missing data, as with streamflow 
#        this may be unnecessary because regular sampling is not completed for WQ data

#Fixme: Search for flow-normalized outliers

# Write water quality data to files----
writeOGR(WQstations_ROI_N, dsn = getwd(), layer = 'NitrogenSites', driver = "ESRI Shapefile")
list.save(x = TN, file = 'TN.yaml', type = "YAML")
writeOGR(WQstations_ROI_P, dsn = getwd(), layer = 'PhosphorusSites', driver = "ESRI Shapefile")
list.save(x = TP, file = 'TP.yaml', type = "YAML")


# BES Water Quality Gauge Data----
setwd("C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\WaterChemistry")
#BES Water quality sample time series
BES_WQ = read.csv('BES-stream-chemistry-data-for-WWW-feb-2018---core-sites-only_JDSprocessed.csv', stringsAsFactors = FALSE)
#BES gauge number reference table
USGS_GaugeMatch = read.csv('Abbreviations_SampleRecordLengths.csv', stringsAsFactors = FALSE)
#Remove spaces from the names of the sites
USGS_GaugeMatch$Abbreviation = gsub(pattern = ' ', replacement = '', x = USGS_GaugeMatch$Abbreviation, fixed = TRUE)
#add leading zeros to gauge numbers that are not NA
USGS_GaugeMatch$USGSGaugeNum[!is.na(USGS_GaugeMatch$USGSGaugeNum)] = paste0("0", USGS_GaugeMatch$USGSGaugeNum[!is.na(USGS_GaugeMatch$USGSGaugeNum)])
BES_WQ$USGSgauge = NA
USGSnums = unique(USGS_GaugeMatch$USGSGaugeNum)[-which(is.na(unique(USGS_GaugeMatch$USGSGaugeNum)) == TRUE)]
#Capitalize all of the site names in the BES and GaugeMatch datasets because there are typos
USGS_GaugeMatch$Abbreviation = toupper(USGS_GaugeMatch$Abbreviation)
BES_WQ$Site = toupper(BES_WQ$Site)
  
#Filter into separate databases by site and add gauge numbers to database
BES_WQ_Sites = list()
uniqueSites = unique(BES_WQ$Site)
for (i in 1:length(uniqueSites)){
  f = BES_WQ[BES_WQ$Site == uniqueSites[i],]
  #Assign gauge number if it has one
  if(f$Site[1] %in% USGS_GaugeMatch$Abbreviation){
    f$USGSgauge = USGS_GaugeMatch$USGSGaugeNum[USGS_GaugeMatch$Abbreviation == f$Site[1]]
  }else{
    f$USGSgauge = NA
  }
  
  f$Dated = NA
  
  #Fix the date format
  #4-digit year - any number greater than current year on computer is going to be assumed 19XX
  curr2DigYear = substring(as.character(as.numeric(strsplit(date(), split = " ", fixed = TRUE)[[1]][5])), first = 3, last = 4)
  for (y in 1:length(f$Date)){
    txtyr = strsplit(f$Date[y], split = '-', fixed = TRUE)[[1]]
    txtyr[3] = ifelse(as.numeric(txtyr[3]) > as.numeric(curr2DigYear), str_c('19', txtyr[3]), str_c('20', txtyr[3]))
    txtyr = as.Date(str_c(txtyr[1], '/', txtyr[2], '/', txtyr[3]), format = '%d/%B/%Y')
    f$Dated[y] = as.character(txtyr)
  }
  BES_WQ_Sites = c(BES_WQ_Sites, list(f))
  names(BES_WQ_Sites) = c(names(BES_WQ_Sites)[-length(names(BES_WQ_Sites))], f$USGSgauge[1])
}
rm(i, f, txtyr)

names(BES_WQ_Sites)

#Plot time series for each of the sites
#Nitrogen
for (i in 1:length(BES_WQ_Sites)){
  if (any(!is.na(BES_WQ_Sites[[i]]$TN..mg.N.L.))){
    png(paste0('BES_N_Timeseries_', BES_WQ_Sites[[i]]$Site[1], '_',  BES_WQ_Sites[[i]]$USGSgauge[1], '_', i, '.png'), res = 300, units = 'in', width = 6, height = 6)
    plot(y = BES_WQ_Sites[[i]]$TN..mg.N.L., x = as.Date(BES_WQ_Sites[[i]]$Dated), type = 'o', pch = 16, cex = 0.3,
         xlab = 'Year', 
         ylab = 'Total Nitrogen (mg N/L)', 
         main = paste0('TN Station #', BES_WQ_Sites[[i]]$USGSgauge[1], ' ', BES_WQ_Sites[[i]]$Site[1]),
         ylim = c(0, max(BES_WQ_Sites[[i]]$TN..mg.N.L., na.rm=TRUE)))
    #xlim = c(as.Date("1980-01-01"), as.Date("2019-06-01")))
    dev.off()
  }
}
rm(i)

#Phosphorus
for (i in 1:length(BES_WQ_Sites)){
  if (any(!is.na(BES_WQ_Sites[[i]]$TP..ugP.L.))){
    png(paste0('BES_P_Timeseries_', BES_WQ_Sites[[i]]$Site[1], '_',  BES_WQ_Sites[[i]]$USGSgauge[1], '_', i, '.png'), res = 300, units = 'in', width = 6, height = 6)
    plot(y = BES_WQ_Sites[[i]]$TP..ugP.L., x = as.Date(BES_WQ_Sites[[i]]$Dated), type = 'o', pch = 16, cex = 0.3,
         xlab = 'Year', 
         ylab = expression(paste('Total Phosphorus (', mu, 'g P/L)')), 
         main = paste0('TP Station #', BES_WQ_Sites[[i]]$USGSgauge[1], ' ', BES_WQ_Sites[[i]]$Site[1]),
         ylim = c(0, max(BES_WQ_Sites[[i]]$TP..ugP.L., na.rm=TRUE)))
    #xlim = c(as.Date("1980-01-01"), as.Date("2019-06-01")))
    dev.off()
  }
}
rm(i)

#Coordinates of BES sites
for (i in 1:length(BES_WQ_Sites)){
  if (!is.na(BES_WQ_Sites[[i]]$USGSgauge[1])){
    if(!exists('BES_WQ_Sites_locs')){
      BES_WQ_Sites_locs = GaugeLocs[GaugeLocs$site_no == BES_WQ_Sites[[i]]$USGSgauge[1],]
    }else{
      BES_WQ_Sites_locs = rbind(BES_WQ_Sites_locs, GaugeLocs[GaugeLocs$site_no == BES_WQ_Sites[[i]]$USGSgauge[1],])
    }
  }
}
rm(i)

#Map the station locations and plot them on a map showing BES data and WQP data
png('TNTP_BES+WQPsites.png', res = 300, units = 'in', width = 6, height = 6)
plot(ROI)
plot(WQstations_ROI_N, pch = 16, add = TRUE, col = 'red')
plot(WQstations_ROI_P, pch = 16, add = TRUE, col = 'blue')
plot(WQstations_ROI_P[WQstations_ROI_N,], pch = 16, add = TRUE, col = 'purple')
plot(WQstations_ROI_N[WQstations_ROI_P,], pch = 16, add = TRUE, col = 'purple')
#BES sites
plot(BES_WQ_Sites_locs, add = TRUE)

# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 365000, yb = 4346000, len = 700, col = 'black', lab = 'N')
legend('topright', title = 'Water Quality Sites', legend = c('T Nitrogen Only', 'T Phosphorus Only', 'Both'), col = c('red', 'blue', 'purple'), pch = 16)
dev.off()

#Plot the streamflow, N, and P data together for each site


#Try plotting the water quality data using R tools----


#Weather Stations----
#Fixme: some AllStations data are climate stations
#       add DEM elevation to dataset and compare
#Load file containing hyperlinks to the climate data
ClimGauges = read.csv("NOAA_HyperlinksToGauges.csv", stringsAsFactors = FALSE)

#Read the hyperlinks and place data into separate text files per gauge
for (i in 1:nrow(ClimGauges)){
  ClimStationData = read.table(ClimGauges[i])
  #Station names are either after NWIS= or AWIC=
  if (length(grep(pattern = 'NWIS=', ClimGauges[i])) == 0){
    StationName = strsplit(x = strsplit(x = ClimGauges[i], split = "NWIS=", fixed = TRUE)[2], split = "%", fixed = TRUE)[1]
  }else{
    StationName = strsplit(x = strsplit(x = ClimGauges[i], split = "AWIC=", fixed = TRUE)[2], split = "%", fixed = TRUE)[1]
  }
  #Fixme: the filename should be the gauge name. Need to find that in the datasets or URL
  write.table(ClimStationData,
              paste0(getwd(), '/Clim_', StationName, ".txt"), 
              sep = "\t")
}
