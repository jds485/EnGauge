#Script to read USGS streamflow and water quality, and NOAA weather station data from gauges

#Author Credits:
# Portions of the function for the streamflow download were provided by:
# Caitline Barber (Caitline.Barber@tufts.edu) and Jonathan Lamontagne (Jonathan.Lamontagne@tufts.edu)
# Modified by Jared Smith (js4yd@virginia.edu) in June, 2019, and started git tracking.
# See git commit history for all later contributions

#Fixme: catch all print statements for each section and write to text file(s)
#Fixme: make each method a separate function or script
#Fixme: add methods for handling detection limits in time series

#Set directory names----
#Region of interest shapefile directory
dir_ROI = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\EnGauge\\EnGauge\\DataForExamples"  
#ColorFunctions.R script directory - from JDS github repo: Geothermal_ESDA
dir_ColFuns = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\USGSGauges"
#EnGauge code repository directory
dir_EnGauge = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\EnGauge\\EnGauge"
#Directory where USGS streamflow gauge data will be downloaded. This directory must already exist.
dir_sfgauges = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\EnGauge\\EnGauge\\DataForExamples"
#Directory where water quality gauge data will be downloaded. This directory must already exist.
dir_wq = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\EnGauge\\EnGauge\\DataForExamples"
#Directory where weather station data will be downloaded. This directory must already exist.
dir_weather = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\EnGauge\\EnGauge\\DataForExamples"

#DEM - specify as a vector of directories if there are multiple tiles to be mosaicked together.
# The directory order has to match the file name order for f_DEM below.
# Fixme: can DEMs be downloaded from a server instead of downloading manually before using this script?
dir_DEM = c("C:\\Users\\js4yd\\OneDrive - University of Virginia\\EnGauge\\EnGauge\\DataForExamples\\grdn40w077_1",
            "C:\\Users\\js4yd\\OneDrive - University of Virginia\\EnGauge\\EnGauge\\DataForExamples\\grdn40w078_1")
#Output directory for processed DEM. This directory must already exist.
dir_DEM_out = 'C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\DEM\\'

#Set input filenames----
#Region of interest shapefile name
f_ROI = "Watershed_GF"

#Streamflow gauges and site coordinates filenames - Only used for Method 2 for Streamflow
f_StreamGaugeData = "BES_USGS_GaugeStations.csv"
#Water quality gauges - Only used for Method 2 for Water Quality
f_WQgauges = "BES_WaterQualityGaugeStations.csv"

#DEM - all separate DEM tiles should be added to this vector (e.g. c("w001001.adf", "w001002.adf") )
f_DEM = c("w001001.adf", "w001001.adf")

#Set output filenames----
#NWIS streamflow gauges in ROI
f_NWIS_ROI_out = 'NWIS_ROI'
#NWIS streamflow gauges in bounding box. Only for Method 2 streamflow.
f_NWIS_bb_out = 'NWIS_bb'
#Streamflow gauge data processing name appendage (_p for processed)
# e.g. 0159384_p
f_sf_processKey = '_p'
#Name for the list of all streamflow gauge timeseries. YAML list type is used below.
f_StreamStationList = 'SF.yaml'

#Name for the list of all TN water quality site timeseries. YAML list type is used below.
f_TNSiteList = 'TN.yaml'
f_TNSiteList_daily = 'TN_d.yaml'
f_TNSiteList_monthly = 'TN_m.yaml'
f_TNSiteList_annual = 'TN_a.yaml'
#Name for the list of all TP water quality site timeseries. YAML list type is used below.
f_TPSiteList = 'TP.yaml'
f_TPSiteList_daily = 'TP_d.yaml'
f_TPSiteList_monthly = 'TP_m.yaml'
f_TPSiteList_annual = 'TP_a.yaml'

#NOAA weather station data locations
f_NOAAstationsROI = 'NOAA_StationLocs'
#Name for the list of all NOAA weather station timeseries. YAML list type is used below.
f_NOAAstationsDataList = 'NOAA_MetStations.yaml'

#DEM - GeoTiff format is the default. You can change in the script.
f_DEM_mosiac = "DEM_mosaic"

#Set project coordinate system----
#This is the coordinate system that all data will be plotted and written in
# It is not the coordinate system of your data (although it could be)
# EPSG codes from: https://spatialreference.org/ref/?page=2
pCRS = '+init=epsg:26918'

#Set plot limits - set to NULL to ignore use----
#Streamflow dates
xlim_dates = c(as.Date("1950-01-01"), as.Date("2020-01-01"))
#Water quality dates
xlim_WQdates = c(as.Date("1980-01-01"), as.Date("2010-01-01"))
#Streamflow y axis limits
ylim_SF = c(0, 7000)
#TN y-axis limits
ylim_TN = c(0, 10)
#TP y-axis limits
ylim_TP = c(0, 0.5)
#Temperature y-axis limits
ylim_Temp = c(-20,50)

#Define a small buffer to use for the ROI, in pCRS map units----
#For streamflow and water quality
ROIbuff = 0
#For weather stations
ROIbuffWeather = 20000

#Set high outlier quantile----
HighOutProb = 0.99

#Load libraries and functions----
#USGS function library - note that a more recent version is available through Github
library(dataRetrieval)
#USGS EGRET library - currently not used
#library(EGRET)
#R libraries
library(stringi)
library(stringr)
library(foreach)
library(parallel)
library(doParallel)
library(sp)
library(rgdal)
library(maptools)
library(GISTools)
library(raster)
library(rlist)
library(lattice)
library(zoo)
library(lubridate)
library(hydroTSM)
library(rnoaa)
library(rappdirs)
#Color functions for plots (R script from Jared Smith's Geothermal_ESDA Github Repo)
setwd(dir_ColFuns)
source('ColorFunctions.R')
#Functions from repository
setwd(dir_EnGauge)
source('fdcDefaultModification.R')
source('missingDates.R')
source('addZerosToGaugeNames.R')
source('processDEM.R')
source('extractWQdata.R')
source('checkDuplicates.R')
source('checkZerosNegs.R')
source('formatMonthlyMatrix.R')
source('matplotDates.R')
source('aggregateTimeseries.R')

#Make Region of Interest (ROI) buffer----
ROI = readOGR(dsn = dir_ROI, layer = f_ROI, stringsAsFactors = FALSE)
ROI = spTransform(ROI, CRS(pCRS))

#buffer
ROI_buff = buffer(ROI, ROIbuff)
ROI_bufferW = buffer(ROI, width = ROIbuffWeather)

#Streamflow----
setwd(dir_sfgauges)
#Create directory to save files
wd_sf = paste0(getwd(), '/Streamflow')
dir.create(path = wd_sf, showWarnings = FALSE)

# Method 1: Get gauges using whatNWISsites()----
AllStations_fn = whatNWISsites(statecode = "MD", parameterCd = '00060')
AllStations_fn = readNWISsite(AllStations_fn$site_no)
#  Make a spatial dataframe: Method 1----
# NOTE: Your data may be all one coordinate system, and therefore not need to split into 2 datasets
#       before joining into 1 dataset.
# NOTE: your coordinate system may be different (epsg code)
# Some of the data are NAD27 projection and others are NAD83 projection. Split the dataset to handle each
#Fixme: function for splitting coordinate systems and returning one same-coordinate system file
#       Would require being able to look up epsg codes.

#vector of unique coordinate systems in station data
StreamUniqueCoords = unique(AllStations_fn$coord_datum_cd)

#Process to the pCRS coordinate system
GaugesLocs_NAD27 = AllStations_fn[which(AllStations_fn$coord_datum_cd == 'NAD27'), ]
coordinates(GaugesLocs_NAD27) = c('dec_long_va', 'dec_lat_va')
proj4string(GaugesLocs_NAD27) = CRS('+init=epsg:4267')
GaugesLocs_NAD83 = AllStations_fn[which(AllStations_fn$coord_datum_cd == 'NAD83'), ]
coordinates(GaugesLocs_NAD83) = c('dec_long_va', 'dec_lat_va')
proj4string(GaugesLocs_NAD83) = CRS('+init=epsg:4269')
#Transform to pCRS
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
NWISstations = AllStations[which(AllStations$Source == 'NWIS'),]

#  Add Gauge locations to that dataset---- 
#   Also contains altitudes of the gauges, which should be crosss-checked with DEM data
#   NOTE: you may have to change the commands to match your file.
GaugesLocs = readNWISsite(NWISstations$GaugeNum)

#  Make spatial dataframe: Method 2----
GaugesLocs_NAD27 = GaugesLocs[which(GaugesLocs$coord_datum_cd == 'NAD27'), ]
coordinates(GaugesLocs_NAD27) = c('dec_long_va', 'dec_lat_va')
proj4string(GaugesLocs_NAD27) = CRS('+init=epsg:4267')
GaugesLocs_NAD83 = GaugesLocs[which(GaugesLocs$coord_datum_cd == 'NAD83'), ]
coordinates(GaugesLocs_NAD83) = c('dec_long_va', 'dec_lat_va')
proj4string(GaugesLocs_NAD83) = CRS('+init=epsg:4269')
#Transform to NAD83 UTM Zone 18N
GaugeLocs_NAD27 = spTransform(GaugesLocs_NAD27, CRS(pCRS))
GaugeLocs_NAD83 = spTransform(GaugesLocs_NAD83, CRS(pCRS))
#Join to one dataset again
GaugeLocs = rbind(GaugeLocs_NAD27, GaugeLocs_NAD83)
#Remove separate datasets
rm(GaugeLocs_NAD27, GaugeLocs_NAD83, GaugesLocs, GaugesLocs_NAD27, GaugesLocs_NAD83)

#  Add coordinates and other data to the NWISstations data:Methodd 2 only----
if (exists(x = "GaugeLocs")){
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
}

# Clip to ROI----
if (exists(x = "GaugeLocs_fn")){
  NWIS_ROI_fn = GaugeLocs_fn[ROI_buff,]
}
if (exists(x = "GaugeLocs")){
  #  Method 2 only: Make NWIS stations a spatial dataframe
  coordinates(NWISstations) = c('dec_long_va', 'dec_lat_va')
  proj4string(NWISstations) = CRS(pCRS)
  NWIS_ROI = NWISstations[ROI_buff,]
}

# Plot locations of gauges----
setwd(wd_sf)
if (exists(x = "GaugeLocs_fn")){
  #  Method 1 only
  png('StremflowGauges_fn.png', res = 300, units = 'in', width = 6, height = 6)
  # ROI
  plot(ROI)
  # All NWIS streamflow gauges in bounding box, in color
  plot(GaugeLocs_fn, pch = 16, add = TRUE)
  # NWIS streamflow gauges in ROI
  plot(NWIS_ROI_fn, pch = 16, col = 'red', add = TRUE)
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
  legend('bottomleft', title = 'Streamflow Stations', legend = c('In ROI', 'Not in ROI'), pch = 16, col = c('red', 'black'))
  dev.off()
}
if (exists(x = "GaugeLocs")){
  #  Method 2 only
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
}

# Statistics to download for each streamflow gauge, if available----
# All codes defined here: https://help.waterdata.usgs.gov/code/stat_cd_nm_query?stat_nm_cd=%25&fmt=html
#Fixme: lookup codes directly from the website
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
Stats = c(Stat.minFlow, Stat.P1, Stat.P5, Stat.P25, Stat.P50, Stat.P75, Stat.P95, Stat.P99, 
          Stat.maxFlow, Stat.avgFlow, Stat.varFlow, Stat.skewFlow)

# Parameters to download for each gauge, if available----
# All codes defined here: https://help.waterdata.usgs.gov/code/parameter_cd_query?fmt=rdb&inline=true&group_cd=%
#Fixme: lookup codes directly from the website
Par.cfsFlow = "00060"
Par.Nflow = "00600"
Par.Lat = "91110"
Par.Long = "91111"

#Fixme: is there a way to collect all variables that begin Par. and collect them into a new vector?
Pars = c(Par.cfsFlow, Par.Nflow, Par.Long, Par.Lat)

# Download the within-ROI stream gauge data in parallel----
#  Use only the unique gauge numbers in the dataset 
#  (repeats numbers occur when multiple variables are available for a gauge)
setwd(wd_sf)
#Get unique numbers
if (exists(x = "NWIS_ROI_fn")){
  #  Method 1
  uniqueNums = unique(NWIS_ROI_fn$site_no)
}else if (exists(x = "NWIS_ROI")){
  #  Method 2
  uniqueNums = unique(NWIS_ROI$site_no)
}

#Download
cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
a = foreach(i = uniqueNums, .packages = 'dataRetrieval') %dopar% {
  # Read all of the Pars data for the provided station number
  #Fixme: unable to test the use of Stats instead of only specifying statistic codes one by one.
  # Need a site that has more than one of these statistic codes reported to test.
  stationData <- readNWISdv(siteNumbers = i, parameterCd = Pars, statCd = Stats)
  write.table(stationData, 
              paste0(getwd(), '/StreamStat_', i, ".txt"), 
              sep = "\t", row.names = FALSE)
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
  f = read.table(Ind_f_StreamStat[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE, 
                 colClasses = c('character', 'character', 'Date', 'numeric', 'character'))
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

#  Check for duplicate observations----
# This step is recommended in USGS publications.
StreamStationList = checkDuplicatesAndRemove(StreamStationList)

#  Check for zeros and negative values in the records----
StreamStationList = checkZerosNegs(StreamStationList)

#   See if any stations have negative or zero values and add to NWIS spatial dataset----
if (exists(x = "NWIS_ROI_fn")){
  #    Method 1
  NWIS_ROI_fn = addNegsToSpatialDataset(StationList = StreamStationList, SpatialDataset = NWIS_ROI_fn, site_D = 'site_no', site_SL = 'site_no')
  NWIS_ROI_fn = addZerosToSpatialDataset(StationList = StreamStationList, SpatialDataset = NWIS_ROI_fn, site_D = 'site_no', site_SL = 'site_no')
}else if (exists(x = "NWIS_ROI")){
  #    Method 2
  NWIS_ROI = addNegsToSpatialDataset(StationList = StreamStationList, SpatialDataset = NWIS_ROI, site_D = 'site_no', site_SL = 'site_no')
  NWIS_ROI = addZerosToSpatialDataset(StationList = StreamStationList, SpatialDataset = NWIS_ROI, site_D = 'site_no', site_SL = 'site_no')
}

#     Plot a map of the stations with zero and negative records----
if (exists(x = "NWIS_ROI_fn")){
  #     Method 1
  png(paste0('Streamflow_ZerosNegsMap_fn.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(ROI)
  #All gauges
  plot(NWIS_ROI_fn, pch = 16, add = TRUE)
  #Gauges with zeros
  plot(NWIS_ROI_fn[!is.na(NWIS_ROI_fn$Zero),], pch = 16, col = 'red', add = TRUE)
  #Gauges with negative values
  plot(NWIS_ROI_fn[!is.na(NWIS_ROI_fn$Neg),], pch = 16, col = 'blue', add = TRUE)
  #Gauges with both
  plot(NWIS_ROI_fn[which(!is.na(NWIS_ROI_fn$Neg) & !is.na(NWIS_ROI_fn$Zero)),], pch = 16, col = 'purple', add = TRUE)
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
  legend('topright', title = 'Streamflow Stations', legend = c('Zeros in Record', 'Negatives in Record', 'Both', 'Neither'), pch = 16, col = c('red', 'blue', 'purple', 'black'), bty = 'n')
  dev.off()
}else if (exists(x = "NWIS_ROI")){
  #     Method 2
  png(paste0('Streamflow_ZerosNegsMap.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(ROI)
  #All gauges
  plot(NWIS_ROI, pch = 16, add = TRUE)
  #Gauges with zeros
  plot(NWIS_ROI[!is.na(NWIS_ROI$Zero),], pch = 16, col = 'red', add = TRUE)
  #Gauges with negative values
  plot(NWIS_ROI[!is.na(NWIS_ROI$Neg),], pch = 16, col = 'blue', add = TRUE)
  #Gauges with both
  plot(NWIS_ROI[which(!is.na(NWIS_ROI$Neg) & !is.na(NWIS_ROI$Zero)),], pch = 16, col = 'purple', add = TRUE)
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
  legend('topright', title = 'Streamflow Stations', legend = c('Zeros in Record', 'Negatives in Record', 'Both', 'Neither'), pch = 16, col = c('red', 'blue', 'purple', 'black'), bty = 'n')
  dev.off()
}

# Identify missing dates and fill them into the timeseries. This also places the timeseries in chronological order----
#Fixme: Missing data fill in with numerical value estimates using prediction in ungauged basins methods for large gaps
if (exists(x = "NWIS_ROI_fn")){
  #  Method 1
  Fills = FillMissingDates_par(Dataset = NWIS_ROI_fn, StationList = StreamStationList, Var = 'X_00060_00003', 
                           Date = 'Date', gapType = 'd', site_no_D = 'site_no', site_no_SL = 'site_no', NoNAcols = 'agency_cd', 
                           NumCores = detectCores()-1)
  NWIS_ROI_fn = Fills$Dataset
  StreamStationList = Fills$StationList
  rm(Fills)
}else if (exists(x = "NWIS_ROI")){
  #  Method 2
  Fills = FillMissingDates_par(Dataset = NWIS_ROI, StationList = StreamStationList, Var = 'X_00060_00003', 
                               Date = 'Date', gapType = 'd', site_no_D = 'site_no', site_no_SL = 'site_no', NoNAcols = 'agency_cd', 
                               NumCores = detectCores()-1)
  NWIS_ROI = Fills$Dataset
  StreamStationList = Fills$StationList
  rm(Fills)
}

# Plot the time series for each gauge, and the eCDF, colored by error code----
#Fixme: add confidence intervals for flow duration curves (fdcu package)
#Fixme: make axes have the same limits on the seasonal plots
for (i in 1:length(StreamStationList)){
  #Timeseries:----
  #Assign colors to the error codes
  colCodes = rainbow(length(ErrCodes))
  
  #Find the indices of the codes that are in this dataset
  codes = which(ErrCodes %in% StreamStationList[[i]]$X_00060_00003_cd)
  
  png(paste0('StreamflowTimeseries_', StreamStationList[[i]]$site_no[1],'.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(x = StreamStationList[[i]]$Date, y = StreamStationList[[i]]$X_00060_00003, type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', ylab = 'Daily Mean Streamflow (cfs)', main = paste0('Station #', StreamStationList[[i]]$site_no[1]),
       ylim = ylim_SF, xlim = xlim_dates)
  #Add colors
  for (cc in 1:length(colCodes)){
    par(new = TRUE)
    plot(x = StreamStationList[[i]]$Date[which(StreamStationList[[i]]$X_00060_00003_cd == ErrCodes[cc])], 
         y = StreamStationList[[i]]$X_00060_00003[which(StreamStationList[[i]]$X_00060_00003_cd == ErrCodes[cc])], 
         pch = 16, cex = 0.3,
         col = colCodes[cc],
         ylim = ylim_SF, xlim = xlim_dates, axes = FALSE, xlab = '', ylab = '')
  }
  legend('topleft', legend = ErrCodes[codes], col = colCodes[codes], pch = 16)
  dev.off()
  
  #Flow-Duration Curve----
  #Make y axis limits for flow duration curve
  at.y <- outer(1:9, 10^(-3:5))
  lab.y <- ifelse(log10(at.y) %% 1 == 0, sapply(log10(at.y),function(k)
    as.expression(bquote(10^ .(k)))), NA)
  
  png(paste0('StreamflowExceedance_', StreamStationList[[i]]$site_no[1],'.png'), res = 300, units = 'in', width = 6, height = 6)
  fdc(x = StreamStationList[[i]]$X_00060_00003, pch = NA, ylab = 'Daily Mean Streamflow (cfs)', 
      main = paste0('Exceedance Probability for Daily Mean Streamflow \n Station #', StreamStationList[[i]]$site_no[1]),
      xlim = c(0,1), lQ.thr = NA, hQ.thr = NA, thr.shw = FALSE, yat = -1000, verbose = FALSE)
  axis(side=2, at=at.y, labels=lab.y, cex.axis=1.2, las=1, tck = -0.01, hadj = 0.7, padj = 0.3)
  axis(side=2, at=10^seq(-3,5,1), labels=FALSE, tck = -0.025)
  dev.off()
  
  #Streamflow exceedance, timeseries, and map----
  png(paste0('StreamflowExceedanceTimeseries_Map_', StreamStationList[[i]]$site_no[1],'.png'), 
      res = 300, units = 'in', width = 10, height = 10)
  layout(rbind(c(1,2), c(3, 2)))
  
  #Exceedance
  fdc(x = StreamStationList[[i]]$X_00060_00003, pch = NA, ylab = 'Daily Mean Streamflow (cfs)', 
      main = paste0('Exceedance Probability for Daily Mean Streamflow \n Station #', StreamStationList[[i]]$site_no[1]),
      xlim = c(0,1), lQ.thr = NA, hQ.thr = NA, thr.shw = FALSE, yat = -1000, verbose = FALSE)
  axis(side=2, at=at.y, labels=lab.y, cex.axis=1.2, las=1, tck = -0.01, hadj = 0.7, padj = 0.3)
  axis(side=2, at=10^seq(-3,5,1), labels=FALSE, tck = -0.025)
  
  #Map of sites with the plotted streamflow exceedance site highlighted in red
  plot(ROI)
  #All gauges
  if (exists(x = "NWIS_ROI_fn")){
    plot(NWIS_ROI_fn, pch = 16, add = TRUE)
    #Gauge selected for exceedance plot
    plot(NWIS_ROI_fn[which(NWIS_ROI_fn$site_no == StreamStationList[[i]]$site_no[1]),], pch = 16, col = 'red', add = TRUE)
    title(main = NWIS_ROI_fn$station_nm[which(NWIS_ROI_fn$site_no == StreamStationList[[i]]$site_no[1])])
  }else if (exists(x = "NWIS_ROI")){
    plot(NWIS_ROI, pch = 16, add = TRUE)
    #Gauge selected for exceedance plot
    plot(NWIS_ROI[which(NWIS_ROI$site_no == StreamStationList[[i]]$site_no[1]),], pch = 16, col = 'red', add = TRUE)
    title(main = NWIS_ROI$station_nm[which(NWIS_ROI$site_no == StreamStationList[[i]]$site_no[1])])
  }
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
  legend('topright', title = 'Streamflow Stations', legend = c('Selected', 'Other'), pch = 16, col = c('red', 'black'), bty = 'n')
  
  #Timeseries
  plot(x = StreamStationList[[i]]$Date, y = StreamStationList[[i]]$X_00060_00003, type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', ylab = 'Daily Mean Streamflow (cfs)', main = paste0('Station #', StreamStationList[[i]]$site_no[1]),
       ylim = c(min(StreamStationList[[i]]$X_00060_00003, na.rm = TRUE), max(StreamStationList[[i]]$X_00060_00003, na.rm = TRUE)), 
       xlim = xlim_dates)
  #Add colors
  for (cc in 1:length(colCodes)){
    par(new = TRUE)
    plot(x = StreamStationList[[i]]$Date[which(StreamStationList[[i]]$X_00060_00003_cd == ErrCodes[cc])], 
         y = StreamStationList[[i]]$X_00060_00003[which(StreamStationList[[i]]$X_00060_00003_cd == ErrCodes[cc])], 
         pch = 16, cex = 0.3, col = colCodes[cc], axes = FALSE, xlab = '', ylab = '',
         ylim = c(min(StreamStationList[[i]]$X_00060_00003, na.rm = TRUE), max(StreamStationList[[i]]$X_00060_00003, na.rm = TRUE)), 
         xlim = xlim_dates)
  }
  legend('topleft', legend = ErrCodes[codes], col = colCodes[codes], pch = 16)
  dev.off()
  
  #hydroTSM plots----
  #daily timeseries
  dts = zoo(StreamStationList[[i]]$X_00060_00003, order.by = StreamStationList[[i]]$Date)
  
  #monthly timeseries for mean daily flows in a month
  mts = daily2monthly(dts, FUN = mean, na.rm = TRUE)
  #plot(mts)
  #monthly total flows boxplot
  #boxplot(coredata(mts) ~ factor(format(time(mts), "%b"), levels=unique(format(time(mts), "%b")), ordered=TRUE))
  #Monthly matrix
  M = formatMonthlyMatrix(mts)
  # Plotting the monthly values as matrixplot 
  png(paste0('StreamflowMonthly_', StreamStationList[[i]]$site_no[1],'.png'), res = 300, units = 'in', width = 6, height = 3)
  print(matrixplot(M, ColorRamp="Precipitation", main="Mean Monthly Streamflow"))
  dev.off()

  #annual timeseries
  #ats = daily2annual(dts, FUN = mean, na.rm = TRUE)
  #plot(ats)
  #annual total flows boxplot
  #boxplot(coredata(ats) ~ factor(format(time(ats), "%b"), levels=unique(format(time(ats), "%b")), ordered=TRUE))
  
  #Summary EDA plots
  png(paste0('StreamflowEDA_', StreamStationList[[i]]$site_no[1],'.png'), res = 300, units = 'in', width = 10, height = 10)
  hydroplot(dts, FUN = mean, col = 'black', var.unit = 'mean cfs')
  dev.off()
  
  #Seasonal flows
  if ((as.numeric((as.Date(max(StreamStationList[[i]]$Date)) - as.Date(min(StreamStationList[[i]]$Date)))) >= 366) 
      & (length(StreamStationList[[i]]$X_00060_00003) >= 366)){
    png(paste0('StreamflowEDA_Seasonal_', StreamStationList[[i]]$site_no[1],'.png'), res = 300, units = 'in', width = 10, height = 10)
    hydroplot(dts, FUN = mean, col = 'black', var.unit = 'mean cfs', pfreq = 'seasonal', ylab = 'Mean Streamflow (cfs)')
    dev.off()
  }
  
  #Fixme: annual timeseries trace plots (overlay of annual hydrographs / shading of 10-25-75-90, plot median sf for water year)
  
  #Fixme: quantitative summary tables like: seasonalfunction(dts, FUN = mean)
  # hydropairs for correlation plot
}
rm(i, cc, colCodes, codes, mts, dts, M, at.y, lab.y)

# Outlier detection----
#Fixme: check for high and low flow outliers in each record, and compare spatially to other gauges on those dates
#Fixme: add temporal outlier detection methods to account for correlation at many lag distances (spectral methods?)
#Fixme: add spatiotemporal outlier detection methods
#Fixme: flags for ST outliers at specified levels

for (i in 1:length(StreamStationList)){
  #Hacky spatiotemporal outlier detection method
  #select all values in the 99th percentile or above
  highs = StreamStationList[[i]][which(StreamStationList[[i]]$X_00060_00003 > quantile(StreamStationList[[i]]$X_00060_00003, probs = HighOutProb, na.rm = TRUE)),c('X_00060_00003', 'Date')]
  #empirical quantile for all of the flows
  qhighs = highs
  qhighs$X_00060_00003 = ecdf(StreamStationList[[i]]$X_00060_00003)(highs$X_00060_00003)
  #Name the streamflow column the station number
  colnames(highs) = c(StreamStationList[[i]]$site_no[1], 'Date')
  colnames(qhighs) = c(StreamStationList[[i]]$site_no[1], 'Date')
  #Flip to have Date be the first column. Makes more sense for reporting and plotting.
  highs = highs[,c(2,1)]
  qhighs = qhighs[,c(2,1)]
  
  #Gather dataframes for the day before and day after the outlier record
  bhighs = ahighs = highs
  bqhighs = aqhighs = qhighs
  #Assign the date
  bhighs$Date = highs$Date - 1
  ahighs$Date = highs$Date + 1
  bqhighs$Date = qhighs$Date - 1
  aqhighs$Date = qhighs$Date + 1
  #Assign NA values to the timeseries
  bhighs[,2] = NA
  ahighs[,2] = NA
  bqhighs[,2] = NA
  aqhighs[,2] = NA
  bInds = which(StreamStationList[[i]]$X_00060_00003 > quantile(StreamStationList[[i]]$X_00060_00003, probs = HighOutProb, na.rm = TRUE))-1
  aInds = which(StreamStationList[[i]]$X_00060_00003 > quantile(StreamStationList[[i]]$X_00060_00003, probs = HighOutProb, na.rm = TRUE))+1
  if (min(bInds) != 0){
    bhighs[,2] = StreamStationList[[i]][bInds,'X_00060_00003']
    bqhighs[,2] = ecdf(StreamStationList[[i]]$X_00060_00003)(bhighs[,2])
  }else{
    zInd = which(bInds == 0)
    #Set to 1 for now, but later will set back to NA
    bInds[zInd] = 1
    bhighs[,2] = StreamStationList[[i]][bInds,'X_00060_00003']
    bqhighs[,2] = ecdf(StreamStationList[[i]]$X_00060_00003)(bhighs[,2])
    bhighs[zInd,2] = NA
    bqhighs[zInd,2] = NA
    rm(zInd)
  }
  if (max(aInds) != (nrow(StreamStationList[[i]])+1)){
    ahighs[,2] = StreamStationList[[i]][aInds,'X_00060_00003']
    aqhighs[,2] = ecdf(StreamStationList[[i]]$X_00060_00003)(ahighs[,2])
  }else{
    zInd = which(aInds == (nrow(StreamStationList[[i]])+1))
    #Set to 1 for now, but later will set back to NA
    aInds[zInd] = 1
    ahighs[,2] = StreamStationList[[i]][aInds,'X_00060_00003']
    aqhighs[,2] = ecdf(StreamStationList[[i]]$X_00060_00003)(ahighs[,2])
    ahighs[zInd,2] = NA
    aqhighs[zInd,2] = NA
    rm(zInd)
  }
  for (j in 1:length(StreamStationList)){
    if (j != i){
      #Search for the dates corresponding to the high flows in other stations
      Inds = which(StreamStationList[[j]]$Date %in% highs$Date)
      if(length(Inds) > 0){
        highsj = StreamStationList[[j]][Inds, c('X_00060_00003', 'Date')]
        #Add new column to the original dataset for this station's flows
        highs[,StreamStationList[[j]]$site_no[1]] = NA
        qhighs[,StreamStationList[[j]]$site_no[1]] = NA
        #Detect the first date that matches
        IndStart = which(highs$Date == min(highsj$Date))
        #Ensure that the dates match in order from the start to the last date that is in the series
        # These should match because timeseries were filled in for missing dates above.
        if(any(highsj$Date[which(highsj$Date %in% highs$Date)] == highs$Date[IndStart:(nrow(highsj)+IndStart-1)]) == FALSE){
          stop('Error: Timeseries dates do not match. Make sure all dates are filled into all stations (no gaps in coverage). NAs are okay.')
        }

        #record the quantile values for that station
        # This could be a problem if the outlier is on the first or last day of the timeseries. That would only affect the first and last values in the timeseries for any one station.
        #empirical quantile for all of the flows
        qhighs[IndStart:(nrow(highsj)+IndStart-1),StreamStationList[[j]]$site_no[1]] = ecdf(StreamStationList[[j]]$X_00060_00003)(highsj$X_00060_00003)
        
        #record the raw values for that station
        highs[IndStart:(nrow(highsj)+IndStart-1),StreamStationList[[j]]$site_no[1]] = highsj$X_00060_00003
        
      }else{
        #Report NAs for that station's values. It doesn't have any dates in common with the station.
        highs[,StreamStationList[[j]]$site_no[1]] = NA
        qhighs[,StreamStationList[[j]]$site_no[1]] = NA
      }
      
      #day before
      #Search for the dates corresponding to the high flows in other stations
      Inds = which(StreamStationList[[j]]$Date %in% bhighs$Date)
      if(length(Inds) > 0){
        highsj = StreamStationList[[j]][Inds, c('X_00060_00003', 'Date')]
        #Add new column to the original dataset for this station's flows
        bhighs[,StreamStationList[[j]]$site_no[1]] = NA
        bqhighs[,StreamStationList[[j]]$site_no[1]] = NA
        #Detect the first date that matches
        IndStart = which(bhighs$Date == min(highsj$Date))
        #Ensure that the dates match in order from the start to the last date that is in the series
        # These should match because timeseries were filled in for missing dates above.
        if(any(highsj$Date[which(highsj$Date %in% bhighs$Date)] == bhighs$Date[IndStart:(nrow(highsj)+IndStart-1)]) == FALSE){
          stop('Error: Timeseries dates do not match. Make sure all dates are filled into all stations (no gaps in coverage). NAs are okay.')
        }
        
        #record the quantile values for that station
        # This could be a problem if the outlier is on the first or last day of the timeseries. That would only affect the first and last values in the timeseries for any one station.
        #empirical quantile for all of the flows
        bqhighs[IndStart:(nrow(highsj)+IndStart-1),StreamStationList[[j]]$site_no[1]] = ecdf(StreamStationList[[j]]$X_00060_00003)(highsj$X_00060_00003)
        
        #record the raw values for that station
        bhighs[IndStart:(nrow(highsj)+IndStart-1),StreamStationList[[j]]$site_no[1]] = highsj$X_00060_00003
        
      }else{
        #Report NAs for that station's values. It doesn't have any dates in common with the station.
        bhighs[,StreamStationList[[j]]$site_no[1]] = NA
        bqhighs[,StreamStationList[[j]]$site_no[1]] = NA
      }
      
      #day after
      #Search for the dates corresponding to the high flows in other stations
      Inds = which(StreamStationList[[j]]$Date %in% ahighs$Date)
      if(length(Inds) > 0){
        highsj = StreamStationList[[j]][Inds, c('X_00060_00003', 'Date')]
        #Add new column to the original dataset for this station's flows
        ahighs[,StreamStationList[[j]]$site_no[1]] = NA
        aqhighs[,StreamStationList[[j]]$site_no[1]] = NA
        #Detect the first date that matches
        IndStart = which(ahighs$Date == min(highsj$Date))
        #Ensure that the dates match in order from the start to the last date that is in the series
        # These should match because timeseries were filled in for missing dates above.
        if(any(highsj$Date[which(highsj$Date %in% ahighs$Date)] == ahighs$Date[IndStart:(nrow(highsj)+IndStart-1)]) == FALSE){
          stop('Error: Timeseries dates do not match. Make sure all dates are filled into all stations (no gaps in coverage). NAs are okay.')
        }
        
        #record the quantile values for that station
        # This could be a problem if the outlier is on the first or last day of the timeseries. That would only affect the first and last values in the timeseries for any one station.
        #empirical quantile for all of the flows
        aqhighs[IndStart:(nrow(highsj)+IndStart-1),StreamStationList[[j]]$site_no[1]] = ecdf(StreamStationList[[j]]$X_00060_00003)(highsj$X_00060_00003)
        
        #record the raw values for that station
        ahighs[IndStart:(nrow(highsj)+IndStart-1),StreamStationList[[j]]$site_no[1]] = highsj$X_00060_00003
        
      }else{
        #Report NAs for that station's values. It doesn't have any dates in common with the station.
        ahighs[,StreamStationList[[j]]$site_no[1]] = NA
        aqhighs[,StreamStationList[[j]]$site_no[1]] = NA
      }
    }
  }
  
  #Save the highs and qhighs data to files
  write.table(highs, 
              paste0(getwd(), '/OutlierCheck_Flows_', colnames(highs)[2], ".txt"), 
              sep = "\t", row.names = FALSE)
  write.table(qhighs, 
              paste0(getwd(), '/OutlierCheck_Quantiles_', colnames(qhighs)[2], ".txt"), 
              sep = "\t", row.names = FALSE)
  write.table(ahighs, 
              paste0(getwd(), '/OutlierCheck_aFlows_', colnames(ahighs)[2], ".txt"), 
              sep = "\t", row.names = FALSE)
  write.table(aqhighs, 
              paste0(getwd(), '/OutlierCheck_aQuantiles_', colnames(aqhighs)[2], ".txt"), 
              sep = "\t", row.names = FALSE)
  write.table(bhighs, 
              paste0(getwd(), '/OutlierCheck_bFlows_', colnames(bhighs)[2], ".txt"), 
              sep = "\t", row.names = FALSE)
  write.table(bqhighs, 
              paste0(getwd(), '/OutlierCheck_bQuantiles_', colnames(bqhighs)[2], ".txt"), 
              sep = "\t", row.names = FALSE)
  
  #Make a lineplot for the quantiles of the outliers across each of the stations in the ROI
  #Order the qhighs and highs dataframes by their column name numbers to ensure that the colors used for each station are the same in each iteration
  qhighs = qhighs[,c(1,1+order(as.numeric(colnames(qhighs[-1]))))]
  highs = highs[,c(1,1+order(as.numeric(colnames(highs[-1]))))]
  bqhighs = bqhighs[,c(1,1+order(as.numeric(colnames(bqhighs[-1]))))]
  bhighs = bhighs[,c(1,1+order(as.numeric(colnames(bhighs[-1]))))]
  aqhighs = aqhighs[,c(1,1+order(as.numeric(colnames(aqhighs[-1]))))]
  ahighs = ahighs[,c(1,1+order(as.numeric(colnames(ahighs[-1]))))]
  
  #also order the stations for plotting with same color scheme
  if (exists(x = "NWIS_ROI_fn")){
    pstations = NWIS_ROI_fn[order(as.numeric(NWIS_ROI_fn$site_no)),]
    #figure titles
    title1 = paste(NWIS_ROI_fn$station_nm[which(NWIS_ROI_fn$site_no == StreamStationList[[i]]$site_no[1])], '\n 99th Percentile Daily Average Flow')
    title2 = paste(NWIS_ROI_fn$station_nm[which(NWIS_ROI_fn$site_no == StreamStationList[[i]]$site_no[1])], '\n Day Before 99th Percentile Daily Average Flow')
    title3 = paste(NWIS_ROI_fn$station_nm[which(NWIS_ROI_fn$site_no == StreamStationList[[i]]$site_no[1])], '\n Day After 99th Percentile Daily Average Flow')
  }else if (exists(x = "NWIS_ROI")){
    pstations = NWIS_ROI[order(as.numeric(NWIS_ROI$site_no)),]
    #figure titles
    title1 = paste(NWIS_ROI$station_nm[which(NWIS_ROI$site_no == StreamStationList[[i]]$site_no[1])], '\n 99th Percentile Daily Average Flow')
    title2 = paste(NWIS_ROI$station_nm[which(NWIS_ROI$site_no == StreamStationList[[i]]$site_no[1])], '\n Day Before 99th Percentile Daily Average Flow')
    title3 = paste(NWIS_ROI$station_nm[which(NWIS_ROI$site_no == StreamStationList[[i]]$site_no[1])], '\n Day After 99th Percentile Daily Average Flow')
  }
  
  png(paste0('Outlier99Quantile_', StreamStationList[[i]]$site_no[1],'.png'), res = 300, units = 'in', width = 12, height = 6)
  layout(rbind(c(1,2)))
  #Plot of outlier timeseries
  cols = rainbow(n = (ncol(qhighs)-1))
  cols[which(colnames(qhighs[-1]) == StreamStationList[[i]]$site_no[1])] = 'black'
  matplotDates(x = qhighs$Date, y = qhighs[,-1], type = 'o', pch = 16, cex = 0.7, col = cols,
               xlab = 'Date', ylab = 'Quantile')
  #Map of sites with the plotted streamflow exceedance site highlighted in red
  plot(ROI)
  #All gauges
  plot(pstations, pch = 16, add = TRUE, col = rainbow(n = (ncol(qhighs)-1)))
  #Gauge selected for 99th percentile flows
  plot(pstations[which(pstations$site_no == StreamStationList[[i]]$site_no[1]),], pch = 16, add = TRUE, col = 'black')
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
  legend('topright', title = 'Streamflow Stations', legend = c('Reference'), pch = 16, col = c('black'), bty = 'n')
  title(main = title1)
  dev.off()
  
  png(paste0('Outlier99bQuantile_', StreamStationList[[i]]$site_no[1],'.png'), res = 300, units = 'in', width = 12, height = 6)
  layout(rbind(c(1,2)))
  #Plot of outlier timeseries
  cols = rainbow(n = (ncol(bqhighs)-1))
  cols[which(colnames(bqhighs[-1]) == StreamStationList[[i]]$site_no[1])] = 'black'
  matplotDates(x = bqhighs$Date, y = bqhighs[,-1], type = 'o', pch = 16, cex = 0.7, col = cols,
               xlab = 'Date', ylab = 'Quantile')
  #Map of sites with the plotted streamflow exceedance site highlighted in red
  plot(ROI)
  #All gauges
  plot(pstations, pch = 16, add = TRUE, col = rainbow(n = (ncol(qhighs)-1)))
  #Gauge selected for 99th percentile flows
  plot(pstations[which(pstations$site_no == StreamStationList[[i]]$site_no[1]),], pch = 16, add = TRUE, col = 'black')
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
  legend('topright', title = 'Streamflow Stations', legend = c('Reference'), pch = 16, col = c('black'), bty = 'n')
  title(main = title2)
  dev.off()
  
  png(paste0('Outlier99aQuantile_', StreamStationList[[i]]$site_no[1],'.png'), res = 300, units = 'in', width = 12, height = 6)
  layout(rbind(c(1,2)))
  #Plot of outlier timeseries
  cols = rainbow(n = (ncol(aqhighs)-1))
  cols[which(colnames(aqhighs[-1]) == StreamStationList[[i]]$site_no[1])] = 'black'
  matplotDates(x = aqhighs$Date, y = aqhighs[,-1], type = 'o', pch = 16, cex = 0.7, col = cols,
               xlab = 'Date', ylab = 'Quantile')
  #Map of sites with the plotted streamflow exceedance site highlighted in red
  plot(ROI)
  #All gauges
  plot(pstations, pch = 16, add = TRUE, col = rainbow(n = (ncol(qhighs)-1)))
  #Gauge selected for 99th percentile flows
  plot(pstations[which(pstations$site_no == StreamStationList[[i]]$site_no[1]),], pch = 16, add = TRUE, col = 'black')
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
  legend('topright', title = 'Streamflow Stations', legend = c('Reference'), pch = 16, col = c('black'), bty = 'n')
  title(main = title3)
  dev.off()
}
rm(i, j, Inds, IndStart, highs, highsj, qhighs, pstations, ahighs, bhighs, aInds, bInds, cols, aqhighs, bqhighs)

# Make a map of gauge locations colored by their record lengths, corrected for the total amount of missing data----
if (exists(x = "NWIS_ROI_fn")){
  NWIS_ROI_fn$RecordLength = NWIS_ROI_fn$RecordLengthMinusGaps = NA
  for (i in 1:length(StreamStationList)){
    NWIS_ROI_fn$RecordLength[which(as.numeric(NWIS_ROI_fn$site_no) == as.numeric(StreamStationList[[i]]$site_no[1]))] = as.numeric((max(StreamStationList[[i]]$Date) - min(StreamStationList[[i]]$Date)))
  }
  rm(i)
  NWIS_ROI_fn$RecordLengthMinusGaps = NWIS_ROI_fn$RecordLength - NWIS_ROI_fn$MissingData_d
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
  legend('right', title = 'Streamflow Station \n Record Lengths \n (years)', legend = seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy), 
         pch = 16, col = colFun(seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy)), bty = 'n')
  dev.off()
}else if (exists(x = "NWIS_ROI")){
  NWIS_ROI$RecordLength = NWIS_ROI$RecordLengthMinusGaps = NA
  for (i in 1:length(StreamStationList)){
    NWIS_ROI$RecordLength[which(as.numeric(NWIS_ROI$site_no) == as.numeric(StreamStationList[[i]]$site_no[1]))] = as.numeric((max(StreamStationList[[i]]$Date) - min(StreamStationList[[i]]$Date)))
  }
  rm(i)
  NWIS_ROI$RecordLengthMinusGaps = NWIS_ROI$RecordLength - NWIS_ROI$MissingData_d
  #in years
  NWIS_ROI$RecordLength = NWIS_ROI$RecordLength/365.25
  NWIS_ROI$RecordLengthMinusGaps = NWIS_ROI$RecordLengthMinusGaps/365.25
  
  #Color by decades
  scaleRange = c(0,70)
  scaleBy = 10
  Pal = rev(rainbow((scaleRange[2] - scaleRange[1])/scaleBy))
  png('StremflowGauges_RecordLengths.png', res = 300, units = 'in', width = 6, height = 6)
  plot(ROI)
  # Gauges colored by their record lengths
  plot(NWIS_ROI, pch = 16, col = colFun(NWIS_ROI$RecordLengthMinusGaps), add = TRUE)
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4346000, len = 700, col = 'black', lab = 'N')
  legend('right', title = 'Streamflow Station \n Record Lengths \n (years)', legend = seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy), 
         pch = 16, col = colFun(seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy)), bty = 'n')
  dev.off()
}

# Write streamflow datasets to files----
setwd(wd_sf)
if (exists(x = "NWIS_ROI_fn")){
  writeOGR(obj = NWIS_ROI_fn, dsn = getwd(), driver = 'ESRI Shapefile', layer = f_NWIS_ROI_out)
}else if (exists(x = "NWIS_ROI")){
  writeOGR(obj = NWIS_ROI, dsn = getwd(), driver = 'ESRI Shapefile', layer = f_NWIS_ROI_out)
  writeOGR(obj = NWISstations, dsn = getwd(), driver = 'ESRI Shapefile', layer = f_NWIS_bb_out)
}
#Stream gauges had missing dates added to the files. Write new files.
for (i in names(StreamStationList)){
  write.table(StreamStationList[[i]], 
              paste0(getwd(), '/StreamStat_', i, f_sf_processKey, ".txt"), 
              sep = "\t", row.names = FALSE)
}
rm(i)
list.save(x = StreamStationList, file = f_StreamStationList, type = "YAML")

# DEM for Streamflow----
#  Add DEM elevation to the gauge datasets----
#  Gather all of the DEM files together and mosaic into one file
DEM = processDEM(dir_DEMs = dir_DEM, f_DEMs = f_DEM, pCRS = pCRS)
setwd(wd_sf)
#Add the DEM elevation to the gauge dataset
#Method 1
if (exists(x = "GaugeLocs_fn")){
  GaugeLocs_fn$ElevDEM = raster::extract(x = DEM, y = GaugeLocs_fn)
}
#Method 2
if (exists(x = "GaugeLocs")){
  GaugeLocs$ElevDEM = raster::extract(x = DEM, y = GaugeLocs)
}

#   Compare the DEM elevation to the listed elevation----
setwd(wd_sf)
#Method 1
if (exists(x = "GaugeLocs_fn")){
  png('CompareStreamflowGaugeElevToDEM_fn.png', res = 300, units = 'in', width = 5, height = 5)
  plot(GaugeLocs_fn$alt_va[which((is.na(GaugeLocs_fn$alt_va) == FALSE) & (is.na(GaugeLocs_fn$ElevDEM) == FALSE))], 
       GaugeLocs_fn$ElevDEM[which((is.na(GaugeLocs_fn$alt_va) == FALSE)  & (is.na(GaugeLocs_fn$ElevDEM) == FALSE))]/.3048,
       xlab = 'USGS Reported Elevation (ft)', ylab = 'DEM Elevation (ft)', main = 'Gauge Elevations')
  lines(c(-100,1100), c(-100,1100), col = 'red')
  dev.off()
}
#Method 2
if (exists(x = "GaugeLocs")){
  png('CompareStreamflowGaugeElevToDEM.png', res = 300, units = 'in', width = 5, height = 5)
  plot(GaugeLocs$alt_va[which((is.na(GaugeLocs$alt_va) == FALSE) & (is.na(GaugeLocs$ElevDEM) == FALSE))], 
       GaugeLocs$ElevDEM[which((is.na(GaugeLocs$alt_va) == FALSE) & (is.na(GaugeLocs$ElevDEM) == FALSE))]/.3048,
       xlab = 'USGS Reported Elevation (ft)', ylab = 'DEM Elevation (ft)', main = 'Gauge Elevations')
  lines(c(-100,1100), c(-100,1100), col = 'red')
  dev.off()
}

#Fixme: Correct the elevations for those gauges that are very different than the DEM
# and have collection codes that suggest lower data quality
# Some of the gauges were likely reported in m instead of in ft in the USGS database
#Identify spatially those gauges that are more than X feet different
#plot(GaugeLocs)
#plot(GaugeLocs[which(abs(GaugeLocs$ElevDEM/.3048 - GaugeLocs$alt_va) >= 50),], col = 'red', add =T)
#plot(GaugeLocs[which(abs(GaugeLocs$ElevDEM/.3048 - GaugeLocs$alt_va) >= 100),], col = 'blue', add =T)
#plot(GaugeLocs[which(abs(GaugeLocs$ElevDEM/.3048 - GaugeLocs$alt_va) >= 200),], col = 'green', add =T)
#plot(GaugeLocs_fn)
#plot(GaugeLocs_fn[which(abs(GaugeLocs_fn$ElevDEM/.3048 - GaugeLocs_fn$alt_va) >= 50),], col = 'red', add =T)
#plot(GaugeLocs_fn[which(abs(GaugeLocs_fn$ElevDEM/.3048 - GaugeLocs_fn$alt_va) >= 100),], col = 'blue', add =T)
#plot(GaugeLocs_fn[which(abs(GaugeLocs_fn$ElevDEM/.3048 - GaugeLocs_fn$alt_va) >= 200),], col = 'green', add =T)

#  Plot reported vs. DEM elevation of gauges within the ROI----
setwd(wd_sf)
if (exists(x = "GaugeLocs_fn")){
  #  Method 1 only
  #assign elevation to the dataset
  for (i in 1:nrow(NWIS_ROI_fn)){
    NWIS_ROI_fn$ElevDEM[i] = GaugeLocs_fn$ElevDEM[GaugeLocs_fn$site_no == NWIS_ROI_fn$site_no[i]]
  }
  
  png('CompareGaugeElevToDEM_ROI_fn.png', res = 300, units = 'in', width = 5, height = 5)
  plot(NWIS_ROI_fn$alt_va, NWIS_ROI_fn$ElevDEM/.3048,
       xlab = 'USGS Reported Elevation (ft)', ylab = 'DEM Elevation (ft)', main = 'Gauge Elevations in ROI')
  lines(c(-100,1100), c(-100,1100), col = 'red')
  dev.off()
}
if (exists(x = "GaugeLocs")){
  #  Method 2 only
  #assign elevation to the dataset
  for (i in 1:nrow(NWIS_ROI)){
    NWIS_ROI$ElevDEM[i] = GaugeLocs$ElevDEM[GaugeLocs$site_no == NWIS_ROI$site_no[i]]
  }
  
  png('CompareGaugeElevToDEM_ROI.png', res = 300, units = 'in', width = 5, height = 5)
  plot(NWIS_ROI$alt_va, NWIS_ROI$ElevDEM/.3048,
       xlab = 'USGS Reported Elevation (ft)', ylab = 'DEM Elevation (ft)', main = 'Gauge Elevations in ROI')
  lines(c(-100,1100), c(-100,1100), col = 'red')
  dev.off()
}

#  Hypsometric plot for DEM in the ROI----
#Clip DEM to the ROI (no buffer) and make a spatial grid dataframe
DEM_ROI = as(mask(DEM, ROI), 'SpatialGridDataFrame')
#Plot curve
setwd(dir_DEM_out)
png('HypsometricCurve.png', res = 300, units = 'in', width = 5, height = 5)
par(mar =c(4,4,1,1))
hypsometric(DEM_ROI, col = 'black', main = 'Gwynns Falls Hypsometric Curve')
dev.off()

#Fixme: add gauges to this plot - requires modifying the function.
# Can shade in the areas above each gauge or select outlet points as horizontal steps
# Show with a map of where the elevation thresholds are located in the DEM
#  This would require plotting the line as points that are colored by their elevation values.

# Write DEM + streamflow datasets to files----
setwd(dir_DEM_out)
writeRaster(x = DEM, filename = f_DEM_mosiac, format = "GTiff")
setwd(wd_sf)
if (exists(x = "NWIS_ROI_fn")){
  writeOGR(obj = NWIS_ROI_fn, dsn = getwd(), driver = 'ESRI Shapefile', layer = paste0(f_NWIS_ROI_out, '_withDEM'))
}else if (exists(x = "NWIS_ROI")){
  writeOGR(obj = NWIS_ROI, dsn = getwd(), driver = 'ESRI Shapefile', layer = paste0(f_NWIS_ROI_out, '_withDEM'))
  writeOGR(obj = NWISstations, dsn = getwd(), driver = 'ESRI Shapefile', layer = paste0(f_NWIS_bb_out, '_withDEM'))
}

#Water Quality----
setwd(dir_wq)
#Fixme: Should have an option to save timeseries as dataframes that have variables as columns and time as rows.
#  This is challenging because the column names would have to be changed automatically in script to something meaningful.
# Method 1: Read water quality station data using the USGS function----
#  Find sites that have any N and P water quality data in a state within the ROI
#Phosphorus
phosSites <- whatWQPsites(statecode="MD", characteristicName="Phosphorus")
#Nitrogen
NitroSites <- whatWQPsites(statecode="MD", characteristicName="Nitrogen")

#  Make a spatial dataframe: Method 1----
#Nitrogen
#vector of unique coordinate systems in station data
NitroUniqueCoords = unique(NitroSites$HorizontalCoordinateReferenceSystemDatumName)

#Process to the pCRS coordinate system
GaugesLocs_NAD27 = NitroSites[which(NitroSites$HorizontalCoordinateReferenceSystemDatumName == 'NAD27'), ]
coordinates(GaugesLocs_NAD27) = c('LongitudeMeasure', 'LatitudeMeasure')
proj4string(GaugesLocs_NAD27) = CRS('+init=epsg:4267')
GaugesLocs_NAD83 = NitroSites[which(NitroSites$HorizontalCoordinateReferenceSystemDatumName == 'NAD83'), ]
coordinates(GaugesLocs_NAD83) = c('LongitudeMeasure', 'LatitudeMeasure')
proj4string(GaugesLocs_NAD83) = CRS('+init=epsg:4269')
GaugesLocs_WGS84 = NitroSites[which(NitroSites$HorizontalCoordinateReferenceSystemDatumName == 'WGS84'), ]
coordinates(GaugesLocs_WGS84) = c('LongitudeMeasure', 'LatitudeMeasure')
proj4string(GaugesLocs_WGS84) = CRS('+init=epsg:4326')
#Assuming that all the unknown coordinate systems are NAD83
GaugesLocs_U = NitroSites[which(NitroSites$HorizontalCoordinateReferenceSystemDatumName == 'UNKWN'), ]
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
#vector of unique coordinate systems in station data
PhosUniqueCoords = unique(phosSites$HorizontalCoordinateReferenceSystemDatumName)

#Process to the pCRS coordinate system
GaugesLocs_NAD27 = phosSites[which(phosSites$HorizontalCoordinateReferenceSystemDatumName == 'NAD27'), ]
coordinates(GaugesLocs_NAD27) = c('LongitudeMeasure', 'LatitudeMeasure')
proj4string(GaugesLocs_NAD27) = CRS('+init=epsg:4267')
GaugesLocs_NAD83 = phosSites[which(phosSites$HorizontalCoordinateReferenceSystemDatumName == 'NAD83'), ]
coordinates(GaugesLocs_NAD83) = c('LongitudeMeasure', 'LatitudeMeasure')
proj4string(GaugesLocs_NAD83) = CRS('+init=epsg:4269')
GaugesLocs_WGS84 = phosSites[which(phosSites$HorizontalCoordinateReferenceSystemDatumName == 'WGS84'), ]
coordinates(GaugesLocs_WGS84) = c('LongitudeMeasure', 'LatitudeMeasure')
proj4string(GaugesLocs_WGS84) = CRS('+init=epsg:4326')
#Assuming that all the unknown coordinate systems are NAD83
GaugesLocs_U = phosSites[which(phosSites$HorizontalCoordinateReferenceSystemDatumName == 'UNKWN'), ]
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

#  Clip to ROI----
WQstations_ROI_N = GaugeLocs_WQN[ROI_buff,]
WQstations_ROI_P = GaugeLocs_WQP[ROI_buff,]

#Fix the dataframe from tibble
WQstations_ROI_N@data = as.data.frame(WQstations_ROI_N@data)
WQstations_ROI_P@data = as.data.frame(WQstations_ROI_P@data)

# Method 2: Read water quality station data from file----
WQstations = read.csv(f_WQgauges, stringsAsFactors = FALSE)
#  Make a spatial dataframe: Method 2----
coordinates(WQstations) = c('LongitudeMeasure', 'LatitudeMeasure')
#Split dataset according to the coordinate reference systems used
#vector of unique coordinate systems in station data
WQStationUniqueCoords = unique(WQstations$HorizontalCoordinateReferenceSystemDatumName)

#Process to the pCRS coordinate system
GaugesLocs_NAD27 = WQstations[which(WQstations$HorizontalCoordinateReferenceSystemDatumName == 'NAD27'), ]
proj4string(GaugesLocs_NAD27) = CRS('+init=epsg:4267')
GaugesLocs_NAD83 = WQstations[which(WQstations$HorizontalCoordinateReferenceSystemDatumName == 'NAD83'), ]
proj4string(GaugesLocs_NAD83) = CRS('+init=epsg:4269')
#Assuming unknown is NAD83
GaugesLocs_U = WQstations[which(WQstations$HorizontalCoordinateReferenceSystemDatumName == 'UNKWN'), ]
proj4string(GaugesLocs_U) = CRS('+init=epsg:4269')
GaugesLocs_WGS84 = WQstations[which(WQstations$HorizontalCoordinateReferenceSystemDatumName == 'WGS84'), ]
proj4string(GaugesLocs_WGS84) = CRS('+init=epsg:4326')
#Transform to NAD83 UTM Zone 18N
GaugeLocs_NAD27 = spTransform(GaugesLocs_NAD27, CRS(pCRS))
GaugeLocs_NAD83 = spTransform(GaugesLocs_NAD83, CRS(pCRS))
GaugeLocs_WGS84 = spTransform(GaugesLocs_WGS84, CRS(pCRS))
GaugeLocs_U = spTransform(GaugesLocs_U, CRS(pCRS))
#Join to one dataset again
WQGaugeLocs = rbind(GaugeLocs_NAD27, GaugeLocs_NAD83, GaugeLocs_WGS84, GaugeLocs_U)
#Remove separate datasets
rm(GaugeLocs_NAD27, GaugeLocs_NAD83, GaugeLocs_WGS84, GaugeLocs_U, GaugesLocs_NAD27, GaugesLocs_NAD83, GaugesLocs_WGS84, GaugesLocs_U)

#  Clip to ROI----
WQstations_ROI = WQGaugeLocs[ROI_buff,]

#Phosphorus sites
phosSites <- whatWQPsites(statecode="MD", characteristicName="Phosphorus")
#Nitrogen sites
NitroSites <- whatWQPsites(statecode="MD", characteristicName="Nitrogen")

#Select only those sites in the ROI that have nitrogen data
WQstations_ROI_N_f = WQstations_ROI[which(WQstations_ROI$MonitoringLocationIdentifier %in% NitroSites$MonitoringLocationIdentifier),]
#Select only those sites in the ROI that have Phosphorus data
WQstations_ROI_P_f = WQstations_ROI[which(WQstations_ROI$MonitoringLocationIdentifier %in% phosSites$MonitoringLocationIdentifier),]

# Plot TN and TP sampling locations on a map----
if (exists(x = "WQstations_ROI_N")){
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
}
if (exists(x = "WQstations_ROI_N_f")){
  png('TNTPsites_f.png', res = 300, units = 'in', width = 6, height = 6)
  plot(ROI)
  plot(WQstations_ROI_N_f, pch = 16, add = TRUE, col = 'red')
  plot(WQstations_ROI_P_f, pch = 16, add = TRUE, col = 'blue')
  plot(WQstations_ROI_P_f[WQstations_ROI_N_f,], pch = 16, add = TRUE, col = 'purple')
  plot(WQstations_ROI_N_f[WQstations_ROI_P_f,], pch = 16, add = TRUE, col = 'purple')
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 365000, yb = 4346000, len = 700, col = 'black', lab = 'N')
  legend('topright', title = 'Water Quality Sites', legend = c('T Nitrogen Only', 'T Phosphorus Only', 'Both'), col = c('red', 'blue', 'purple'), pch = 16)
  dev.off()
}

# Download data for those sites in parallel----
#   Use only the unique gauge numbers in the dataset 
#   (repeats occur when multiple variables are available for a gauge)
uniqueWQNums_N = unique(WQstations_ROI_N$MonitoringLocationIdentifier)
uniqueWQNums_P = unique(WQstations_ROI_P$MonitoringLocationIdentifier)

#Make directories for storing Nitrogen and Phosphorous data
wd_N = paste0(getwd(), '/Nitrogen')
wd_P = paste0(getwd(), '/Phosphorus')
dir.create(path = wd_N, showWarnings = FALSE)
dir.create(path = wd_P, showWarnings = FALSE)

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
              sep = "\t", row.names = FALSE)
}
p = foreach(i = uniqueWQNums_P, .packages = 'dataRetrieval') %dopar% {
  setwd(wd_P)
  # Read all of the data for the provided station number
  stationData <- readWQPqw(siteNumbers = i, parameterCd = "")
  write.table(stationData, 
              paste0(getwd(), '/Phosphorus_', i, ".txt"), 
              sep = "\t", row.names = FALSE)
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
#  Gather the records for each gauge into a list of dataframes----
NitroStationList = makeWQStationList(pattern = 'Nitrogen_', wd = wd_N, Sites = WQstations_ROI_N)

#  Extract timeseries for each variable for each site----
#Cycle through all of the unique combinations of ResultSampleFractionText and CharacteristicName 
# and write separate text files for each variable
extractWQdata(NitroStationList, fName = "Nitrogen")

#Fixme: should have a function for users to specify for which of the extracted variables they want to have detailed plots.
# Template is processing of nitrogen and phosphorus data. Phosphorus has detection limits components.

#   Process Total Nitrogen Data----
#List of all Total Nitrogen gauge datasets
TN = selectWQDataType(wd = wd_N, charName = '_cnNitrogen', resName = 'Total')

#    Check for duplicate observations----
#Strict check on the specified column names only, disregardinig all other columns. 
TN = checkDuplicatesAndRemove(TN, colNames = c('ResultMeasureValue', 'SortDateTime'))

#    Check for zeros and negative values in the records----
TN = checkZerosNegs(TN, Var = 'ResultMeasureValue')

#See if any stations have negative values and add to WQstations_ROI_N spatial dataset
WQstations_ROI_N = addNegsToSpatialDataset(StationList = TN, SpatialDataset = WQstations_ROI_N, site_D = 'MonitoringLocationIdentifier', site_SL = 'MonitoringLocationIdentifier')

#See if any stations have zero values and add to spatial dataset
WQstations_ROI_N = addZerosToSpatialDataset(StationList = TN, SpatialDataset = WQstations_ROI_N, site_D = 'MonitoringLocationIdentifier', site_SL = 'MonitoringLocationIdentifier')

#     Plot a map of the stations with zero and negative records----
png(paste0('TN_WQP_ZerosNegsMap.png'), res = 300, units = 'in', width = 6, height = 6)
plot(ROI)
#All gauges
plot(WQstations_ROI_N, pch = 16, add = TRUE)
#Gauges with zeros
plot(WQstations_ROI_N[!is.na(WQstations_ROI_N$Zero),], pch = 16, col = 'red', add = TRUE)
#Gauges with negative values
plot(WQstations_ROI_N[!is.na(WQstations_ROI_N$Neg),], pch = 16, col = 'blue', add = TRUE)
#Gauges with both
plot(WQstations_ROI_N[which(!is.na(WQstations_ROI_N$Neg) & !is.na(WQstations_ROI_N$Zero)),], pch = 16, col = 'purple', add = TRUE)
# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
legend('topright', title = 'Streamflow Stations', legend = c('Zeros in Record', 'Negatives in Record', 'Both', 'Neither'), pch = 16, col = c('red', 'blue', 'purple', 'black'), bty = 'n')
dev.off()

#    Aggregate into average concentrations, if desired---- 
TN_agg = aggregateTimesteps(StationList = TN, aggVal = c('d', 'm', 'a'), aggVar = 'ResultMeasureValue', date = 'SortDate', site = 'MonitoringLocationIdentifier', fun = 'mean')
TN_d = TN_agg$daily
TN_m = TN_agg$mthyr
TN_a = TN_agg$ann
rm(TN_agg)

#     Handle missing data in the daily, monthly, and annual aggregated timeseries----
TN_d2 = FillMissingDates_par(Dataset = WQstations_ROI_N, StationList = TN_d, Var = 'ResultMeasureValue', 
                        Date = 'SortDate', gapType = 'd', site_no_D = 'MonitoringLocationIdentifier', 
                        site_no_SL = 'MonitoringLocationIdentifier', NoNAcols = 'MonitoringLocationIdentifier', 
                        NumCores = detectCores()-1)
WQstations_ROI_N = TN_d2$Dataset
TN_d = TN_d2$StationList
rm(TN_d2)

TN_m2 = FillMissingDates_par(Dataset = WQstations_ROI_N, StationList = TN_m, Var = 'ResultMeasureValue', 
                         Date = 'YrMthDy', gapType = 'm', site_no_D = 'MonitoringLocationIdentifier', 
                         site_no_SL = 'MonitoringLocationIdentifier', NoNAcols = c('MonitoringLocationIdentifier'), 
                         NumCores = detectCores()-1)
WQstations_ROI_N = TN_m2$Dataset
TN_m = TN_m2$StationList
rm(TN_m2)

TN_a2 = FillMissingDates(Dataset = WQstations_ROI_N, StationList = TN_a, Var = 'ResultMeasureValue', 
                         Date = 'YrMthDy', gapType = 'a', site_no_D = 'MonitoringLocationIdentifier', 
                         site_no_SL = 'MonitoringLocationIdentifier', NoNAcols = c('MonitoringLocationIdentifier'))
WQstations_ROI_N = TN_a2$Dataset
TN_a = TN_a2$StationList
rm(TN_a2)
  
#   Plot TN timeseries----
#Fixme: allow for plotting instantaneous time using POSIX format for date and time. Currently plots all by date, instead of an average by day.
for (i in 1:length(TN)){
  png(paste0('TN_Timeseries_', TN[[i]]$MonitoringLocationIdentifier[1],'.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = TN[[i]]$ResultMeasureValue, x = as.Date(TN[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = paste0(TN[[i]]$ResultSampleFractionText[1], " ", TN[[i]]$CharacteristicName[1], " (", TN[[i]]$ResultMeasure.MeasureUnitCode[1], ")"), 
       main = paste0('Station #', TN[[i]]$MonitoringLocationIdentifier[1]),
       ylim = ylim_TN, xlim = xlim_WQdates)
  
  #Check for and add detection limits
  dl = unique(TN[[i]]$DetectionQuantitationLimitMeasure.MeasureValue)
  #remove NAs
  dl = dl[!is.na(dl)]
  if(length(dl) != 0){
    #Add detection limits
    #Lower detection limit
    par(new = TRUE)
    plot(y = TN[[i]]$DetectionQuantitationLimitMeasure.MeasureValue[grep(pattern = "Low", x = TN[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)], 
         x = as.Date(TN[[i]]$SortDate[grep(pattern = "Low", x = TN[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)]), pch = 16, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = ylim_TN, xlim = xlim_WQdates, axes = FALSE,
         col = 'blue')
    #Upper detection limit
    par(new = TRUE)
    plot(y = TN[[i]]$DetectionQuantitationLimitMeasure.MeasureValue[grep(pattern = "Up", x = TN[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)], 
         x = as.Date(TN[[i]]$SortDate[grep(pattern = "Up", x = TN[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)]), pch = 16, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = ylim_TN, xlim = xlim_WQdates, axes = FALSE,
         col = 'red')
    
    legend('topright', title = 'Detection Limits', legend = c('Upper', 'Lower'), col = c('red', 'blue'), pch = 16)
  }
  if(length(as.Date(TN[[i]]$SortDate[which(is.na(TN[[i]]$ResultMeasureValue) & is.na(TN[[i]]$DetectionQuantitationLimitMeasure.MeasureValue))])) > 0){
    #Plot where there are NA days that do not correspond to detection limits
    par(new = TRUE)
    plot(y = rep(0, length(as.Date(TN[[i]]$SortDate[which(is.na(TN[[i]]$ResultMeasureValue) & is.na(TN[[i]]$DetectionQuantitationLimitMeasure.MeasureValue))]))), 
         x = as.Date(TN[[i]]$SortDate[which(is.na(TN[[i]]$ResultMeasureValue) & is.na(TN[[i]]$DetectionQuantitationLimitMeasure.MeasureValue))]), pch = 4, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = ylim_TN, xlim = xlim_WQdates, axes = FALSE,
         col = 'grey')
    legend('right', legend = 'NA values', col = 'grey', pch = 4)
  }
  dev.off()
  
  #With map and timeseries
  png(paste0('TN_Timeseries_Map_', TN[[i]]$MonitoringLocationIdentifier[1],'.png'), res = 300, units = 'in', width = 10, height = 10)
  layout(rbind(c(1,2)))
  
  plot(y = TN[[i]]$ResultMeasureValue, x = as.Date(TN[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = paste0(TN[[i]]$ResultSampleFractionText[1], " ", TN[[i]]$CharacteristicName[1], " (", TN[[i]]$ResultMeasure.MeasureUnitCode[1], ")"), 
       main = paste0('Station #', TN[[i]]$MonitoringLocationIdentifier[1]),
       ylim = ylim_TN, xlim = xlim_WQdates)
  if(length(dl) != 0){
    #Add detection limits
    #Lower detection limit
    par(new = TRUE)
    plot(y = TN[[i]]$DetectionQuantitationLimitMeasure.MeasureValue[grep(pattern = "Low", x = TN[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)], 
         x = as.Date(TN[[i]]$SortDate[grep(pattern = "Low", x = TN[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)]), pch = 16, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = ylim_TN, xlim = xlim_WQdates, axes = FALSE,
         col = 'blue')
    #Upper detection limit
    par(new = TRUE)
    plot(y = TN[[i]]$DetectionQuantitationLimitMeasure.MeasureValue[grep(pattern = "Up", x = TN[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)], 
         x = as.Date(TN[[i]]$SortDate[grep(pattern = "Up", x = TN[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)]), pch = 16, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = ylim_TN, xlim = xlim_WQdates, axes = FALSE,
         col = 'red')
    
    legend('topright', title = 'Detection Limits', legend = c('Upper', 'Lower'), col = c('red', 'blue'), pch = 16)
  }
  if(length(as.Date(TN[[i]]$SortDate[which(is.na(TN[[i]]$ResultMeasureValue) & is.na(TN[[i]]$DetectionQuantitationLimitMeasure.MeasureValue))])) > 0){
    #Plot where there are NA days that do not correspond to detection limits
    par(new = TRUE)
    plot(y = rep(0, length(as.Date(TN[[i]]$SortDate[which(is.na(TN[[i]]$ResultMeasureValue) & is.na(TN[[i]]$DetectionQuantitationLimitMeasure.MeasureValue))]))), 
         x = as.Date(TN[[i]]$SortDate[which(is.na(TN[[i]]$ResultMeasureValue) & is.na(TN[[i]]$DetectionQuantitationLimitMeasure.MeasureValue))]), pch = 4, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = ylim_TN, xlim = xlim_WQdates, axes = FALSE,
         col = 'grey')
    legend('right', legend = 'NA values', col = 'grey', pch = 4)
  }
  
  #map
  plot(ROI)
  #All gauges
  plot(WQstations_ROI_N, pch = 16, add = TRUE)
  #Gauge selected for timeseries plot
  plot(WQstations_ROI_N[which(WQstations_ROI_N$MonitoringLocationIdentifier == TN[[i]]$MonitoringLocationIdentifier[1]),], pch = 16, col = 'red', add = TRUE)
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
  legend('topright', title = 'Nitrogen Stations', legend = c('Selected', 'Other'), pch = 16, col = c('red', 'black'), bty = 'n')
  title(main = paste0('Station #', TN[[i]]$MonitoringLocationIdentifier[1]))
  dev.off()
  
  #Fixme: 20 is temporary - the whole if statement should be made more sophisticated.
  if (length(which(!is.na(TN_d[[i]]$ResultMeasureValue)))>20){
    #hydroTSM plots----
    #daily timeseries
    dts = zoo(TN_d[[i]]$ResultMeasureValue, order.by = TN_d[[i]]$SortDate)
    
    #monthly timeseries for mean daily flows in a month
    mts = zoo(TN_m[[i]]$ResultMeasureValue, order.by = as.Date(TN_m[[i]]$YrMthDy))
    #Monthly matrix
    M = formatMonthlyMatrix(mts)
    # Plotting the monthly values as matrixplot 
    png(paste0('TNMonthly_', TN[[i]]$MonitoringLocationIdentifier[1],'.png'), res = 300, units = 'in', width = 6, height = 3)
    print(matrixplot(M, ColorRamp="Precipitation", main="Mean Monthly Nitrogen"))
    dev.off()
    
    #Summary EDA plots
    png(paste0('TN_EDA_', TN[[i]]$MonitoringLocationIdentifier[1],'.png'), res = 300, units = 'in', width = 10, height = 10)
    hydroplot(dts, FUN = mean, col = 'black', var.unit = 'mean mg/L as N')
    dev.off()
    
    #Seasonal flows
    #Fixme: check that every season is represented before using this.
    png(paste0('TN_EDA_Seasonal_', TN[[i]]$MonitoringLocationIdentifier[1],'.png'), res = 300, units = 'in', width = 10, height = 10)
    hydroplot(dts, FUN = mean, col = 'black', var.unit = 'mean mg/L as N', pfreq = 'seasonal', ylab = 'Mean Nitrogen (mg/L)')
    dev.off()
  }
}
rm(i, dts, mts, M, dl)

#Fixme: Search for flow-normalized outliers

#   Plot histograms of data with detection limits----
for (i in 1:length(TN)){
  png(paste0('TN_hist_', TN[[i]]$MonitoringLocationIdentifier, '.png'), res = 300, units = 'in', width = 5, height = 5)
  hist(c(log10(TN[[i]]$ResultMeasureValue), log10(TN[[i]]$DetectionQuantitationLimitMeasure.MeasureValue)),
       xlab = 'log10(Total Nitrogen [mg/L as N])', main = TN[[i]]$MonitoringLocationIdentifier[1])
  dev.off()
}
rm(i)

#Fixme: note any changes in rounding that occur from measurement device precision (USGS recommended)

# Process Phosphorus Data----
setwd(wd_P)
#  Gather the records for each gauge into a list of dataframes----
PhosStationList = makeWQStationList(pattern = 'Phosphorus_', wd = wd_P, Sites = WQstations_ROI_P)

#  Extract timeseries for each variable for each site----
#Cycle through all of the unique combinations of ResultSampleFractionText and CharacteristicName and return separate text files for each variable
extractWQdata(StationList = PhosStationList, fName = "Phosphorus")

#   Process Total Phosphorus Data----
#List of all Total Phosphorus gauge datasets
TP = selectWQDataType(wd = wd_P, charName = '_cnPhosphorus', resName = 'Total')

#Some stations have detection limits

#    Check for duplicate observations----
#Strict check on the specified column names only, disregardinig all other columns.
TP = checkDuplicatesAndRemove(TP, colNames = c('ResultMeasureValue', 'SortDateTime'))
#    Check for zeros and negative values in the records----
TP = checkZerosNegs(TP, Var = 'ResultMeasureValue')

#See if any stations have negative values and add to WQstations_ROI_P spatial dataset
WQstations_ROI_P = addNegsToSpatialDataset(StationList = TP, SpatialDataset = WQstations_ROI_P, site_D = 'MonitoringLocationIdentifier', site_SL = 'MonitoringLocationIdentifier')

#See if any stations have zero values and add to spatial dataset
WQstations_ROI_P = addZerosToSpatialDataset(StationList = TP, SpatialDataset = WQstations_ROI_P, site_D = 'MonitoringLocationIdentifier', site_SL = 'MonitoringLocationIdentifier')

#     Plot a map of the stations with zero and negative records----
png(paste0('TP_WQP_ZerosNegsMap.png'), res = 300, units = 'in', width = 6, height = 6)
plot(ROI)
#All gauges
plot(WQstations_ROI_P, pch = 16, add = TRUE)
#Gauges with zeros
plot(WQstations_ROI_P[!is.na(WQstations_ROI_P$Zero),], pch = 16, col = 'red', add = TRUE)
#Gauges with negative values
plot(WQstations_ROI_P[!is.na(WQstations_ROI_P$Neg),], pch = 16, col = 'blue', add = TRUE)
#Gauges with both
plot(WQstations_ROI_P[which(!is.na(WQstations_ROI_P$Neg) & !is.na(WQstations_ROI_P$Zero)),], pch = 16, col = 'purple', add = TRUE)
# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
legend('topright', title = 'Streamflow Stations', legend = c('Zeros in Record', 'Negatives in Record', 'Both', 'Neither'), pch = 16, col = c('red', 'blue', 'purple', 'black'), bty = 'n')
dev.off()

#    Aggregate into average concentrations, if desired---- 
TP_agg = aggregateTimesteps(StationList = TP, aggVal = c('d', 'm', 'a'), aggVar = 'ResultMeasureValue', date = 'SortDate', site = 'MonitoringLocationIdentifier', fun = 'mean')
TP_d = TP_agg$daily
TP_m = TP_agg$mthyr
TP_a = TP_agg$ann
rm(TP_agg)

#     Handle missing data in the daily, monthly, and annual aggregated timeseries----
TP_d2 = FillMissingDates_par(Dataset = WQstations_ROI_P, StationList = TP_d, Var = 'ResultMeasureValue', 
                         Date = 'SortDate', gapType = 'd', site_no_D = 'MonitoringLocationIdentifier', 
                         site_no_SL = 'MonitoringLocationIdentifier', NoNAcols = 'MonitoringLocationIdentifier', 
                         NumCores = detectCores()-1)
WQstations_ROI_P = TP_d2$Dataset
TP_d = TP_d2$StationList
rm(TP_d2)

TP_m2 = FillMissingDates_par(Dataset = WQstations_ROI_P, StationList = TP_m, Var = 'ResultMeasureValue', 
                         Date = 'YrMthDy', gapType = 'm', site_no_D = 'MonitoringLocationIdentifier', 
                         site_no_SL = 'MonitoringLocationIdentifier', NoNAcols = c('MonitoringLocationIdentifier'), 
                         NumCores = detectCores()-1)
WQstations_ROI_P = TP_m2$Dataset
TP_m = TP_m2$StationList
rm(TP_m2)

TP_a2 = FillMissingDates(Dataset = WQstations_ROI_P, StationList = TP_a, Var = 'ResultMeasureValue', 
                         Date = 'YrMthDy', gapType = 'a', site_no_D = 'MonitoringLocationIdentifier', 
                         site_no_SL = 'MonitoringLocationIdentifier', NoNAcols = c('MonitoringLocationIdentifier'))
WQstations_ROI_P = TP_a2$Dataset
TP_a = TP_a2$StationList
rm(TP_a2)

#   Plot TP timeseries----
for (i in 1:length(TP)){
  png(paste0('TP_Timeseries_', TP[[i]]$MonitoringLocationIdentifier[1],'.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = TP[[i]]$ResultMeasureValue, x = as.Date(TP[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = paste0(TP[[i]]$ResultSampleFractionText[1], " ", TP[[i]]$CharacteristicName[1], " (", TP[[i]]$ResultMeasure.MeasureUnitCode[1], ")"), 
       main = paste0('Station #', TP[[i]]$MonitoringLocationIdentifier[1]),
       ylim = ylim_TP, xlim = xlim_WQdates)
  
  #Check for and add detection limits
  dl = unique(TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue)
  #remove NAs
  dl = dl[!is.na(dl)]
  if(length(dl) != 0){
    #Add detection limits
    #Lower detection limit
    par(new = TRUE)
    plot(y = TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue[grep(pattern = "Low", x = TP[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)], 
         x = as.Date(TP[[i]]$SortDate[grep(pattern = "Low", x = TP[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)]), pch = 16, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = ylim_TP, xlim = xlim_WQdates, axes = FALSE,
         col = 'blue')
    #Upper detection limit
    par(new = TRUE)
    plot(y = TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue[grep(pattern = "Up", x = TP[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)], 
         x = as.Date(TP[[i]]$SortDate[grep(pattern = "Up", x = TP[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)]), pch = 16, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = ylim_TP, xlim = xlim_WQdates, axes = FALSE,
         col = 'red')
    
    legend('topright', title = 'Detection Limits', legend = c('Upper', 'Lower'), col = c('red', 'blue'), pch = 16)
  }
  if(length(as.Date(TP[[i]]$SortDate[which(is.na(TP[[i]]$ResultMeasureValue) & is.na(TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue))])) > 0){
    #Plot where there are NA days that do not correspond to detection limits
    par(new = TRUE)
    plot(y = rep(0, length(as.Date(TP[[i]]$SortDate[which(is.na(TP[[i]]$ResultMeasureValue) & is.na(TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue))]))), 
         x = as.Date(TP[[i]]$SortDate[which(is.na(TP[[i]]$ResultMeasureValue) & is.na(TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue))]), pch = 4, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = ylim_TP, xlim = xlim_WQdates, axes = FALSE,
         col = 'grey')
    legend('right', legend = 'NA values', col = 'grey', pch = 4)
  }
  dev.off()
  
  #With map and timeseries
  png(paste0('TP_Timeseries_Map_', TP[[i]]$MonitoringLocationIdentifier[1],'.png'), res = 300, units = 'in', width = 10, height = 10)
  layout(rbind(c(1,2)))
  
  plot(y = TP[[i]]$ResultMeasureValue, x = as.Date(TP[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = paste0(TP[[i]]$ResultSampleFractionText[1], " ", TP[[i]]$CharacteristicName[1], " (", TP[[i]]$ResultMeasure.MeasureUnitCode[1], ")"), 
       main = paste0('Station #', TP[[i]]$MonitoringLocationIdentifier[1]),
       ylim = ylim_TP, xlim = xlim_WQdates)
  
  #Check for and add detection limits
  dl = unique(TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue)
  #remove NAs
  dl = dl[!is.na(dl)]
  if(length(dl) != 0){
    #Add detection limits
    #Lower detection limit
    par(new = TRUE)
    plot(y = TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue[grep(pattern = "Low", x = TP[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)], 
         x = as.Date(TP[[i]]$SortDate[grep(pattern = "Low", x = TP[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)]), pch = 16, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = ylim_TP, xlim = xlim_WQdates, axes = FALSE,
         col = 'blue')
    #Upper detection limit
    par(new = TRUE)
    plot(y = TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue[grep(pattern = "Up", x = TP[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)], 
         x = as.Date(TP[[i]]$SortDate[grep(pattern = "Up", x = TP[[i]]$DetectionQuantitationLimitTypeName, ignore.case = TRUE)]), pch = 16, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = ylim_TP, xlim = xlim_WQdates, axes = FALSE,
         col = 'red')
    
    legend('topright', title = 'Detection Limits', legend = c('Upper', 'Lower'), col = c('red', 'blue'), pch = 16)
  }
  if(length(as.Date(TP[[i]]$SortDate[which(is.na(TP[[i]]$ResultMeasureValue) & is.na(TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue))])) > 0){
    #Plot where there are NA days that do not correspond to detection limits
    par(new = TRUE)
    plot(y = rep(0, length(as.Date(TP[[i]]$SortDate[which(is.na(TP[[i]]$ResultMeasureValue) & is.na(TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue))]))), 
         x = as.Date(TP[[i]]$SortDate[which(is.na(TP[[i]]$ResultMeasureValue) & is.na(TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue))]), pch = 4, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = ylim_TP, xlim = xlim_WQdates, axes = FALSE,
         col = 'grey')
    legend('right', legend = 'NA values', col = 'grey', pch = 4)
  }
  
  plot(ROI)
  #All gauges
  plot(WQstations_ROI_P, pch = 16, add = TRUE)
  #Gauge selected for timeseries plot
  plot(WQstations_ROI_P[which(WQstations_ROI_P$MonitoringLocationIdentifier == TP[[i]]$MonitoringLocationIdentifier[1]),], pch = 16, col = 'red', add = TRUE)
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
  legend('topright', title = 'Phosphorus Stations', legend = c('Selected', 'Other'), pch = 16, col = c('red', 'black'), bty = 'n')
  title(main = paste0('Station #', TP[[i]]$MonitoringLocationIdentifier[1]))
  dev.off()
  
  #Fixme: 20 is temporary - the whole if statement should be made more sophisticated.
  if (length(which(!is.na(TP_d[[i]]$ResultMeasureValue)))>20){
    #hydroTSM plots----
    #daily timeseries
    dts = zoo(TP_d[[i]]$ResultMeasureValue, order.by = TP_d[[i]]$SortDate)
    
    #monthly timeseries for mean daily flows in a month
    mts = zoo(TP_m[[i]]$ResultMeasureValue, order.by = as.Date(TP_m[[i]]$YrMthDy))
    #Monthly matrix
    M = formatMonthlyMatrix(mts)
    # Plotting the monthly values as matrixplot 
    png(paste0('TPMonthly_', TP[[i]]$MonitoringLocationIdentifier[1],'.png'), res = 300, units = 'in', width = 6, height = 3)
    print(matrixplot(M, ColorRamp="Precipitation", main="Mean Monthly Phosphorus"))
    dev.off()
    
    #Summary EDA plots
    png(paste0('TP_EDA_', TP[[i]]$MonitoringLocationIdentifier[1],'.png'), res = 300, units = 'in', width = 10, height = 10)
    hydroplot(dts, FUN = mean, col = 'black', var.unit = 'mean ug/L as P')
    dev.off()
    
    #Seasonal flows
    #Fixme: check that every season is represented before using this.
    png(paste0('TP_EDA_Seasonal_', TP[[i]]$MonitoringLocationIdentifier[1],'.png'), res = 300, units = 'in', width = 10, height = 10)
    hydroplot(dts, FUN = mean, col = 'black', var.unit = 'mean ug/L as P', pfreq = 'seasonal', ylab = 'Mean Phosphorus (ug/L)')
    dev.off()
  }
}
rm(i, dts, mts, M, dl)

#Fixme: Search for flow-normalized outliers

#   Plot histograms of data with detection limits----
for (i in 1:length(TP)){
  png(paste0('TP_hist_', TP[[i]]$MonitoringLocationIdentifier, '.png'), res = 300, units = 'in', width = 5, height = 5)
  hist(c(log10(TP[[i]]$ResultMeasureValue), log10(TP[[i]]$DetectionQuantitationLimitMeasure.MeasureValue)),
       xlab = 'log10(Total Phosphorus [mg/L as P])', main = TP[[i]]$MonitoringLocationIdentifier[1])
  dev.off()
}
rm(i)

# Write water quality data to files----
setwd(dir = wd_N)
writeOGR(WQstations_ROI_N, dsn = getwd(), layer = 'NitrogenSites', driver = "ESRI Shapefile")
list.save(x = TN, file = f_TNSiteList, type = "YAML")
list.save(x = TN_d, file = f_TNSiteList_daily, type = "YAML")
list.save(x = TN_m, file = f_TNSiteList_monthly, type = "YAML")
list.save(x = TN_a, file = f_TNSiteList_annual, type = "YAML")
setwd(dir = wd_P)
writeOGR(WQstations_ROI_P, dsn = getwd(), layer = 'PhosphorusSites', driver = "ESRI Shapefile")
list.save(x = TP, file = f_TPSiteList, type = "YAML")
list.save(x = TP_d, file = f_TPSiteList_daily, type = "YAML")
list.save(x = TP_m, file = f_TPSiteList_monthly, type = "YAML")
list.save(x = TP_a, file = f_TPSiteList_annual, type = "YAML")

# Fixme: Try plotting the water quality data using R tools like EGRET----
#NOAA weather station data----
setwd(dir_weather)
#Gather all ghcnd stations
AllNOAAstations = ghcnd_stations()
#Make a spatial dataframe, and clip to ROI
coordinates(AllNOAAstations) = c('longitude', 'latitude')
proj4string(AllNOAAstations) = CRS('+init=epsg:4326')
#buffer the ROI a bit - assumes a UTM coordinate system. Will not work with other systems.
ROI_bufferW_WGS = spTransform(ROI_bufferW, CRS('+init=epsg:4326'))
NOAAstations_locs = AllNOAAstations[ROI_bufferW_WGS,]
NOAAstations_locs = spTransform(NOAAstations_locs, CRS(pCRS))

# Download data and store in a list----
MetStations = list()
#Fixme: downloading in parallel throws errors.
#cl = makeCluster(detectCores() - 1)
#registerDoParallel(cl)
#foreach (i = 1:length(unique(NOAAstations_locs$id)), .packages = 'rnoaa') %dopar%{
#  ghcnd(stationid = unique(NOAAstations_locs$id)[i], refresh = TRUE)
#  m = meteo_tidy_ghcnd(stationid = unique(NOAAstations_locs$id)[i], var = 'all', keep_flags = TRUE)
#  MetStations = c(MetStations, list(m))
#}
#stopCluster(cl)
#rm(cl)

#Serial
for (i in 1:length(unique(NOAAstations_locs$id))){
  ghcnd(stationid = unique(NOAAstations_locs$id)[i], refresh = FALSE)
  m = meteo_tidy_ghcnd(stationid = unique(NOAAstations_locs$id)[i], var = 'all', keep_flags = TRUE)
  MetStations = c(MetStations, list(m))
}
names(MetStations) = unique(NOAAstations_locs$id)
rm(m)

#Move files from the cache directory to a permanent directory
wd_NOAA = paste0(dir_weather, '\\NOAA')
dir.create(path = wd_NOAA)
file.copy(from = user_cache_dir(appname = 'rnoaa', version = NULL, opinion = TRUE, expand = TRUE), to = wd_NOAA, recursive = TRUE)
setwd(wd_NOAA)

# Map of the type of data at each station----
for (i in 1:length(unique(NOAAstations_locs$element))){
  #map for each type of data recorded
  png(paste0('NOAA_Datatype_', unique(NOAAstations_locs$element)[i],'.png'), width = 5, height = 5, units = 'in', res = 300)
  plot(NOAAstations_locs, col = 'white')
  plot(ROI, add = T)
  plot(NOAAstations_locs, add = T)
  plot(NOAAstations_locs[NOAAstations_locs$element == unique(NOAAstations_locs$element)[i],], col = 'blue', add = T)
  
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4343000, len = 700, col = 'black', lab = 'N')
  title(main = paste0(unique(NOAAstations_locs$element)[i], ' Stations'))
  dev.off()
}
rm(i)

# Fill in missing dates----
#make dataframe of the spatial dataset
NOAAstations_locs@data = as.data.frame(NOAAstations_locs@data)

#  Precipitation----
MetStations_Precip = FillMissingDates_par(Dataset = NOAAstations_locs[NOAAstations_locs$element == 'PRCP',], 
                                          StationList = MetStations[names(MetStations) %in% NOAAstations_locs[NOAAstations_locs$element == 'PRCP',]$id], 
                                          Var = 'prcp', Date = 'date', gapType = 'd', site_no_D = 'id', site_no_SL = 'id', NoNAcols = 'id',
                                          NumCores = detectCores()-1)
#Extract data from the function return
NOAAstations_locs_Precip = MetStations_Precip$Dataset
MetStations_Precip = MetStations_Precip$StationList

#  Maximum Temperature----
MetStations_TMAX = FillMissingDates_par(Dataset = NOAAstations_locs[NOAAstations_locs$element == 'TMAX',], 
                                        StationList = MetStations[names(MetStations) %in% NOAAstations_locs[NOAAstations_locs$element == 'TMAX',]$id], 
                                        Var = 'tmax', Date = 'date', gapType = 'd', site_no_D = 'id', site_no_SL = 'id', NoNAcols = 'id',
                                        NumCores = detectCores()-1)
#Extract data from the function return
NOAAstations_locs_TMAX = MetStations_TMAX$Dataset
MetStations_TMAX = MetStations_TMAX$StationList

#  Minimum Temperature----
MetStations_TMIN = FillMissingDates_par(Dataset = NOAAstations_locs[NOAAstations_locs$element == 'TMIN',], 
                                        StationList = MetStations[names(MetStations) %in% NOAAstations_locs[NOAAstations_locs$element == 'TMIN',]$id], 
                                        Var = 'tmin', Date = 'date', gapType = 'd', site_no_D = 'id', site_no_SL = 'id', NoNAcols = 'id',
                                        NumCores = detectCores()-1)
#Extract data from the function return
NOAAstations_locs_TMIN = MetStations_TMIN$Dataset
MetStations_TMIN = MetStations_TMIN$StationList

# Plot the met station timeseries----
#  Precipitation----
for (i in 1:length(MetStations_Precip)){
  png(paste0('Precip_Timeseries_', MetStations_Precip[[i]]$id[1], '.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = MetStations_Precip[[i]]$prcp/10, x = as.Date(MetStations_Precip[[i]]$date), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = 'Precipitation (mm)', 
       main = paste0('Station #', MetStations_Precip[[i]]$id[1]), cex.lab = 1.5, cex.axis = 1.5)
  dev.off()
  
  #With map and timeseries
  png(paste0('Precip_Timeseries_Map', MetStations_Precip[[i]]$id[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
  layout(rbind(c(1,2)))
  
  plot(y = MetStations_Precip[[i]]$prcp/10, x = as.Date(MetStations_Precip[[i]]$date), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = 'Precipitation (mm)', 
       main = paste0('Station #', MetStations_Precip[[i]]$id[1]), cex.lab = 1.5, cex.axis = 1.5)
  
  #map
  plot(NOAAstations_locs, col = 'white')
  plot(ROI, add = T)
  #All temperature gauges
  plot(NOAAstations_locs_Precip, pch = 16, col = 'black', add = TRUE)
  #Gauge selected for timeseries plot
  plot(NOAAstations_locs_Precip[i,], pch = 16, col = 'blue', add = TRUE)
  
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4343000, len = 700, col = 'black', lab = 'N')
  legend('bottomleft', title = 'NOAA Precipitation Stations', legend = c('Selected', 'Other'), pch = 16, col = c('blue', 'black'), bty = 'n')
  title(main = paste0('Station #', MetStations_Precip[[i]]$id[1]))
  dev.off()
}
rm(i)

#  Max and Min Temperature----
for (i in 1:length(MetStations_TMAX)){
  png(paste0('Temp_Timeseries_', MetStations_TMAX[[i]]$id[1], '.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = MetStations_TMAX[[i]]$tmax/10, x = as.Date(MetStations_TMAX[[i]]$date), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = expression(paste('Temperature (', degree, 'C)')), 
       main = paste0('Station #', MetStations_TMAX[[i]]$id[1]), col = 'red', ylim = ylim_Temp, cex.lab = 1.5, cex.axis = 1.5)
  par(new = TRUE)
  plot(y = MetStations_TMIN[[i]]$tmin/10, x = as.Date(MetStations_TMIN[[i]]$date), type = 'o', pch = 16, cex = 0.3,
       xlab = '', ylab = '', axes = FALSE, 
       main = paste0('Station #', MetStations_TMIN[[i]]$id[1]), col = 'blue', ylim = ylim_Temp)
  dev.off()
  
  #With map and timeseries
  png(paste0('Temp_Timeseries_Map', MetStations_TMAX[[i]]$id[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
  layout(rbind(c(1,2)))
  
  plot(y = MetStations_TMAX[[i]]$tmax/10, x = as.Date(MetStations_TMAX[[i]]$date), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = expression(paste('Temperature (', degree, 'C)')), 
       main = paste0('Station #', MetStations_TMAX[[i]]$id[1]), col = 'red', ylim = ylim_Temp, cex.lab = 1.5, cex.axis = 1.5)
  par(new = TRUE)
  plot(y = MetStations_TMIN[[i]]$tmin/10, x = as.Date(MetStations_TMIN[[i]]$date), type = 'o', pch = 16, cex = 0.3,
       xlab = '', ylab = '', axes = FALSE, 
       main = paste0('Station #', MetStations_TMIN[[i]]$id[1]), col = 'blue', ylim = ylim_Temp)
  
  #map
  plot(NOAAstations_locs, col = 'white')
  plot(ROI, add = T)
  #All temperature gauges
  plot(NOAAstations_locs_TMAX, pch = 16, col = 'black', add = TRUE)
  #Gauge selected for timeseries plot
  plot(NOAAstations_locs_TMAX[i,], pch = 16, col = 'purple', add = TRUE)
  
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4343000, len = 700, col = 'black', lab = 'N')
  legend('bottomleft', title = 'Temperature Stations', legend = c('Selected', 'Other'), pch = 16, col = c('purple', 'black'), bty = 'n')
  title(main = paste0('Station #', MetStations_TMAX[[i]]$id[1]))
  dev.off()
}
rm(i)

# Fixme: ACF for precip and temp datasets----
# Fixme: evaluate outliers for precip and temp (possible multivariate outliers) - make scatterplots of precip and temp----
# Fixme: Spatial predicton of precip and temperature----
# Write met station data to files----
setwd(dir = wd_NOAA)
writeOGR(NOAAstations_locs, dsn = getwd(), layer = f_NOAAstationsROI, driver = "ESRI Shapefile")
list.save(x = MetStations, file = f_NOAAstationsDataList, type = "YAML")

# DEM for Weather Gauges----
#  Plot the DEM of the station vs. the DEM of the grid cell----
DEM = processDEM(dir_DEMs = dir_DEM, f_DEMs = f_DEM, pCRS = pCRS)

NOAAstations_locs$DEMelev = raster::extract(DEM, y = NOAAstations_locs)
png('CompareNOAAGaugeElevToDEM.png', res = 300, units = 'in', width = 5, height = 5)
plot(NOAAstations_locs$elevation, NOAAstations_locs$DEMelev,
     xlab = 'Reported Elevation (m)', ylab = 'DEM Elevation (m)', main = 'NOAA Gauge Elevations')
lines(c(-100,1100), c(-100,1100), col = 'red')
dev.off()

# Write DEM + met station data to files----
setwd(dir = wd_NOAA)
writeOGR(NOAAstations_locs, dsn = getwd(), layer = paste0(f_NOAAstationsROI, '_withDEM'), driver = "ESRI Shapefile")

# #Other Weather Station Data Processing - In Development----
# setwd(dir_sfgauges)
# #Load file containing hyperlinks to the climate data
# ClimGauges = read.csv("NOAA_HyperlinksToGauges.csv", stringsAsFactors = FALSE, header = FALSE)
# 
# #Make new directory for weather data
# wd_clim = paste0(getwd(), '\\Weather')
# dir.create(wd_clim)
# setwd(wd_clim)
# 
# #Fixme: some AllStations data are climate stations when loaded from the csv file from the website specified in the readme.
# #Fixme: add DEM elevation to dataset and compare
# 
# #Load from NWIS
# TempCStations = whatNWISsites(statecode = "MD", parameterCd = '00020')
# TempFStations = whatNWISsites(statecode = "MD", parameterCd = '00021')
# TotalPrecipStations = whatNWISsites(statecode = "MD", parameterCd = '00045')
# TempCStations = readNWISsite(TempCStations$site_no)
# TempFStations = readNWISsite(TempFStations$site_no)
# TotalPrecipStations = readNWISsite(TotalPrecipStations$site_no)
# 
# #Make spatial data
# GaugesLocs_NAD27 = TempCStations[which(TempCStations$coord_datum_cd == 'NAD27'), ]
# coordinates(GaugesLocs_NAD27) = c('dec_long_va', 'dec_lat_va')
# proj4string(GaugesLocs_NAD27) = CRS('+init=epsg:4267')
# GaugesLocs_NAD83 = TempCStations[which(TempCStations$coord_datum_cd == 'NAD83'), ]
# coordinates(GaugesLocs_NAD83) = c('dec_long_va', 'dec_lat_va')
# proj4string(GaugesLocs_NAD83) = CRS('+init=epsg:4269')
# #Transform to NAD83 UTM Zone 18N
# GaugeLocs_NAD27 = spTransform(GaugesLocs_NAD27, CRS(pCRS))
# GaugeLocs_NAD83 = spTransform(GaugesLocs_NAD83, CRS(pCRS))
# #Join to one dataset again
# TempCStations = rbind(GaugeLocs_NAD27, GaugeLocs_NAD83)
# #Remove separate datasets
# rm(GaugeLocs_NAD27, GaugeLocs_NAD83, GaugesLocs_NAD27, GaugesLocs_NAD83)
# 
# GaugesLocs_NAD83 = TempFStations[which(TempFStations$coord_datum_cd == 'NAD83'), ]
# coordinates(GaugesLocs_NAD83) = c('dec_long_va', 'dec_lat_va')
# proj4string(GaugesLocs_NAD83) = CRS('+init=epsg:4269')
# #Transform to NAD83 UTM Zone 18N
# GaugeLocs_NAD83 = spTransform(GaugesLocs_NAD83, CRS(pCRS))
# #Join to one dataset again
# TempFStations = GaugeLocs_NAD83
# #Remove separate datasets
# rm(GaugeLocs_NAD83, GaugesLocs_NAD83)
# 
# GaugesLocs_NAD27 = TotalPrecipStations[which(TotalPrecipStations$coord_datum_cd == 'NAD27'), ]
# coordinates(GaugesLocs_NAD27) = c('dec_long_va', 'dec_lat_va')
# proj4string(GaugesLocs_NAD27) = CRS('+init=epsg:4267')
# GaugesLocs_NAD83 = TotalPrecipStations[which(TotalPrecipStations$coord_datum_cd == 'NAD83'), ]
# coordinates(GaugesLocs_NAD83) = c('dec_long_va', 'dec_lat_va')
# proj4string(GaugesLocs_NAD83) = CRS('+init=epsg:4269')
# #Transform to NAD83 UTM Zone 18N
# GaugeLocs_NAD27 = spTransform(GaugesLocs_NAD27, CRS(pCRS))
# GaugeLocs_NAD83 = spTransform(GaugesLocs_NAD83, CRS(pCRS))
# #Join to one dataset again
# TotalPrecipStations = rbind(GaugeLocs_NAD27, GaugeLocs_NAD83)
# #Remove separate datasets
# rm(GaugeLocs_NAD27, GaugeLocs_NAD83, GaugesLocs_NAD27, GaugesLocs_NAD83)
# 
# #Read the hyperlinks and place data into separate text files per gauge
# #List of climate stations named by the station number
# ClimateStationList = list()
# for (i in 1:nrow(ClimGauges)){
#   #First save the file as a *.txt
#   if (is.na(strsplit(strsplit(ClimGauges[i,], split = 'NWIS=')[[1]][2], split = '%', fixed = TRUE)[[1]][1])){
#     #ACIS stations
#     fname = paste0(strsplit(strsplit(ClimGauges[i,], split = 'ACIS=')[[1]][2], split = '%', fixed = TRUE)[[1]][1], '.txt')
#   }else{
#     #NWIS stations
#     fname = paste0(strsplit(strsplit(ClimGauges[i,], split = 'NWIS=')[[1]][2], split = '%', fixed = TRUE)[[1]][1], '.txt')
#   }
#   
#   download.file(url = ClimGauges[i,], destfile = fname, mode = 'wb', method = "curl")
#   ClimStationData = read.table(fname, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
#   #Reformat the time to standard format
#   ClimStationData$time = as.POSIXct(ClimStationData$time, format = '%m/%d/%Y %H:%M')
#   
#   #Place the measurements in chronological order by the sort date and time
#   ClimStationData = ClimStationData[order(as.POSIXct(ClimStationData$time)),]
#   
#   #Add to list
#   ClimateStationList = c(ClimateStationList, list(ClimStationData))
#   names(ClimateStationList)[length(ClimateStationList)] = strsplit(fname, split = '.txt', fixed = TRUE)[[1]][1]
# }
# rm(i, ClimStationData, fname)

#Get the site information for each of those gauges

#Another alternative - not recommended, and not tested.
#Get the FIPS code for Maryland
#MDfips = rnoaa::fipscodes[which(rnoaa::fipscodes$state == 'Maryland'),]
#Download GHCND data - requires access code requested from: https://www.ncdc.noaa.gov/cdo-web/token
# This code takes a while to be emailed to you. AND one can only download a year at a time - not recommended
#NOAAstations = ncdc(datasetid='GHCND', locationid = paste0('FIPS:', MDfips$fips_state), startdate = '1998-01-01', enddate = '2018-12-31')
#For metadata
#rnoaa::homr(state = 'MD', county = 'Baltimore')
#Python alternative: https://k3.cicsnc.org/jared/GHCNpy