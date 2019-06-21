#Script to read USGS streamflow and water quality data from gauges
# Function help at: https://github.com/USGS-R/dataRetrieval 
#             Blog: https://owi.usgs.gov/R/training-curriculum/usgs-packages/
#        Slideshow: https://owi.usgs.gov/R/dataRetrieval.html#1
# Also has a component for downloading weather station data stored on the USGS database (which includes NOAA ACIS data)

#You can find USGS gauges of any type in your region of interest on this website:
# https://cida.usgs.gov/enddat/dataDiscovery.jsp
# Record the coordinates of the bounding box for your region of interest! They could be important for the next step. 
# Select the data you want and copy the resulting table of values into a txt file
#  txt is important to retain leading zeros on gauge numbers
# If you inhereted a csv or txt file without leading zeros, a script below can add them to NWIS gauges.
#  Other gauge datasets seem to not require leading zeros.

#Coordinates of gauges are not reported in that download :( 
# The USGS function readNWISsite() can look up the coordinates, given the gauge numbers.
# You can also find coordinates for your gauges on this site if you provide a bounding box of coordinates:
# https://waterdata.usgs.gov/nwis/inventory?search_criteria=lat_long_bounding_box&submitted_form=introduction
# Obtained coordinates may include altitude of the gauge, which can be important vs. DEM elevation 
#  (e.g. if there's a cliff at the gauge vs. DEM mean elevation of the pixel)
#  Select the altitude features that you want from the scroll list on that website to download them
#  A plot of DEM vs. Gauge reported elevation is made below to visually detect discrepancies
#   (e.g. USGS data reported in m instead of in ft)
#  Ensure that the units and the vertical datum are the same for your DEM and all gauge altitudes
#   DEMs in the US tend to be NAVD88 in m, whereas gauges tend to be referenced to NGVD29 in ft
#   Differences tend to be minor in the US, except in the West: https://www.ngs.noaa.gov/TOOLS/Vertcon/vertcon.html
#  Altitude datum codes: https://help.waterdata.usgs.gov/code/alt_datum_cd_query?fmt=html
#   collection method codes: https://help.waterdata.usgs.gov/code/alt_meth_cd_query?fmt=html

#You can find water quality gauges on this website:
# https://www.waterqualitydata.us/portal/
# Download "site data only" to receive a csv file with gauge/site information with coordinates.
# The MonitoringLocationIdentifier field is used to download data for each gauge/site in the script below.

#You can find data quality codes for USGS datasets here:
# https://help.waterdata.usgs.gov/codes-and-parameters/codes#discharge_cd
# You should always look at these quality codes, and process your data accordingly.
# This script colors streamflow time series by error code, but doesn't process further than that.
# codes for streamflow: https://help.waterdata.usgs.gov/codes-and-parameters/daily-value-qualification-code-dv_rmk_cd
#                  AND: https://help.waterdata.usgs.gov/codes-and-parameters/instantaneous-and-daily-value-status-codes
#                  Yes, there are 2 separate reference schemes for streamflow data.

#Portions of the function for the streamflow download were provided by:
# Caitline Barber (Caitline.Barber@tufts.edu) and Jonathan Lamontagne (Jonathan.Lamontagne@tufts.edu)
# Modified by Jared Smith (js4yd@virginia.edu) in June, 2019, and started git tracking.
# See git commit history for all later contributions

#Set directory names----
#Region of interest shapefile
dir_ROI = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\BES-Watersheds-Land-Cover-Analysis"  
#Color functions
dir_ColFuns = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\USGSGauges"
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
f_StreamGaugeSites = "USGS_GaugeSites.txt"
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
#Color functions for plots (R script from Jared Smith's Geothermal_ESDA Github Repo)
setwd(dir_ColFuns)
source('ColorFunctions.R')

#Streamflow----
setwd(dir_sfgauges)
# Read USGS station data from .csv file----
AllStations <- read.csv(f_StreamGaugeData, stringsAsFactors = FALSE)

#Add a leading 0 to the NWIS gauges to look up their values on the server
# Some of the gauges are not numbers so only add 0 to number gauges
# NOTE: This step may not be necessary for your dataset.
# NOTE: suppressing warnings for NAs introduced by coercion, which is intended.
#       Users should check that other warnings are not also being suppressed.
#Fixme: Function for leading zeros to NWIS gauges
StationStart = substr(AllStations$GaugeNum, start = 1, stop = 1)
for (i in 1:nrow(AllStations)){
  if(AllStations$Source[i] == 'NWIS'){
    if (suppressWarnings(is.na(as.numeric(StationStart[i])))){
      #This station starts with a character. Retain original name
      AllStations$GaugeNum[i] <- AllStations$GaugeNum[i]
    }else{
      #Add a leading 0
      AllStations$GaugeNum[i] <- paste0("0", AllStations$GaugeNum[i])
    }
  }
}
rm(StationStart, i)

# Statistics to report for each streamflow gauge, if available----
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

# Parameters to report for each gauge, if available----
# All codes defined here: https://help.waterdata.usgs.gov/code/parameter_cd_query?fmt=rdb&inline=true&group_cd=%
#Fixme: lookup codes from the website
Par.cfsFlow = "00060"
Par.Nflow = "00600"
Par.Lat = "91110"
Par.Long = "91111"

#Fixme: is there a way to collect all variables that begin Par. and collect them into a new vector?
Pars = c(Par.cfsFlow, Par.Nflow, Par.Long, Par.Lat)

# Collect data for each NWIS gauge----
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
# NOTE: Your data may be all one coordinate system, and therefore not need to split into 2 datasets
#       before joining into 1 dataset.
# NOTE: your coordinate system may be different (epsg code)
# Some of the data are NAD27 projection and others are NAD83 projection. Split the dataset to handle each
#Fixme: function for splitting coordinate systems and returning one same-coordinate system file
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

#  Add DEM elevation in resolution of your choice for modeling----
#   Gather all of the DEMs together and mosaic into one file
#Fixme: function for DEM mosaikinng from separate tiles
for (d in 1:length(dir_DEM)){
  if (d > 1){
    DEM2 = raster(x = paste0(dir_DEM[d], "/", f_DEM[d]))
    DEM = mosaic(DEM, DEM2, fun = mean)
  }else{
    DEM = raster(x = paste0(dir_DEM[d], "/", f_DEM[d]))
  }
}
rm(d, DEM2)
#Project to project coordinate system
DEM = projectRaster(DEM, crs = CRS(pCRS))
#Add the DEM elevation to the gauge dataset
elev = extract(x = DEM, y = GaugeLocs)
GaugeLocs$ElevDEM = elev
rm(elev)

#   Compare the DEM elevation to the listed elevation----
png('CompareGaugeElev.png', res = 300, units = 'in', width = 5, height = 5)
plot(GaugeLocs$alt_va[is.na(GaugeLocs$alt_va) == FALSE], GaugeLocs$ElevDEM[is.na(GaugeLocs$alt_va) == FALSE]/.3048,
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

#  Add coordinates and elevation data to the NWISstations data----
for (i in 1:nrow(NWISstations)){
  #NOTE: if your data are not both character, errors stating the following will appear:
  # Error in data.frame(..., check.names = FALSE) : 
  #  arguments imply differing number of rows: 1, 0
  Ind = which(GaugeLocs$site_no == NWISstations$GaugeNum[i])
  if (length(Ind) > 1){
    print('More than one gauge number matches the uniqueNum gauge ', i, '. Using only first match.')
  }
  test = cbind(NWISstations[i,], GaugeLocs@data[Ind,], GaugeLocs@coords[Ind,][1], GaugeLocs@coords[Ind,][2])
  colnames(test) = c(colnames(NWISstations), colnames(GaugeLocs@data), colnames(GaugeLocs@coords))
  if (i == 1){ 
    NewData = test
  }else{
    NewData = rbind(NewData, test)
  }
}
NWISstations = NewData
rm(i, Ind, NewData, test)

#Make NWIS stations a spatial dataframe
coordinates(NWISstations) = c('dec_long_va', 'dec_lat_va')
proj4string(NWISstations) = CRS(pCRS)

#Clip the NWIS gauges to the region of interest (shapefile) 
ROI = readOGR(dsn = dir_ROI, layer = f_ROI, stringsAsFactors = FALSE)
ROI = spTransform(ROI, CRS(pCRS))
#Clip the streamflow gauges to the region of interest
NWIS_ROI = NWISstations[ROI,]

# Plot locations of NWIS gauges----
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

# Compare reported vs. DEM elevation of gauges within the ROI----
png('CompareGaugeElev.png', res = 300, units = 'in', width = 5, height = 5)
plot(NWIS_ROI$alt_va, NWIS_ROI$ElevDEM/.3048,
     xlab = 'USGS Reported Elevation (ft)', ylab = 'DEM Elevation (ft)', main = 'Gauge Elevations')
lines(c(-100,1100), c(-100,1100), col = 'red')
dev.off()

#One of these gauge elevations is a lot lower than DEM. Likely that the gauge was reported in m in USGS database
#identify(NWIS_ROI$alt_va, NWIS_ROI$ElevDEM/.3048)

# Download the within-ROI stream gauge data in parallel----
#  Use only the unique gauge numbers in the dataset 
#  (repeats occur when multiple variables are available for a gauge)
uniqueNums = unique(NWIS_ROI$GaugeNum)
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

# Gather the records for each gauge into a list of dataframes----
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

# Identify missing data and add the total number of missing days to the NWIS_ROI file----
#Fixme: make this a function
NWIS_ROI$MissingData = NA
for (i in 1:length(StreamStationList)){
  #Missing data that are reported as NA cells
  NWIS_ROI$MissingData[which(as.numeric(NWIS_ROI$GaugeNum) == as.numeric(StreamStationList[[i]]$site_no[1]))] = length(which(is.na(StreamStationList[[i]]$X_00060_00003)))
  #Add to those the missing data resulting from gaps > 1 day in the record
  # NOTE: this assumes data are daily streamflow
  gaps = c(with(data = StreamStationList[[i]], as.numeric(Date[-1]) - as.numeric(Date[-nrow(StreamStationList[[i]])])))
  
  NWIS_ROI$MissingData[which(as.numeric(NWIS_ROI$GaugeNum) == as.numeric(StreamStationList[[i]]$site_no[1]))] = NWIS_ROI$MissingData[which(as.numeric(NWIS_ROI$GaugeNum) == as.numeric(StreamStationList[[i]]$site_no[1]))] + sum(gaps[which(gaps > 1)])
  
  #Fill in the missing data dates with NA values to have a complete time series for all records
  Inds = which(gaps > 1)
  if (length(Inds) > 0){
    for (j in 1:length(Inds)){
      NumNAs = as.numeric(StreamStationList[[i]][(Inds[j]+1),]$Date - StreamStationList[[i]][Inds[j],]$Date) - 1
      DateNAs = seq(StreamStationList[[i]][Inds[j],]$Date+1, StreamStationList[[i]][Inds[j]+1,]$Date-1, 1)
      #Assign NumNAs new rows to this dataframe with NA streamflow
      for (k in 1:NumNAs){
        r = cbind(StreamStationList[[i]][1,1:2], DateNAs[k], NA, NA)
        colnames(r) = colnames(StreamStationList[[i]]) 
        StreamStationList[[i]] = rbind(StreamStationList[[i]], r)
      }
    }
  }
  #Sort the streamflow series by date
  StreamStationList[[i]] = StreamStationList[[i]][order(StreamStationList[[i]]$Date),]
}
rm(gaps, Inds, NumNAs, DateNAs, r, i, j, k)

#Fixme: Missing data fill in with numerical value estimates using prediction in ungauged basins methods for large gaps
#Fixme: check for high and low flow outliers in each record, and compare spatially to other gauges on those dates
# Can include both FFA and daily flow outlier analysis

# Make a map of points colored by their record lengths, corrected for the total amount of missing data----
NWIS_ROI$RecordLength = NWIS_ROI$RecordLengthMinusGaps = NA
for (i in 1:length(StreamStationList)){
  NWIS_ROI$RecordLength[which(as.numeric(NWIS_ROI$GaugeNum) == as.numeric(StreamStationList[[i]]$site_no[1]))] = as.numeric((max(StreamStationList[[i]]$Date) - min(StreamStationList[[i]]$Date)))
}
rm(i)
NWIS_ROI$RecordLengthMinusGaps = NWIS_ROI$RecordLength - NWIS_ROI$MissingData
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
legend('right', title = 'Streamflow Station \n Record Lengths \n (years)', legend = seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy), pch = 16, col = colFun(seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy)), bty = 'n')
dev.off()

# Write streamflow datasets to files----
setwd(dir_DEM_out)
writeRaster(x = DEM, filename = f_DEM_mosiac, format = "GTiff")
setwd(dir_sfgauges)
writeOGR(obj = NWIS_ROI, dsn = getwd(), driver = 'ESRI Shapefile', layer = f_NWIS_ROI_out)
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

#Find sites that have any N and P water quality data in Maryland
#Phosphorous
phosSites <- whatWQPsites(statecode="MD", characteristicName="Phosphorus")
#Nitrogen
NitroSites <- whatWQPsites(statecode="MD", characteristicName="Nitrogen")

#Select only those sites that have nitrogen data
WQstations_ROI_N = WQstations_ROI[WQstations_ROI$MonitoringLocationIdentifier %in% NitroSites$MonitoringLocationIdentifier,]
#Select only those sites that have phosphorous data
WQstations_ROI_P = WQstations_ROI[WQstations_ROI$MonitoringLocationIdentifier %in% phosSites$MonitoringLocationIdentifier,]

#Use only the unique gauge numbers in the dataset (repeats occur when multiple variables are available for a gauge)
uniqueWQNums_N = unique(WQstations_ROI_N$MonitoringLocationIdentifier)
uniqueWQNums_P = unique(WQstations_ROI_P$MonitoringLocationIdentifier)

#Run the downloads in parallel.
# NOTE: These downloads occasionally fail when run in parallel and return internal server errors.
# If that happens to you, try running in serial and see if you still get the errors.
# I'm not sure what to do if you still get them. Running in serial has worked for me.
#Fixme: which function is better? This also exists: readNWISqw
cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
n = foreach(i = uniqueWQNums_N, .packages = 'dataRetrieval') %dopar% {
  # Read all of the data for the provided station number
  stationData <- readWQPqw(siteNumbers = i, parameterCd = "")
  write.table(stationData, 
              paste0(getwd(), '/Nitrogen_', i, ".txt"), 
              sep = "\t")
}
p = foreach(i = uniqueWQNums_P, .packages = 'dataRetrieval') %dopar% {
  # Read all of the data for the provided station number
  stationData <- readWQPqw(siteNumbers = i, parameterCd = "")
  write.table(stationData, 
              paste0(getwd(), '/Phosphorous_', i, ".txt"), 
              sep = "\t")
}
stopCluster(cl)
rm(cl, stationData)
#Check that the run was successful
if (any(!is.null(unlist(n)))){
  print('NITROGEN DOWNLOAD UNSUCCESSFUL')
}else{
  print('Nitrogen gauge data download complete!')
  rm(n)
}
if (any(!is.null(unlist(p)))){
  print('PHOSPHOROUS DOWNLOAD UNSUCCESSFUL')
}else{
  print('Phosphorous gauge data download complete!')
  rm(p)
}

# Process Nitrogen data----
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

#Filter the CharacteristicName for nitrogen measurements only. Then look at the ResultSampleFractionText
#For now, taking only total nitrogen, but many sites have other nitrogen measurements
#Fixme: ideally, this would cycle through all of the unique combinations of ResultSampleFractionText and CharacteristicName and return separate lists for each variable
for (i in 1:length(NitroStationList)){
  #Find all of the "Total" or "total" in the ResultSampleFractionText field
  T_Ind = grep(x = NitroStationList[[i]][,"ResultSampleFractionText"], pattern = 'Total', ignore.case = TRUE)
  #There are some datasets that have multiple total nitrogen CharacteristcName.
  #Find all and make separate time series for each
  N_names = unique(NitroStationList[[i]][T_Ind,][grep(x = NitroStationList[[i]][T_Ind,"CharacteristicName"], pattern = 'nitrogen', ignore.case = TRUE), "CharacteristicName"])
  for (n in N_names){
    N_Ind = which(NitroStationList[[i]][T_Ind,"CharacteristicName"] == n)
    write.table(NitroStationList[[i]][T_Ind,][N_Ind,], 
                paste0(getwd(), '/Nitrogen_', NitroStationList[[i]]$MonitoringLocationIdentifier[1], "_", gsub(pattern = " ", replacement = "", x = n, fixed = TRUE),".txt"), 
                sep = "\t")
  }
}

#Need to check that the units are all the same for each measurement type in each record and across records
#ResultMeasure.MeasureUnitCode
#ResultMeasureValue

#Check any detection limits in DetectionQuantitationLimitMeasure.MeasureValue

#Plot the time series for each gauge, and the eCDF, colored by error code
# Fixme: make the legend only include the error codes that are in the station's dataset
for (i in 1:length(StreamStationList)){
  #Assign colors to the error codes
  colCodes = rainbow(length(ErrCodes))
  
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
  legend('topleft', legend = ErrCodes, col = colCodes, pch = 16)
  dev.off()
  
  png(paste0('StreamflowExceedance_', StreamStationList[[i]]$site_no[1],'.png'), res = 300, units = 'in', width = 6, height = 6)
  qqnorm(StreamStationList[[i]]$X_00060_00003, pch = 1, 
         ylab = 'Daily Mean Streamflow (cfs)', main = paste0('Non-Exceedance Probability for Daily Mean Streamflow \n Station #', StreamStationList[[i]]$site_no[1]))
  dev.off()
}
rm(i, c, colCodes)

#Missing data

#Evaluate nondetects in result_cd

#Search for flow-normalized outliers

# Process Phosphorous Data----

# BES Water Quality Gauge Data----
setwd("C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\WaterChemistry")
#BES Water quality sample time series
BES_WQ = read.csv('BES-stream-chemistry-data-for-WWW-feb-2018---core-sites-only_JDSprocessed.csv', stringsAsFactors = FALSE)
#BES gauge number reference table
USGS_GaugeMatch = read.csv('Abbreviations_SampleRecordLengths.csv', stringsAsFactors = FALSE)
#Remove spaces from the names of the sites
USGS_GaugeMatch$Abbreviation = gsub(pattern = ' ', replacement = '', x = USGS_GaugeMatch$Abbreviation, fixed = TRUE)
#add leading zeros to gauge numbers
USGS_GaugeMatch$USGSGaugeNum = paste0("0", USGS_GaugeMatch$USGSGaugeNum)
BES_WQ$USGSgauge = NA
USGSnums = unique(USGS_GaugeMatch$USGSGaugeNum[-grep(pattern = 'NA', x = USGS_GaugeMatch$USGSGaugeNum)])

#Filter into separate databases by site and add gauge numbers to database
BES_WQ_Sites = list()
uniqueSites = unique(BES_WQ$Site)
for (i in 1:length(uniqueSites)){
  f = BES_WQ[BES_WQ$Site == uniqueSites[i],]
  if(f$Site[1] %in% USGS_GaugeMatch$Abbreviation){
    f$USGSgauge = USGS_GaugeMatch$USGSGaugeNum[USGS_GaugeMatch$Abbreviation == f$Site[1]]
  }
  BES_WQ_Sites = c(BES_WQ_Sites, list(f))
}

#Plot time series for each of the sites
#Nitrogen
png('NitrogenBESgauges.png', res = 300, height = 6, width = 6, units = 'in')
plot(BES_WQ_Sites[[1]]$Date, BES_WQ_Sites[[1]]$TN..mg.N.L.)
dev.off()

#Associate each of the sites to a station, if available

#Map the station locations and plot them on a map showing BES data and WQP data


#Weather Stations----
#Fixme: some AllStations data are climate stations
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