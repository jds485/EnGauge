#Loading data from the output of Method1Example to use with additional datasets acquired for a specific site

#Fixme: add satellite background images to maps

#Set directory names----
#EnGauge repository
dir_EnGauge = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\EnGauge\\EnGauge"
#Color functions - from JDS github repo: Geothermal_ESDA
dir_ColFuns = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\USGSGauges"
#Region of interest shapefile
dir_ROI = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\BES-Watersheds-Land-Cover-Analysis"  
#Watersheds shapefiles
dir_Sheds = "C:\\Users\\js4yd\\Documents\\DEMtest"
#Color functions - from JDS github repo: Geothermal_ESDA
dir_ColFuns = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\USGSGauges"
#DEM
dir_DEM_out = 'C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\DEM\\'
#Streamflow gauges
wd_sf = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\USGSGauges\\Streamflow"
#Nitrogen
wd_N = 'C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\USGSGauges\\Nitrogen'
#Phosphorus
wd_P = 'C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\USGSGauges\\Phosphorus'
#Water Chemistry datasets
dir_WChem = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\WaterChemistry"
#Synoptic Water Chemistry
dir_SynWChem_Kenworth = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\WaterChemistry\\BARN_Synoptic\\WaterChemical_Kenworth_01-02"
dir_SynWChem_Smith = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\WaterChemistry\\BARN_Synoptic\\WaterChemical_Smith_06-07"
#WRTDS Output
wd_WRTDS_PondBranch = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\RHESSysFilePreparation\\obs\\POBR\\"
#Precipitation
dir_precip = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\Precipitation"
dir_Nexrad = paste0(dir_precip, '/', "BES_Nexrad")
#Streams
dir_streams = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\NHD_H_Maryland_State_Shape\\Shape"

#Set input filenames----
#Region of interest shapefile name
f_ROI = "BaismanGFMerge2"
#Watersheds
f_DR = 'DRFrank9'
f_BR = 'BaisRun5'
f_SL = 'ScottsLevel7'
#DEM
f_DEM_mosiac = 'DEM_mosaic.tif'
#NWIS streamflow gauges in ROI
f_NWIS_ROI = 'NWIS_ROI'
#NWIS Streamflow timeseries list
f_sflist = 'SF.yaml'
#Nitrogen timeseries list
f_Nlist = 'TN.yaml'
f_Nalist = 'TN_a.yaml'
f_Nmlist = 'TN_m.yaml'
f_Ndlist = 'TN_d.yaml'
#Nitrogen sites shapefile
f_WQNsites = 'NitrogenSites'
#Phosphorus timeseries list
f_Plist = 'TP.yaml'
f_Palist = 'TP_a.yaml'
f_Pmlist = 'TP_m.yaml'
f_Pdlist = 'TP_d.yaml'
#Phosphorus sites shapefile
f_WQPsites = 'PhosphorusSites'
#Water Chemistry from BES stations
f_WChem = 'BES-stream-chemistry-data-for-WWW-feb-2018---core-sites-only_JDSprocessed.csv'
#USGS stream gauge matches
f_USGS_GaugeMatch = 'Abbreviations_SampleRecordLengths.csv'
#Precipitation Measurements
f_Precip = 'BES_precipitation.csv'
f_Precip_locs = 'BES_GaugeLocs.csv'
#Stream shapefile
f_streams = 'NHDFlowline'

#Set output filenames----
#all chemical data as lists by site
f_BESallchem = 'BES_AllChem.yaml'
#BES N, P
# Locations
f_Nsites = 'BES_NitrogenSites'
f_Psites = 'BES_PhosphorusSites'
# all chemical data at sites with TN data
f_BES_TNallchem = 'BES_TN_AllChem.yaml'
# BES TN aggregated to daily, monthly, yearly info
f_BES_TNd = 'BES_TN_d.yaml'
f_BES_TNm = 'BES_TN_m.yaml'
f_BES_TNa = 'BES_TN_a.yaml'
# all chemical data at sites with TP data
f_BES_TPallchem = 'BES_TP_AllChem.yaml'
# BES TP aggregated to daily, monthly, yearly info
f_BES_TPd = 'BES_TP_d.yaml'
f_BES_TPm = 'BES_TP_m.yaml'
f_BES_TPa = 'BES_TP_a.yaml'

#Precipitation from BES
f_BESPrecipList = 'BES_PrecipList_AllGauges.yaml'
# Aggregated to annual, monthly, daily
f_BESPrecipAggYear = 'BES_Precip_a_AllGauges.yaml'
f_BESPrecipAggMonth = 'BES_Precip_m_AllGauges.yaml'
f_BESPrecipAggDay = 'BES_Precip_d_AllGauges.yaml'
# Daily + site average of 2 gauges
f_BESPrecipAggDay_SiteAverage = 'BES_Precip_d_Sites.yaml'
# Locations
f_BESprecipSites = 'BES_PrecipSites'
#NOAA station data - all
f_NOAAstationsROI = 'NOAA_StationLocs_20km'
f_NOAAstationsDataList = 'NOAA_MetStations.yaml'
#Nexrad Data
f_NexradDataList = 'NexradYearList.yaml'

#Fixme: provide output filenames for these weather station gauge datasets
# Precipitation

# Temperature

#Set project coordinate system----
#This is the coordinate system that all data will be plotted and written in
# It is not the coordinate system of your data (although it could be)
# EPSG codes from: https://spatialreference.org/ref/?page=2
pCRS = '+init=epsg:26918'

#Define a buffer to use for the ROI to download weather gauges, in map units
ROIbuffPrecip = 20000

#Load libraries----
library(dataRetrieval)
library(sp)
library(raster)
library(rgdal)
library(rlist)
#library(nhdR)
library(stringr)
library(stringi)
library(GISTools)
library(foreach)
library(doParallel)
library(lubridate)
library(zoo)
library(hydroTSM)
library(rnoaa)
library(rappdirs)
library(EGRET)
library(survival)
library(pracma)
library(psych)
library(car)
#Functions from repository
setwd(dir_EnGauge)
source('missingDates.R')
source('checkDuplicates.R')
source('checkZerosNegs.R')
source('formatMonthlyMatrix.R')
source('matplotDates.R')
source('aggregateTimeseries.R')
source('scatterHistCols.R')
#Color functions for plots (R script from Jared Smith's Geothermal_ESDA Github Repo)
setwd(dir_ColFuns)
source('ColorFunctions.R')
#Load modified WRTDS functions and the interpolation tables for regression parameters
setwd('C:\\Users\\js4yd\\OneDrive - University of Virginia\\RHESSys_ParameterSA')
source('WRTDS_modifiedFunctions.R')

#Load information from EnGauge downloads and place into the project coordinate system----
DEM = raster(x = paste0(dir_DEM_out, f_DEM_mosiac))
DEM = projectRaster(DEM, crs = CRS(pCRS))
#Region of interest
ROI = readOGR(dsn = dir_ROI, layer = f_ROI, stringsAsFactors = FALSE)
ROI = spTransform(ROI, CRS(pCRS))
#Streamflow stations in the ROI from NWIS
NWIS_ROI_fn = readOGR(dsn = wd_sf, layer = f_NWIS_ROI, stringsAsFactors = FALSE)
NWIS_ROI_fn = spTransform(NWIS_ROI_fn, CRS(pCRS))
StreamStationList = list.load(file = paste0(wd_sf, '/', f_sflist))
#Total nitrogen stations in ROI from the WQP
WQstations_ROI_N = readOGR(dsn = wd_N, layer = f_WQNsites, stringsAsFactors = FALSE)
WQstations_ROI_N = spTransform(WQstations_ROI_N, CRS(pCRS))
TN = list.load(file = paste0(wd_N, '\\', f_Nlist))
TN_a = list.load(file = paste0(wd_N, '\\', f_Nalist))
TN_d = list.load(file = paste0(wd_N, '\\', f_Ndlist))
TN_m = list.load(file = paste0(wd_N, '\\', f_Nmlist))
#Total phosphorus stations in ROI from the WQP
WQstations_ROI_P = readOGR(dsn = wd_P, layer = f_WQPsites, stringsAsFactors = FALSE)
WQstations_ROI_P = spTransform(WQstations_ROI_P, CRS(pCRS))
TP = list.load(file = paste0(wd_P, '\\', f_Plist))
TP_a = list.load(file = paste0(wd_P, '\\', f_Palist))
TP_d = list.load(file = paste0(wd_P, '\\', f_Pdlist))
TP_m = list.load(file = paste0(wd_P, '\\', f_Pmlist))

#Load watersheds----
setwd(dir_Sheds)
DeadRun_Outlet = readOGR(dsn = getwd(), layer = f_DR, stringsAsFactors = FALSE)
BaisRun_Outlet = readOGR(dsn = getwd(), layer = f_BR, stringsAsFactors = FALSE)
ScottsLevel = readOGR(dsn = getwd(), layer = f_SL, stringsAsFactors = FALSE)
DeadRun_Outlet = spTransform(DeadRun_Outlet, CRS(pCRS))
BaisRun_Outlet = spTransform(BaisRun_Outlet, CRS(pCRS))
ScottsLevel = spTransform(ScottsLevel, CRS(pCRS))

#Load streams----
#Fixme: These files are not loading as spatial files when they're downloaded
#Download streams for the state of Maryland
#nhd_get(state = 'MD')
#Load streams into session
#MDstreams = nhd_load(state = 'MD', dsn = 'NHDFlowlineVAA', file_ext = 'shp')
setwd(dir_streams)
MDstreams = readOGR(dsn = getwd(), layer = f_streams, stringsAsFactors = FALSE)
MDstreams = spTransform(MDstreams, CRS(pCRS))
#Clip streams to the ROI
MDstreams = MDstreams[ROI,]

#BES Water Quality Gauge Data----
setwd(dir_WChem)
#BES Water quality time series
BES_WQ = read.csv(f_WChem, stringsAsFactors = FALSE)
#BES gauge number reference table
USGS_GaugeMatch = read.csv(f_USGS_GaugeMatch, stringsAsFactors = FALSE)
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
#corresponding gauge numbers
BES_WQ_GaugeNums = vector('character', length = length(uniqueSites))
for (i in 1:length(uniqueSites)){
  f = BES_WQ[which(BES_WQ$Site == uniqueSites[i]),]
  #Assign gauge number if it has one
  if(f$Site[1] %in% USGS_GaugeMatch$Abbreviation){
    f$USGSgauge = USGS_GaugeMatch$USGSGaugeNum[which(USGS_GaugeMatch$Abbreviation == f$Site[1])]
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
  rm(y)
  
  #Remove any NA dates. These entries are useless.
  f = f[!is.na(f$Dated),]
  
  #Make a POSIX date-time entry
  #Check for errors in the : portion of the time field
  if (length(grep(f$time, pattern = ';', fixed = TRUE)) > 0){
    #semicolon needs to be replaced with colon
    f$time[grep(f$time, pattern = ';', fixed = TRUE)] = str_replace(string = f$time[grep(f$time, pattern = ';', fixed = TRUE)], pattern = ';', replacement = ':')
  }
  if (length(grep(f$time, pattern = '.', fixed = TRUE)) > 0){
    #period needs to be replaced with colon
    f$time[grep(f$time, pattern = '.', fixed = TRUE)] = str_replace(string = f$time[grep(f$time, pattern = '.', fixed = TRUE)], pattern = fixed("."), replacement = ':')
  }
  
  f$DateTime = NA
  f$time[f$time == ''] = NA
  if(any(!is.na(f$time))){
    f$DateTime[!is.na(f$time)] = paste(as.POSIXct(paste0(f$Dated[!is.na(f$time)], ' ', f$time[!is.na(f$time)])))
  }
  
  #Make a new column for the sort date.
  f$SortDate = as.Date(f$Dated)
  f$SortDateTime = as.POSIXct(f$DateTime)
  
  #Place the measurements in chronological order by the sort date and time
  # First order by time (which may have some NA values), then order by date (no NA values)
  f = f[order(as.POSIXct(f$SortDateTime)),][order(f[order(as.POSIXct(f$SortDateTime)),]$SortDate),]
  
  BES_WQ_Sites = c(BES_WQ_Sites, list(f))
  names(BES_WQ_Sites) = c(names(BES_WQ_Sites)[-length(names(BES_WQ_Sites))], f$Site[1])
  BES_WQ_GaugeNums[i] = f$USGSgauge[1]
}
rm(i, f, txtyr, curr2DigYear, BES_WQ)

#  Check for duplicate records----
BES_WQ_Sites = checkDuplicatesAndRemove(StationList = BES_WQ_Sites, colNames = c('SortDate', 'SortDateTime', 'TN..mg.N.L.', 'TP..ugP.L.'))

#  Check for zeros and negative values----
#TN
BES_WQ_Sites_TN = checkZerosNegs(StationList = BES_WQ_Sites, Var = 'TN..mg.N.L.', NegReplace = NA, ZeroReplace = 0)

#TP
BES_WQ_Sites_TP = checkZerosNegs(StationList = BES_WQ_Sites, Var = 'TP..ugP.L.', NegReplace = NA, ZeroReplace = 0)

#  Check for detection limits----
#NA values are gaps in sampling.
#Real detection limits must be discovered for each site by plotting the data
#Check for detection limits in nitrate, phosphorus, TN, TP. Detection limits in NO3 and PO4 would result in the TN and TP values being censored.

#Looking that the timeseries and histograms for each gauge, it is possible to see likely detection limits.
#Only interested in correcting sites that have a large proportion of the data as possibly censored. 
#Not able to tell if there are limits for small amounts of possibly censored data.
#Timeseries useful if limits changed over time.
i = 5
hist(log10(BES_WQ_Sites_TN[[i]]$NO3..mg.N.L.), breaks = 10000)
plot(BES_WQ_Sites_TN[[i]]$SortDateTime, BES_WQ_Sites_TN[[i]]$NO3..mg.N.L., log = 'y')
hist(log10(BES_WQ_Sites_TN[[i]]$TN..mg.N.L.), breaks = 10000)
plot(BES_WQ_Sites_TN[[i]]$SortDateTime, BES_WQ_Sites_TN[[i]]$TN..mg.N.L., log = 'y', ylim = c(0.01,5), xlim = c(as.POSIXct('2000-01-01'), as.POSIXct('2020-01-01')))

#NoTNDL = c(1, 2, 3, 4, 6, 7, 8)
#TNDL = c(5)

#Add detection limit info to Pond Branch
BES_WQ_Sites_TN$POBR$DetectionLimit = NA
#All that are 0.01 are detection limits
BES_WQ_Sites_TN$POBR$DetectionLimit[BES_WQ_Sites_TN$POBR$TN..mg.N.L. == 0.01] = 1
#All that are 0.05 after 2012 are detection limits
BES_WQ_Sites_TN$POBR$DetectionLimit[which((BES_WQ_Sites_TN$POBR$TN..mg.N.L. == 0.05) & (BES_WQ_Sites_TN$POBR$SortDate > as.Date('2012-01-01')))] = 2

i = 5
plot(BES_WQ_Sites_TN[[i]]$SortDateTime, BES_WQ_Sites_TN[[i]]$TN..mg.N.L., log = 'y', ylim = c(0.01,5), xlim = c(as.POSIXct('2000-01-01'), as.POSIXct('2020-01-01')))
par(new = T)
plot(BES_WQ_Sites_TN[[i]]$SortDateTime[is.na(BES_WQ_Sites_TN[[i]]$DetectionLimit) == FALSE], BES_WQ_Sites_TN[[i]]$TN..mg.N.L.[is.na(BES_WQ_Sites_TN[[i]]$DetectionLimit) == FALSE], log = 'y', ylim = c(0.01,5), xlim = c(as.POSIXct('2000-01-01'), as.POSIXct('2020-01-01')), col = 'red')

# i = 7
# hist(log10(BES_WQ_Sites_TP[[i]]$PO4..ug.P.L.), breaks = 10000)
# plot(BES_WQ_Sites_TP[[i]]$SortDateTime, BES_WQ_Sites_TP[[i]]$PO4..ug.P.L., log = 'y')
# hist(log10(BES_WQ_Sites_TP[[i]]$TP..ugP.L.), breaks = 10000)
# plot(BES_WQ_Sites_TP[[i]]$SortDateTime, BES_WQ_Sites_TP[[i]]$TP..ugP.L., log = 'y')

#Looks like detection limits changed to significantly higher values in the late 00's, early 10's
#Then, after 2015 all values are notably higher than previous set - calibration issues??
#DL = c(1, 2, 3, 4, 5, 6, 7, 8)


#  Coordinates of BES sites----
for (i in 1:length(BES_WQ_Sites)){
  if (!is.na(BES_WQ_Sites[[i]]$USGSgauge[1])){
    if(!exists('BES_WQ_Sites_locs')){
      BES_WQ_Sites_locs = readNWISsite(BES_WQ_Sites[[i]]$USGSgauge[1])
      #Add a uniqueID from the original dataset. Multiple WQ samples are taken at each USGSgauge site
      BES_WQ_Sites_locs = cbind(BES_WQ_Sites_locs, Site = BES_WQ_Sites[[i]]$Site[1], stringsAsFactors = FALSE)
    }else{
      BES_WQ_Sites_locs = rbind(BES_WQ_Sites_locs, cbind(readNWISsite(BES_WQ_Sites[[i]]$USGSgauge[1]), Site = BES_WQ_Sites[[i]]$Site[1], stringsAsFactors = FALSE))
    }
  }
}
rm(i)

#Make spatial dataframe
coordinates(BES_WQ_Sites_locs) = c('dec_long_va', 'dec_lat_va')
proj4string(BES_WQ_Sites_locs) = CRS('+init=epsg:4269')
BES_WQ_Sites_locs = spTransform(BES_WQ_Sites_locs, CRS(pCRS))

#Add zeros and negatives to spatial dataset
BES_WQ_Sites_locs_TN = addZerosToSpatialDataset(StationList = BES_WQ_Sites_TN, SpatialDataset = BES_WQ_Sites_locs, site_D = 'Site', site_SL = 'Site')
BES_WQ_Sites_locs_TN = addNegsToSpatialDataset(StationList = BES_WQ_Sites_TN, SpatialDataset = BES_WQ_Sites_locs_TN, site_D = 'Site', site_SL = 'Site')
BES_WQ_Sites_locs_TP = addZerosToSpatialDataset(StationList = BES_WQ_Sites_TP, SpatialDataset = BES_WQ_Sites_locs, site_D = 'Site', site_SL = 'Site')
BES_WQ_Sites_locs_TP = addNegsToSpatialDataset(StationList = BES_WQ_Sites_TP, SpatialDataset = BES_WQ_Sites_locs_TP, site_D = 'Site', site_SL = 'Site')

#   Map the station locations and plot them on a map showing BES data and WQP data----
png('TNTP_BES+WQPsites.png', res = 300, units = 'in', width = 6, height = 6)
plot(ROI, lwd = 2)
plot(MDstreams, add = TRUE, col = 'skyblue')
plot(DeadRun_Outlet, add = TRUE)
plot(ScottsLevel, add = TRUE)
plot(BaisRun_Outlet, add = TRUE)
plot(WQstations_ROI_N, pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(WQstations_ROI_P, pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(WQstations_ROI_P[WQstations_ROI_N,], pch = 16, add = TRUE, col = 'purple', cex = 0.7)
plot(WQstations_ROI_N[WQstations_ROI_P,], pch = 16, add = TRUE, col = 'purple', cex = 0.7)
#BES sites
plot(BES_WQ_Sites_locs, add = TRUE, pch = 4, col = 'purple')

# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 365000, yb = 4348000, len = 700, col = 'black', lab = 'N')
legend('bottomleft', title = 'Water Quality Types', legend = c('T Nitrogen Only', 'T Phosphorus Only', 'Both'), col = c('red', 'blue', 'purple'), pch = c(16,16,16))
legend(x = 334000, y = 4358000, title = 'Water Quality Sites', legend = c('At Stream Gauge', 'Other Location'), col = 'black', pch = c(4, 16))
dev.off()

#  Plot time series for each of the sites----
wd_BESN = paste0(dir_WChem, '/Nitrogen')
dir.create(path = wd_BESN, showWarnings = FALSE)
setwd(wd_BESN)
#Nitrogen
for (i in 1:length(BES_WQ_Sites_TN)){
  if (any(!is.na(BES_WQ_Sites_TN[[i]]$TN..mg.N.L.))){
    png(paste0('BES_N_Timeseries_', BES_WQ_Sites_TN[[i]]$Site[1], '_',  BES_WQ_Sites_TN[[i]]$USGSgauge[1], '_', i, '.png'), res = 300, units = 'in', width = 6, height = 6)
    plot(y = BES_WQ_Sites_TN[[i]]$TN..mg.N.L., x = as.Date(BES_WQ_Sites_TN[[i]]$Dated), type = 'o', pch = 16, cex = 0.3,
         xlab = 'Year', 
         ylab = 'Total Nitrogen (mg N/L)', 
         main = paste0('TN Station #', BES_WQ_Sites_TN[[i]]$USGSgauge[1], ' ', BES_WQ_Sites_TN[[i]]$Site[1]),
         ylim = c(0, max(BES_WQ_Sites_TN[[i]]$TN..mg.N.L., na.rm=TRUE)))
    #xlim = c(as.Date("1980-01-01"), as.Date("2019-06-01")))
    dev.off()
  }
}
rm(i)

#Phosphorus
wd_BESP = paste0(dir_WChem, '/Phosphorus')
dir.create(path = wd_BESP, showWarnings = FALSE)
setwd(wd_BESP)
for (i in 1:length(BES_WQ_Sites_TP)){
  if (any(!is.na(BES_WQ_Sites_TP[[i]]$TP..ugP.L.))){
    png(paste0('BES_P_Timeseries_', BES_WQ_Sites_TP[[i]]$Site[1], '_',  BES_WQ_Sites_TP[[i]]$USGSgauge[1], '_', i, '.png'), res = 300, units = 'in', width = 6, height = 6)
    plot(y = BES_WQ_Sites_TP[[i]]$TP..ugP.L., x = as.Date(BES_WQ_Sites_TP[[i]]$Dated), type = 'o', pch = 16, cex = 0.3,
         xlab = 'Year', 
         ylab = expression(paste('Total Phosphorus (', mu, 'g P/L)')), 
         main = paste0('TP Station #', BES_WQ_Sites_TP[[i]]$USGSgauge[1], ' ', BES_WQ_Sites_TP[[i]]$Site[1]),
         ylim = c(0, max(BES_WQ_Sites_TP[[i]]$TP..ugP.L., na.rm=TRUE)))
    #xlim = c(as.Date("1980-01-01"), as.Date("2019-06-01")))
    dev.off()
  }
}
rm(i)

# For Water Quality at stream gauge sites----
setwd(dir_WChem)
#  Fill in missing dates----
BES_WQ_Sites_TN_Gauged = BES_WQ_Sites_TN[!is.na(BES_WQ_GaugeNums)]
BES_WQ_Sites_TP_Gauged = BES_WQ_Sites_TP[!is.na(BES_WQ_GaugeNums)]
BES_TN = FillMissingDates_par(Dataset = BES_WQ_Sites_locs_TN, StationList = BES_WQ_Sites_TN_Gauged, Var = 'TN..mg.N.L.', Date = 'SortDate', gapType = 'd', site_no_D = 'Site', site_no_SL = 'Site', NoNAcols = c('USGSgauge', 'site_no'))
BES_TP = FillMissingDates_par(Dataset = BES_WQ_Sites_locs_TP, StationList = BES_WQ_Sites_TP_Gauged, Var = 'TP..ugP.L.', Date = 'SortDate', gapType = 'd', site_no_D = 'Site', site_no_SL = 'Site', NoNAcols = c('USGSgauge', 'site_no'))
#Extract data from the function return
BES_WQ_Sites_locs_TN = BES_TN$Dataset
BES_WQ_Sites_locs_TP = BES_TP$Dataset
BES_TN = BES_TN$StationList
BES_TP = BES_TP$StationList
rm(BES_WQ_Sites_TN_Gauged, BES_WQ_Sites_TP_Gauged)

# Plot a map of the stations with zero and negative records----
png(paste0('BES_TNTP_ZerosNegsMap.png'), res = 300, units = 'in', width = 6, height = 6)
plot(ROI)
plot(MDstreams, col = 'skyblue', add = TRUE)
#All gauges
plot(BES_WQ_Sites_locs, pch = 16, add = TRUE)
#Gauges with zeros
plot(BES_WQ_Sites_locs_TN[!is.na(BES_WQ_Sites_locs_TN$Zero),], pch = 16, col = 'red', add = TRUE)
plot(BES_WQ_Sites_locs_TP[!is.na(BES_WQ_Sites_locs_TP$Zero),], pch = 16, col = 'red', add = TRUE)
#Gauges with negative values
plot(BES_WQ_Sites_locs_TN[!is.na(BES_WQ_Sites_locs_TN$Neg),], pch = 16, col = 'red', add = TRUE)
plot(BES_WQ_Sites_locs_TP[!is.na(BES_WQ_Sites_locs_TP$Neg),], pch = 16, col = 'red', add = TRUE)
#Gauges with both
plot(BES_WQ_Sites_locs_TN[which(!is.na(BES_WQ_Sites_locs_TN$Neg) & !is.na(BES_WQ_Sites_locs_TN$Zero)),], pch = 16, col = 'purple', add = TRUE)
plot(BES_WQ_Sites_locs_TP[which(!is.na(BES_WQ_Sites_locs_TP$Neg) & !is.na(BES_WQ_Sites_locs_TP$Zero)),], pch = 16, col = 'purple', add = TRUE)

# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
legend('bottomleft', title = 'Water Quality Stations', legend = c('Zeros in Record', 'Negatives in Record', 'Both', 'Neither'), pch = 16, col = c('red', 'blue', 'purple', 'black'), bty = 'n')
dev.off()

# Make a map of gauge locations colored by their record lengths, corrected for the total amount of missing data----
BES_WQ_Sites_locs_TN$RecordLength = BES_WQ_Sites_locs_TN$RecordLengthMinusGaps = NA
BES_WQ_Sites_locs_TP$RecordLength = BES_WQ_Sites_locs_TP$RecordLengthMinusGaps = NA
for (i in 1:length(BES_TN)){
  BES_WQ_Sites_locs_TN$RecordLength[which(BES_WQ_Sites_locs_TN$Site == BES_TN[[i]]$Site[1])] = as.numeric((max(BES_TN[[i]]$SortDate) - min(BES_TN[[i]]$SortDate)))
}
rm(i)
for (i in 1:length(BES_TP)){
  BES_WQ_Sites_locs_TP$RecordLength[which(BES_WQ_Sites_locs_TN$Site == BES_TN[[i]]$Site[1])] = as.numeric((max(BES_TP[[i]]$SortDate) - min(BES_TP[[i]]$SortDate)))
}
rm(i)
BES_WQ_Sites_locs_TN$RecordLengthMinusGaps = BES_WQ_Sites_locs_TN$RecordLength - BES_WQ_Sites_locs_TN$MissingData_d
BES_WQ_Sites_locs_TP$RecordLengthMinusGaps = BES_WQ_Sites_locs_TP$RecordLength - BES_WQ_Sites_locs_TP$MissingData_d
#in years
BES_WQ_Sites_locs_TN$RecordLength = BES_WQ_Sites_locs_TN$RecordLength/365.25
BES_WQ_Sites_locs_TN$RecordLengthMinusGaps = BES_WQ_Sites_locs_TN$RecordLengthMinusGaps/365.25
BES_WQ_Sites_locs_TP$RecordLength = BES_WQ_Sites_locs_TP$RecordLength/365.25
BES_WQ_Sites_locs_TP$RecordLengthMinusGaps = BES_WQ_Sites_locs_TP$RecordLengthMinusGaps/365.25

#Color by decades
scaleRange = c(0,20)
scaleBy = 5
Pal = rev(rainbow((scaleRange[2] - scaleRange[1])/scaleBy))

setwd(wd_BESN)
png('TN_RecordLengths.png', res = 300, units = 'in', width = 6, height = 6)
plot(ROI)
plot(MDstreams, col = 'skyblue', add = TRUE)
# Gauges colored by their record lengths
plot(BES_WQ_Sites_locs_TN, pch = 16, col = colFun(BES_WQ_Sites_locs_TN$RecordLength), add = TRUE)
# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 370000, yb = 4346000, len = 700, col = 'black', lab = 'N')
legend('right', title = 'Nitrogen Station \n Record Lengths \n (years)', legend = seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy), pch = 16, col = colFun(seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy)), bty = 'n')
dev.off()

setwd(wd_BESP)
png('TP_RecordLengths.png', res = 300, units = 'in', width = 6, height = 6)
plot(ROI)
plot(MDstreams, col = 'skyblue', add = TRUE)
# Gauges colored by their record lengths
plot(BES_WQ_Sites_locs_TP, pch = 16, col = colFun(BES_WQ_Sites_locs_TP$RecordLength), add = TRUE)
# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 370000, yb = 4346000, len = 700, col = 'black', lab = 'N')
legend('right', title = 'Phosphorus Station \n Record Lengths \n (years)', legend = seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy), pch = 16, col = colFun(seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy)), bty = 'n')
dev.off()

scaleRange = c(0,10)
scaleBy = 1
Pal = rev(rainbow((scaleRange[2] - scaleRange[1])/scaleBy))

setwd(wd_BESN)
png('TN_RecordLengthsMinusGaps.png', res = 300, units = 'in', width = 6, height = 6)
plot(ROI)
plot(MDstreams, col = 'skyblue', add = TRUE)
# Gauges colored by their record lengths
plot(BES_WQ_Sites_locs_TN, pch = 16, col = colFun(BES_WQ_Sites_locs_TN$RecordLengthMinusGaps), add = TRUE)
# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 370000, yb = 4346000, len = 700, col = 'black', lab = 'N')
legend('right', title = 'Nitrogen Station \n Record Lengths \n (years)', legend = seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy), pch = 16, col = colFun(seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy)), bty = 'n')
dev.off()

setwd(wd_BESP)
png('TP_RecordLengthsMinusGaps.png', res = 300, units = 'in', width = 6, height = 6)
plot(ROI)
plot(MDstreams, col = 'skyblue', add = TRUE)
# Gauges colored by their record lengths
plot(BES_WQ_Sites_locs_TP, pch = 16, col = colFun(BES_WQ_Sites_locs_TP$RecordLengthMinusGaps), add = TRUE)
# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 370000, yb = 4346000, len = 700, col = 'black', lab = 'N')
legend('right', title = 'Phosphorus Station \n Record Lengths \n (years)', legend = seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy), pch = 16, col = colFun(seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy)), bty = 'n')
dev.off()

# Aggregate into average concentrations, if desired---- 
setwd(dir_WChem)
BES_TN_agg = aggregateTimesteps(StationList = BES_TN, aggVal = c('d', 'm', 'a'), aggVar = 'TN..mg.N.L.', date = 'SortDate', site = 'Site', fun = 'mean')
BES_TP_agg = aggregateTimesteps(StationList = BES_TP, aggVal = c('d', 'm', 'a'), aggVar = 'TP..ugP.L.', date = 'SortDate', site = 'Site', fun = 'mean')
BES_TN_d = BES_TN_agg$daily
BES_TN_m = BES_TN_agg$mthyr
BES_TN_a = BES_TN_agg$ann
BES_TP_d = BES_TP_agg$daily
BES_TP_m = BES_TP_agg$mthyr
BES_TP_a = BES_TP_agg$ann
rm(BES_TN_agg, BES_TP_agg)

#  Handle missing data in the daily, monthly, and annual aggregated timeseries----
BES_TN_d2 = FillMissingDates_par(Dataset = BES_WQ_Sites_locs_TN, StationList = BES_TN_d, Var = 'TN..mg.N.L.', 
                         Date = 'SortDate', gapType = 'd', site_no_D = 'Site', 
                         site_no_SL = 'Site', NoNAcols = c('USGSgauge', 'site_no'))
BES_WQ_Sites_locs_TN = BES_TN_d2$Dataset
BES_TN_d = BES_TN_d2$StationList
rm(BES_TN_d2)

BES_TN_m2 = FillMissingDates_par(Dataset = BES_WQ_Sites_locs_TN, StationList = BES_TN_m, Var = 'TN..mg.N.L.', 
                             Date = 'YrMthDy', gapType = 'm', site_no_D = 'Site', 
                             site_no_SL = 'Site', NoNAcols = c('USGSgauge', 'site_no'))
BES_WQ_Sites_locs_TN = BES_TN_m2$Dataset
BES_TN_m = BES_TN_m2$StationList
rm(BES_TN_m2)

BES_TN_a2 = FillMissingDates(Dataset = BES_WQ_Sites_locs_TN, StationList = BES_TN_a, Var = 'TN..mg.N.L.', 
                             Date = 'YrMthDy', gapType = 'a', site_no_D = 'Site', 
                             site_no_SL = 'Site', NoNAcols = c('USGSgauge', 'site_no'))
BES_WQ_Sites_locs_TN = BES_TN_a2$Dataset
BES_TN_a = BES_TN_a2$StationList
rm(BES_TN_a2)

#Phosphorus
BES_TP_d2 = FillMissingDates_par(Dataset = BES_WQ_Sites_locs_TP, StationList = BES_TP_d, Var = 'TP..ugP.L.', 
                             Date = 'SortDate', gapType = 'd', site_no_D = 'Site', 
                             site_no_SL = 'Site', NoNAcols = c('USGSgauge', 'site_no'))
BES_WQ_Sites_locs_TP = BES_TP_d2$Dataset
BES_TP_d = BES_TP_d2$StationList
rm(BES_TP_d2)

BES_TP_m2 = FillMissingDates_par(Dataset = BES_WQ_Sites_locs_TP, StationList = BES_TP_m, Var = 'TP..ugP.L.', 
                             Date = 'YrMthDy', gapType = 'm', site_no_D = 'Site', 
                             site_no_SL = 'Site', NoNAcols = c('USGSgauge', 'site_no'))
BES_WQ_Sites_locs_TP = BES_TP_m2$Dataset
BES_TP_m = BES_TP_m2$StationList
rm(BES_TP_m2)

BES_TP_a2 = FillMissingDates(Dataset = BES_WQ_Sites_locs_TP, StationList = BES_TP_a, Var = 'TP..ugP.L.', 
                             Date = 'YrMthDy', gapType = 'a', site_no_D = 'Site', 
                             site_no_SL = 'Site', NoNAcols = c('USGSgauge', 'site_no'))
BES_WQ_Sites_locs_TP = BES_TP_a2$Dataset
BES_TP_a = BES_TP_a2$StationList
rm(BES_TP_a2)

# Add streamflow data to the timeseries----
#Fixme: write monthly and annual average streamflows in output from USGS function
for(i in 1:length(BES_TN_d)){
  BES_TN_d[[i]]$Flow = NA
  BES_TN_d[[i]]$FlowCode = NA
  BES_TN_d[[i]]$Zero = NA
  BES_TN_d[[i]]$Neg = NA
  station = as.data.frame(StreamStationList[names(StreamStationList) == BES_TN[names(BES_TN) == BES_TN_d[[i]]$Site[1]][[1]]$USGSgauge[1]][[1]])
  #Find the first date in BES_TN_d
  Ind1 = which(as.Date(station$Date) == BES_TN_d[[i]]$SortDate[1])
  #Find the last date in BES_TN_d
  Indl = which(as.Date(station$Date) == BES_TN_d[[i]]$SortDate[length(BES_TN_d[[i]]$SortDate)])
  #Check that the length is the same
  if (((Indl+1) - Ind1) != nrow(BES_TN_d[[i]])){
    stop('Streamflow and nitrogen dates are different from nitrogen start to end dates')
  }
  BES_TN_d[[i]]$Flow = station[Ind1:Indl,'X_00060_00003']
  BES_TN_d[[i]]$FlowCode = as.character(station[Ind1:Indl,'X_00060_00003_cd'])
  BES_TN_d[[i]]$Zero = station[Ind1:Indl,'Zero']
  BES_TN_d[[i]]$Neg = station[Ind1:Indl,'Neg']
}
rm(i, station, Ind1, Indl)
for(i in 1:length(BES_TP_d)){
  BES_TP_d[[i]]$Flow = NA
  BES_TP_d[[i]]$FlowCode = NA
  BES_TP_d[[i]]$Zero = NA
  BES_TP_d[[i]]$Neg = NA
  station = as.data.frame(StreamStationList[names(StreamStationList) == BES_TP[names(BES_TP) == BES_TP_d[[i]]$Site[1]][[1]]$USGSgauge[1]][[1]])
  #Find the first date in BES_TP_d
  Ind1 = which(as.Date(station$Date) == BES_TP_d[[i]]$SortDate[1])
  #Find the last date in BES_TP_d
  Indl = which(as.Date(station$Date) == BES_TP_d[[i]]$SortDate[length(BES_TP_d[[i]]$SortDate)])
  #Check that the length is the same
  if (((Indl+1) - Ind1) != nrow(BES_TP_d[[i]])){
    stop('Streamflow and phosphorus dates are different from phosphorus start to end dates')
  }
  BES_TP_d[[i]]$Flow = station[Ind1:Indl,'X_00060_00003']
  BES_TP_d[[i]]$FlowCode = as.character(station[Ind1:Indl,'X_00060_00003_cd'])
  BES_TP_d[[i]]$Zero = station[Ind1:Indl,'Zero']
  BES_TP_d[[i]]$Neg = station[Ind1:Indl,'Neg']
}
rm(i, station, Ind1, Indl)

# Plot TN timeseries----
setwd(wd_BESN)
#Fixme: parallelize these plots
for (i in 1:length(BES_TN)){
  png(paste0('TN_Timeseries_', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = BES_TN[[i]]$TN..mg.N.L., x = as.Date(BES_TN[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = 'Total Nitrogen (mg/L as N)', 
       main = paste0('Station #', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1]),
       ylim = c(0, 10), xlim = c(as.Date("1995-01-01"), as.Date("2010-01-01")))
  if(length(which(is.na(BES_TN[[i]]$TN..mg.N.L.))) > 0){
    #Plot where there are NA days that do not correspond to detection limits
    par(new = TRUE)
    plot(y = rep(0, length(which(is.na(BES_TN[[i]]$TN..mg.N.L.)))), x = as.Date(BES_TN[[i]]$SortDate[which(is.na(BES_TN[[i]]$TN..mg.N.L.))]), pch = 4, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = c(0, 0.5), xlim = c(as.Date("1995-01-01"), as.Date("2010-01-01")), axes = FALSE,
         col = 'grey')
    legend('topright', legend = 'NA values', col = 'grey', pch = 4)
  }
  dev.off()
  
  #With map and timeseries
  png(paste0('TN_Timeseries_Map_', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
  layout(rbind(c(1,2)))
  
  plot(y = BES_TN[[i]]$TN..mg.N.L., x = as.Date(BES_TN[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = 'Total Nitrogen (mg/L as N)', 
       main = paste0('Station #', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1]),
       ylim = c(0, 10), xlim = c(as.Date("1995-01-01"), as.Date("2010-01-01")))
  if(length(which(is.na(BES_TN[[i]]$TN..mg.N.L.))) > 0){
    #Plot where there are NA days that do not correspond to detection limits
    par(new = TRUE)
    plot(y = rep(0, length(which(is.na(BES_TN[[i]]$TN..mg.N.L.)))), x = as.Date(BES_TN[[i]]$SortDate[which(is.na(BES_TN[[i]]$TN..mg.N.L.))]), pch = 4, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = c(0, 0.5), xlim = c(as.Date("1995-01-01"), as.Date("2010-01-01")), axes = FALSE,
         col = 'grey')
    legend('topright', legend = 'NA values', col = 'grey', pch = 4)
  }
  
  #map
  plot(ROI)
  plot(MDstreams, col = 'skyblue', add = TRUE)
  #All gauges
  plot(BES_WQ_Sites_locs_TN, pch = 16, add = TRUE)
  #Gauge selected for timeseries plot
  plot(BES_WQ_Sites_locs_TN[which(BES_WQ_Sites_locs_TN$Site == BES_TN[[i]]$Site[1]),], pch = 16, col = 'red', add = TRUE)
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
  legend('bottomleft', title = 'Nitrogen Stations', legend = c('Selected', 'Other'), pch = 16, col = c('red', 'black'), bty = 'n')
  title(main = paste0('Station #', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1]))
  dev.off()
  
  #Fixme: 20 is temporary - the whole if statement condition should be improved.
  if (length(which(!is.na(BES_TN_d[[i]]$TN..mg.N.L.)))>20){
    #hydroTSM plots----
    #daily timeseries
    dts = zoo(BES_TN_d[[i]]$TN..mg.N.L., order.by = BES_TN_d[[i]]$SortDate)
    
    #monthly timeseries for mean daily flows in a month
    mts = zoo(BES_TN_m[[i]]$TN..mg.N.L., order.by = as.Date(BES_TN_m[[i]]$YrMthDy))
    #Monthly matrix
    M = formatMonthlyMatrix(mts)
    # Plotting the monthly values as matrixplot 
    png(paste0('TNMonthly_', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 6, height = 3)
    print(matrixplot(M, ColorRamp="Precipitation", main="Mean Monthly Nitrogen"))
    dev.off()
    
    #Summary EDA plots
    png(paste0('TN_EDA_', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
    hydroplot(dts, FUN = mean, col = 'black', var.unit = 'mean mg/L as N')
    dev.off()
    
    #Seasonal flows
    #Fixme: check that every season is represented before using this.
    png(paste0('TN_EDA_Seasonal_', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
    hydroplot(dts, FUN = mean, col = 'black', var.unit = 'mean mg/L as N', pfreq = 'seasonal', ylab = 'Mean Nitrogen (mg/L)')
    dev.off()
  }
  
  #N/streamflow
  png(paste0('TN_Flow_Timeseries_', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = BES_TN_d[[i]]$TN..mg.N.L./BES_TN_d[[i]]$Flow, x = as.Date(BES_TN_d[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = 'Total Nitrogen (mg/L as N) / Streamflow (cfs)', 
       main = paste0('Station #', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1]))
  dev.off()
  
  #Concentration vs. Flow
  at.x <- outer(1:9, 10^(-3:5))
  lab.x <- ifelse(log10(at.x) %% 1 == 0, sapply(log10(at.x),function(k)
    as.expression(bquote(10^ .(k)))), NA)
  #Color by season
  BES_TN_d[[i]]$cols = NA
  BES_TN_d[[i]]$cols[which(month(as.Date(BES_TN_d[[i]]$SortDate)) %in% c(1,2,3))] = 'blue'
  BES_TN_d[[i]]$cols[which(month(as.Date(BES_TN_d[[i]]$SortDate)) %in% c(4,5,6))] = 'green'
  BES_TN_d[[i]]$cols[which(month(as.Date(BES_TN_d[[i]]$SortDate)) %in% c(7,8,9))] = 'red'
  BES_TN_d[[i]]$cols[which(month(as.Date(BES_TN_d[[i]]$SortDate)) %in% c(10,11,12))] = 'yellow'
  
  png(paste0('TNvsFlow_', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = BES_TN_d[[i]]$TN..mg.N.L., x = BES_TN_d[[i]]$Flow, pch = 16, cex = 0.3,
       xlab = expression(paste('log'[10] ~ '(Flow [cfs])')), 
       ylab = 'Total Nitrogen (mg/L as N)', 
       main = paste0('Station #', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1]), log = 'x', axes = FALSE, col = BES_TN_d[[i]]$cols)
  axis(side=1, at=at.x, labels=lab.x, cex.axis=1.2, las=1, tck = -0.01)
  axis(side=1, at=10^seq(-3,5,1), labels=FALSE, tck = -0.025)
  axis(side = 2)
  box()
  legend('topright', legend = c('Winter', 'Spring', 'Summer', 'Fall'), col = c('blue', 'green', 'red', 'yellow'), pch = 16)
  dev.off()
  
  #Streamflow and nitrogen vertical
  #N/streamflow
  png(paste0('TN_Flow_Both_Timeseries_', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 6, height = 6)
  layout(c(1,2))
  plot(y = BES_TN_d[[i]]$Flow, x = as.Date(BES_TN_d[[i]]$SortDate), pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = 'Streamflow (cfs)', 
       main = paste0('Station #', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1]))
  plot(y = BES_TN_d[[i]]$TN..mg.N.L., x = as.Date(BES_TN_d[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = 'Total Nitrogen (mg/L as N)', 
       main = paste0('Station #', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1]))
  dev.off()
  
  #Map, Q-C, Streamflow, Nitrogen
  png(paste0('TN_EDA_Map', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
  layout(rbind(c(1,2), c(1,2), c(3,2), c(4,2)))
  
  plot(y = BES_TN_d[[i]]$TN..mg.N.L., x = BES_TN_d[[i]]$Flow, pch = 16, cex = 0.3,
       xlab = expression(paste('log'[10] ~ '(Flow [cfs])')), 
       ylab = 'Total Nitrogen (mg/L as N)', 
       main = paste0('Station #', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1]), log = 'x', axes = FALSE, col = BES_TN_d[[i]]$cols)
  axis(side=1, at=at.x, labels=lab.x, cex.axis=1.2, las=1, tck = -0.01)
  axis(side=1, at=10^seq(-3,5,1), labels=FALSE, tck = -0.025)
  axis(side = 2)
  box()
  legend('topright', legend = c('Winter', 'Spring', 'Summer', 'Fall'), col = c('blue', 'green', 'red', 'yellow'), pch = 16)
  
  #map
  plot(ROI)
  plot(MDstreams, col = 'skyblue', add = TRUE)
  #All gauges
  plot(BES_WQ_Sites_locs_TN, pch = 16, add = TRUE)
  #Gauge selected for timeseries plot
  plot(BES_WQ_Sites_locs_TN[which(BES_WQ_Sites_locs_TN$Site == BES_TN[[i]]$Site[1]),], pch = 16, col = 'red', add = TRUE)
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
  legend('bottomleft', title = 'Nitrogen Stations', legend = c('Selected', 'Other'), pch = 16, col = c('red', 'black'), bty = 'n')
  title(main = paste0('Station #', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1]))
  
  plot(y = BES_TN_d[[i]]$Flow, x = as.Date(BES_TN_d[[i]]$SortDate), pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = 'Streamflow (cfs)', 
       main = paste0('Station #', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1]))
  
  plot(y = BES_TN_d[[i]]$TN..mg.N.L., x = as.Date(BES_TN_d[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = 'Total Nitrogen (mg/L as N)', 
       main = paste0('Station #', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1]))
  dev.off()
  
  png(paste0('TN_ScatterHist', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
  scatterHist_mod(x = log10(BES_TN_d[[i]]$Flow[BES_TN_d[[i]]$Flow > 0]), y = BES_TN_d[[i]]$TN..mg.N.L.[BES_TN_d[[i]]$Flow > 0], density = FALSE, freq = FALSE, x.breaks = 100, y.breaks = 100, correl = FALSE, xlab = expression(paste('log'[10] ~ '(Flow [cfs])')), ylab = 'Total Nitrogen (mg/L as N)', cex.lab = 1.5, cex.axis = 1.5, title = '', pch = 16, ellipse = FALSE, smooth = FALSE, col = BES_TN_d[[i]]$cols[BES_TN_d[[i]]$Flow > 0])
  dev.off()
}
rm(i, dts, mts, M, at.x, lab.x)

# Plot TP timeseries----
setwd(wd_BESP)
for (i in 1:length(BES_TP)){
  png(paste0('TP_Timeseries_', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = BES_TP[[i]]$TP..ugP.L., x = as.Date(BES_TP[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = expression(paste('Total Phosphorus ( ', mu, 'g/L as P)')), 
       main = paste0('Station #', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1]),
       ylim = c(0, 10), xlim = c(as.Date("1995-01-01"), as.Date("2010-01-01")))
  if(length(which(is.na(BES_TP[[i]]$TP..ugP.L.))) > 0){
    #Plot where there are NA days that do not correspond to detection limits
    par(new = TRUE)
    plot(y = rep(0, length(which(is.na(BES_TP[[i]]$TP..ugP.L.)))), x = as.Date(BES_TP[[i]]$SortDate[which(is.na(BES_TP[[i]]$TP..ugP.L.))]), pch = 4, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = c(0, 0.5), xlim = c(as.Date("1995-01-01"), as.Date("2010-01-01")), axes = FALSE,
         col = 'grey')
    legend('topright', legend = 'NA values', col = 'grey', pch = 4)
  }
  dev.off()
  
  #With map and timeseries
  png(paste0('TP_Timeseries_Map_', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
  layout(rbind(c(1,2)))
  
  plot(y = BES_TP[[i]]$TP..ugP.L., x = as.Date(BES_TP[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = expression(paste('Total Phosphorus ( ', mu, 'g/L as P)')), 
       main = paste0('Station #', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1]),
       ylim = c(0, 10), xlim = c(as.Date("1995-01-01"), as.Date("2010-01-01")))
  if(length(which(is.na(BES_TP[[i]]$TP..ugP.L.))) > 0){
    #Plot where there are NA days that do not correspond to detection limits
    par(new = TRUE)
    plot(y = rep(0, length(which(is.na(BES_TP[[i]]$TP..ugP.L.)))), x = as.Date(BES_TP[[i]]$SortDate[which(is.na(BES_TP[[i]]$TP..ugP.L.))]), pch = 4, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = c(0, 0.5), xlim = c(as.Date("1995-01-01"), as.Date("2010-01-01")), axes = FALSE,
         col = 'grey')
    legend('topright', legend = 'NA values', col = 'grey', pch = 4)
  }
  
  #map
  plot(ROI)
  plot(MDstreams, col = 'skyblue', add = TRUE)
  #All gauges
  plot(BES_WQ_Sites_locs_TP, pch = 16, add = TRUE)
  #Gauge selected for timeseries plot
  plot(BES_WQ_Sites_locs_TP[which(BES_WQ_Sites_locs_TP$Site == BES_TP[[i]]$Site[1]),], pch = 16, col = 'red', add = TRUE)
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
  legend('bottomleft', title = 'Phosphorus Stations', legend = c('Selected', 'Other'), pch = 16, col = c('red', 'black'), bty = 'n')
  title(main = paste0('Station #', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1]))
  dev.off()
  
  #Fixme: 20 is temporary - see above for nitrogen
  if (length(which(!is.na(BES_TP_d[[i]]$TP..ugP.L.)))>20){
    #hydroTSM plots----
    #daily timeseries
    dts = zoo(BES_TP_d[[i]]$TP..ugP.L., order.by = BES_TP_d[[i]]$SortDate)
    
    #monthly timeseries for mean daily flows in a month
    mts = zoo(BES_TP_m[[i]]$TP..ugP.L., order.by = as.Date(BES_TP_m[[i]]$YrMthDy))
    #Monthly matrix
    M = formatMonthlyMatrix(mts)
    # Plotting the monthly values as matrixplot 
    png(paste0('TPMonthly_', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 6, height = 3)
    print(matrixplot(M, ColorRamp="Precipitation", main="Mean Monthly Phosphorus"))
    dev.off()
    
    #Summary EDA plots
    png(paste0('TP_EDA_', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
    hydroplot(dts, FUN = mean, col = 'black', var.unit = 'mean ug/L as P')
    dev.off()
    
    #Seasonal flows
    #Fixme: check that every season is represented before using this.
    png(paste0('TP_EDA_Seasonal_', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
    hydroplot(dts, FUN = mean, col = 'black', var.unit = 'mean ug/L as P', pfreq = 'seasonal', ylab = expression(paste('Mean Phosphorus (', mu, 'g/L)')))
    dev.off()
  }
  
  #P/streamflow
  png(paste0('TP_Flow_Timeseries_', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = BES_TP_d[[i]]$TP..ugP.L./BES_TP_d[[i]]$Flow, x = as.Date(BES_TP_d[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = expression(paste('Total Phosphorus ( ', mu, 'g/L as P) / Streamflow (cfs)')), 
       main = paste0('Station #', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1]))
  dev.off()
  
  #Concentration vs. Flow
  at.x <- outer(1:9, 10^(-3:5))
  lab.x <- ifelse(log10(at.x) %% 1 == 0, sapply(log10(at.x),function(k)
    as.expression(bquote(10^ .(k)))), NA)
  #Color by season
  BES_TP_d[[i]]$cols = NA
  BES_TP_d[[i]]$cols[which(month(as.Date(BES_TP_d[[i]]$SortDate)) %in% c(1,2,3))] = 'blue'
  BES_TP_d[[i]]$cols[which(month(as.Date(BES_TP_d[[i]]$SortDate)) %in% c(4,5,6))] = 'green'
  BES_TP_d[[i]]$cols[which(month(as.Date(BES_TP_d[[i]]$SortDate)) %in% c(7,8,9))] = 'red'
  BES_TP_d[[i]]$cols[which(month(as.Date(BES_TP_d[[i]]$SortDate)) %in% c(10,11,12))] = 'yellow'
  
  png(paste0('TPvsFlow_', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = BES_TP_d[[i]]$TP..ugP.L., x = BES_TP_d[[i]]$Flow, pch = 16, cex = 0.3,
       xlab = expression(paste('log'[10] ~ '(Flow [cfs])')), 
       ylab = expression(paste('Total Phosphorus ( ', mu, 'g/L as P)')), 
       main = paste0('Station #', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1]), log = 'x', axes = FALSE, col = BES_TP_d[[i]]$cols)
  axis(side=1, at=at.x, labels=lab.x, cex.axis=1.2, las=1, tck = -0.01)
  axis(side=1, at=10^seq(-3,5,1), labels=FALSE, tck = -0.025)
  axis(side = 2)
  box()
  legend('topright', legend = c('Winter', 'Spring', 'Summer', 'Fall'), col = c('blue', 'green', 'red', 'yellow'), pch = 16)
  dev.off()
  
  #Streamflow and phosphorus vertical
  #P/streamflow
  png(paste0('TP_Flow_Both_Timeseries_', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 6, height = 6)
  layout(c(1,2))
  plot(y = BES_TP_d[[i]]$Flow, x = as.Date(BES_TP_d[[i]]$SortDate), pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = 'Streamflow (cfs)', 
       main = paste0('Station #', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1]))
  plot(y = BES_TP_d[[i]]$TP..ugP.L., x = as.Date(BES_TP_d[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = expression(paste('Total Phosphorus ( ', mu, 'g/L as P)')),
       main = paste0('Station #', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1]))
  dev.off()
  
  #Map, Q-C, Streamflow, Phosphorus
  png(paste0('TP_EDA_Map', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
  layout(rbind(c(1,2), c(1,2), c(3,2), c(4,2)))
  
  plot(y = BES_TP_d[[i]]$TP..ugP.L., x = BES_TP_d[[i]]$Flow, pch = 16, cex = 0.3,
       xlab = expression(paste('log'[10] ~ '(Flow [cfs])')), 
       ylab = expression(paste('Total Phosphorus ( ',mu,'g/L as P)')), 
       main = paste0('Station #', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1]), log = 'x', axes = FALSE, col = BES_TP_d[[i]]$cols)
  axis(side=1, at=at.x, labels=lab.x, cex.axis=1.2, las=1, tck = -0.01)
  axis(side=1, at=10^seq(-3,5,1), labels=FALSE, tck = -0.025)
  axis(side = 2)
  box()
  legend('topright', legend = c('Winter', 'Spring', 'Summer', 'Fall'), col = c('blue', 'green', 'red', 'yellow'), pch = 16)
  
  #map
  plot(ROI)
  plot(MDstreams, col = 'skyblue', add = TRUE)
  #All gauges
  plot(BES_WQ_Sites_locs_TP, pch = 16, add = TRUE)
  #Gauge selected for timeseries plot
  plot(BES_WQ_Sites_locs_TP[which(BES_WQ_Sites_locs_TP$Site == BES_TP[[i]]$Site[1]),], pch = 16, col = 'red', add = TRUE)
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
  legend('bottomleft', title = 'Phosphorus Stations', legend = c('Selected', 'Other'), pch = 16, col = c('red', 'black'), bty = 'n')
  title(main = paste0('Station #', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1]))
  
  plot(y = BES_TP_d[[i]]$Flow, x = as.Date(BES_TP_d[[i]]$SortDate), pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = 'Streamflow (cfs)', 
       main = paste0('Station #', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1]))
  
  plot(y = BES_TP_d[[i]]$TP..ugP.L., x = as.Date(BES_TP_d[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = expression(paste('Total Phosphorus ( ',mu,'g/L as P)')), 
       main = paste0('Station #', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1]))
  dev.off()
  
  png(paste0('TP_ScatterHist', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
  scatterHist_mod(x = log10(BES_TP_d[[i]]$Flow[BES_TP_d[[i]]$Flow > 0]), y = BES_TP_d[[i]]$TP..ugP.L.[BES_TP_d[[i]]$Flow > 0], density = FALSE, freq = FALSE, x.breaks = 100, y.breaks = 100, correl = FALSE, xlab = expression(paste('log'[10] ~ '(Flow [cfs])')), ylab = expression(paste('Total Phosphorus ( ', mu, 'g/L as P)')), cex.lab = 1.5, cex.axis = 1.5, title = '', pch = 16, ellipse = FALSE, smooth = FALSE, col = BES_TP_d[[i]]$cols[BES_TP_d[[i]]$Flow > 0])
  dev.off()
}
rm(i, dts, mts, M, at.x, lab.x)

# Plot histograms of data----
setwd(wd_BESN)
for (i in 1:length(BES_TN)){
  png(paste0('TN_hist_', BES_TN[[i]]$USGSgauge[1], ' ', BES_TN[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 5, height = 5)
  hist(log10(BES_TN[[i]]$TN..mg.N.L.),
       xlab = expression(paste('log'[10] ~ '(Total Nitrogen [mg/L as N])')), main = paste('Station #', BES_TN[[i]]$USGSgauge[1], BES_TN[[i]]$Site[1]))
  dev.off()
}
rm(i)
setwd(wd_BESP)
for (i in 1:length(BES_TP)){
  png(paste0('TP_hist_', BES_TP[[i]]$USGSgauge[1], ' ', BES_TP[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 5, height = 5)
  hist(log10(BES_TP[[i]]$TP..ugP.L.),
       xlab = expression(paste('log'[10] ~ '(Total Phosphorus [',mu,'g/L as P])')), main = paste('Station #', BES_TN[[i]]$USGSgauge[1], BES_TN[[i]]$Site[1]))
  dev.off()
}
rm(i)

#Fixme: Search for flow-normalized outliers----
# #Select only the non-NA values
# t = BES_TN_d[[i]][which((is.na(BES_TN_d[[i]]$Flow) == FALSE) & (is.na(BES_TN_d[[i]]$TN..mg.N.L.) == FALSE)),c('Flow', 'TN..mg.N.L.')]
# mu = apply(X = t, MARGIN = 2, FUN = mean)
# #Compute the covariance matrix for only the non-NA values
# s = cov(t)
# mahalanobis(x = t, center = FALSE, cov = s)

#Fixme: Check correlation between N, P at all sites to see if it's worthwhile to use both for calibration of model----
# FOR NOW, NOT NEEDED BECAUSE RHESSYS DOESN'T PROVIDE P OUTPUT!
# Color by time of year
# for (i in 1:length(BES_TN_d)){
#   #Gather all of the dates that match in both timeseries
#   
# }
# plot(BES_TN_d$POBR$TN..mg.N.L.[362:nrow(BES_TN_d$POBR)][which((is.na(BES_TN_d$POBR$TN..mg.N.L.[362:nrow(BES_TN_d$POBR)]) == FALSE) & (is.na(BES_TP_d$POBR$TP..ugP.L.) == FALSE))], BES_TP_d$POBR$TP..ugP.L.[which((is.na(BES_TN_d$POBR$TN..mg.N.L.[362:nrow(BES_TN_d$POBR)]) == FALSE) & (is.na(BES_TP_d$POBR$TP..ugP.L.) == FALSE))], log = 'xy')
# plot(BES_TN_d$BARN$TN..mg.N.L.[which((is.na(BES_TN_d$BARN$TN..mg.N.L.) == FALSE) & (is.na(BES_TP_d$BARN$TP..ugP.L.) == FALSE))], BES_TP_d$BARN$TP..ugP.L.[which((is.na(BES_TN_d$BARN$TN..mg.N.L.) == FALSE) & (is.na(BES_TP_d$BARN$TP..ugP.L.) == FALSE))], log = 'xy')

#Fixme: Check for gaps in the TN/NO3 series that can be filled in with TP correlation with TN
# Hypsometric plot for DEM in the ROI----
#Fixme: Can shade in the areas above each gauge or select outlet points as horizontal steps
# Allocation of BMPs could be in certain elevation zones
# Couple with map next to it of the elevations above certain thresholds

#Clip DEM to the ROI and make a spatial grid dataframe
DEM_DR = as(mask(DEM, DeadRun_Outlet), 'SpatialGridDataFrame')
DEM_BR = as(mask(DEM, BaisRun_Outlet), 'SpatialGridDataFrame')
DEM_SL = as(mask(DEM, ScottsLevel), 'SpatialGridDataFrame')
#Plot curve
setwd(dir_Sheds)
png('DR_HypsometricCurve.png', res = 300, units = 'in', width = 5, height = 5)
par(mar =c(4,4,1,1))
hypsometric(DEM_DR, col = 'black', main = 'Dead Run Hypsometric Curve')
dev.off()
png('BR_HypsometricCurve.png', res = 300, units = 'in', width = 5, height = 5)
par(mar =c(4,4,1,1))
hypsometric(DEM_BR, col = 'black', main = 'Baisman Run Hypsometric Curve')
dev.off()
png('SL_HypsometricCurve.png', res = 300, units = 'in', width = 5, height = 5)
par(mar =c(4,4,1,1))
hypsometric(DEM_SL, col = 'black', main = "Scott's Level Hypsometric Curve")
dev.off()

# Write water quality data to files----
setwd(dir = dir_WChem)
list.save(x = BES_WQ_Sites, file = f_BESallchem, type = "YAML")
setwd(wd_BESN)
writeOGR(BES_WQ_Sites_locs_TN, dsn = getwd(), layer = f_Nsites, driver = "ESRI Shapefile")
list.save(x = BES_TN, file = f_BES_TNallchem, type = "YAML")
list.save(x = BES_TN_d, file = f_BES_TNd, type = "YAML")
list.save(x = BES_TN_m, file = f_BES_TNm, type = "YAML")
list.save(x = BES_TN_a, file = f_BES_TNa, type = "YAML")
setwd(wd_BESP)
writeOGR(BES_WQ_Sites_locs_TP, dsn = getwd(), layer = f_Psites, driver = "ESRI Shapefile")
list.save(x = BES_TP, file = f_BES_TPallchem, type = "YAML")
list.save(x = BES_TP_d, file = f_BES_TPd, type = "YAML")
list.save(x = BES_TP_m, file = f_BES_TPm, type = "YAML")
list.save(x = BES_TP_a, file = f_BES_TPa, type = "YAML")

#Weather Station Data----
setwd(dir_precip)
# BES Precipitation Stations----
BES_Precip = read.csv(file = f_Precip, stringsAsFactors = FALSE)

#Make a list
BES_PrecipList = list()
for(i in 1:length(unique(BES_Precip$Rain_Gauge_ID))){
  #get all the entries from the same station
  f = BES_Precip[BES_Precip$Rain_Gauge_ID == unique(BES_Precip$Rain_Gauge_ID)[i],]
  
  #Process the dates to be standard format
  #Some of the gauges use 2-digit year numbers, others use 4-digit year numbers.
  f$Date = NA
  f$DateTime = NA
  if(nchar(strsplit(strsplit(x = f$Date_Time_.EST.[1], split = '/', fixed = TRUE)[[1]][3], split = ' ', fixed = TRUE)[[1]][1]) == 2){
    f$Date = as.Date(f$Date_Time_.EST., format = '%m/%d/%y')
    f$DateTime = paste(as.POSIXct(f$Date_Time_.EST., format = '%m/%d/%y %H:%M'))
    #Fixme: There are NA DateTimes as a result of daylight savings time not automatically adjusting for these gauges. 
    # Add an hour to spring times, and subtract an hour from fall times.
    # For now, adding 1 hour to the times because they're all the same date (3/14/2010)
    if(length(which(f$DateTime == 'NA')) > 0){
      tvec = NA
      for (t in 1:length(which(f$DateTime == 'NA'))){
        tvec[t] = paste(as.POSIXct(paste(strsplit(f$Date_Time_.EST.[which(f$DateTime == 'NA')][t], split = ' ', fixed = TRUE)[[1]][1], strsplit(as.character(as.POSIXct(strsplit(f$Date_Time_.EST.[which(f$DateTime == 'NA')][t], split = ' ', fixed = TRUE)[[1]][2], format = '%H:%M') + 1*3600), split = ' ', fixed = TRUE)[[1]][2]), format = '%m/%d/%y %H:%M'))
      }
      f$DateTime[which(f$DateTime == 'NA')] = tvec
    }
  }else{
    f$Date = as.Date(f$Date_Time_.EST., format = '%m/%d/%Y')
    f$DateTime = paste(as.POSIXct(f$Date_Time_.EST., format = '%m/%d/%Y %H:%M'))
    #Fixme: There are NA DateTimes as a result of daylight savings time not automatically adjusting for these gauges. 
    # Add an hour to spring times, and subtract an hour from fall times.
    # For now, adding 1 hour to the times because they're all the same date (3/14/2010)
    if(length(which(f$DateTime == 'NA')) > 0){
      tvec = NA
      for (t in 1:length(which(f$DateTime == 'NA'))){
        tvec[t] = paste(as.POSIXct(paste(strsplit(f$Date_Time_.EST.[which(f$DateTime == 'NA')][t], split = ' ', fixed = TRUE)[[1]][1], strsplit(as.character(as.POSIXct(strsplit(f$Date_Time_.EST.[which(f$DateTime == 'NA')][t], split = ' ', fixed = TRUE)[[1]][2], format = '%H:%M') + 1*3600), split = ' ', fixed = TRUE)[[1]][2]), format = '%m/%d/%Y %H:%M'))
      }
      f$DateTime[which(f$DateTime == 'NA')] = tvec
    }
  }
  #Remove any NA dates. These entries are useless.
  f = f[!is.na(f$Date),]
  
  #Make a new column for the sort date.
  f$SortDate = as.Date(f$Date)
  f$SortDateTime = as.POSIXct(f$DateTime)
  
  #Place the measurements in chronological order by the sort date and time
  # First order by time (which may have some NA values), then order by date (no NA values)
  f = f[order(as.POSIXct(f$SortDateTime)),][order(f[order(as.POSIXct(f$SortDateTime)),]$SortDate),]
  
  #Add to list
  BES_PrecipList = c(BES_PrecipList, list(f))
  names(BES_PrecipList)[length(BES_PrecipList)] = f$Rain_Gauge_ID[1]
}
rm(f,i,t,tvec)

#Add coordinates
BES_Precip_locs = read.csv(file = f_Precip_locs, stringsAsFactors = FALSE)
coordinates(BES_Precip_locs) = c('Long', 'Lat')
proj4string(BES_Precip_locs) = CRS('+init=epsg:4269')
BES_Precip_locs = spTransform(BES_Precip_locs, CRS(pCRS))

#Plot locations of precip gauges
png('Precip_BES.png', res = 300, units = 'in', width = 6, height = 6)
plot(ROI, lwd = 2)
plot(DeadRun_Outlet, add = TRUE)
plot(ScottsLevel, add = TRUE)
plot(BaisRun_Outlet, add = TRUE)
#BES sites
plot(BES_Precip_locs, add = TRUE, pch = 16, col = 'blue')

# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 365000, yb = 4348000, len = 700, col = 'black', lab = 'N')
legend('bottomleft', title = 'Precipitation Sites', legend = c('BES'), col = c('blue'), pch = c(16))
dev.off()

#  Check for duplicate records----
BES_PrecipList = checkDuplicatesAndRemove(StationList = BES_PrecipList, colNames = c('Date_Time_.EST.', 'SortDate', 'SortDateTime', 'Precipitation_.mm.'))

#  Check for zeros and negative values----
BES_PrecipList = checkZerosNegs(StationList = BES_PrecipList, Var = 'Precipitation_.mm.', NegReplace = NA, ZeroReplace = 0)

#Add zeros and negatives to spatial dataset
BES_Precip_locs = addZerosToSpatialDataset(StationList = BES_PrecipList, SpatialDataset = BES_Precip_locs, site_D = 'Gauge', site_SL = 'Rain_Gauge_ID')
BES_Precip_locs = addNegsToSpatialDataset(StationList = BES_PrecipList, SpatialDataset = BES_Precip_locs, site_D = 'Gauge', site_SL = 'Rain_Gauge_ID')

#  Compute rainfall rates----
# Fixme: compute rates and compare for each gauge. Use for error identification

#  Aggregate into total precip per day, month, year---- 
BES_Precip_agg = aggregateTimesteps(StationList = BES_PrecipList, aggVal = c('d', 'm', 'a'), aggVar = 'Precipitation_.mm.', date = 'SortDate', site = 'Rain_Gauge_ID', fun = 'sum')
BES_Precip_d = BES_Precip_agg$daily
BES_Precip_m = BES_Precip_agg$mthyr
BES_Precip_a = BES_Precip_agg$ann
rm(BES_Precip_agg)

#  Handle missing data in the daily, monthly, and annual aggregated timeseries----
BES_Precip_d2 = FillMissingDates_par(Dataset = BES_Precip_locs, StationList = BES_Precip_d, Var = 'Precipitation_.mm.', 
                             Date = 'SortDate', gapType = 'd', site_no_D = 'Gauge', 
                             site_no_SL = 'Rain_Gauge_ID', NoNAcols = c('Site','Rain_Gauge_ID'))
BES_Precip_locs = BES_Precip_d2$Dataset
BES_Precip_d = BES_Precip_d2$StationList
rm(BES_Precip_d2)

BES_Precip_m2 = FillMissingDates(Dataset = BES_Precip_locs, StationList = BES_Precip_m, Var = 'Precipitation_.mm.', 
                                 Date = 'YrMthDy', gapType = 'm', site_no_D = 'Gauge', 
                                 site_no_SL = 'Rain_Gauge_ID', NoNAcols = c('Site','Rain_Gauge_ID'))
BES_Precip_locs = BES_Precip_m2$Dataset
BES_Precip_m = BES_Precip_m2$StationList
rm(BES_Precip_m2)

BES_Precip_a2 = FillMissingDates(Dataset = BES_Precip_locs, StationList = BES_Precip_a, Var = 'Precipitation_.mm.', 
                                 Date = 'YrMthDy', gapType = 'a', site_no_D = 'Gauge', 
                                 site_no_SL = 'Rain_Gauge_ID', NoNAcols = c('Site','Rain_Gauge_ID'))
BES_Precip_locs = BES_Precip_a2$Dataset
BES_Precip_a = BES_Precip_a2$StationList
rm(BES_Precip_a2)

#  Plot total precip timeseries----
for (i in 1:length(BES_PrecipList)){
  png(paste0('Precip_Timeseries_', BES_PrecipList[[i]]$Rain_Gauge_ID[1], ' ', BES_PrecipList[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = BES_PrecipList[[i]]$Precipitation_.mm., x = as.Date(BES_PrecipList[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = 'Total Precipitation (mm)', 
       main = paste0('Station #', BES_PrecipList[[i]]$Rain_Gauge_ID[1], '\n', BES_PrecipList[[i]]$Site[1]))
  dev.off()
  
  #With map and timeseries
  png(paste0('Precip_Timeseries_Map', BES_PrecipList[[i]]$Rain_Gauge_ID[1], ' ', BES_PrecipList[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
  layout(rbind(c(1,2)))
  
  plot(y = BES_PrecipList[[i]]$Precipitation_.mm., x = as.Date(BES_PrecipList[[i]]$SortDate), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = 'Total Precipitation (mm)', 
       main = paste0('Station #', BES_PrecipList[[i]]$Rain_Gauge_ID[1], '\n', BES_PrecipList[[i]]$Site[1]))
  
  #map
  plot(ROI)
  #All gauges
  plot(BES_Precip_locs, pch = 16, add = TRUE)
  #Gauge selected for timeseries plot
  plot(BES_Precip_locs[which(BES_Precip_locs$Site == BES_PrecipList[[i]]$Site[1]),], pch = 16, col = 'red', add = TRUE)
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4347000, len = 700, col = 'black', lab = 'N')
  legend('bottomleft', title = 'Precipitation Stations', legend = c('Selected', 'Other'), pch = 16, col = c('red', 'black'), bty = 'n')
  title(main = paste0('Station #', BES_PrecipList[[i]]$Rain_Gauge_ID[1], '\n', BES_PrecipList[[i]]$Site[1]))
  dev.off()
  
  #Fixme: 20 is temporary - same as above
  if (length(which(!is.na(BES_Precip_d[[i]]$Precipitation_.mm.)))>20){
    #hydroTSM plots----
    #daily timeseries
    dts = zoo(BES_Precip_d[[i]]$Precipitation_.mm., order.by = BES_Precip_d[[i]]$SortDate)
    
    #monthly timeseries for mean daily flows in a month
    mts = zoo(BES_Precip_m[[i]]$Precipitation_.mm., order.by = as.Date(BES_Precip_m[[i]]$YrMthDy))
    #Monthly matrix
    M = formatMonthlyMatrix(mts)
    # Plotting the monthly values as matrixplot 
    png(paste0('PrecipMonthly_', BES_PrecipList[[i]]$Rain_Gauge_ID[1], ' ', BES_PrecipList[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 6, height = 3)
    print(matrixplot(M, ColorRamp="Precipitation", main="Mean Monthly Precipitation"))
    dev.off()
    
    #Summary EDA plots
    png(paste0('Precip_EDA_', BES_PrecipList[[i]]$Rain_Gauge_ID[1], ' ', BES_PrecipList[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
    hydroplot(dts, FUN = sum, col = 'black', var.unit = 'mm')
    dev.off()
    
    #Seasonal flows
    #Fixme: check that every season is represented before using this.
    png(paste0('Precip_EDA_Seasonal_', BES_PrecipList[[i]]$Rain_Gauge_ID[1], ' ', BES_PrecipList[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
    hydroplot(dts, FUN = sum, col = 'black', var.unit = 'mm', pfreq = 'seasonal', ylab = 'Precipitationn (mm)')
    dev.off()
  }
}
rm(i, dts, mts, M)

#  Spatially aggregate the gauges by averaging the 2 gauges per site----
#Working on daily precip
# sites are in pairs of 2 in the list
BES_Precip_Avg_d = list()
for(i in seq(1,(length(BES_Precip_d)-1),2)){
  #Check for the start date at each site
  Inds1 = BES_Precip_d[[i]]$SortDate[1]
  Inds2 = BES_Precip_d[[i+1]]$SortDate[1]
  #If necessary, extend the length of the shorter site with NAs
  if (Inds1 != Inds2){
    dif = as.numeric(Inds1 - Inds2)
    if (dif > 0){
      #Add dif number of days to start of second dataset, using the data from the start of the first dataset
      BES_Precip_d[[i+1]] = rbind(BES_Precip_d[[i]][1:dif,], BES_Precip_d[[i+1]])
      #Set the precip to NA and the gauge ID to that of the second dataset
      BES_Precip_d[[i+1]][1:dif,'Precipitation_.mm.'] = NA
      BES_Precip_d[[i+1]][1:dif,'Rain_Gauge_ID'] = BES_Precip_d[[i+1]]$Rain_Gauge_ID[dif+1]
    }else{
      #Add dif number of days to start of second dataset, using the data from the start of the first dataset
      BES_Precip_d[[i]] = rbind(BES_Precip_d[[i+1]][1:dif,], BES_Precip_d[[i]])
      #Set the precip to NA and the gauge ID to that of the second dataset
      BES_Precip_d[[i]][1:dif,'Precipitation_.mm.'] = NA
      BES_Precip_d[[i]][1:dif,'Rain_Gauge_ID'] = BES_Precip_d[[i]]$Rain_Gauge_ID[dif+1]
    }
  }
  #Check for the end date at each site
  Inde1 = BES_Precip_d[[i]]$SortDate[nrow(BES_Precip_d[[i]])]
  Inde2 = BES_Precip_d[[i+1]]$SortDate[nrow(BES_Precip_d[[i+1]])]
  #If necessary, extend the length of the shorter site with NAs
  if (Inde1 != Inde2){
    dif = as.numeric(Inde1 - Inde2)
    if (dif > 0){
      #Add dif number of days to end of second dataset, using the data from the end of the first dataset
      BES_Precip_d[[i+1]] = rbind(BES_Precip_d[[i+1]], BES_Precip_d[[i]][(nrow(BES_Precip_d[[i]])-(dif-1)):nrow(BES_Precip_d[[i]]),])
      #Set the precip to NA and the gauge ID to that of the second dataset
      BES_Precip_d[[i+1]][(nrow(BES_Precip_d[[i+1]])-(dif-1)):nrow(BES_Precip_d[[i+1]]),'Precipitation_.mm.'] = NA
      BES_Precip_d[[i+1]][(nrow(BES_Precip_d[[i+1]])-(dif-1)):nrow(BES_Precip_d[[i+1]]),'Rain_Gauge_ID'] = BES_Precip_d[[i+1]]$Rain_Gauge_ID[1]
    }else{
      #Add dif number of days to end of first dataset, using the data from the end of the second dataset
      BES_Precip_d[[i]] = rbind(BES_Precip_d[[i]], BES_Precip_d[[i+1]][(nrow(BES_Precip_d[[i+1]])-(dif-1)):nrow(BES_Precip_d[[i+1]]),])
      #Set the precip to NA and the gauge ID to that of the second dataset
      BES_Precip_d[[i]][(nrow(BES_Precip_d[[i]])-(dif-1)):nrow(BES_Precip_d[[i]]),'Precipitation_.mm.'] = NA
      BES_Precip_d[[i]][(nrow(BES_Precip_d[[i]])-(dif-1)):nrow(BES_Precip_d[[i]]),'Rain_Gauge_ID'] = BES_Precip_d[[i]]$Rain_Gauge_ID[1]
    }
  }
  #Combine into dataframe of both precipitation vectors and the date vector
  sitePrecip = cbind(BES_Precip_d[[i]][,1:2], BES_Precip_d[[i+1]]$Precipitation_.mm.)
  colnames(sitePrecip) = c('SortDate', paste0('P_', BES_Precip_d[[i]]$Rain_Gauge_ID[1]), paste0('P_', BES_Precip_d[[i+1]]$Rain_Gauge_ID[1]))
  
  #Any time only one site has an NA, use the other site for that day
  #Compute the mean value for each day
  sitePrecip$mean = apply(X = sitePrecip[,2:3], MARGIN = 1, FUN = mean, na.rm = TRUE)
  #Compute the difference from the mean value
  sitePrecip$meanDiff1 = sitePrecip[,2] - sitePrecip$mean
  sitePrecip$meanDiff2 = sitePrecip[,3] - sitePrecip$mean
  #Compute the total difference between the two gauges
  sitePrecip$GaugeDiff = sitePrecip[,2] - sitePrecip[,3]
  #return the dataframe
  BES_Precip_Avg_d = c(BES_Precip_Avg_d, list(sitePrecip))
  names(BES_Precip_Avg_d)[length(BES_Precip_Avg_d)] = BES_Precip_locs$Site[i]
}
rm(sitePrecip, dif, i, Inde1, Inde2, Inds1, Inds2)

#   compare site daily measurements----
# make plots comparing values of one gauge vs other at a site
for (i in seq(1,(length(BES_PrecipList)-1),2)){
  colPale = colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))
  png(paste0('Precip_SiteGaugeCompare_', BES_PrecipList[[i]]$Site[1], '.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = BES_Precip_d[[i+1]]$Precipitation_.mm., x = BES_Precip_d[[i]]$Precipitation_.mm., pch = 1, cex = .5,
       xlab = 'Rain Gauge 1 (mm)', 
       ylab = 'Rain Gauge 2 (mm)', 
       main = paste0('Station ', BES_PrecipList[[i]]$Site[1]), col = colPale(nrow(BES_Precip_d[[i]])))
  lines(c(-100,300),c(-100,300), col = 'red')
  dev.off()
}

# Fixme: Do some diagnostics on the sites to come up with a single timeseries for each site----
# Fixme: check if the high bias days are also biased on the previous and following day----
# Fixme: color the bias by the rainfall rate (> 1.5"/hr is prone to errors)
# NOAA weather station data----
#Get the FIPS code for Maryland
MDfips = rnoaa::fipscodes[which(rnoaa::fipscodes$state == 'Maryland'),]
#Download GHCND data - requires access code requested from: https://www.ncdc.noaa.gov/cdo-web/token
# This code takes a while to be emailed to you. AND one can only download a year at a time - not recommended
#NOAAstations = ncdc(datasetid='GHCND', locationid = paste0('FIPS:', MDfips$fips_state), startdate = '1998-01-01', enddate = '2018-12-31')
#For metadata
#rnoaa::homr(state = 'MD', county = 'Baltimore')
#Python alternative: https://k3.cicsnc.org/jared/GHCNpy

#Alternative method without need for authentication:
AllNOAAstations = ghcnd_stations()
#Make a spatial dataframe, and clip to ROI
coordinates(AllNOAAstations) = c('longitude', 'latitude')
proj4string(AllNOAAstations) = CRS('+init=epsg:4326')
#buffer the ROI a bit - assumes a UTM coordinate system. Will not work with other systems.
ROI_buffer = buffer(ROI, width = ROIbuffPrecip)
ROI_buffer_WGS = spTransform(ROI_buffer, CRS('+init=epsg:4326'))
NOAAstations_locs = AllNOAAstations[ROI_buffer_WGS,]
NOAAstations_locs = spTransform(NOAAstations_locs, CRS(pCRS))

#Download data and store in a list
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

for (i in 1:length(unique(NOAAstations_locs$id))){
  ghcnd(stationid = unique(NOAAstations_locs$id)[i], refresh = FALSE)
  m = meteo_tidy_ghcnd(stationid = unique(NOAAstations_locs$id)[i], var = 'all', keep_flags = TRUE)
  MetStations = c(MetStations, list(m))
}
names(MetStations) = unique(NOAAstations_locs$id)
rm(m)

#Move files from the cache to a permanent directory
wd_NOAA = paste0(dir_precip, '\\NOAA2')
dir.create(path = wd_NOAA)
file.copy(from = user_cache_dir(appname = 'rnoaa', version = NULL, opinion = TRUE, expand = TRUE), to = wd_NOAA, recursive = TRUE)
setwd(wd_NOAA)

#Select only the met stations in a 5 km buffer for now
ROI_buffer_5km = buffer(ROI, width = 5000)
ROI_buffer_5km_WGS = spTransform(ROI_buffer_5km, CRS('+init=epsg:4326'))
NOAAstations_5km_locs = AllNOAAstations[ROI_buffer_5km_WGS,]
NOAAstations_5km_locs = spTransform(NOAAstations_5km_locs, CRS(pCRS))

MetStations_5km = MetStations[(names(MetStations) %in% NOAAstations_5km_locs$id)]

#  Map of the type of data at each station----
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

#  Plot the DEM of the station vs. the DEM of the grid cell----
NOAAstations_locs$DEMelev = raster::extract(DEM, y = NOAAstations_locs)
png('CompareNOAAGaugeElevToDEM.png', res = 300, units = 'in', width = 5, height = 5)
plot(NOAAstations_locs$elevation, NOAAstations_locs$DEMelev,
     xlab = 'Reported Elevation (m)', ylab = 'DEM Elevation (m)', main = 'NOAA Gauge Elevations')
lines(c(-100,1100), c(-100,1100), col = 'red')
dev.off()

#  Fill in missing dates----
NOAAstations_5km_locs@data = as.data.frame(NOAAstations_5km_locs@data)
MetStations_5km_Precip = FillMissingDates_par(Dataset = NOAAstations_5km_locs[NOAAstations_5km_locs$element == 'PRCP',], StationList = MetStations_5km[names(MetStations_5km) %in% NOAAstations_5km_locs[NOAAstations_5km_locs$element == 'PRCP',]$id], Var = 'prcp', Date = 'date', gapType = 'd', site_no_D = 'id', site_no_SL = 'id', NoNAcols = 'id')
#Extract data from the function return
NOAAstations_5km_locs_Precip = MetStations_5km_Precip$Dataset
MetStations_5km_Precip = MetStations_5km_Precip$StationList

MetStations_5km_TMAX = FillMissingDates_par(Dataset = NOAAstations_5km_locs[NOAAstations_5km_locs$element == 'TMAX',], StationList = MetStations_5km[names(MetStations_5km) %in% NOAAstations_5km_locs[NOAAstations_5km_locs$element == 'TMAX',]$id], Var = 'tmax', Date = 'date', gapType = 'd', site_no_D = 'id', site_no_SL = 'id', NoNAcols = 'id')
#Extract data from the function return
NOAAstations_5km_locs_TMAX = MetStations_5km_TMAX$Dataset
MetStations_5km_TMAX = MetStations_5km_TMAX$StationList

MetStations_5km_TMIN = FillMissingDates_par(Dataset = NOAAstations_5km_locs[NOAAstations_5km_locs$element == 'TMIN',], StationList = MetStations_5km[names(MetStations_5km) %in% NOAAstations_5km_locs[NOAAstations_5km_locs$element == 'TMIN',]$id], Var = 'tmin', Date = 'date', gapType = 'd', site_no_D = 'id', site_no_SL = 'id', NoNAcols = 'id')
#Extract data from the function return
NOAAstations_5km_locs_TMIN = MetStations_5km_TMIN$Dataset
MetStations_5km_TMIN = MetStations_5km_TMIN$StationList

# Plot the met station timeseries----
#  Precip----
for (i in 1:length(MetStations_5km_Precip)){
  png(paste0('Precip_Timeseries_', MetStations_5km_Precip[[i]]$id[1], '.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = MetStations_5km_Precip[[i]]$prcp/10, x = as.Date(MetStations_5km_Precip[[i]]$date), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = 'Precipitation (mm)', 
       main = paste0('Station #', MetStations_5km_Precip[[i]]$id[1]), cex.lab = 1.5, cex.axis = 1.5)
  dev.off()
  
  #With map and timeseries
  png(paste0('Precip_Timeseries_Map', MetStations_5km_Precip[[i]]$id[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
  layout(rbind(c(1,2)))
  
  plot(y = MetStations_5km_Precip[[i]]$prcp/10, x = as.Date(MetStations_5km_Precip[[i]]$date), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = 'Precipitation (mm)', 
       main = paste0('Station #', MetStations_5km_Precip[[i]]$id[1]), cex.lab = 1.5, cex.axis = 1.5)
  
  #map
  plot(NOAAstations_locs, col = 'white')
  plot(ROI, add = T)
  #All temperature gauges
  plot(NOAAstations_5km_locs_Precip, pch = 16, col = 'black', add = TRUE)
  #Gauge selected for timeseries plot
  plot(NOAAstations_5km_locs_Precip[i,], pch = 16, col = 'blue', add = TRUE)
  
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4343000, len = 700, col = 'black', lab = 'N')
  legend('bottomleft', title = 'NOAA Precipitation Stations', legend = c('Selected', 'Other'), pch = 16, col = c('blue', 'black'), bty = 'n')
  title(main = paste0('Station #', MetStations_5km_Precip[[i]]$id[1]))
  dev.off()
}
rm(i)

#  Max and Min Temperature----
for (i in 1:length(MetStations_5km_TMAX)){
  png(paste0('Temp_Timeseries_', MetStations_5km_TMAX[[i]]$id[1], '.png'), res = 300, units = 'in', width = 6, height = 6)
  plot(y = MetStations_5km_TMAX[[i]]$tmax/10, x = as.Date(MetStations_5km_TMAX[[i]]$date), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = expression(paste('Temperature (', degree, 'C)')), 
       main = paste0('Station #', MetStations_5km_TMAX[[i]]$id[1]), col = 'red', ylim = c(-20,50), cex.lab = 1.5, cex.axis = 1.5)
  par(new = TRUE)
  plot(y = MetStations_5km_TMIN[[i]]$tmin/10, x = as.Date(MetStations_5km_TMIN[[i]]$date), type = 'o', pch = 16, cex = 0.3,
       xlab = '', ylab = '', axes = FALSE, 
       main = paste0('Station #', MetStations_5km_TMIN[[i]]$id[1]), col = 'blue', ylim = c(-20,50))
  dev.off()
  
  #With map and timeseries
  png(paste0('Temp_Timeseries_Map', MetStations_5km_TMAX[[i]]$id[1], '.png'), res = 300, units = 'in', width = 10, height = 10)
  layout(rbind(c(1,2)))
  
  plot(y = MetStations_5km_TMAX[[i]]$tmax/10, x = as.Date(MetStations_5km_TMAX[[i]]$date), type = 'o', pch = 16, cex = 0.3,
       xlab = 'Year', 
       ylab = expression(paste('Temperature (', degree, 'C)')), 
       main = paste0('Station #', MetStations_5km_TMAX[[i]]$id[1]), col = 'red', ylim = c(-20,50), cex.lab = 1.5, cex.axis = 1.5)
  par(new = TRUE)
  plot(y = MetStations_5km_TMIN[[i]]$tmin/10, x = as.Date(MetStations_5km_TMIN[[i]]$date), type = 'o', pch = 16, cex = 0.3,
       xlab = '', ylab = '', axes = FALSE, 
       main = paste0('Station #', MetStations_5km_TMIN[[i]]$id[1]), col = 'blue', ylim = c(-20,50))
  
  #map
  plot(NOAAstations_locs, col = 'white')
  plot(ROI, add = T)
  #All temperature gauges
  plot(NOAAstations_5km_locs_TMAX, pch = 16, col = 'black', add = TRUE)
  #Gauge selected for timeseries plot
  plot(NOAAstations_5km_locs_TMAX[i,], pch = 16, col = 'purple', add = TRUE)
  
  # Add coordinates
  axis(side = 1)
  axis(side = 2)
  box()
  north.arrow(xb = 370000, yb = 4343000, len = 700, col = 'black', lab = 'N')
  legend('bottomleft', title = 'Temperature Stations', legend = c('Selected', 'Other'), pch = 16, col = c('purple', 'black'), bty = 'n')
  title(main = paste0('Station #', MetStations_5km_TMAX[[i]]$id[1]))
  dev.off()
}
rm(i)

# Map of record lengths for precip and temperature stations----
BES_Precip_locs$RecordLength = BES_Precip_locs$RecordLengthMinusGaps = NA
for (i in 1:length(BES_PrecipList)){
  BES_Precip_locs$RecordLength[which(BES_Precip_locs$Gauge == BES_PrecipList[[i]]$Rain_Gauge_ID[1])] = as.numeric((max(BES_PrecipList[[i]]$SortDate) - min(BES_PrecipList[[i]]$SortDate)))
}
rm(i)
BES_Precip_locs$RecordLengthMinusGaps = BES_Precip_locs$RecordLength - BES_Precip_locs$MissingData_d
#in years
BES_Precip_locs$RecordLength = BES_Precip_locs$RecordLength/365.25
BES_Precip_locs$RecordLengthMinusGaps = BES_Precip_locs$RecordLengthMinusGaps/365.25

#Color by decades
scaleRange = c(0,10)
scaleBy = 2
Pal = rev(rainbow((scaleRange[2] - scaleRange[1])/scaleBy))

setwd(dir_precip)
png('BESPrecip_RecordLengths.png', res = 300, units = 'in', width = 6, height = 6)
plot(ROI)
# Gauges colored by their record lengths
plot(BES_Precip_locs, pch = 16, col = colFun(BES_Precip_locs$RecordLength), add = TRUE)
# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 370000, yb = 4346000, len = 700, col = 'black', lab = 'N')
legend('right', title = 'Precipitation Station \n Record Lengths \n (years)', legend = seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy), pch = 16, col = colFun(seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy)), bty = 'n')
dev.off()

scaleRange = c(0,10)
scaleBy = 1
Pal = rev(rainbow((scaleRange[2] - scaleRange[1])/scaleBy))

png('BESPrecip_RecordLengthsMinusGaps.png', res = 300, units = 'in', width = 6, height = 6)
plot(ROI)
# Gauges colored by their record lengths
plot(BES_Precip_locs, pch = 16, col = colFun(BES_Precip_locs$RecordLengthMinusGaps), add = TRUE)
# Add coordinates
axis(side = 1)
axis(side = 2)
box()
north.arrow(xb = 370000, yb = 4346000, len = 700, col = 'black', lab = 'N')
legend('right', title = 'Precipitation Station \n Record Lengths \n (years)', legend = seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy), pch = 16, col = colFun(seq(scaleRange[1], scaleRange[2]-scaleBy, scaleBy)), bty = 'n')
dev.off()

#Fixme: ACF for precip and temp datasets----
#Fixme: evaluate outliers for precip and temp (possible multivariate outliers) - make scatterplots of precip and temp----
#Write met station data to files----
#BES data
setwd(dir = dir_precip)
#Fixme: writing the whole list throws an error.
list.save(x = BES_PrecipList, file = f_BESPrecipList, type = "YAML")
#Saving as an RData file for now
list.save(x = BES_PrecipList, file = paste0(f_BESPrecipList, '.RData'), type = "RData")
list.save(x = BES_Precip_a, file = f_BESPrecipAggYear, type = "YAML")
list.save(x = BES_Precip_m, file = f_BESPrecipAggMonth, type = "YAML")
list.save(x = BES_Precip_d, file = f_BESPrecipAggDay, type = "YAML")
list.save(x = BES_Precip_Avg_d, file = f_BESPrecipAggDay_SiteAverage, type = "YAML")
#This YAML does not work now for some reason, so saving as RData file
list.save(x = BES_Precip_Avg_d, file = paste0(f_BESPrecipAggDay_SiteAverage, '.RData'), type = "RData")
writeOGR(BES_Precip_locs, dsn = getwd(), layer = f_BESprecipSites, driver = "ESRI Shapefile")
#NOAA data
setwd(dir = wd_NOAA)
writeOGR(NOAAstations_locs, dsn = getwd(), layer = f_NOAAstationsROI, driver = "ESRI Shapefile")
list.save(x = MetStations, file = f_NOAAstationsDataList, type = "YAML")

#Fixme: Spatial predicton of precip and temperature----
#Fixme: correlation plots of streamflow and N across watershed for different sites----
#       need to have matching timeseries to do this. Hydropairs function may work.
#       can help the spatial prediction
#Process hydro nexrad radar info----
#Initial test with year 2000
setwd(paste0(dir_Nexrad, "/Balto2000"))
Nexrad = read.table('200004040830.txt', header = FALSE)
Nexrad2 = read.table('200009040830.txt', header = FALSE)
colnames(Nexrad) = c('lat', 'long', 'precip')
colnames(Nexrad2) = c('lat', 'long', 'precip')
coordinates(Nexrad) = c('long', 'lat')
proj4string(Nexrad) = CRS('+init=epsg:4326')
Nexrad = spTransform(Nexrad, CRSobj = pCRS)
coordinates(Nexrad2) = c('long', 'lat')
proj4string(Nexrad2) = CRS('+init=epsg:4326')
Nexrad2 = spTransform(Nexrad2, CRSobj = pCRS)

#Check if the samples are in the same locations
any((Nexrad@coords == Nexrad2@coords) == FALSE)

#Plot the raster point (cell center) locations
plot(ROI)
plot(Nexrad, col = 'red', pch = 15, cex = 0.3, add = T)

BaispCRS = spTransform(BaisRun_Outlet, CRSobj = pCRS)
plot(BaispCRS)
plot(Nexrad, col = 'red', pch = 15, cex = 0.5, add = T)

#Extract the Baisman Run pixels
BaisPix = Nexrad[BaispCRS,]
BaisPix$Pix = c(1,2,3,4)

# Map of which pixel is which----
plot(BaispCRS)
plot(Nexrad, pch = 15, add = TRUE)
plot(BaisPix, col = c('red', 'orange', 'green', 'blue'), pch = 15, add = TRUE)

# Loop over the available NEXRAD txt files and extract these 4 pixels for Baisman Run----
#  Use 2007 as a test year----
#Store as a dataframe with the pixels as columns and the dates/times as rows
#setwd(paste0(dir_Nexrad, "/Balto2007"))
#fs2007 = list.files()

#Serial
#BaisNexMat = as.data.frame(matrix(NA, nrow = length(fs2007), ncol = 5))
#colnames(BaisNexMat) = c('Date', 'Pix1', 'Pix2', 'Pix3', 'Pix4')
# for (i in 1:length(fs2007)[1:1000]){
#   #Open the file and project into spatial dataframe
#   f = read.table(fs2007[i], header = FALSE)
#   colnames(f) = c('lat', 'long', 'precip')
#   coordinates(f) = c('long', 'lat')
#   proj4string(f) = CRS('+init=epsg:4326')
#   f = spTransform(f, CRSobj = pCRS)
#   
#   #Extract the Baisman Run pixels
#   f_BaisPix = f[BaispCRS,]
#   rm(f)
#   
#   #Extract the date and time information
#   Date = paste(paste(substr(fs2007[i], start = 1, stop = 4), substr(fs2007[i], start = 5, stop = 6), substr(fs2007[i], start = 7, stop = 8), sep = '-'), paste(substr(fs2007[i], start = 9, stop = 10), substr(fs2007[i], start = 11, stop = 12), sep = ':'), sep = ' ') 
#   
#   #Store the date and time as Posix
#   BaisNexMat$Date[i] = as.character(as.POSIXct(Date))
#   
#   #Store the precip info (mm/ha)
#   BaisNexMat[i, 2:5] = f_BaisPix$precip  
# }
# 
# #Parallel
# cl = makeCluster(detectCores() - 1)
# registerDoParallel(cl)
# BaisNexMat = foreach (i = 1:length(fs2007), .combine = rbind, .inorder = FALSE, .packages = 'sp') %dopar% {
#   #Open the file and project into spatial dataframe
#   f = read.table(fs2007[i], header = FALSE)
#   colnames(f) = c('lat', 'long', 'precip')
#   coordinates(f) = c('long', 'lat')
#   proj4string(f) = CRS('+init=epsg:4326')
#   f = spTransform(f, CRSobj = pCRS)
#   
#   #Extract the Baisman Run pixels
#   f_BaisPix = f[BaispCRS,]
#   rm(f)
#   
#   #Extract the date and time information
#   Date = paste(paste(substr(fs2007[i], start = 1, stop = 4), substr(fs2007[i], start = 5, stop = 6), substr(fs2007[i], start = 7, stop = 8), sep = '-'), paste(substr(fs2007[i], start = 9, stop = 10), substr(fs2007[i], start = 11, stop = 12), sep = ':'), sep = ' ') 
#   
#   #Store the date and time as Posix
#   Date = as.character(as.POSIXct(Date))
#   
#   #Store the precip info (mm/ha) and return as vector
#   retvals = c(Date, f_BaisPix$precip)
# }
# stopCluster(cl)
# rm(fs2007)
# 
# rownames(BaisNexMat) = NULL
# colnames(BaisNexMat) = c('Date', 'Pix1', 'Pix2', 'Pix3', 'Pix4')
# 
# #Fixme: why are there negative values in this processed product?
# #Give all negative values NAs
# BaisNexMat[BaisNexMat < 0] = NA

##Plot of the 4 Baisman pixels over time
#matplotDates(x = as.POSIXct(BaisNexMat[,1]), y = BaisNexMat[,2:5], type = 'l', col = c('red', 'orange', 'green', 'blue'), 
#        xlab = 'Time', ylab = 'Precipitation (mm/ha)')
##xlim = c(as.POSIXct(x = '2007-01-01', origin = '1970-01-01'), as.POSIXct('2007-12-31', origin = '1970-01-01'))

##Scatterplot Matrix
#plot(as.data.frame(BaisNexMat[,2:5]))

#  Parallel over all available Nexrad years----
setwd(dir_Nexrad)
folsNexrad = grep(list.files(), pattern = 'Balto', value = TRUE)
cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
NexradList = list()
for (j in 1:length(folsNexrad)){
  fs = list.files(paste0(getwd(), '/', folsNexrad[j]))
  #Get all of the nexrad data for this year.
  BaisNexMat = foreach (i = 1:length(fs), .combine = rbind, .inorder = FALSE, .packages = 'sp') %dopar% {
    #Open the file and project into spatial dataframe
    f = read.table(paste0(getwd(), '/', folsNexrad[j], '/', fs[i]), header = FALSE)
    colnames(f) = c('lat', 'long', 'precip')
    coordinates(f) = c('long', 'lat')
    proj4string(f) = CRS('+init=epsg:4326')
    f = spTransform(f, CRSobj = pCRS)
    
    #Extract the Baisman Run pixels
    f_BaisPix = f[BaispCRS,]
    rm(f)
    
    #Extract the date and time information
    Date = paste(paste(substr(fs[i], start = 1, stop = 4), substr(fs[i], start = 5, stop = 6), substr(fs[i], start = 7, stop = 8), sep = '-'), paste(substr(fs[i], start = 9, stop = 10), substr(fs[i], start = 11, stop = 12), sep = ':'), sep = ' ') 
    
    #Store the date and time as Posix
    Date = as.character(as.POSIXct(Date))
    
    #Store the precip info (mm/ha) and return as vector
    retvals = c(Date, f_BaisPix$precip)
  }
  
  rownames(BaisNexMat) = NULL
  colnames(BaisNexMat) = c('Date', 'Pix1', 'Pix2', 'Pix3', 'Pix4')
  
  #Save to list
  NexradList = c(NexradList, list(BaisNexMat))
  names(NexradList) = c(names(NexradList)[1:(j-1)], folsNexrad[j])
}
stopCluster(cl)

#Save the Nexrad product list
list.save(x = NexradList, file = f_NexradDataList, type = 'YAML')

#Go from a list to a matrix / data.frame
NexradMat = as.data.frame(NexradList[[1]])
NexradMat$Date = as.character(NexradMat$Date)
NexradMat[,2] = as.numeric(matrix(NexradMat[,2]))
NexradMat[,3] = as.numeric(matrix(NexradMat[,3]))
NexradMat[,4] = as.numeric(matrix(NexradMat[,4]))
NexradMat[,5] = as.numeric(matrix(NexradMat[,5]))

for(i in 2:length(NexradList)){
  NexradMat = rbind(NexradMat, NexradList[[i]])
}
NexradMat[,2] = as.numeric(matrix(NexradMat[,2]))
NexradMat[,3] = as.numeric(matrix(NexradMat[,3]))
NexradMat[,4] = as.numeric(matrix(NexradMat[,4]))
NexradMat[,5] = as.numeric(matrix(NexradMat[,5]))

rm(NexradList)

#Reassigned dates using this method
# NexradMat$Date[1:26600] = as.character(test[[1]][1:26600])
# NexradMat$Date[26601:(26600+25732)] = as.character(test[[2]][1:25732])
# NexradMat$Date[(26600+25732+1):(26600+25732+27020)] = as.character(test[[3]][1:27020])
# NexradMat$Date[(26600+25732+27020+1):(26600+25732+27020+30621)] = as.character(test[[4]][1:30621])
# NexradMat$Date[(26600+25732+27020+30621+1):(26600+25732+27020+30621+30648)] = as.character(test[[5]][1:30648])
# NexradMat$Date[(26600+25732+27020+30621+30648+1):(26600+25732+27020+30621+30648+30717)] = as.character(test[[6]][1:30717])
# NexradMat$Date[(26600+25732+27020+30621+30648+30717+1):(26600+25732+27020+30621+30648+30717+32732)] = as.character(test[[7]][1:32732])
# NexradMat$Date[(26600+25732+27020+30621+30648+30717+32732+1):(26600+25732+27020+30621+30648+30717+32732+33077)] = as.character(test[[8]][1:33077])
# NexradMat$Date[(26600+25732+27020+30621+30648+30717+32732+33077+1):(26600+25732+27020+30621+30648+30717+32732+33077+33727)] = as.character(test[[9]][1:33727])
# NexradMat$Date[(26600+25732+27020+30621+30648+30717+32732+33077+33727+1):(26600+25732+27020+30621+30648+30717+32732+33077+33727+34166)] = as.character(test[[10]][1:34166])

#   Check for duplicate records----
NexradDuplicates = which(duplicated(NexradMat))
#All of these are for the timestamp 2:00 - 2:45 AM. 
#I'm not sure why only a few dates have this problem. It does not affect further analysis, so leaving as is.

#   Check for negative values----
NexradMat$Pix1NegRm = NexradMat$Pix1
NexradMat$Pix2NegRm = NexradMat$Pix2
NexradMat$Pix3NegRm = NexradMat$Pix3
NexradMat$Pix4NegRm = NexradMat$Pix4

NexradMat$Pix1NegRm[NexradMat$Pix1NegRm < 0] = NA
NexradMat$Pix2NegRm[NexradMat$Pix2NegRm < 0] = NA
NexradMat$Pix3NegRm[NexradMat$Pix3NegRm < 0] = NA
NexradMat$Pix4NegRm[NexradMat$Pix4NegRm < 0] = NA

#   Aggregate to daily timeseries----
#Make a list of the pixels to use the aggregate Timeseries function
BaisPixList = list(NexradMat[,c(1,6)], NexradMat[,c(1,7)], NexradMat[,c(1,8)], NexradMat[,c(1,9)])
for(i in 1:length(BaisPixList)){
  BaisPixList[[i]]$Pix = i
  BaisPixList[[i]]$DateTime = BaisPixList[[i]]$Date
  colnames(BaisPixList[[i]]) = c("Date", 'Precip', 'Pix', "DateTime")
  BaisPixList[[i]]$Date = as.Date(BaisPixList[[i]]$Date)
}
BaisNex_Precip_agg = aggregateTimesteps(StationList = BaisPixList, aggVal = 'd', aggVar = 'Precip', date = 'Date', site = 'Pix', fun = 'sum')
BaisNex_Precip_d = BaisNex_Precip_agg$daily
rm(BaisNex_Precip_agg)

#   Fill in missing dates in the aggregated timeseries----
BaisNex_Precip_d2 = FillMissingDates_par(Dataset = BaisPix, StationList = BaisNex_Precip_d, Var = 'Precip', 
                                     Date = 'Date', gapType = 'd', site_no_D = 'Pix', 
                                     site_no_SL = 'Pix', NoNAcols = c('Pix'))
BaisPix_Processed = BaisNex_Precip_d2$Dataset
BaisNex_Precip_d = BaisNex_Precip_d2$StationList
rm(BaisNex_Precip_d2)

#   Convert units to mm/km^2 from mm/ha----
#1 ha = .01 km^2
for(i in 1:length(BaisNex_Precip_d)){
  BaisNex_Precip_d[[i]]$Precip/.01
}
rm(i)

# Compare the pixels with each other----
#make a matrix of the data to get a scatterplot matrix
BaisNexPrecipMat = matrix(NA, nrow = nrow(BaisNex_Precip_d[[1]]), ncol = 5)
BaisNexPrecipMat[,1] = BaisNex_Precip_d[[1]]$Date
BaisNexPrecipMat[,2] = BaisNex_Precip_d[[1]]$Precip
BaisNexPrecipMat[,3] = BaisNex_Precip_d[[2]]$Precip
BaisNexPrecipMat[,4] = BaisNex_Precip_d[[3]]$Precip
BaisNexPrecipMat[,5] = BaisNex_Precip_d[[4]]$Precip
BaisNexPrecipMat = as.data.frame(BaisNexPrecipMat)
colnames(BaisNexPrecipMat) = c('Date', paste0('Pix', seq(1,4,1)))

png('BaismanNexradPixPrecipTimeseries.png', res = 300, height = 5, width = 5, units = 'in')
par(mar = c(5,5,2,2))
matplotDates(as.Date(BaisNexPrecipMat[,1]), BaisNexPrecipMat[,-1], type = 'l', col = c('red', 'orange', 'green', 'blue'), 
             xlab = 'Time', ylab = expression(paste('Precipitation (mm/km'^2,')')))
dev.off()

png('BaismanNexradPixPrecipScatterPlotMatrix.png', res = 300, height = 8, width = 8, units = 'in')
pairs.panels(x = log10(BaisNexPrecipMat[,-1]+.01), scale = FALSE, density = FALSE, ellipses = FALSE, smooth = FALSE, 
             digits = 3, lm = TRUE, jiggle = FALSE, rug = FALSE, cex.cor = .7, method = 'spearman')
dev.off()

# Compare to the Oregon Ridge rain gauge----
min(as.Date(BES_Precip_d$WXORDG_RG1$SortDate))
max(as.Date(BaisNex_Precip_d[[1]]$Date))
IndStart = which(as.Date(BaisNex_Precip_d[[1]]$Date) == min(as.Date(BES_Precip_d$WXORDG_RG1$SortDate)))

#Plot the rain gauge vs. this gauge
BaisNex_Precip_d[[1]][IndStart:nrow(BaisNex_Precip_d[[1]]),]

BaisNexPrecipMat$ORGauge = NA 
BaisNexPrecipMat$ORGauge[IndStart:nrow(BaisNexPrecipMat)] = BES_Precip_Avg_d$`Oregon Ridge Park`$mean[as.Date(BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate) <= max(as.Date(BaisNex_Precip_d[[1]]$Date))]

png('BaismanNexradPixPrecipScatterPlotMatrix_WithRainGauge.png', res = 300, height = 8, width = 8, units = 'in')
pairs.panels(x = log10(BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),-1][which((BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),-1]$Pix1 > 0) & !is.nan(BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),-1]$ORGauge)),]), scale = FALSE, density = FALSE, ellipses = FALSE, smooth = FALSE, 
             digits = 3, lm = TRUE, jiggle = FALSE, rug = FALSE, cex.cor = .7, method = 'spearman', xlim = c(-2.5,2.5), ylim = c(-2.5, 2.5))
dev.off()


# Fixme: Fill in rasters with gauge values on days without radar information----
# Fixme: Convert to a raster for use in RHESSys patch precip assignment----
gridded(Nexrad) = TRUE
gridded(Nexrad2) = TRUE
NexradRas = raster(Nexrad)
NexradRas = projectRaster(NexradRas, crs = CRS(pCRS))

NexradRas2 = raster(Nexrad2)
NexradRas2 = projectRaster(NexradRas2, crs = CRS(pCRS))

#Fixme: Download climate indices (ENSO, PDO, etc.) and evaluate timeseries relative to those----
#Create Observation timeseries for RHESSys calibration----
#Baisman USGS Streamflow and TN----
setwd("C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\RHESSysFilePreparation\\obs")

#Streamflow and Nitrogen are both available in the BES_TN_d list
BES_TN_d_load = list.load('C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\WaterChemistry\\Nitrogen\\BES_TN_d.yaml', type = "YAML")

#Earliest possible start date for calibration
min(as.Date(BES_TN_d_load$BARN$SortDate, origin="1970-01-01"))

#Latest possible end date for validation
max(as.Date(BES_TN_d_load$BARN$SortDate, origin="1970-01-01"))

#Total number of samples
#Fixme: edit the hard coded numbers
length(which(!is.na(BES_TN_d_load$BARN$TN..mg.N.L.)))
#Want to have at least 15-20% of record for validation
876*.85
which(!is.na(BES_TN_d_load$BARN$TN..mg.N.L.))[743]
as.Date(BES_TN_d_load$BARN$SortDate[5365], origin="1970-01-01")

876*.80
which(!is.na(BES_TN_d_load$BARN$TN..mg.N.L.))[700]
as.Date(BES_TN_d_load$BARN$SortDate[5050], origin="1970-01-01")

#Calendar year end - about 18.4% of observations
as.Date(BES_TN_d_load$BARN$SortDate[5162], origin="1970-01-01")
length(which((!is.na(BES_TN_d_load$BARN$TN..mg.N.L.)) & (as.Date(BES_TN_d_load$BARN$SortDate, origin="1970-01-01") >= '2014-01-01')))/length(which((!is.na(BES_TN_d_load$BARN$TN..mg.N.L.))))

#Flow % about 18.7%
length(which(as.Date(BES_TN_d_load$BARN$SortDate, origin="1970-01-01") < '2014-01-01'))/length(BES_TN_d_load$BARN$SortDate)

#Water year peak on April 1 - about 17% of observations
as.Date(BES_TN_d_load$BARN$SortDate[5252], origin="1970-01-01")
length(which((!is.na(BES_TN_d_load$BARN$TN..mg.N.L.)) & (as.Date(BES_TN_d_load$BARN$SortDate, origin="1970-01-01") >= '2014-03-31')))/length(which((!is.na(BES_TN_d_load$BARN$TN..mg.N.L.))))

#Flow % about 17.25% of observations
1 - length(which(as.Date(BES_TN_d_load$BARN$SortDate, origin="1970-01-01") < '2014-03-31'))/length(BES_TN_d_load$BARN$SortDate)

#Streamflow - Calibration
StreamCal = data.frame(Date = as.character(as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01')[as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01') < "2013-10-01"]), Flow = BES_TN_d_load$BARN$Flow[as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01') < "2013-10-01"], stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = StreamCal, file = 'BaismanStreamflow_Cal.txt', sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)
#Streamflow - Validation
StreamVal = data.frame(Date = as.character(as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01')[as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01') >= "2008-11-15"]), Flow = BES_TN_d_load$BARN$Flow[as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01') >= "2008-11-15"], stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = StreamVal, file = 'BaismanStreamflow_Val.txt', sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)

#Streamflow - Calibration - Pond Branch
StreamCal_POBR = data.frame(Date = as.character(as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01')[as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01') < "2013-10-01"]), Flow = BES_TN_d_load$POBR$Flow[as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01') < "2013-10-01"], stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = StreamCal_POBR, file = 'BaismanStreamflow_POBR_Cal.txt', sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)
#Streamflow - Validation - Pond Branch
StreamVal_POBR = data.frame(Date = as.character(as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01')[as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01') >= "2008-11-15"]), Flow = BES_TN_d_load$POBR$Flow[as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01') >= "2008-11-15"], stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = StreamVal_POBR, file = 'BaismanStreamflow_POBR_Val.txt', sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)

#TN - Calibration
TNCal = data.frame(Date = as.character(as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01')[as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01') < "2013-10-01"]), TN = BES_TN_d_load$BARN$TN..mg.N.L.[as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01') < "2013-10-01"], stringsAsFactors = FALSE)
TNCal_WRTDS = data.frame(Date = as.character(as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01')[as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01') < "2013-10-01"]), remark = rep('', length(BES_TN_d_load$BARN$TN..mg.N.L.[as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01') < "2013-10-01"])), TN = BES_TN_d_load$BARN$TN..mg.N.L.[as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01') < "2013-10-01"], stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = TNCal, file = 'TN_Cal.txt', sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
write.table(x = TNCal_WRTDS, file = 'TN_Cal_WRTDS.txt', sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)
#TN - Validation
TNVal = data.frame(Date = as.character(as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01')[as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01') >= "2008-11-15"]), TN = BES_TN_d_load$BARN$TN..mg.N.L.[as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01') >= "2008-11-15"], stringsAsFactors = FALSE)
TNVal_WRTDS = data.frame(Date = as.character(as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01')[as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01') >= "2008-11-15"]), remark = rep('', length(BES_TN_d_load$BARN$TN..mg.N.L.[as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01') >= "2008-11-15"])), TN = BES_TN_d_load$BARN$TN..mg.N.L.[as.Date(BES_TN_d_load$BARN$SortDate, origin='1970-01-01') >= "2008-11-15"], stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = TNVal, file = 'TN_Val.txt', sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
write.table(x = TNVal_WRTDS, file = 'TN_Val_WRTDS.txt', sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)

#TN - Calibration - Pond Branch
TNCal_POBR = data.frame(Date = as.character(as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01')[as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01') < "2013-10-01"]), TN = BES_TN_d_load$POBR$TN..mg.N.L.[as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01') < "2013-10-01"], stringsAsFactors = FALSE)
TNCal_POBR_WRTDS = data.frame(Date = as.character(as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01')[as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01') < "2013-10-01"]), remark = rep('', length(BES_TN_d_load$POBR$TN..mg.N.L.[as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01') < "2013-10-01"])), TN = BES_TN_d_load$POBR$TN..mg.N.L.[as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01') < "2013-10-01"], stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = TNCal_POBR, file = 'TN_POBR_Cal.txt', sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
write.table(x = TNCal_POBR_WRTDS, file = 'TN_POBR_Cal_WRTDS.txt', sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)
#TN - Validation - Pond Branch
TNVal_POBR = data.frame(Date = as.character(as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01')[as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01') >= "2008-11-15"]), TN = BES_TN_d_load$POBR$TN..mg.N.L.[as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01') >= "2008-11-15"], stringsAsFactors = FALSE)
TNVal_POBR_WRTDS = data.frame(Date = as.character(as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01')[as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01') >= "2008-11-15"]), remark = rep('', length(BES_TN_d_load$POBR$TN..mg.N.L.[as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01') >= "2008-11-15"])), TN = BES_TN_d_load$POBR$TN..mg.N.L.[as.Date(BES_TN_d_load$POBR$SortDate, origin='1970-01-01') >= "2008-11-15"], stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = TNVal_POBR, file = 'TN_POBR_Val.txt', sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
write.table(x = TNVal_POBR_WRTDS, file = 'TN_POBR_Val_WRTDS.txt', sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)

#Create climate station timeseries for RHESSys----
#Rainfall for Oregon Ridge----
#Laurence's dataset
OR_Lin = read.table(file = 'C:\\Users\\js4yd\\Documents\\rhessys30m_Pond\\clim\\Oregon.rain', sep = '\t', stringsAsFactors = FALSE)

x = as.Date('2006-01-01')
x = x + seq(0, nrow(OR_Lin)-2, 1)
plot(x, OR_Lin[-1,], type ='l')


#Baisman Oregon Ridge vs. MD Science Center precip
plot(MetStations_5km$USW00093784[3938:7118,]$prcp/10, BES_Precip_Avg_d$`Oregon Ridge Park`$mean, log = 'xy')
lines(c(1,2000), c(1,2000), col = 'red')

#Compare to BES OR precip
plot(BES_Precip_Avg_d$`Oregon Ridge Park`$mean/1000, as.numeric(OR_Lin[-1,][1197:4377]), log = 'xy')
#Laurence used BES gauge 1, did not average the two gauges
plot(BES_Precip_d$WXORDG_RG1$Precipitation_.mm./1000, as.numeric(OR_Lin[-1,][1197:4377]), log = 'xy')

#Compare to MD Science Center USW00093784
plot(MetStations_5km$USW00093784[2742:3937,]$prcp/10000, as.numeric(OR_Lin[-1,][1:1196]), log = 'xy')

#Compare to BALTIMORE WASH INTL AP USW00093721
BWI = meteo_tidy_ghcnd(stationid = 'USW00093721', var = 'all', keep_flags = TRUE)
plot(BWI[24292:25487,]$prcp/10000, as.numeric(OR_Lin[-1,][1:1196]), log = 'xy')

#Compare BWI and MD Sci Center
plot(BWI[24292:25487,]$prcp/10000, MetStations_5km$USW00093784[2742:3937,]$prcp/10000, log = 'xy')

#Compare BES and BWI
plot(BES_Precip_d$WXORDG_RG1$Precipitation_.mm./1000, BWI[25488:(25488+3180),]$prcp/10000, log = 'xy')
#Compare BWI and MD Sci Center
plot(MetStations_5km$USW00093784[3938:(3938+3180),]$prcp/10000, BWI[25488:(25488+3180),]$prcp/10, log = 'xy', cex = 0.2, pch = 16)
plot(BWI[(22053):(22053+6697),]$prcp/10, MetStations_5km$USW00093784[503:(503+6697),]$prcp/10, pch = 16, cex = 0.2, log='xy')

#More correlated with the MD Science Center - Use that station to fill in before OR Ridge data are available.
RainTimeseries = MetStations_5km$USW00093784[503:3937,]$prcp/10

#Fill in NAs in the MD Science Center data with the BWI data. Most are zeros or a small amount of rain
RainTimeseries[which(is.na(MetStations_5km$USW00093784[503:3937,]$prcp))] = BWI$prcp[which(BWI$date %in% MetStations_5km$USW00093784[503:3937,]$date[which(is.na(MetStations_5km$USW00093784[503:3937,]$prcp))])]/10

#Add the corrected/fixed Oregon Ridge data to this timeseries
#Correct Oregon Ridge dataset
plot(BES_Precip_d$WXORDG_RG1$Precipitation_.mm., BES_Precip_d$WXORDG_RG2$Precipitation_.mm., ylim = c(0,150), xlim = c(0,150), xlab = 'Gauge 1', ylab = 'Gauge 2')
par(new = TRUE)
plot(BES_Precip_Avg_d$`Oregon Ridge Park`$P_WXORDG_RG1[which(abs(BES_Precip_Avg_d$`Oregon Ridge Park`$GaugeDiff) > 15)], BES_Precip_Avg_d$`Oregon Ridge Park`$P_WXORDG_RG2[which(abs(BES_Precip_Avg_d$`Oregon Ridge Park`$GaugeDiff) > 15)], col = 'red', ylim = c(0,150), xlim = c(0,150), xlab = '', ylab = '', axes = FALSE)
par(new = TRUE)
plot(BES_Precip_Avg_d$`Oregon Ridge Park`$P_WXORDG_RG1[which(abs(BES_Precip_Avg_d$`Oregon Ridge Park`$GaugeDiff) > 20)], BES_Precip_Avg_d$`Oregon Ridge Park`$P_WXORDG_RG2[which(abs(BES_Precip_Avg_d$`Oregon Ridge Park`$GaugeDiff) > 20)], col = 'blue', ylim = c(0,150), xlim = c(0,150), xlab = '', ylab = '', axes = FALSE)
#Cloud cover lines
lines(c(.0254*1000,.0254*1000), c(0,150))
lines(y = c(.0254*1000,.0254*1000), x = c(0,150))

#There are two instances of very high precip at one gauge and no precip at the other. 
# Several others differ by more than 20 mm. Use a buffer of "If greater than 20 mm different, take the larger instead of the mean"
BES_Precip_Avg_d$`Oregon Ridge Park`$Selected = BES_Precip_Avg_d$`Oregon Ridge Park`$mean
BES_Precip_Avg_d$`Oregon Ridge Park`$Selected[which(abs(BES_Precip_Avg_d$`Oregon Ridge Park`$GaugeDiff) > 20)] = apply(X = rbind(BES_Precip_Avg_d$`Oregon Ridge Park`$P_WXORDG_RG1[which(abs(BES_Precip_Avg_d$`Oregon Ridge Park`$GaugeDiff) > 20)], BES_Precip_Avg_d$`Oregon Ridge Park`$P_WXORDG_RG2[which(abs(BES_Precip_Avg_d$`Oregon Ridge Park`$GaugeDiff) > 20)]), MARGIN = 2, FUN = max)

#Check the NA rain dates with the MD Science Center precip. It's possible that some of these NA dates were days that the gauges did not record.
#BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[which(is.na(BES_Precip_Avg_d$`Oregon Ridge Park`$Selected))]
plot(MetStations_5km$USW00093784[3938:(3938+3180),]$date[MetStations_5km$USW00093784[3938:(3938+3180),]$date %in% BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[which(is.na(BES_Precip_Avg_d$`Oregon Ridge Park`$Selected))]],
     MetStations_5km$USW00093784[3938:(3938+3180),]$prcp[MetStations_5km$USW00093784[3938:(3938+3180),]$date %in% BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[which(is.na(BES_Precip_Avg_d$`Oregon Ridge Park`$Selected))]]/10)
#Any storm larger than 1 in of rain will be considered falling in OR Ridge, too.
FillInDates = MetStations_5km$USW00093784[3938:(3938+3180),]$date[MetStations_5km$USW00093784[3938:(3938+3180),]$date %in% BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[which(is.na(BES_Precip_Avg_d$`Oregon Ridge Park`$Selected))]][which(MetStations_5km$USW00093784[3938:(3938+3180),]$prcp[MetStations_5km$USW00093784[3938:(3938+3180),]$date %in% BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[which(is.na(BES_Precip_Avg_d$`Oregon Ridge Park`$Selected))]]/10 >= 25.4)]
FillInPrecips = MetStations_5km$USW00093784[3938:(3938+3180),]$prcp[MetStations_5km$USW00093784[3938:(3938+3180),]$date %in% BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[which(is.na(BES_Precip_Avg_d$`Oregon Ridge Park`$Selected))]][which(MetStations_5km$USW00093784[3938:(3938+3180),]$prcp[MetStations_5km$USW00093784[3938:(3938+3180),]$date %in% BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[which(is.na(BES_Precip_Avg_d$`Oregon Ridge Park`$Selected))]]/10 >= 25.4)]
BES_Precip_Avg_d$`Oregon Ridge Park`$Selected[which(BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate %in% FillInDates)] = FillInPrecips/10

#Make all other NA values equal to 0
BES_Precip_Avg_d$`Oregon Ridge Park`$Selected[is.na(BES_Precip_Avg_d$`Oregon Ridge Park`$Selected)] = 0

RainTimeseries = c(RainTimeseries, BES_Precip_Avg_d$`Oregon Ridge Park`$Selected)
RainTimeseries = RainTimeseries/1000
RainDates = seq(as.Date(MetStations_5km$USW00093784[503,]$date), max(as.Date(BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate)), 1) 

#Save text file for RHESSys input - in m instead of mm
#FullTimeseries
options(scipen = 999)
write.table(x = c("1999 11 15 1", RainTimeseries), file = 'AllTN.rain', sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote = FALSE)
options(scipen = 0)

#SA
#1999-11-15 to 2010-09-30
options(scipen = 999)
write.table(x = c("1999 11 15 1", RainTimeseries[1:which(RainDates == '2010-09-30')]), file = 'SA.rain', sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote = FALSE)
options(scipen = 0)

#Calibration
#1999-11-15 to 2013-09-30
options(scipen = 999)
write.table(x = c("1999 11 15 1", RainTimeseries[1:which(RainDates == '2013-09-30')]), file = 'Cal.rain', sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote = FALSE)
options(scipen = 0)

#Validation
#2008-11-15 to 2017-04-01
options(scipen = 999)
write.table(x = c("2008 11 15 1", RainTimeseries[which(RainDates == '2008-11-15'):which(RainDates == '2017-04-01')]), file = 'Val.rain', sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote = FALSE)
options(scipen = 0)

#Temperature min and max from MD Science Center and BWI Airport ----
#Laurence's dataset
OR_Lin_tmin = read.table(file = 'C:\\Users\\js4yd\\Documents\\rhessys30m_Pond\\clim\\Oregon.tmin', sep = '\t', stringsAsFactors = FALSE)
OR_Lin_tmax = read.table(file = 'C:\\Users\\js4yd\\Documents\\rhessys30m_Pond\\clim\\Oregon.tmax', sep = '\t', stringsAsFactors = FALSE)

x = as.Date('1996-04-01')
x = x + seq(0, nrow(OR_Lin_tmin)-2, 1)
plot(x, as.numeric(OR_Lin_tmin[-1,]), type ='l', col = 'blue', ylim = c(-20,50), xlab = 'Years', ylab = expression(paste('Temperature (',degree,'C)')))
par(new=TRUE)
plot(x, as.numeric(OR_Lin_tmax[-1,]), type ='l', col = 'red', ylim = c(-20,50), axes=FALSE, xlab = '', ylab = '')

#Compare to MD Science Center USW00093784
plot(MetStations_5km$USW00093784[503:(503+6697),]$tmin/10, as.numeric(OR_Lin_tmin[-1,][1324:(nrow(OR_Lin_tmin)-1)]))
plot(MetStations_5km$USW00093784[503:(503+6697),]$tmax/10, as.numeric(OR_Lin_tmax[-1,][1324:(nrow(OR_Lin_tmin)-1)]))

#Compare to BALTIMORE WASH INTL AP USW00093721
plot(BWI[20730:(20730+8020),]$tmin/10, as.numeric(OR_Lin_tmin[-1,]))
plot(BWI[(22053):(22053+6697),]$tmin/10, as.numeric(OR_Lin_tmin[-1,][1324:(nrow(OR_Lin_tmin)-1)]))

#Laurence used BWI
plot(BWI[20730:(20730+8020),]$tmax/10, as.numeric(OR_Lin_tmax[-1,]))
plot(BWI[(22053):(22053+6697),]$tmax/10, as.numeric(OR_Lin_tmax[-1,][1324:(nrow(OR_Lin_tmax)-1)]))

#Compare BWI and MD Sci Center
plot(BWI[(22053):(22053+6697),]$tmin/10, MetStations_5km$USW00093784[503:(503+6697),]$tmin/10, pch = 16, cex = 0.2)
plot(BWI[(22053):(22053+6697),]$tmax/10, MetStations_5km$USW00093784[503:(503+6697),]$tmax/10, pch = 16, cex = 0.2)

#Make temperature timeseries for MD Sci Center
TminTimeseries = MetStations_5km$USW00093784[503:(503+6615),]$tmin/10
TmaxTimeseries = MetStations_5km$USW00093784[503:(503+6615),]$tmax/10
TempDates = MetStations_5km$USW00093784[503:(503+6615),]$date

#Fill in the NA dates with the BWI airport temperatures
TminTimeseries[which(is.na(TminTimeseries))] = BWI$tmin[which(BWI$date %in% TempDates[which(is.na(TminTimeseries))])]/10
TmaxTimeseries[which(is.na(TmaxTimeseries))] = BWI$tmax[which(BWI$date %in% TempDates[which(is.na(TmaxTimeseries))])]/10

#Check that the min is less than the max
TmaxTimeseries[which((TmaxTimeseries - TminTimeseries) <0)]
TminTimeseries[which((TmaxTimeseries - TminTimeseries) <0)]
TempDates[which((TmaxTimeseries - TminTimeseries) <0)]
#Fill in the date of disagreement with BWI temperature
TminTimeseries[which((TmaxTimeseries - TminTimeseries) <0)] = BWI$tmin[which(BWI$date %in% TempDates[which((TmaxTimeseries - TminTimeseries) <0)])]/10
TmaxTimeseries[which((TmaxTimeseries - TminTimeseries) <0)] = BWI$tmax[which(BWI$date %in% TempDates[which((TmaxTimeseries - TminTimeseries) <0)])]/10

#Full Timeseries
options(scipen = 999)
write.table(x = c("1999 11 15 1", TminTimeseries), file = 'AllTN.tmin', sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
write.table(x = c("1999 11 15 1", TmaxTimeseries), file = 'AllTN.tmax', sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
options(scipen = 0)

#SA
#1999-11-15 to 2010-09-30
options(scipen = 999)
write.table(x = c("1999 11 15 1", TminTimeseries[1:which(TempDates == '2010-09-30')]), file = 'SA.tmin', sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
write.table(x = c("1999 11 15 1", TmaxTimeseries[1:which(TempDates == '2010-09-30')]), file = 'SA.tmax', sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
options(scipen = 0)

#Calibration
#1999-11-15 to 2013-09-30
options(scipen = 999)
write.table(x = c("1999 11 15 1", TminTimeseries[1:which(TempDates == '2013-09-30')]), file = 'Cal.tmin', sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
write.table(x = c("1999 11 15 1", TmaxTimeseries[1:which(TempDates == '2013-09-30')]), file = 'Cal.tmax', sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
options(scipen = 0)

#Validation
#2008-11-15 to 2017-04-01
options(scipen = 999)
write.table(x = c("2008 11 15 1", TminTimeseries[which(TempDates == '2008-11-15'):which(TempDates == '2017-04-01')]), file = 'Val.tmin', sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
write.table(x = c("2008 11 15 1", TmaxTimeseries[which(TempDates == '2008-11-15'):which(TempDates == '2017-04-01')]), file = 'Val.tmax', sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
options(scipen = 0)

#Elevation of BWI and MD Sci Center----
#No lapse rate correction made for these stations. 
#Fixme: hard coded numbers
#MD Sci Center
NOAAstations_locs@data[874,]
#BWI
NOAAstations_locs@data[871,]
#BES Oregon Ridge: 178.22 according to Laurence. 60 m screen height, which has been changed to 6 m.
#Using these elevations in calculations.


#CO2 from Mauna Loa - just for the general trend ----

#Nitrogen deposition from CAST sources----

#Baisman rainfall-runoff----
setwd(wd_sf)
#Find streamflow date that matches the precip date
Ind17 = which(StreamStationList$`01583570`$Date == BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[1])
Ind27 = which(StreamStationList$`01583570`$Date == BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[nrow(BES_Precip_Avg_d$`Oregon Ridge Park`)])

# Precip
png('RainfallRunoff_POBR.png', res = 300, units = 'in', width = 6, height = 6)
par(xaxs="i", yaxs="i", mar=c(5,5,5,5))
plot(as.Date(BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate), BES_Precip_Avg_d$`Oregon Ridge Park`$mean, type="h", ylim=c(max(BES_Precip_Avg_d$`Oregon Ridge Park`$mean, na.rm=TRUE)*1.5,0),
     axes=FALSE, xlab=NA, ylab=NA, col="blue",
     lwd=2, lend="square")
axis(4)
mtext("Precipitation (mm)", side=4, line=3)

# Streamflow
par(new=TRUE)
plot(as.Date(StreamStationList$`01583570`$Date[Ind17:Ind27]), StreamStationList$`01583570`$X_00060_00003[Ind1:Ind2], type="l", lwd=1, ylim=c(0, max(StreamStationList$`01583570`$X_00060_00003[Ind1:Ind2])*1.2),
     ylab = 'Streamflow (cfs)', xlab = 'Year', main = 'Pond Branch')
dev.off()

Ind1 = which(StreamStationList$`01583580`$Date == BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[1])
Ind2 = which(StreamStationList$`01583580`$Date == BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[nrow(BES_Precip_Avg_d$`Oregon Ridge Park`)])

# Precip
png('RainfallRunoff_BARN.png', res = 300, units = 'in', width = 6, height = 6)
par(xaxs="i", yaxs="i", mar=c(5,5,5,5))
plot(as.Date(BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate), BES_Precip_Avg_d$`Oregon Ridge Park`$mean, type="h", ylim=c(max(BES_Precip_Avg_d$`Oregon Ridge Park`$mean, na.rm=TRUE)*1.5,0),
     axes=FALSE, xlab=NA, ylab=NA, col="blue",
     lwd=2, lend="square")
axis(4)
mtext("Precipitation (mm)", side=4, line=3)

# Streamflow
par(new=TRUE)
plot(as.Date(StreamStationList$`01583580`$Date[Ind1:Ind2]), StreamStationList$`01583580`$X_00060_00003[Ind1:Ind2], type="l", lwd=1, ylim=c(0, max(StreamStationList$`01583580`$X_00060_00003[Ind1:Ind2])*1.2), 
     ylab = 'Streamflow (cfs)', xlab = 'Year', main = 'Baisman Run Outlet')
dev.off()

#Pond branch and outlet streamflow correlation
png('Streamflow_BARN&POBR.png', res = 300, units = 'in', width = 6, height = 6)
plot(StreamStationList$`01583580`$X_00060_00003[Ind1:Ind2], StreamStationList$`01583570`$X_00060_00003[Ind17:Ind27], 
     ylab = 'Streamflow Pond Branch (cfs)', xlab = 'Streamflow Baisman Run Outlet (cfs)')
text(x = 3, y = 2.5, paste('R^2 = ', round(cor(StreamStationList$`01583580`$X_00060_00003[Ind1:Ind2], StreamStationList$`01583570`$X_00060_00003[Ind17:Ind27])^2, 2)))
dev.off()

cor(StreamStationList$`01583580`$X_00060_00003[Ind1:Ind2], StreamStationList$`01583570`$X_00060_00003[Ind17:Ind27])

a = lm(StreamStationList$`01583580`$X_00060_00003[Ind1:Ind2] ~ StreamStationList$`01583570`$X_00060_00003[Ind17:Ind27])


#Estimate WRTDS interpolation tables----
#Fixme: Was the regression made with cms instead of cfs? Does it matter?
#Load the streamflow data into WRTDS format
Daily = readUserDaily(filePath = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\RHESSysFilePreparation\\obs", fileName = 'BaismanStreamflow_Cal.txt', hasHeader = TRUE, separator = '\t', qUnit = 1, verbose = FALSE)
#Read the TN data
Sample = readUserSample(filePath = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\RHESSysFilePreparation\\obs", fileName = 'TN_Cal_WRTDS.txt', hasHeader = TRUE, separator = '\t', verbose = FALSE)
#Set the required information
INFO = readUserInfo(filePath = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\RHESSysFilePreparation\\obs", fileName = 'WRTDS_INFO.csv', interactive = FALSE)
eList = mergeReport(INFO = INFO, Daily = Daily, Sample = Sample)
saveResults("C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\RHESSysFilePreparation\\obs\\", eList)

# Default WRTDS parameters----
WRTDSmod = modelEstimation(eList = eList, windowY = 7, windowQ = 2, windowS = .5, minNumObs = 100, minNumUncen = 50, edgeAdjust = TRUE)

setwd("C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\RHESSysFilePreparation\\obs\\")
png('ConcFluxTime.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod)
plotFluxTimeDaily(WRTDSmod)
dev.off()

png('ConcFluxTimeErrs.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod)
plotFluxPred(WRTDSmod)
dev.off()

plotResidPred(WRTDSmod)
plotResidQ(WRTDSmod)
plotResidTime(WRTDSmod)
boxResidMonth(WRTDSmod)
boxConcThree(WRTDSmod)
boxQTwice(WRTDSmod)
plotFluxHist(WRTDSmod)
plotConcHist(WRTDSmod)

png('BiasPlot.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,2.6,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

# Model#2-5 WRTDS----
WRTDSmod2 = modelEstimation(eList = eList, windowY = 7, windowQ = 2, windowS = .25, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE)
WRTDSmod3 = modelEstimation(eList = eList, windowY = 4, windowQ = 2, windowS = .25, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE)
WRTDSmod4 = modelEstimation(eList = eList, windowY = 2, windowQ = 2, windowS = .25, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE)
WRTDSmod5 = modelEstimation(eList = eList, windowY = 1.5, windowQ = 2, windowS = .25, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE)

png('ConcFluxTime_mod2.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod2)
plotFluxTimeDaily(WRTDSmod2)
dev.off()

png('ConcFluxTimeErrs_mod2.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod2)
plotFluxPred(WRTDSmod2)
dev.off()

plotResidPred(WRTDSmod2)
plotResidQ(WRTDSmod2)
plotResidTime(WRTDSmod2)
boxResidMonth(WRTDSmod2)
boxConcThree(WRTDSmod2)
boxQTwice(WRTDSmod2)
plotFluxHist(WRTDSmod2)
plotConcHist(WRTDSmod2)

png('BiasPlot_mod2.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod2, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_mod2.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod2, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,2.6,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr_mod2.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod2, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

png('ConcFluxTime_mod3.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod3)
plotFluxTimeDaily(WRTDSmod3)
dev.off()

png('ConcFluxTimeErrs_mod3.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod3)
plotFluxPred(WRTDSmod3)
dev.off()

plotResidPred(WRTDSmod3)
plotResidQ(WRTDSmod3)
plotResidTime(WRTDSmod3)
boxResidMonth(WRTDSmod3)
boxConcThree(WRTDSmod3)
boxQTwice(WRTDSmod3)
plotFluxHist(WRTDSmod3)
plotConcHist(WRTDSmod3)

png('BiasPlot_mod3.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod3, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_mod3.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod3, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,2.6,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr_mod3.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod3, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

png('ConcFluxTime_mod4.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod4)
plotFluxTimeDaily(WRTDSmod4)
dev.off()

png('ConcFluxTimeErrs_mod4.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod4)
plotFluxPred(WRTDSmod4)
dev.off()

plotResidPred(WRTDSmod4)
plotResidQ(WRTDSmod4)
plotResidTime(WRTDSmod4)
boxResidMonth(WRTDSmod4)
boxConcThree(WRTDSmod4)
boxQTwice(WRTDSmod4)
plotFluxHist(WRTDSmod4)
plotConcHist(WRTDSmod4)

png('BiasPlot_mod4.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod4, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_mod4.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod4, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,2.6,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr_mod4.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod4, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

png('ConcFluxTime_mod5.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod5)
plotFluxTimeDaily(WRTDSmod5)
dev.off()

png('ConcFluxTimeErrs_mod5.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod5)
plotFluxPred(WRTDSmod5)
dev.off()

plotResidPred(WRTDSmod5)
plotResidQ(WRTDSmod5)
plotResidTime(WRTDSmod5)
boxResidMonth(WRTDSmod5)
boxConcThree(WRTDSmod5)
boxQTwice(WRTDSmod5)
plotFluxHist(WRTDSmod5)
plotConcHist(WRTDSmod5)

png('BiasPlot_mod5.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod5, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_mod5.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod5, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,2.6,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr_mod5.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod5, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

#LOOCV SE for sample and discharge---
sum(WRTDSmod$Sample$SE^2)
sum(WRTDSmod2$Sample$SE^2)
sum(WRTDSmod3$Sample$SE^2)
sum(WRTDSmod4$Sample$SE^2)
sum(WRTDSmod5$Sample$SE^2)

sum(WRTDSmod$Daily$SE^2)
sum(WRTDSmod2$Daily$SE^2)
sum(WRTDSmod3$Daily$SE^2)
sum(WRTDSmod4$Daily$SE^2)
sum(WRTDSmod5$Daily$SE^2)

# MSE for sample----
mean((WRTDSmod$Sample$ConcHat - WRTDSmod$Sample$ConcAve))^2 + sum((WRTDSmod$Sample$ConcHat - WRTDSmod$Sample$ConcAve)^2)/(nrow(Sample)-1)
mean((WRTDSmod2$Sample$ConcHat - WRTDSmod2$Sample$ConcAve))^2 + sum((WRTDSmod2$Sample$ConcHat - WRTDSmod2$Sample$ConcAve)^2)/(nrow(Sample)-1)
mean((WRTDSmod3$Sample$ConcHat - WRTDSmod3$Sample$ConcAve))^2 + sum((WRTDSmod3$Sample$ConcHat - WRTDSmod3$Sample$ConcAve)^2)/(nrow(Sample)-1)
mean((WRTDSmod4$Sample$ConcHat - WRTDSmod4$Sample$ConcAve))^2 + sum((WRTDSmod4$Sample$ConcHat - WRTDSmod4$Sample$ConcAve)^2)/(nrow(Sample)-1)
mean((WRTDSmod5$Sample$ConcHat - WRTDSmod5$Sample$ConcAve))^2 + sum((WRTDSmod5$Sample$ConcHat - WRTDSmod5$Sample$ConcAve)^2)/(nrow(Sample)-1)


# Concentration and Discharge - shift after 2005?----
png('ConcDischarge.png', res = 300, units ='in', width = 5, height = 5)
par(mar=c(4,4,.5,.5))
layout(c(1,2))
plot(Daily$Date, log(Daily$Q), type = 'l', ylab = 'log(Flow [cms])', xlab = 'Year')
plot(Sample$Date, Sample$ConcAve, type = 'l', ylab = 'Total Nitrogen (mg/L)', xlab = 'Year')
dev.off()


# Running selected model with the modified functions that report the parameters of the surfaces----
setwd('C:\\Users\\js4yd\\OneDrive - University of Virginia\\RHESSys_ParameterSA')
WRTDSmod4m = modelEstimation(eList = eList, windowY = 2, windowQ = 2, windowS = .25, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE, numTsteps = 50, numQsteps = 100)

#Make interpolation tables from the surfaces information and save to files that can be loaded in
#Error in location 2, then locations 4-8 are Intercept, DecYear, LogQ, SinDY, CosDY
TabInt = WRTDSmod4m$surfaces[,,4]
TabYear = WRTDSmod4m$surfaces[,,5]
TabLogQ = WRTDSmod4m$surfaces[,,6]
TabSinYear = WRTDSmod4m$surfaces[,,7]
TabCosYear = WRTDSmod4m$surfaces[,,8]
TabLogErr = WRTDSmod4m$surfaces[,,2]

rownames(TabInt) = rownames(TabYear) = rownames(TabLogQ) = rownames(TabSinYear) = rownames(TabCosYear) = rownames(TabLogErr) = attr(WRTDSmod4m$surfaces, which = 'LogQ')
colnames(TabInt) = colnames(TabYear) = colnames(TabLogQ) = colnames(TabSinYear) = colnames(TabCosYear) = colnames(TabLogErr) = attr(WRTDSmod4m$surfaces, which = 'Year')

#There are 3 cells that have NA values that need to be adjusted using parameter interpolation from the tables.
IndReplace = which(is.na(TabInt))
#6 is one for each of the parameters from WRTDS
RepVal = matrix(0, nrow=length(IndReplace), ncol = 6)
for (ir in 1:length(IndReplace)){
  #Make a new temporary interpolation table that does not include the missing cell in it
  TempColInd = ceiling(IndReplace[ir]/nrow(TabInt))
  TempRowInd = which(rownames(TabInt) == names(which(is.na(TabInt[,TempColInd]))))
  
  #Check if this is the first column
  if (TempColInd == 1){
    #Use the parameters from the next column without NA values for both interpolation columns
    count = 1
    while (is.na(TabInt[TempRowInd,TempColInd+count])){
      count = count + 1
      if ((TempColInd+count) == (ncol(TabInt)+1)){
        stop(paste('No columns in row', TempRowInd, 'of the interpolation table have non-NA values.'))
      }
    }
    Col1 = Col2 = TempColInd+count
  }else if (TempColInd == ncol(TabInt)){
    #Used the second to last column's parameters for this column as well
    count = 1
    while (is.na(TabInt[TempRowInd,TempColInd-count])){
      count = count - 1
      if ((TempColInd-count) == 0){
        stop(paste('No columns in row', TempRowInd, 'of the interpolation table have non-NA values.'))
      }
    }
    Col1 = Col2 = TempColInd-count
  }else{
    #Use the columns on either side of the TempColInd in the new interpolation table.
    
    #Col1 first
    count = 1
    #Check on if column 2 exists
    Col2Exists = 0
    while (is.na(TabInt[TempRowInd,TempColInd-count])){
      count = count - 1
      if ((TempColInd-count) == 0){
        print(paste('Columns in row', TempRowInd, 'of the interpolation table are all NA values for columns less than the target column to be filled in. Checking the greater columns.'))
        count2 = 1
        while (is.na(TabInt[TempRowInd,TempColInd+count2])){
          count2 = count2 + 1
          if ((TempColInd+count2) == (ncol(TabInt)+1)){
            stop(paste('No columns in row', TempRowInd, 'of the interpolation table have non-NA values.'))
          }
        }
        #Found a value to use. Assign that column to both Col1 and Col2
        Col1 = Col2 = TempColInd+count2
        Col2Exists = 1
        rm(count2)
      }
      if ((TempColInd-count) == 0){
        #If it makes it here, it found a value. break out of the while
        break
      }
    }
    
    if (Col2Exists == 0){
      Col1 = TempColInd-count
      
      #Col2
      count = 1
      while (is.na(TabInt[TempRowInd,TempColInd+count])){
        count = count + 1
        if ((TempColInd+count) == (ncol(TabInt)+1)){
          print(paste('Columns in row', TempRowInd, 'of the interpolation table are all NA values for columns greater than the target column to be filled in. Using the column less than the target.'))
          Col2 = Col1
          break
        }
      }
      Col2 = TempColInd+count
    }
  }
  
  #Get the rows
  #Check if this is the first row
  if (TempRowInd == 1){
    #Use the parameters from the next row without NA values for both interpolation columns
    count = 1
    while (is.na(TabInt[TempRowInd+count, TempColInd])){
      count = count + 1
      if ((TempRowInd+count) == (nrow(TabInt)+1)){
        stop(paste('No rows in column', TempColInd, 'of the interpolation table have non-NA values.'))
      }
    }
    Row1 = Row2 = TempRowInd+count
  }else if (TempRowInd == nrow(TabInt)){
    #Used the second to last row's parameters for this row as well
    count = 1
    while (is.na(TabInt[TempRowInd-count,TempColInd])){
      count = count - 1
      if ((TempRowInd-count) == 0){
        stop(paste('No rows in column', TempColInd, 'of the interpolation table have non-NA values.'))
      }
    }
    Row1 = Row2 = TempRowInd-count
  }else{
    #Use the rows on either side of the TempRowInd in the new interpolation table.
    
    #Row1 first
    count = 1
    #Check on if Row2 exists
    Row2Exists = 0
    while (is.na(TabInt[TempRowInd-count,TempColInd])){
      count = count - 1
      if ((TempRowInd-count) == 0){
        print(paste('Rows in column', TempColInd, 'of the interpolation table are all NA values for rows less than the target row to be filled in. Checking the greater rows.'))
        count2 = 1
        while (is.na(TabInt[TempRowInd+count2,TempColInd])){
          count2 = count2 + 1
          if ((TempRowInd+count2) == (nrow(TabInt)+1)){
            stop(paste('No rows in column', TempColInd, 'of the interpolation table have non-NA values.'))
          }
        }
        #Found a value to use. Assign that column to both Col1 and Col2
        Row1 = Row2 = TempRowInd+count2
        Row2Exists = 1
        rm(count2)
      }
      if ((TempRowInd-count) == 0){
        #If it makes it here, it found a value. break out of the while
        break
      }
    }
    
    if (Row2Exists == 0){
      Row1 = TempRowInd-count
      
      #Col2
      count = 1
      while (is.na(TabInt[TempRowInd+count,TempColInd])){
        count = count + 1
        if ((TempRowInd+count) == (nrow(TabInt)+1)){
          print(paste('Rows in column', TempColInd, 'of the interpolation table are all NA values for rows greater than the target row to be filled in. Using the row less than the target.'))
          Row2 = Row1
          break
        }
      }
      Row2 = TempRowInd+count
    }
  }
  
  #Get new table using the corners of the tables that bound the NA values
  TempTabInt = TabInt[unique(c(Row1:TempRowInd, TempRowInd:Row2)), unique(c(Col1:TempColInd, TempColInd:Col2))]
  TempTabYear = TabYear[unique(c(Row1:TempRowInd, TempRowInd:Row2)), unique(c(Col1:TempColInd, TempColInd:Col2))]
  TempTabLogQ = TabLogQ[unique(c(Row1:TempRowInd, TempRowInd:Row2)), unique(c(Col1:TempColInd, TempColInd:Col2))]
  TempTabSinYear = TabSinYear[unique(c(Row1:TempRowInd, TempRowInd:Row2)), unique(c(Col1:TempColInd, TempColInd:Col2))]
  TempTabCosYear = TabCosYear[unique(c(Row1:TempRowInd, TempRowInd:Row2)), unique(c(Col1:TempColInd, TempColInd:Col2))]
  TempTabLogErr = TabLogErr[unique(c(Row1:TempRowInd, TempRowInd:Row2)), unique(c(Col1:TempColInd, TempColInd:Col2))]
  
  #Extract the new table values for the NA location using the available data in the tables that are not NA
  #Find the column and row index corresponding to the point
  
  RepVal[ir,] = FillTableNAs(DateInd = which(colnames(TempTabInt) == colnames(TabInt)[TempColInd]), FlowInd = which(rownames(TempTabInt) == rownames(TabInt)[TempRowInd]))
}
#Replace the table values
TabInt[IndReplace] = RepVal[,1]
TabYear[IndReplace] = RepVal[,2]
TabLogQ[IndReplace] = RepVal[,3]
TabSinYear[IndReplace] = RepVal[,4]
TabCosYear[IndReplace] = RepVal[,5]
TabLogErr[IndReplace] = RepVal[,6]

rm(ir, TempTabInt, TempTabCosYear, TempColInd, TempRowInd, TempTabLogErr, TempTabLogQ, TempTabSinYear, TempTabYear)
rm(Row2, Row1, Col1, Col2, count, Col2Exists, Row2Exists, RepVal, IndReplace)

#Write the tables
setwd('C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\WRTDS')
options(scipen = 999)
write.table(round(TabInt,5), file = 'TabIntMod4_p5.txt', sep = '\t', row.names = attr(WRTDSmod4m$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod4m$surfaces, which = 'Year'))
write.table(signif(TabYear,4), file = 'TabYearMod4_p4.txt', sep = '\t', row.names = attr(WRTDSmod4m$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod4m$surfaces, which = 'Year'))
write.table(signif(TabLogQ,4), file = 'TabLogQMod4_p4.txt', sep = '\t', row.names = attr(WRTDSmod4m$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod4m$surfaces, which = 'Year'))
write.table(signif(TabSinYear,4), file = 'TabSinYearMod4_p4.txt', sep = '\t', row.names = attr(WRTDSmod4m$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod4m$surfaces, which = 'Year'))
write.table(signif(TabCosYear,4), file = 'TabCosYearMod4_p4.txt', sep = '\t', row.names = attr(WRTDSmod4m$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod4m$surfaces, which = 'Year'))
write.table(round(TabLogErr,5), file = 'TabLogErrMod4_p5.txt', sep = '\t', row.names = attr(WRTDSmod4m$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod4m$surfaces, which = 'Year'))
options(scipen = 0)

#Contour plots of the parameters
png('TabInt.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4m, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2500,2100,100),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4m, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1.2,1.4,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4m, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1.8,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4m, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1.6,1.2,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4m, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1.4,1.8,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4m, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()
# Load WRTDS Interpolation Tables----
setwd("C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\WRTDS")
#Baisman
TabInt = as.matrix(read.table(file = 'TabIntMod4_p5.txt', sep = '\t', header = TRUE, check.names = FALSE))
TabYear = as.matrix(read.table(file = 'TabYearMod4_p4.txt', sep = '\t', header = TRUE, check.names = FALSE))
TabLogQ = as.matrix(read.table(file = 'TabLogQMod4_p4.txt', sep = '\t', header = TRUE, check.names = FALSE))
TabSinYear = as.matrix(read.table(file = 'TabSinYearMod4_p4.txt', sep = '\t', header = TRUE, check.names = FALSE))
TabCosYear = as.matrix(read.table(file = 'TabCosYearMod4_p4.txt', sep = '\t', header = TRUE, check.names = FALSE))
TabLogErr = as.matrix(read.table(file = 'TabLogErrMod4_p5.txt', sep = '\t', header = TRUE, check.names = FALSE))

#WRTDS for hillslopes----
# Data: Kenworth 2001-2002----
setwd(dir_SynWChem_Kenworth)
K_Q = read.csv(file = "Kenworth_Q.csv",stringsAsFactors = FALSE)
#Convert to m^3/s
K_Q[,-1] = K_Q[,-1]/1000
K_TN = read.csv(file = "Kenworth_TN.csv",stringsAsFactors = FALSE)
K_Sites = readOGR(dsn = getwd(), layer = "Kenworth_Chem_Sites", stringsAsFactors = FALSE)

# Data: Smith 2006-2007----
setwd(dir_SynWChem_Smith)
S_Q = read.csv(file = "Smith_discharge.csv",stringsAsFactors = FALSE)
S_TN = read.csv(file = "Smith_TN.csv",stringsAsFactors = FALSE)
S_Sites = readOGR(dsn = getwd(), layer = "Smith_Chem_Sites", stringsAsFactors = FALSE)

# Make regressions for the hillslopes that were sampled by both Kenworth and Smith----
#  BR3----
BR3_Q = K_Q[K_Q$Site == 'BR3',]
BR3_Q = cbind(BR3_Q, S_Q[S_Q$Site == 'BA3',-1])

BR3_TN = K_TN[K_TN$Site == 'BR3',]
BR3_TN = cbind(BR3_TN, S_TN[S_TN$Site == 'BA3',-1])

#   Relation using only the dates that match the BARN outlet sampling----
Dates_BR3 = colnames(BR3_TN)[which(!is.na(BR3_TN))][-1]
for(i in 1:length(Dates_BR3)){
  Dates_BR3[i] = strcat(strsplit(x = Dates_BR3[i], split = '.', fixed = TRUE)[[1]][-1], collapse = '-')
}

BR3_NSamps = which(as.Date(Dates_BR3) %in% as.Date(Sample$Date))

#   Relation using all available data for BR3----
Dates_BR3 = colnames(BR3_TN)[which(!is.na(BR3_TN))][-1]
for(i in 1:length(Dates_BR3)){
  Dates_BR3[i] = strcat(strsplit(x = Dates_BR3[i], split = '.', fixed = TRUE)[[1]][-1], collapse = '-')
}
Flows_BR3 = as.numeric(BR3_Q[which(!is.na(BR3_TN))][-1])

BR3_PredTN = matrix(NA, ncol = 3, nrow = length(Flows_BR3))
for(i in 1:length(Flows_BR3)){
  BR3_PredTN[i,] = predictWRTDS(Date = as.character(as.Date(Dates_BR3))[i], Flow = Flows_BR3[i], rowt = as.numeric(rownames(TabInt)), colt = as.numeric(colnames(TabInt)), TabInt = TabInt, TabYear = TabYear, TabLogQ = TabLogQ, TabSinYear = TabSinYear, TabCosYear = TabCosYear, TabLogErr = TabLogErr)
}
rm(i)
colnames(BR3_PredTN) = c('05', 'Med', '95')

BR3_PredTN = as.data.frame(BR3_PredTN)
BR3_PredTN$True = as.numeric(BR3_TN[which(!is.na(BR3_TN))][-1])
BR3_PredTN$DiffMed = BR3_PredTN$True - BR3_PredTN$Med

setwd(wd_BESN)
png('BR3_TrueTNvsWRTDS.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR3_PredTN$True, BR3_PredTN$Med, ylim = c(0,9), xlim = c(0,9), xlab = 'True BR3 TN (mg N/L)', ylab = 'WRTDS Predicted BR3 TN (mg N/L)')
lines(c(-10,10), c(-10,10))
dev.off()

plot(BR3_PredTN$True, BR3_PredTN$`05`)
plot(BR3_PredTN$True, BR3_PredTN$`95`)

#Date numbers in a leap year
DateNumMat = matrix(NA, nrow = 366, ncol = 2)
DateNumMat[,1] = seq(as.Date("2000-01-01"), as.Date("2000-12-31"), 1)
DateNumMat[,2] = seq(1, 366, 1)

DatesNums_BR3 = vector('numeric', length=length(Dates_BR3))
for(i in 1:length(Dates_BR3)){
  DatesNums_BR3[i] = DateNumMat[which((month(as.Date(DateNumMat[,1])) == month(as.Date(Dates_BR3[i]))) & (day(as.Date(DateNumMat[,1])) == day(as.Date(Dates_BR3[i])))), 2]
}

lm_BR3 = lm(BR3_PredTN$DiffMed ~ BR3_PredTN$Med)
lm_BR3_Dates = lm(BR3_PredTN$DiffMed ~ BR3_PredTN$Med + sin(2*pi*DatesNums_BR3/366) + cos(2*pi*DatesNums_BR3/366))

summary(lm_BR3)
plot(lm_BR3)

summary(lm_BR3_Dates)
plot(lm_BR3_Dates)

#plot of regression y vs. x
plot(y = BR3_PredTN$DiffMed, x = BR3_PredTN$Med)
plot(y = BR3_PredTN$DiffMed, x = Flows_BR3)
plot(y = BR3_PredTN$DiffMed, x = sin(2*pi*DatesNums_BR3/366))
plot(y = BR3_PredTN$DiffMed, x = cos(2*pi*DatesNums_BR3/366))

#   Relation dropping suspeccted outlier----
lm_BR3 = lm(BR3_PredTN$DiffMed[BR3_PredTN$True > 1] ~ BR3_PredTN$Med[BR3_PredTN$True > 1])
lm_BR3_Flows = lm(BR3_PredTN$DiffMed[BR3_PredTN$True > 1] ~ BR3_PredTN$Med[BR3_PredTN$True > 1]*Flows_BR3[BR3_PredTN$True > 1])
lm_BR3_Dates = lm(BR3_PredTN$DiffMed[BR3_PredTN$True > 1] ~ BR3_PredTN$Med[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))

summary(lm_BR3)
plot(lm_BR3)

summary(lm_BR3_Flows)
plot(lm_BR3_Flows)

summary(lm_BR3_Dates)
plot(lm_BR3_Dates)

#plot of regression y vs. x
plot(y = BR3_PredTN$DiffMed[BR3_PredTN$True > 1], x = BR3_PredTN$Med[BR3_PredTN$True > 1])
plot(y = BR3_PredTN$DiffMed[BR3_PredTN$True > 1], x = Flows_BR3[BR3_PredTN$True > 1])
plot(y = BR3_PredTN$DiffMed[BR3_PredTN$True > 1], x = sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))
plot(y = BR3_PredTN$DiffMed[BR3_PredTN$True > 1], x = cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))

#  BR5----
BR5_Q = K_Q[K_Q$Site == 'BR5A',]
BR5_Q = cbind(BR5_Q, S_Q[S_Q$Site == 'BA5AJC2',-1])

BR5_TN = K_TN[K_TN$Site == 'BR5A',]
BR5_TN = cbind(BR5_TN, S_TN[S_TN$Site == 'BA5AJC2',-1])

#   Relation using all of the dates that match the BARN outlet sampling----
Dates_BR5 = colnames(BR5_TN)[which(!is.na(BR5_TN))][-1]
for(i in 1:length(Dates_BR5)){
  Dates_BR5[i] = strcat(strsplit(x = Dates_BR5[i], split = '.', fixed = TRUE)[[1]][-1], collapse = '-')
}

BR5_NSamps = which(as.Date(Dates_BR5) %in% as.Date(Sample$Date))

#   Relation using all available data for BR5----
Dates_BR5 = colnames(BR5_TN)[which(!is.na(BR5_TN))][-1]
for(i in 1:length(Dates_BR5)){
  Dates_BR5[i] = strcat(strsplit(x = Dates_BR5[i], split = '.', fixed = TRUE)[[1]][-1], collapse = '-')
}
Flows_BR5 = as.numeric(BR5_Q[which(!is.na(BR5_TN))][-1])
#Drop the NA flow date
Dates_BR5 = Dates_BR5[-length(Dates_BR5)]
Flows_BR5 = Flows_BR5[-length(Flows_BR5)]

BR5_PredTN = matrix(NA, ncol = 3, nrow = length(Flows_BR5))
for(i in 1:length(Flows_BR5)){
  BR5_PredTN[i,] = predictWRTDS(Date = as.character(as.Date(Dates_BR5))[i], Flow = Flows_BR5[i], rowt = as.numeric(rownames(TabInt)), colt = as.numeric(colnames(TabInt)), TabInt = TabInt, TabYear = TabYear, TabLogQ = TabLogQ, TabSinYear = TabSinYear, TabCosYear = TabCosYear, TabLogErr = TabLogErr)
}
rm(i)
colnames(BR5_PredTN) = c('05', 'Med', '95')

BR5_PredTN = as.data.frame(BR5_PredTN)
BR5_PredTN$True = as.numeric(BR5_TN[which(!is.na(BR5_TN))][-1])[-ncol(BR5_TN[which(!is.na(BR5_TN))][-1])]
BR5_PredTN$DiffMed = BR5_PredTN$True - BR5_PredTN$Med

png('BR5_TrueTNvsWRTDS.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR5_PredTN$True, BR5_PredTN$Med, ylim = c(0,9), xlim = c(0,9), xlab = 'True BR5 TN (mg N/L)', ylab = 'WRTDS Predicted BR5 TN (mg N/L)')
lines(c(-10,10), c(-10,10))
dev.off()

plot(BR5_PredTN$True, BR5_PredTN$`05`)
plot(BR5_PredTN$True, BR5_PredTN$`95`)

#Date numbers in the year
DateNumMat = matrix(NA, nrow = 366, ncol = 2)
DateNumMat[,1] = seq(as.Date("2000-01-01"), as.Date("2000-12-31"), 1)
DateNumMat[,2] = seq(1, 366, 1)

DatesNums_BR5 = vector('numeric', length=length(Dates_BR5))
for(i in 1:length(Dates_BR5)){
  DatesNums_BR5[i] = DateNumMat[which((month(as.Date(DateNumMat[,1])) == month(as.Date(Dates_BR5[i]))) & (day(as.Date(DateNumMat[,1])) == day(as.Date(Dates_BR5[i])))), 2]
}

lm_BR5 = lm(BR5_PredTN$DiffMed ~ BR5_PredTN$Med)
lm_BR5_Dates = lm(BR5_PredTN$DiffMed ~ BR5_PredTN$Med + sin(2*pi*DatesNums_BR5/366) + cos(2*pi*DatesNums_BR5/366))

summary(lm_BR5)
plot(lm_BR5)

summary(lm_BR5_Dates)
plot(lm_BR5_Dates)

#plot of regression y vs. x
plot(y = BR5_PredTN$DiffMed, x = BR5_PredTN$Med)
plot(y = BR5_PredTN$DiffMed, x = sin(2*pi*DatesNums_BR5/366))
plot(y = BR5_PredTN$DiffMed, x = cos(2*pi*DatesNums_BR5/366))

#   Relation dropping suspeccted outlier----
lm_BR5 = lm(BR5_PredTN$DiffMed[BR5_PredTN$True > 1] ~ BR5_PredTN$Med[BR5_PredTN$True > 1])
lm_BR5_Flows = lm(BR5_PredTN$DiffMed[BR5_PredTN$True > 1] ~ BR5_PredTN$Med[BR5_PredTN$True > 1] + Flows_BR5[BR5_PredTN$True > 1])
lm_BR5_Dates = lm(BR5_PredTN$DiffMed[BR5_PredTN$True > 1] ~ BR5_PredTN$Med[BR5_PredTN$True > 1] + sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + Flows_BR5[BR5_PredTN$True > 1])

summary(lm_BR5)
plot(lm_BR5)

summary(lm_BR5_Flows)
plot(lm_BR5_Flows)

summary(lm_BR5_Dates)
plot(lm_BR5_Dates)

#plot of regression y vs. x
plot(y = BR5_PredTN$DiffMed[BR5_PredTN$True > 1], x = BR5_PredTN$Med[BR5_PredTN$True > 1], xlab = 'Predicted Mean TN (mg N/L)', ylab = 'Difference from True TN (mg N/L)')
plot(y = BR5_PredTN$DiffMed[BR5_PredTN$True > 1], x = Flows_BR5[BR5_PredTN$True > 1])
plot(y = BR5_PredTN$DiffMed[BR5_PredTN$True > 1], x = sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))
plot(y = BR5_PredTN$DiffMed[BR5_PredTN$True > 1], x = cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))

#  Pond Branch----
setwd(wd_BESN)
BES_TN_d = list.load(file = "BES_TN_d.yaml", type = 'YAML')

#   Relation using only the dates that match the BARN outlet sampling----
#length(which(as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]) %in% as.Date(BES_TN_d$BARN$SortDate[!is.na(BES_TN_d$BARN$TN..mg.N.L.)])))
#which(as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]) %in% as.Date(BES_TN_d$BARN$SortDate[!is.na(BES_TN_d$BARN$TN..mg.N.L.)]))
#which(as.Date(BES_TN_d$BARN$SortDate[!is.na(BES_TN_d$BARN$TN..mg.N.L.)]) %in% as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]))

#Plot of TN in POBR vs. BARN
plot(x = BES_TN_d$BARN$TN..mg.N.L.[!is.na(BES_TN_d$BARN$TN..mg.N.L.)][which(as.Date(BES_TN_d$BARN$SortDate[!is.na(BES_TN_d$BARN$TN..mg.N.L.)]) %in% as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]))],
     y = BES_TN_d$POBR$TN..mg.N.L.[!is.na(BES_TN_d$POBR$TN..mg.N.L.)][which(as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]) %in% as.Date(BES_TN_d$BARN$SortDate[!is.na(BES_TN_d$BARN$TN..mg.N.L.)]))],
     ylab = 'Pond Branch TN (mg N/L)', xlab = 'Baisman Run TN (mg N/L)', )

#Now use BARN WRTDS with POBR flows to predict POBR TN
#Get the POBR streamflows that match those dates. Only use the data up until 2014 because interplation tables stop at 2014.
Flows_POBR = BES_TN_d$POBR$Flow[!is.na(BES_TN_d$POBR$TN..mg.N.L.)][which((as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]) %in% as.Date(BES_TN_d$BARN$SortDate[!is.na(BES_TN_d$BARN$TN..mg.N.L.)])) & (as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]) <= '2014-01-01'))]
Dates_POBR = BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)][which((as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]) %in% as.Date(BES_TN_d$BARN$SortDate[!is.na(BES_TN_d$BARN$TN..mg.N.L.)])) & (as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]) <= '2014-01-01'))]

POBR_PredTN = matrix(NA, ncol = 3, nrow = length(Flows_POBR))
for(i in 1:length(Flows_POBR)){
  POBR_PredTN[i,] = predictWRTDS(Date = as.character(as.Date(Dates_POBR))[i], Flow = Flows_POBR[i]*12^3*2.54^3/100^3, rowt = as.numeric(rownames(TabInt)), colt = as.numeric(colnames(TabInt)), TabInt = TabInt, TabYear = TabYear, TabLogQ = TabLogQ, TabSinYear = TabSinYear, TabCosYear = TabCosYear, TabLogErr = TabLogErr)
}
rm(i)
colnames(POBR_PredTN) = c('05', 'Med', '95')

POBR_PredTN = as.data.frame(POBR_PredTN)
POBR_PredTN$True = BES_TN_d$POBR$TN..mg.N.L.[!is.na(BES_TN_d$POBR$TN..mg.N.L.)][which((as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]) %in% as.Date(BES_TN_d$BARN$SortDate[!is.na(BES_TN_d$BARN$TN..mg.N.L.)])) & (as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]) <= '2014-01-01'))]
POBR_PredTN$DiffMed = POBR_PredTN$True - POBR_PredTN$Med

plot(POBR_PredTN$True, POBR_PredTN$Med)
lines(c(-10,10), c(-10,10))
plot(POBR_PredTN$True, POBR_PredTN$`05`)
plot(POBR_PredTN$True, POBR_PredTN$`95`)

#Assign number values to the dates to be used in regression. 1 - 366
DateNumMat = matrix(NA, nrow = 366, ncol = 2)
DateNumMat[,1] = seq(as.Date("2000-01-01"), as.Date("2000-12-31"), 1)
DateNumMat[,2] = seq(1, 366, 1)

DatesNums_POBR = vector('numeric', length=length(Dates_POBR))
for(i in 1:length(Dates_POBR)){
  DatesNums_POBR[i] = DateNumMat[which((month(as.Date(DateNumMat[,1])) == month(as.Date(Dates_POBR[i]))) & (day(as.Date(DateNumMat[,1])) == day(as.Date(Dates_POBR[i])))), 2]
}

lm_POBR = lm(POBR_PredTN$DiffMed ~ I(Flows_POBR*12^3*2.54^3/100^3) + POBR_PredTN$Med + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))

summary(lm_POBR)
plot(lm_POBR)

#   Relation using all available data for Pond Branch----
Flows_POBR = BES_TN_d$POBR$Flow[which((!is.na(BES_TN_d$POBR$TN..mg.N.L.)) & (as.Date(BES_TN_d$POBR$SortDate) <= '2014-01-01'))]
Dates_POBR = as.Date(BES_TN_d$POBR$SortDate[which((!is.na(BES_TN_d$POBR$TN..mg.N.L.)) & (as.Date(BES_TN_d$POBR$SortDate) <= '2014-01-01'))])

POBR_PredTN = matrix(NA, ncol = 3, nrow = length(Flows_POBR))
for(i in 1:length(Flows_POBR)){
  POBR_PredTN[i,] = predictWRTDS(Date = as.character(as.Date(Dates_POBR))[i], Flow = Flows_POBR[i]*12^3*2.54^3/100^3, rowt = as.numeric(rownames(TabInt)), colt = as.numeric(colnames(TabInt)), TabInt = TabInt, TabYear = TabYear, TabCosYear = TabCosYear, TabSinYear = TabSinYear, TabLogQ = TabLogQ, TabLogErr = TabLogErr)
}
rm(i)
colnames(POBR_PredTN) = c('05', 'Med', '95')

POBR_PredTN = as.data.frame(POBR_PredTN)
POBR_PredTN$True = BES_TN_d$POBR$TN..mg.N.L.[which((!is.na(BES_TN_d$POBR$TN..mg.N.L.)) & (as.Date(BES_TN_d$POBR$SortDate) <= '2014-01-01'))]
POBR_PredTN$DiffMed = POBR_PredTN$True - POBR_PredTN$Med

png('POBR_TrueTNvsWRTDS.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN$True, POBR_PredTN$Med, ylim = c(0,9), xlim = c(0,9), xlab = 'True Pond Branch TN (mg N/L)', ylab = 'WRTDS Predicted Pond Branch TN (mg N/L)')
lines(c(-10,10), c(-10,10))
dev.off()

plot(POBR_PredTN$True, POBR_PredTN$`05`)
plot(POBR_PredTN$True, POBR_PredTN$`95`)

#Assign number values to the dates to be used in regression. 1 - 366
DateNumMat = matrix(NA, nrow = 366, ncol = 2)
DateNumMat[,1] = seq(as.Date("2000-01-01"), as.Date("2000-12-31"), 1)
DateNumMat[,2] = seq(1, 366, 1)

DatesNums_POBR = vector('numeric', length=length(Dates_POBR))
for(i in 1:length(Dates_POBR)){
  DatesNums_POBR[i] = DateNumMat[which((month(as.Date(DateNumMat[,1])) == month(as.Date(Dates_POBR[i]))) & (day(as.Date(DateNumMat[,1])) == day(as.Date(Dates_POBR[i])))), 2]
}

#Fixme - this needs to be a censored regression because of the y censoring. Thresholds change over time. Great.
lm_POBR = lm(POBR_PredTN$DiffMed ~ I(Flows_POBR*12^3*2.54^3/100^3) + POBR_PredTN$Med)
lm_POBR_Dates = lm(POBR_PredTN$DiffMed ~ I(Flows_POBR*12^3*2.54^3/100^3) + POBR_PredTN$Med + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))
lm_POBR_log = lm(log10(POBR_PredTN$DiffMed + abs(min(POBR_PredTN$DiffMed)) + .001) ~ POBR_PredTN$Med + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))

summary(lm_POBR)
plot(lm_POBR)

summary(lm_POBR_Dates)
plot(lm_POBR_Dates)

summary(lm_POBR_log)

#plot of regression y vs. x
plot(y = POBR_PredTN$DiffMed, x = POBR_PredTN$Med)
plot(y = POBR_PredTN$DiffMed, x = sin(2*pi*DatesNums_POBR/366))
plot(y = POBR_PredTN$DiffMed, x = cos(2*pi*DatesNums_POBR/366))

# Make regressions based on the TN load instead of the TN concentration----
#  Load at Basin Outlet----
#This load must result from the loadings at all other sites, minus the uptake by plants and organisms along the river. Ignoring uptake for now.
BES_TN_d$BARN$Load = BES_TN_d$BARN$TN..mg.N.L.*1000*BES_TN_d$BARN$Flow*12^3*2.54^3/100^3

#  Load at Pond Branch----
BES_TN_d$POBR$Load = BES_TN_d$POBR$TN..mg.N.L.*1000*BES_TN_d$POBR$Flow*12^3*2.54^3/100^3

#Make indicators for likely detection limits
BES_TN_d$POBR$LowTNLim[BES_TN_d$POBR$TN..mg.N.L. == min(BES_TN_d$POBR$TN..mg.N.L.[1:5000], na.rm = TRUE)] = 1
BES_TN_d$POBR$LowTNLim[BES_TN_d$POBR$TN..mg.N.L. == min(BES_TN_d$POBR$TN..mg.N.L.[5001:length(BES_TN_d$POBR$TN..mg.N.L.)], na.rm = TRUE)] = 2
POBR_PredTN$LowTNLim = NA
POBR_PredTN$LowTNLim[POBR_PredTN$True == min(POBR_PredTN$True[1:650])] = 1
POBR_PredTN$LowTNLim[POBR_PredTN$True == min(POBR_PredTN$True[651:length(POBR_PredTN$True)])] = 2

#Use WRTDS regression for predicted concentration and transform to loads. Compare predicted to true load.
POBR_PredTN$TrueLoad = POBR_PredTN$True*1000*Flows_POBR*12^3*2.54^3/100^3
POBR_PredTN$Load05 = POBR_PredTN$`05`*1000*Flows_POBR*12^3*2.54^3/100^3
POBR_PredTN$MedLoad = POBR_PredTN$Med*1000*Flows_POBR*12^3*2.54^3/100^3
POBR_PredTN$Load95 = POBR_PredTN$`95`*1000*Flows_POBR*12^3*2.54^3/100^3

png('POBR_TrueTNLoadvsWRTDSLoad.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN$TrueLoad, POBR_PredTN$MedLoad, ylim=range(c(POBR_PredTN$MedLoad-POBR_PredTN$Load05, POBR_PredTN$MedLoad+POBR_PredTN$Load95)), xlim = c(0, 50), xlab = 'True Pond Branch TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of Pond Branch TN (mg N/s)')
arrows(POBR_PredTN$TrueLoad, POBR_PredTN$MedLoad-POBR_PredTN$Load05, POBR_PredTN$TrueLoad, POBR_PredTN$MedLoad+POBR_PredTN$Load95, length=0.05, angle=90, code=3)
lines(c(-10,50), c(-10,50))
dev.off()

png('POBR_TrueTNLoadvsWRTDSLoad_log.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN$TrueLoad, POBR_PredTN$MedLoad, ylim=range(c(POBR_PredTN$MedLoad-POBR_PredTN$Load05, POBR_PredTN$MedLoad+POBR_PredTN$Load95)), xlab = 'True Pond Branch TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of Pond Branch TN (mg N/s)', log = 'xy')
arrows(POBR_PredTN$TrueLoad, POBR_PredTN$MedLoad-POBR_PredTN$Load05, POBR_PredTN$TrueLoad, POBR_PredTN$MedLoad+POBR_PredTN$Load95, length=0.05, angle=90, code=3)
lines(c(1e-7,1e7), c(1e-7,1e7))
dev.off()

png('POBR_TrueTNLoadvsWRTDSLoad_log_colors.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN$TrueLoad, POBR_PredTN$MedLoad, ylim=range(c(POBR_PredTN$MedLoad-POBR_PredTN$Load05, POBR_PredTN$MedLoad+POBR_PredTN$Load95)),
     xlim=range(POBR_PredTN$TrueLoad),
     xlab = 'True Pond Branch TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of Pond Branch TN (mg N/s)', log = 'xy')
par(new = TRUE)
plot(POBR_PredTN$TrueLoad[POBR_PredTN$LowTNLim == 1], POBR_PredTN$MedLoad[POBR_PredTN$LowTNLim == 1], 
     xlim=range(POBR_PredTN$TrueLoad),
     ylim=range(c(POBR_PredTN$MedLoad-POBR_PredTN$Load05, POBR_PredTN$MedLoad+POBR_PredTN$Load95)), xlab = '', ylab = '', 
     log = 'xy', col = 'blue', axes = FALSE)
par(new = TRUE)
plot(POBR_PredTN$TrueLoad[POBR_PredTN$LowTNLim == 2], POBR_PredTN$MedLoad[POBR_PredTN$LowTNLim == 2], 
     xlim=range(POBR_PredTN$TrueLoad),
     ylim=range(c(POBR_PredTN$MedLoad-POBR_PredTN$Load05, POBR_PredTN$MedLoad+POBR_PredTN$Load95)), xlab = '', ylab = '', 
     log = 'xy', col = 'red', axes = FALSE)
#arrows(POBR_PredTN$TrueLoad, POBR_PredTN$MedLoad-POBR_PredTN$Load05, POBR_PredTN$TrueLoad, POBR_PredTN$MedLoad+POBR_PredTN$Load95, length=0.05, angle=90, code=3)
lines(c(1e-7,1e7), c(1e-7,1e7))
legend('bottomright', title = 'Detection Limits (mg N/L)', legend = c('0.01', '0.05'), pch = 1, col = c('blue', 'red'))
dev.off()

#Add colors for seasons to this dataset
POBR_PredTN$cols = BES_TN_d$POBR$cols[as.Date(BES_TN_d$POBR$SortDate) %in% as.Date(Dates_POBR)]

png('POBR_TN_ScatterHist.png', res = 300, units = 'in', width = 10, height = 10)
scatterHist_mod(x = log10(POBR_PredTN$TrueLoad), y = log10(POBR_PredTN$MedLoad), density = FALSE, freq = FALSE, x.breaks = 100, y.breaks = 100, correl = FALSE, xlab = expression(paste('log'[10] ~ '(True TN Load [mg N/s])')), ylab = expression(paste('log'[10] ~ 'Predicted TN Load (mg N/s)')), cex.lab = 1.5, cex.axis = 1.5, title = '', pch = 16, ellipse = FALSE, smooth = FALSE, col = POBR_PredTN$cols)
dev.off()

POBR_PredTN$DiffMedLoad = POBR_PredTN$TrueLoad - POBR_PredTN$MedLoad

#plot of regression y vs. x
plot(y = POBR_PredTN$DiffMedLoad, x = POBR_PredTN$MedLoad)
plot(y = POBR_PredTN$DiffMedLoad, x = Flows_POBR*12^3*2.54^3/100^3)
plot(y = POBR_PredTN$DiffMedLoad, x = sin(2*pi*DatesNums_POBR/366))
plot(y = POBR_PredTN$DiffMedLoad, x = cos(2*pi*DatesNums_POBR/366))

#Fixme - this needs to be a censored regression because of the y censoring. Thresholds change over time. Great.
lm_POBRLoad = lm(POBR_PredTN$DiffMedLoad ~ I(Flows_POBR*12^3*2.54^3/100^3) + POBR_PredTN$MedLoad)
lm_POBRLoad_Dates = lm(POBR_PredTN$DiffMedLoad ~ I(Flows_POBR*12^3*2.54^3/100^3) + POBR_PredTN$MedLoad + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))
lm_POBRLoad_Dates_Interaction = lm(POBR_PredTN$DiffMedLoad ~ I(Flows_POBR*12^3*2.54^3/100^3)*POBR_PredTN$Med + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))
lm_POBRLoad_log = lm(log10(POBR_PredTN$DiffMedLoad + abs(min(POBR_PredTN$DiffMedLoad)) + .001) ~ POBR_PredTN$MedLoad + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))

summary(lm_POBRLoad)
plot(lm_POBRLoad)

summary(lm_POBRLoad_Dates)
plot(lm_POBRLoad_Dates)

summary(lm_POBRLoad_Dates_Interaction)
plot(lm_POBRLoad_Dates_Interaction)

summary(lm_POBRLoad_log)
plot(lm_POBRLoad_log)

#From the POBRLoad_Dates, use Box-Cox transform to attempt normality
BoxCoxtest = boxcox(I(POBR_PredTN$DiffMedLoad + 40) ~ I(Flows_POBR*12^3*2.54^3/100^3)*POBR_PredTN$Med + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366), lambda = seq(-2,2, length = 1000))
lambda = BoxCoxtest$x[which(BoxCoxtest$y == max(BoxCoxtest$y))]
POBR_PredTN$DiffMedLoad_BC = bcPower(U = POBR_PredTN$DiffMedLoad+40, gamma = 0, lambda = lambda, jacobian.adjusted = FALSE)
lm_POBRLoad_Dates_BC = lm((((POBR_PredTN$DiffMedLoad + 40)^lambda - 1)/lambda) ~ I(Flows_POBR*12^3*2.54^3/100^3)*POBR_PredTN$Med + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))

summary(lm_POBRLoad_Dates_BC)
plot(lm_POBRLoad_Dates_BC)

#plot of regression y vs. x
plot(y = POBR_PredTN$DiffMedLoad_BC, x = POBR_PredTN$MedLoad)
plot(y = POBR_PredTN$DiffMedLoad_BC, x = I(Flows_POBR*12^3*2.54^3/100^3))
plot(y = POBR_PredTN$DiffMedLoad_BC, x = sin(2*pi*DatesNums_POBR/366))
plot(y = POBR_PredTN$DiffMedLoad_BC, x = cos(2*pi*DatesNums_POBR/366))

#Color outlier
plot(POBR_PredTN$DiffMedLoad_BC, POBR_PredTN$MedLoad, xlim = c(0, 40), ylim = c(0,40))
par(new = T)
plot(POBR_PredTN$DiffMedLoad_BC[719], POBR_PredTN$MedLoad[719], col = 'red', xlim = c(0, 40), ylim = c(0,40), axes = FALSE, xlab ='', ylab = '')

#  Load at BR3----
#mg/s
BR3_Load = BR3_Q[-1]*BR3_TN[-1]*1000

#Use WRTDS regression for predicted concentration and transform to loads. Compare predicted to true load.
BR3_PredTN$TrueLoad = BR3_PredTN$True*1000*Flows_BR3*12^3*2.54^3/100^3
BR3_PredTN$Load05 = BR3_PredTN$`05`*1000*Flows_BR3*12^3*2.54^3/100^3
BR3_PredTN$MedLoad = BR3_PredTN$Med*1000*Flows_BR3*12^3*2.54^3/100^3
BR3_PredTN$Load95 = BR3_PredTN$`95`*1000*Flows_BR3*12^3*2.54^3/100^3

png('BR3_TrueTNLoadvsWRTDSLoad.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$MedLoad[BR3_PredTN$True > 1], ylim=range(c(BR3_PredTN$MedLoad[BR3_PredTN$True > 1]-BR3_PredTN$Load05[BR3_PredTN$True > 1], BR3_PredTN$MedLoad[BR3_PredTN$True > 1]+BR3_PredTN$Load95[BR3_PredTN$True > 1])), xlab = 'True BR3 TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of BR3 TN (mg N/s)')
arrows(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$MedLoad[BR3_PredTN$True > 1]-BR3_PredTN$Load05[BR3_PredTN$True > 1], BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$MedLoad[BR3_PredTN$True > 1]+BR3_PredTN$Load95[BR3_PredTN$True > 1], length=0.05, angle=90, code=3)
lines(c(-10,10), c(-10,10))
dev.off()

png('BR3_TrueTNLoadvsWRTDSLoad_log.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$MedLoad[BR3_PredTN$True > 1], ylim=range(c(BR3_PredTN$MedLoad[BR3_PredTN$True > 1]-BR3_PredTN$Load05[BR3_PredTN$True > 1], BR3_PredTN$MedLoad[BR3_PredTN$True > 1]+BR3_PredTN$Load95[BR3_PredTN$True > 1])), xlab = expression(paste('log'[10] ~ 'True BR3 TN Load (mg N/s)')), ylab = expression(paste('log'[10] ~ 'WRTDS Predicted Load of BR3 TN (mg N/s)')), log='xy')
arrows(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$MedLoad[BR3_PredTN$True > 1]-BR3_PredTN$Load05[BR3_PredTN$True > 1], BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$MedLoad[BR3_PredTN$True > 1]+BR3_PredTN$Load95[BR3_PredTN$True > 1], length=0.05, angle=90, code=3)
lines(c(1e-7,10), c(1e-7,10))
dev.off()

BR3_PredTN$DiffMedLoad = BR3_PredTN$TrueLoad - BR3_PredTN$MedLoad

#plot of regression y vs. x
plot(y = BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1], x = BR3_PredTN$MedLoad[BR3_PredTN$True > 1])
plot(y = BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1], x = I(Flows_BR3[BR3_PredTN$True > 1]*12^3*2.54^3/100^3))
plot(y = BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1], x = sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))
plot(y = BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1], x = cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))

#Regression w/o suspected outlier
lm_BR3Load = lm(BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]*12^3*2.54^3/100^3) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1])
lm_BR3Load_Dates = lm(BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]*12^3*2.54^3/100^3) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))
lm_BR3Load_Dates_Interaction = lm(BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]*12^3*2.54^3/100^3)*BR3_PredTN$Med[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))
lm_BR3Load_Interaction = lm(BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]*12^3*2.54^3/100^3)*BR3_PredTN$Med[BR3_PredTN$True > 1])
lm_BR3Load_log = lm(log10(BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1] + abs(min(BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1])) + .001) ~ BR3_PredTN$MedLoad[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))

summary(lm_BR3Load)
plot(lm_BR3Load)

summary(lm_BR3Load_Dates)
plot(lm_BR3Load_Dates)

summary(lm_BR3Load_Dates_Interaction)
plot(lm_BR3Load_Dates_Interaction)

summary(lm_BR3Load_Interaction)
plot(lm_BR3Load_Interaction)

summary(lm_BR3Load_log)
plot(lm_BR3Load_log)


#From the BR3Load_Dates, use Box-Cox transform to attempt normality
BoxCoxtest_BR3 = boxcox(BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]*12^3*2.54^3/100^3) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366), lambda = seq(0,1, length = 1000))
lambda_BR3 = BoxCoxtest_BR3$x[which(BoxCoxtest_BR3$y == max(BoxCoxtest_BR3$y))]

BR3_PredTN$DiffMedLoad_BC[BR3_PredTN$True > 1] = bcPower(U = BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1], gamma = 0, lambda = lambda_BR3, jacobian.adjusted = FALSE)
lm_BR3Load_Dates_BC = lm(BR3_PredTN$DiffMedLoad_BC[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]*12^3*2.54^3/100^3) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))

summary(lm_BR3Load_Dates_BC)
plot(lm_BR3Load_Dates_BC)

#plot of regression y vs. x
plot(y = BR3_PredTN$DiffMedLoad_BC, x = BR3_PredTN$MedLoad)
plot(y = BR3_PredTN$DiffMedLoad_BC, x = I(Flows_BR3*12^3*2.54^3/100^3))
plot(y = BR3_PredTN$DiffMedLoad_BC, x = sin(2*pi*DatesNums_BR3/366))
plot(y = BR3_PredTN$DiffMedLoad_BC, x = cos(2*pi*DatesNums_BR3/366))

#  Load at BR5----
#mg/s
BR5_Load = BR5_Q[-1]*BR5_TN[-1]*1000

#Use WRTDS regression for predicted concentration and transform to loads. Compare predicted to true load.
BR5_PredTN$TrueLoad = BR5_PredTN$True*1000*Flows_BR5*12^3*2.54^3/100^3
BR5_PredTN$Load05 = BR5_PredTN$`05`*1000*Flows_BR5*12^3*2.54^3/100^3
BR5_PredTN$MedLoad = BR5_PredTN$Med*1000*Flows_BR5*12^3*2.54^3/100^3
BR5_PredTN$Load95 = BR5_PredTN$`95`*1000*Flows_BR5*12^3*2.54^3/100^3

png('BR5_TrueTNLoadvsWRTDSLoad.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR5_PredTN$TrueLoad, BR5_PredTN$MedLoad, ylim=range(c(BR5_PredTN$MedLoad-BR5_PredTN$Load05, BR5_PredTN$MedLoad+BR5_PredTN$Load95)), xlab = 'True BR5 TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of BR5 TN (mg N/s)')
arrows(BR5_PredTN$TrueLoad, BR5_PredTN$MedLoad-BR5_PredTN$Load05, BR5_PredTN$TrueLoad, BR5_PredTN$MedLoad+BR5_PredTN$Load95, length=0.05, angle=90, code=3)
lines(c(-10,10), c(-10,10))
dev.off()

png('BR5_TrueTNLoadvsWRTDSLoad_log.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], BR5_PredTN$MedLoad[BR5_PredTN$True > 1], ylim=range(c(BR5_PredTN$MedLoad[BR5_PredTN$True > 1]-BR5_PredTN$Load05[BR5_PredTN$True > 1], BR5_PredTN$MedLoad[BR5_PredTN$True > 1]+BR5_PredTN$Load95[BR5_PredTN$True > 1])), xlab = expression(paste('log'[10] ~ 'True BR5 TN Load (mg N/s)')), ylab = expression(paste('log'[10] ~ 'WRTDS Predicted Load of BR5 TN (mg N/s)')), log='xy')
arrows(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], BR5_PredTN$MedLoad[BR5_PredTN$True > 1]-BR5_PredTN$Load05[BR5_PredTN$True > 1], BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], BR5_PredTN$MedLoad[BR5_PredTN$True > 1]+BR5_PredTN$Load95[BR5_PredTN$True > 1], length=0.05, angle=90, code=3)
lines(c(1e-7,10), c(1e-7,10))
dev.off()

BR5_PredTN$DiffMedLoad = BR5_PredTN$TrueLoad - BR5_PredTN$MedLoad

#plot of regression y vs. x
plot(y = BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1], x = BR5_PredTN$MedLoad[BR5_PredTN$True > 1])
plot(y = BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1], x = I(Flows_BR5[BR5_PredTN$True > 1]*12^3*2.54^3/100^3))
plot(y = BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1], x = sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))
plot(y = BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1], x = cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))

#Regression w/o suspected outlier
lm_BR5Load = lm(BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1]*12^3*2.54^3/100^3) + BR5_PredTN$MedLoad[BR5_PredTN$True > 1])
lm_BR5Load_Interaction = lm(BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1]*12^3*2.54^3/100^3)*BR5_PredTN$MedLoad[BR5_PredTN$True > 1])
lm_BR5Load_Dates = lm(BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1]*12^3*2.54^3/100^3) + BR5_PredTN$MedLoad[BR5_PredTN$True > 1] + sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))
lm_BR5Load_Dates_Interaction = lm(BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1]*12^3*2.54^3/100^3)*BR5_PredTN$Med[BR5_PredTN$True > 1] + sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))
lm_BR5Load_log = lm(log10(BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1] + abs(min(BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1])) + .001) ~ BR5_PredTN$MedLoad[BR5_PredTN$True > 1] + sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))

summary(lm_BR5Load)
plot(lm_BR5Load)

summary(lm_BR5Load_Interaction)
plot(lm_BR5Load_Interaction)

summary(lm_BR5Load_Dates)
plot(lm_BR5Load_Dates)

summary(lm_BR5Load_Dates_Interaction)
plot(lm_BR5Load_Dates_Interaction)

summary(lm_BR5Load_log)
plot(lm_BR5Load_log)


# Fraction of Baisman Outlet Load----
ProportionPond = BES_TN_d$POBR$Load[which(as.Date(BES_TN_d$POBR$SortDate) %in% as.Date(BES_TN_d$BARN$SortDate))]/BES_TN_d$BARN$Load[which(as.Date(BES_TN_d$BARN$SortDate) %in% as.Date(BES_TN_d$POBR$SortDate))]
boxplot(ProportionPond)

ProportionBR3 = BR3_Load[which(as.Date(Dates_BR3) %in% as.Date(BES_TN_d$BARN$SortDate))]/BES_TN_d$BARN$Load[which(as.Date(BES_TN_d$BARN$SortDate) %in% as.Date(Dates_BR3))]
plot(x = rep(1, length(as.numeric(ProportionBR3))), y = as.numeric(ProportionBR3))

ProportionBR5 = BR5_Load[which(as.Date(Dates_BR5) %in% as.Date(BES_TN_d$BARN$SortDate))]/BES_TN_d$BARN$Load[which(as.Date(BES_TN_d$BARN$SortDate) %in% as.Date(Dates_BR5))]
plot(x = rep(1, length(as.numeric(ProportionBR5))), y = as.numeric(ProportionBR5))

#Sum of the proportions for the matching dates
ProportionSumPond = BES_TN_d$POBR$Load[which((as.Date(BES_TN_d$POBR$SortDate) %in% as.Date(BES_TN_d$BARN$SortDate)) & (as.Date(BES_TN_d$POBR$SortDate) %in% as.Date(Dates_BR5[c(1,2,3,4)])))]/BES_TN_d$BARN$Load[which((as.Date(BES_TN_d$BARN$SortDate) %in% as.Date(BES_TN_d$POBR$SortDate)) & (as.Date(BES_TN_d$BARN$SortDate) %in% as.Date(Dates_BR5[c(1,2,3,4)])))]

ProportionBR5[c(2,3,4,5)] + ProportionBR3[c(2,3,4,5)] + ProportionSumPond


# Regressions to predict load directly instead of predict the correction----
#Fixme: ensure that regression loads are non-negative for all positive flows
#  Pond Branch - not as good as above----
#plot of regression y vs. x
plot(y = POBR_PredTN$TrueLoad, x = POBR_PredTN$MedLoad)
plot(y = POBR_PredTN$TrueLoad, x = Flows_POBR*12^3*2.54^3/100^3)
plot(y = POBR_PredTN$TrueLoad, x = sin(2*pi*DatesNums_POBR/366))
plot(y = POBR_PredTN$TrueLoad, x = cos(2*pi*DatesNums_POBR/366))

#Fixme - this needs to be a censored regression because of the y censoring. Thresholds change over time. Great.
lm_POBRLoad = lm(POBR_PredTN$TrueLoad ~ I(Flows_POBR*12^3*2.54^3/100^3) + POBR_PredTN$MedLoad)
lm_POBRLoad_Dates = lm(POBR_PredTN$TrueLoad ~ I(Flows_POBR*12^3*2.54^3/100^3) + POBR_PredTN$MedLoad + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))
lm_POBRLoad_Dates_Interaction = lm(POBR_PredTN$TrueLoad ~ I(Flows_POBR*12^3*2.54^3/100^3)*POBR_PredTN$Med + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))
lm_POBRLoad_log = lm(log10(POBR_PredTN$TrueLoad) ~ POBR_PredTN$MedLoad + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))

summary(lm_POBRLoad)
plot(lm_POBRLoad)

summary(lm_POBRLoad_Dates)
plot(lm_POBRLoad_Dates)

summary(lm_POBRLoad_Dates_Interaction)
plot(lm_POBRLoad_Dates_Interaction)

summary(lm_POBRLoad_log)
plot(lm_POBRLoad_log)

#  BR3 - a little better than above----
#plot of regression y vs. x
plot(y = BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], x = BR3_PredTN$MedLoad[BR3_PredTN$True > 1])
plot(y = BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], x = Flows_BR3[BR3_PredTN$True > 1]*12^3*2.54^3/100^3)
plot(y = BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], x = sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))
plot(y = BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], x = cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))

#Regression w/o suspected outlier
lm_BR3Load = lm(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]*12^3*2.54^3/100^3) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1])
lm_BR3Load_Interaction = lm(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]*12^3*2.54^3/100^3)*BR3_PredTN$MedLoad[BR3_PredTN$True > 1])
lm_BR3Load_Dates = lm(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]*12^3*2.54^3/100^3) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))
lm_BR3Load_Dates_Interaction = lm(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]*12^3*2.54^3/100^3)*BR3_PredTN$Med[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))
lm_BR3Load_log = lm(log10(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1]) ~ BR3_PredTN$MedLoad[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))

summary(lm_BR3Load)
plot(lm_BR3Load)

summary(lm_BR3Load_Interaction)
plot(lm_BR3Load_Interaction)

summary(lm_BR3Load_Dates)
plot(lm_BR3Load_Dates)

summary(lm_BR3Load_Dates_Interaction)
plot(lm_BR3Load_Dates_Interaction)

summary(lm_BR3Load_log)
plot(lm_BR3Load_log)


#From the BR3Load, use Box-Cox transform to attempt normality
BoxCoxtest_BR3_True = boxcox(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]*12^3*2.54^3/100^3) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1], lambda = seq(0,1, length = 1000))
lambda_BR3_True = BoxCoxtest_BR3_True$x[which(BoxCoxtest_BR3_True$y == max(BoxCoxtest_BR3_True$y))]

BR3_PredTN$DiffMedLoad_BC_True[BR3_PredTN$True > 1] = bcPower(U = BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], gamma = 0, lambda = lambda_BR3_True, jacobian.adjusted = FALSE)
lm_BR3Load_Dates_BC_True = lm(BR3_PredTN$DiffMedLoad_BC_True[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]*12^3*2.54^3/100^3) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1])

summary(lm_BR3Load_Dates_BC_True)
plot(lm_BR3Load_Dates_BC)

#plot of regression y vs. x
plot(y = BR3_PredTN$DiffMedLoad_BC_True, x = BR3_PredTN$MedLoad)
plot(y = BR3_PredTN$DiffMedLoad_BC_True, x = I(Flows_BR3*12^3*2.54^3/100^3))
plot(y = BR3_PredTN$DiffMedLoad_BC_True, x = sin(2*pi*DatesNums_BR3/366))
plot(y = BR3_PredTN$DiffMedLoad_BC_True, x = cos(2*pi*DatesNums_BR3/366))

#  BR5 - essentially same as above----
#plot of regression y vs. x
plot(y = BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], x = BR5_PredTN$MedLoad[BR5_PredTN$True > 1])
plot(y = BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], x = Flows_BR5[BR5_PredTN$True > 1]*12^3*2.54^3/100^3)
plot(y = BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], x = sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))
plot(y = BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], x = cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))

#Regression w/o suspected outlier
lm_BR5Load = lm(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1]*12^3*2.54^3/100^3) + BR5_PredTN$MedLoad[BR5_PredTN$True > 1])
lm_BR5Load_Interaction = lm(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1]*12^3*2.54^3/100^3)*BR5_PredTN$MedLoad[BR5_PredTN$True > 1])
lm_BR5Load_Dates = lm(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1]*12^3*2.54^3/100^3) + BR5_PredTN$MedLoad[BR5_PredTN$True > 1] + sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))
lm_BR5Load_Dates_Interaction = lm(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1]*12^3*2.54^3/100^3)*BR5_PredTN$Med[BR5_PredTN$True > 1] + sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))
lm_BR5Load_log = lm(log10(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1]) ~ BR5_PredTN$MedLoad[BR5_PredTN$True > 1] + sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))

summary(lm_BR5Load)
plot(lm_BR5Load)

summary(lm_BR5Load_Interaction)
plot(lm_BR5Load_Interaction)

summary(lm_BR5Load_Dates)
plot(lm_BR5Load_Dates)

summary(lm_BR5Load_Dates_Interaction)
plot(lm_BR5Load_Dates_Interaction)

summary(lm_BR5Load_log)
plot(lm_BR5Load_log)

#From the BR5Load, use Box-Cox transform to attempt normality
BoxCoxtest_BR5_True = boxcox(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1]*12^3*2.54^3/100^3)*BR5_PredTN$MedLoad[BR5_PredTN$True > 1], lambda = seq(0,2, length = 1000))
lambda_BR5_True = BoxCoxtest_BR5_True$x[which(BoxCoxtest_BR5_True$y == max(BoxCoxtest_BR5_True$y))]

BR5_PredTN$DiffMedLoad_BC_True[BR5_PredTN$True > 1] = bcPower(U = BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], gamma = 0, lambda = lambda_BR5_True, jacobian.adjusted = FALSE)
lm_BR5Load_Dates_BC_True = lm(BR5_PredTN$DiffMedLoad_BC_True[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1]*12^3*2.54^3/100^3)*BR5_PredTN$MedLoad[BR5_PredTN$True > 1])

summary(lm_BR5Load_Dates_BC_True)
plot(lm_BR5Load_Dates_BC_True)

#plot of regression y vs. x
plot(y = BR5_PredTN$DiffMedLoad_BC_True, x = BR5_PredTN$MedLoad)
plot(y = BR5_PredTN$DiffMedLoad_BC_True, x = I(Flows_BR5*12^3*2.54^3/100^3))
plot(y = BR5_PredTN$DiffMedLoad_BC_True, x = sin(2*pi*DatesNums_BR5/366))
plot(y = BR5_PredTN$DiffMedLoad_BC_True, x = cos(2*pi*DatesNums_BR5/366))

plot(y = BR5_PredTN$MedLoad, x = I(Flows_BR5*12^3*2.54^3/100^3))
cor(y = BR5_PredTN$MedLoad, x = I(Flows_BR5*12^3*2.54^3/100^3))

# WRTDS Interpolation Tables for Pond Branch----
#Load the streamflow data into WRTDS format
Daily_POBR = readUserDaily(filePath = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\RHESSysFilePreparation\\obs", fileName = 'BaismanStreamflow_POBR_Cal.txt', hasHeader = TRUE, separator = '\t', qUnit = 2, verbose = FALSE)
#Read the TN data
Sample_POBR = readUserSample(filePath = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\RHESSysFilePreparation\\obs", fileName = 'TN_POBR_Cal_WRTDS.txt', hasHeader = TRUE, separator = '\t', verbose = FALSE)
#Set the required information
INFO_POBR = readUserInfo(filePath = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\RHESSysFilePreparation\\obs", fileName = 'WRTDS_INFO_POBR.csv', interactive = FALSE)
eList_POBR = mergeReport(INFO = INFO_POBR, Daily = Daily_POBR, Sample = Sample_POBR)

#Make edits to include left-censored data. Nore two different thresholds.
Sample_POBR$ConcLow = ifelse((as.Date(Sample_POBR$Date) >= '2012-03-15') & (Sample_POBR$ConcLow == 0.05), NA, Sample_POBR$ConcLow)
Sample_POBR$ConcLow = ifelse((as.Date(Sample_POBR$Date) < '2012-03-15') & (Sample_POBR$ConcLow == 0.01), NA, Sample_POBR$ConcLow)
Sample_POBR$ConcHigh = ifelse((as.Date(Sample_POBR$Date) >= '2012-03-15') & (Sample_POBR$ConcHigh == 0.05), 0.05, Sample_POBR$ConcHigh)
Sample_POBR$ConcHigh = ifelse((as.Date(Sample_POBR$Date) < '2012-03-15') & (Sample_POBR$ConcHigh == 0.01), 0.01, Sample_POBR$ConcHigh)

#Fix the eList
eList_POBR$Sample = Sample_POBR
eList_POBR = fixSampleFrame(eList_POBR)

saveResults(savePath = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\RHESSysFilePreparation\\obs\\", eList = eList_POBR)

#  Default WRTDS parameters----
WRTDSmod_POBR = modelEstimation(eList = eList_POBR, windowY = 7, windowQ = 2, windowS = .5, minNumObs = 100, minNumUncen = 50, edgeAdjust = TRUE)

setwd(wd_WRTDS_PondBranch)
png('ConcFluxTime_POBR.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod_POBR)
plotFluxTimeDaily(WRTDSmod_POBR)
dev.off()

png('ConcFluxTimeErrs_POBR.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod_POBR)
plotFluxPred(WRTDSmod_POBR)
dev.off()

plotResidPred(WRTDSmod_POBR)
plotResidQ(WRTDSmod_POBR)
plotResidTime(WRTDSmod_POBR)
boxResidMonth(WRTDSmod_POBR)
boxConcThree(WRTDSmod_POBR)
boxQTwice(WRTDSmod_POBR)
plotFluxHist(WRTDSmod_POBR)
plotConcHist(WRTDSmod_POBR)

png('BiasPlot_POBR.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod_POBR, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_POBR.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,4,0.2),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 1)
dev.off()

png('ContourPlotErr_POBR.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,2,0.2),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 2)
dev.off()

#  Model#2-5 WRTDS----
WRTDSmod2_POBR = modelEstimation(eList = eList_POBR, windowY = 7, windowQ = 2, windowS = .25, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE)
WRTDSmod3_POBR = modelEstimation(eList = eList_POBR, windowY = 4, windowQ = 3, windowS = .5, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE)
WRTDSmod4_POBR = modelEstimation(eList = eList_POBR, windowY = 2, windowQ = 3, windowS = .5, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE)
#Almost no trend with flow, so making weights not matter as much.
WRTDSmod5_POBR = modelEstimation(eList = eList_POBR, windowY = 4, windowQ = 5, windowS = .5, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE)
WRTDSmod6_POBR = modelEstimation(eList = eList_POBR, windowY = 4, windowQ = 5, windowS = .25, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE)

png('ConcFluxTime_POBR_mod2.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod2_POBR)
plotFluxTimeDaily(WRTDSmod2_POBR)
dev.off()

png('ConcFluxTimeErrs_POBR_mod2.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod2_POBR)
plotFluxPred(WRTDSmod2_POBR)
dev.off()

plotResidPred(WRTDSmod2_POBR)
plotResidQ(WRTDSmod2_POBR)
plotResidTime(WRTDSmod2_POBR)
boxResidMonth(WRTDSmod2_POBR)
boxConcThree(WRTDSmod2_POBR)
boxQTwice(WRTDSmod2_POBR)
plotFluxHist(WRTDSmod2_POBR)
plotConcHist(WRTDSmod2_POBR)

png('BiasPlot_POBR_mod2.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod2_POBR, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_POBR_mod2.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod2_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,6,0.2),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 1)
dev.off()

png('ContourPlotErr_POBR_mod2.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod2_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,2.4,0.2),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 2)
dev.off()

png('ConcFluxTime_POBR_mod3.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod3_POBR)
plotFluxTimeDaily(WRTDSmod3_POBR)
dev.off()

png('ConcFluxTimeErrs_POBR_mod3.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod3_POBR)
plotFluxPred(WRTDSmod3_POBR)
dev.off()

plotResidPred(WRTDSmod3_POBR)
plotResidQ(WRTDSmod3_POBR)
plotResidTime(WRTDSmod3_POBR)
boxResidMonth(WRTDSmod3_POBR)
boxConcThree(WRTDSmod3_POBR)
boxQTwice(WRTDSmod3_POBR)
plotFluxHist(WRTDSmod3_POBR)
plotConcHist(WRTDSmod3_POBR)

png('BiasPlot_POBR_mod3.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod3_POBR, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_POBR_mod3.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod3_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,6,0.5),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 1)
dev.off()

png('ContourPlotErr_POBR_mod3.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod3_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,2.4,0.2),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 2)
dev.off()

png('ConcFluxTime_POBR_mod4.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod4_POBR)
plotFluxTimeDaily(WRTDSmod4_POBR)
dev.off()

png('ConcFluxTimeErrs_POBR_mod4.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod4_POBR)
plotFluxPred(WRTDSmod4_POBR)
dev.off()

plotResidPred(WRTDSmod4_POBR)
plotResidQ(WRTDSmod4_POBR)
plotResidTime(WRTDSmod4_POBR)
boxResidMonth(WRTDSmod4_POBR)
boxConcThree(WRTDSmod4_POBR)
boxQTwice(WRTDSmod4_POBR)
plotFluxHist(WRTDSmod4_POBR)
plotConcHist(WRTDSmod4_POBR)

png('BiasPlot_POBR_mod4.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod4_POBR, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_POBR_mod4.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod4_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,10,0.5),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 1)
dev.off()

png('ContourPlotErr_POBR_mod4.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod4_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,3,0.5),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 2)
dev.off()

png('ConcFluxTime_POBR_mod5.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod5_POBR)
plotFluxTimeDaily(WRTDSmod5_POBR)
dev.off()

png('ConcFluxTimeErrs_POBR_mod5.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod5_POBR)
plotFluxPred(WRTDSmod5_POBR)
dev.off()

plotResidPred(WRTDSmod5_POBR)
plotResidQ(WRTDSmod5_POBR)
plotResidTime(WRTDSmod5_POBR)
boxResidMonth(WRTDSmod5_POBR)
boxConcThree(WRTDSmod5_POBR)
boxQTwice(WRTDSmod5_POBR)
plotFluxHist(WRTDSmod5_POBR)
plotConcHist(WRTDSmod5_POBR)

png('BiasPlot_POBR_mod5.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod5_POBR, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_POBR_mod5.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod5_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,6,0.2),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 1)
dev.off()

png('ContourPlotErr_POBR_mod5.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod5_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,3,0.5),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 2)
dev.off()

png('ConcFluxTime_POBR_mod6.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod6_POBR)
plotFluxTimeDaily(WRTDSmod6_POBR)
dev.off()

png('ConcFluxTimeErrs_POBR_mod6.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod6_POBR)
plotFluxPred(WRTDSmod6_POBR)
dev.off()

plotResidPred(WRTDSmod6_POBR)
plotResidQ(WRTDSmod6_POBR)
plotResidTime(WRTDSmod6_POBR)
boxResidMonth(WRTDSmod6_POBR)
boxConcThree(WRTDSmod6_POBR)
boxQTwice(WRTDSmod6_POBR)
plotFluxHist(WRTDSmod6_POBR)
plotConcHist(WRTDSmod6_POBR)

png('BiasPlot_POBR_mod6.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod6_POBR, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_POBR_mod6.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod6_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,6,0.2),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 1)
dev.off()

png('ContourPlotErr_POBR_mod6.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod6_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,3,0.5),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 2)
dev.off()


#LOOCV SE for sample and discharge---
sum(WRTDSmod_POBR$Sample$SE^2)
sum(WRTDSmod2_POBR$Sample$SE^2)
sum(WRTDSmod3_POBR$Sample$SE^2)
sum(WRTDSmod4_POBR$Sample$SE^2)
sum(WRTDSmod5_POBR$Sample$SE^2)
sum(WRTDSmod6_POBR$Sample$SE^2)

sum(WRTDSmod_POBR$Daily$SE^2)
sum(WRTDSmod2_POBR$Daily$SE^2)
sum(WRTDSmod3_POBR$Daily$SE^2)
sum(WRTDSmod4_POBR$Daily$SE^2)
sum(WRTDSmod5_POBR$Daily$SE^2)
sum(WRTDSmod6_POBR$Daily$SE^2)

#  MSE for sample----
mean((WRTDSmod_POBR$Sample$ConcHat - WRTDSmod_POBR$Sample$ConcAve))^2 + sum((WRTDSmod_POBR$Sample$ConcHat - WRTDSmod_POBR$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)
mean((WRTDSmod2_POBR$Sample$ConcHat - WRTDSmod2_POBR$Sample$ConcAve))^2 + sum((WRTDSmod2_POBR$Sample$ConcHat - WRTDSmod2_POBR$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)
mean((WRTDSmod3_POBR$Sample$ConcHat - WRTDSmod3_POBR$Sample$ConcAve))^2 + sum((WRTDSmod3_POBR$Sample$ConcHat - WRTDSmod3_POBR$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)
mean((WRTDSmod4_POBR$Sample$ConcHat - WRTDSmod4_POBR$Sample$ConcAve))^2 + sum((WRTDSmod4_POBR$Sample$ConcHat - WRTDSmod4_POBR$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)
mean((WRTDSmod5_POBR$Sample$ConcHat - WRTDSmod5_POBR$Sample$ConcAve))^2 + sum((WRTDSmod5_POBR$Sample$ConcHat - WRTDSmod5_POBR$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)
mean((WRTDSmod6_POBR$Sample$ConcHat - WRTDSmod6_POBR$Sample$ConcAve))^2 + sum((WRTDSmod6_POBR$Sample$ConcHat - WRTDSmod6_POBR$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)

#  Make tables for selected model----
setwd('C:\\Users\\js4yd\\OneDrive - University of Virginia\\RHESSys_ParameterSA')
WRTDSmod5m_POBR = modelEstimation(eList = eList_POBR, windowY = 4, windowQ = 5, windowS = .5, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE, numTsteps = 50, numQsteps = 100)

#Make interpolation tables from the surfaces information and save to files that can be loaded in
#Error in location 2, then locations 4-8 are Intercept, DecYear, LogQ, SinDY, CosDY
TabInt_POBR = WRTDSmod5m_POBR$surfaces[,,4]
TabYear_POBR = WRTDSmod5m_POBR$surfaces[,,5]
TabLogQ_POBR = WRTDSmod5m_POBR$surfaces[,,6]
TabSinYear_POBR = WRTDSmod5m_POBR$surfaces[,,7]
TabCosYear_POBR = WRTDSmod5m_POBR$surfaces[,,8]
TabLogErr_POBR = WRTDSmod5m_POBR$surfaces[,,2]

rownames(TabInt_POBR) = rownames(TabYear_POBR) = rownames(TabLogQ_POBR) = rownames(TabSinYear_POBR) = rownames(TabCosYear_POBR) = rownames(TabLogErr_POBR) = attr(WRTDSmod5m_POBR$surfaces, which = 'LogQ')
colnames(TabInt_POBR) = colnames(TabYear_POBR) = colnames(TabLogQ_POBR) = colnames(TabSinYear_POBR) = colnames(TabCosYear_POBR) = colnames(TabLogErr_POBR) = attr(WRTDSmod5m_POBR$surfaces, which = 'Year')

#Contour plots of the parameters
png('TabInt_POBR.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5m_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1000,4000,100),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 4)
dev.off()
png('TabYear_POBR.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5m_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 5)
dev.off()
png('TabLogFlow_POBR.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5m_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 6)
dev.off()
png('TabSinYear_POBR.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5m_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2,2,0.2),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 7)
dev.off()
png('TabCosYear_POBR.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5m_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-3,1,0.2),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 8)
dev.off()
png('TabLogErr_POBR.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5m_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,2.5,0.1),qUnit=1, qBottom = 0.0001, qTop = 5, whatSurface = 2)
dev.off()

#Write the tables
setwd('C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\WRTDS')
options(scipen = 999)
write.table(round(TabInt_POBR,5), file = 'TabInt_POBRMod5_p5.txt', sep = '\t', row.names = attr(WRTDSmod5m_POBR$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_POBR$surfaces, which = 'Year'))
write.table(signif(TabYear_POBR,4), file = 'TabYear_POBRMod5_p4.txt', sep = '\t', row.names = attr(WRTDSmod5m_POBR$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_POBR$surfaces, which = 'Year'))
write.table(signif(TabLogQ_POBR,4), file = 'TabLogQ_POBRMod5_p4.txt', sep = '\t', row.names = attr(WRTDSmod5m_POBR$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_POBR$surfaces, which = 'Year'))
write.table(signif(TabSinYear_POBR,4), file = 'TabSinYear_POBRMod5_p4.txt', sep = '\t', row.names = attr(WRTDSmod5m_POBR$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_POBR$surfaces, which = 'Year'))
write.table(signif(TabCosYear_POBR,4), file = 'TabCosYear_POBRMod5_p4.txt', sep = '\t', row.names = attr(WRTDSmod5m_POBR$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_POBR$surfaces, which = 'Year'))
write.table(round(TabLogErr_POBR,5), file = 'TabLogErr_POBRMod5_p5.txt', sep = '\t', row.names = attr(WRTDSmod5m_POBR$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_POBR$surfaces, which = 'Year'))
options(scipen = 0)

# Load WRTDS Interpolation Tables----
setwd("C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\WRTDS")
#Pond Branch
TabInt_POBR = as.matrix(read.table(file = 'TabInt_POBRMod5_p5.txt', sep = '\t', header = TRUE, check.names = FALSE))
TabYear_POBR = as.matrix(read.table(file = 'TabYear_POBRMod5_p4.txt', sep = '\t', header = TRUE, check.names = FALSE))
TabLogQ_POBR = as.matrix(read.table(file = 'TabLogQ_POBRMod5_p4.txt', sep = '\t', header = TRUE, check.names = FALSE))
TabSinYear_POBR = as.matrix(read.table(file = 'TabSinYear_POBRMod5_p4.txt', sep = '\t', header = TRUE, check.names = FALSE))
TabCosYear_POBR = as.matrix(read.table(file = 'TabCosYear_POBRMod5_p4.txt', sep = '\t', header = TRUE, check.names = FALSE))
TabLogErr_POBR = as.matrix(read.table(file = 'TabLogErr_POBRMod5_p5.txt', sep = '\t', header = TRUE, check.names = FALSE))

#  Compare predictions of Pond Branch WRTDS and the other regression to true value----
POBR_PredTN_POBRWRTDS = matrix(NA, ncol = 3, nrow = length(Flows_POBR))
for(i in 1:length(Flows_POBR)){
  POBR_PredTN_POBRWRTDS[i,] = predictWRTDS(Date = as.character(as.Date(Dates_POBR))[i], Flow = Flows_POBR[i], rowt = as.numeric(rownames(TabInt_POBR)), colt = as.numeric(colnames(TabInt_POBR)), TabInt = TabInt_POBR, TabYear = TabYear_POBR, TabLogQ = TabLogQ_POBR, TabSinYear = TabSinYear_POBR, TabCosYear = TabCosYear_POBR, TabLogErr = TabLogErr_POBR)
}
rm(i)
colnames(POBR_PredTN_POBRWRTDS) = c('05', 'Med', '95')

POBR_PredTN_POBRWRTDS = as.data.frame(POBR_PredTN_POBRWRTDS)
POBR_PredTN_POBRWRTDS$True = BES_TN_d$POBR$TN..mg.N.L.[which((!is.na(BES_TN_d$POBR$TN..mg.N.L.)) & (as.Date(BES_TN_d$POBR$SortDate) < '2014-01-01'))]
POBR_PredTN_POBRWRTDS$DiffMed = POBR_PredTN_POBRWRTDS$True - POBR_PredTN_POBRWRTDS$Med

png('POBR_TrueTNvsPOBRWRTDS.png', res = 300, units = 'in', height = 5, width = 10)
layout(rbind(c(1,2)))
plot(POBR_PredTN$True, POBR_PredTN$Med, ylim = c(0,9), xlim = c(0,9), xlab = 'True Pond Branch TN (mg N/L)', ylab = 'WRTDS Predicted Pond Branch TN (mg N/L)', main = 'Basin WRTDS Model')
arrows(POBR_PredTN$True, POBR_PredTN$Med - POBR_PredTN$`05`, POBR_PredTN$True, POBR_PredTN$Med + POBR_PredTN$`95`, length=0.05, angle=90, code=3)
lines(c(-10,10), c(-10,10))

plot(POBR_PredTN_POBRWRTDS$True, POBR_PredTN_POBRWRTDS$Med, ylim = c(0,9), xlim = c(0,9), xlab = 'True Pond Branch TN (mg N/L)', ylab = 'WRTDS Predicted Pond Branch TN (mg N/L)', main = 'Pond Branch WRTDS Model')
arrows(POBR_PredTN_POBRWRTDS$True, POBR_PredTN_POBRWRTDS$Med - POBR_PredTN_POBRWRTDS$`05`, POBR_PredTN_POBRWRTDS$True, POBR_PredTN_POBRWRTDS$Med + POBR_PredTN_POBRWRTDS$`95`, length=0.05, angle=90, code=3)
lines(c(-10,10), c(-10,10))
dev.off()

png('POBR_TrueTNvsPOBRWRTDS_noArrows.png', res = 300, units = 'in', height = 5, width = 10)
layout(rbind(c(1,2)))
plot(POBR_PredTN$True, POBR_PredTN$Med, ylim = c(0,9), xlim = c(0,9), xlab = 'True Pond Branch TN (mg N/L)', ylab = 'WRTDS Predicted Pond Branch TN (mg N/L)', main = 'Basin WRTDS Model')
lines(c(-10,10), c(-10,10))

plot(POBR_PredTN_POBRWRTDS$True, POBR_PredTN_POBRWRTDS$Med, ylim = c(0,9), xlim = c(0,9), xlab = 'True Pond Branch TN (mg N/L)', ylab = 'WRTDS Predicted Pond Branch TN (mg N/L)', main = 'Pond Branch WRTDS Model')
lines(c(-10,10), c(-10,10))
dev.off()

#Load
POBR_PredTN_POBRWRTDS$LowTNLim = NA
POBR_PredTN_POBRWRTDS$LowTNLim[POBR_PredTN_POBRWRTDS$True == min(POBR_PredTN_POBRWRTDS$True[1:650])] = 1
POBR_PredTN_POBRWRTDS$LowTNLim[POBR_PredTN_POBRWRTDS$True == min(POBR_PredTN_POBRWRTDS$True[651:length(POBR_PredTN_POBRWRTDS$True)])] = 2

#Use WRTDS regression for predicted concentration and transform to loads. Compare predicted to true load.
POBR_PredTN_POBRWRTDS$TrueLoad = POBR_PredTN_POBRWRTDS$True*Flows_POBR*12^3*2.54^3/100^3
POBR_PredTN_POBRWRTDS$Load05 = POBR_PredTN_POBRWRTDS$`05`*Flows_POBR*12^3*2.54^3/100^3
POBR_PredTN_POBRWRTDS$MedLoad = POBR_PredTN_POBRWRTDS$Med*Flows_POBR*12^3*2.54^3/100^3
POBR_PredTN_POBRWRTDS$Load95 = POBR_PredTN_POBRWRTDS$`95`*Flows_POBR*12^3*2.54^3/100^3

png('POBR_TrueTNLoadvsPOBRWRTDSLoad.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$MedLoad, ylim=range(c(POBR_PredTN_POBRWRTDS$MedLoad-POBR_PredTN_POBRWRTDS$Load05, POBR_PredTN_POBRWRTDS$MedLoad+POBR_PredTN_POBRWRTDS$Load95)), xlim = c(0,.1), xlab = 'True Pond Branch TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of Pond Branch TN (mg N/s)')
arrows(POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$MedLoad-POBR_PredTN_POBRWRTDS$Load05, POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$MedLoad+POBR_PredTN_POBRWRTDS$Load95, length=0.05, angle=90, code=3)
lines(c(-10,10), c(-10,10))
dev.off()

png('POBR_TrueTNLoadvsPOBRWRTDSLoad_log.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$MedLoad, ylim=range(c(POBR_PredTN_POBRWRTDS$MedLoad-POBR_PredTN_POBRWRTDS$Load05, POBR_PredTN_POBRWRTDS$MedLoad+POBR_PredTN_POBRWRTDS$Load95)), xlab = 'True Pond Branch TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of Pond Branch TN (mg N/s)', log = 'xy')
arrows(POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$MedLoad-POBR_PredTN_POBRWRTDS$Load05, POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$MedLoad+POBR_PredTN_POBRWRTDS$Load95, length=0.05, angle=90, code=3)
lines(c(1e-7,1e7), c(1e-7,1e7))
dev.off()

png('POBR_TrueTNLoadvsPOBRWRTDSLoad_log_colors.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$MedLoad, ylim=range(c(POBR_PredTN_POBRWRTDS$MedLoad-POBR_PredTN_POBRWRTDS$Load05, POBR_PredTN_POBRWRTDS$MedLoad+POBR_PredTN_POBRWRTDS$Load95)),
     xlim=range(POBR_PredTN_POBRWRTDS$TrueLoad),
     xlab = 'True Pond Branch TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of Pond Branch TN (mg N/s)', log = 'xy')
par(new = TRUE)
plot(POBR_PredTN_POBRWRTDS$TrueLoad[POBR_PredTN_POBRWRTDS$LowTNLim == 1], POBR_PredTN_POBRWRTDS$MedLoad[POBR_PredTN_POBRWRTDS$LowTNLim == 1], 
     xlim=range(POBR_PredTN_POBRWRTDS$TrueLoad),
     ylim=range(c(POBR_PredTN_POBRWRTDS$MedLoad-POBR_PredTN_POBRWRTDS$Load05, POBR_PredTN_POBRWRTDS$MedLoad+POBR_PredTN_POBRWRTDS$Load95)), xlab = '', ylab = '', 
     log = 'xy', col = 'blue', axes = FALSE)
par(new = TRUE)
plot(POBR_PredTN_POBRWRTDS$TrueLoad[POBR_PredTN_POBRWRTDS$LowTNLim == 2], POBR_PredTN_POBRWRTDS$MedLoad[POBR_PredTN_POBRWRTDS$LowTNLim == 2], 
     xlim=range(POBR_PredTN_POBRWRTDS$TrueLoad),
     ylim=range(c(POBR_PredTN_POBRWRTDS$MedLoad-POBR_PredTN_POBRWRTDS$Load05, POBR_PredTN_POBRWRTDS$MedLoad+POBR_PredTN_POBRWRTDS$Load95)), xlab = '', ylab = '', 
     log = 'xy', col = 'red', axes = FALSE)
#arrows(POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$MedLoad-POBR_PredTN_POBRWRTDS$Load05, POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$MedLoad+POBR_PredTN_POBRWRTDS$Load95, length=0.05, angle=90, code=3)
lines(c(1e-7,1e7), c(1e-7,1e7))
legend('topleft', title = 'Detection Limits (mg N/L)', legend = c('0.01', '0.05'), pch = 1, col = c('blue', 'red'))
dev.off()

#Add colors for seasons to this dataset
POBR_PredTN_POBRWRTDS$cols = BES_TN_d$POBR$cols[as.Date(BES_TN_d$POBR$SortDate) %in% as.Date(Dates_POBR)]

png('POBR_TN_ScatterHist_POBRWRTDS.png', res = 300, units = 'in', width = 10, height = 10)
scatterHist_mod(x = log10(POBR_PredTN_POBRWRTDS$TrueLoad), y = log10(POBR_PredTN_POBRWRTDS$MedLoad), density = FALSE, freq = FALSE, x.breaks = 100, y.breaks = 100, correl = FALSE, xlab = expression(paste('log'[10] ~ '(True TN Load [mg N/s])')), ylab = expression(paste('log'[10] ~ 'Predicted TN Load (mg N/s)')), cex.lab = 1.5, cex.axis = 1.5, title = '', pch = 16, ellipse = FALSE, smooth = FALSE, col = POBR_PredTN_POBRWRTDS$cols)
dev.off()

POBR_PredTN_POBRWRTDS$DiffMedLoad = POBR_PredTN_POBRWRTDS$TrueLoad - POBR_PredTN_POBRWRTDS$MedLoad

#TN Load as modeled from regression
POBR_PredTN$ModeledTNLoad = POBR_PredTN$MedLoad + lm_POBRLoad_Dates_Interaction$fitted.values

png('POBR_TrueTNLoadvsWRTDSLoad+CorrModel_log_colors.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN$TrueLoad, POBR_PredTN$ModeledTNLoad, ylim = c(1e-5, 1e-1), xlim = c(1e-5, 1e-1),
     xlab = 'True Pond Branch TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of Pond Branch TN (mg N/s)', log = 'xy')
par(new = TRUE)
plot(POBR_PredTN$TrueLoad[POBR_PredTN$LowTNLim == 1], POBR_PredTN$ModeledTNLoad[POBR_PredTN$LowTNLim == 1], 
     ylim = c(1e-5, 1e-1), xlim = c(1e-5, 1e-1), xlab = '', ylab = '', 
     log = 'xy', col = 'blue', axes = FALSE)
par(new = TRUE)
plot(POBR_PredTN$TrueLoad[POBR_PredTN$LowTNLim == 2], POBR_PredTN$ModeledTNLoad[POBR_PredTN$LowTNLim == 2], 
     ylim = c(1e-5, 1e-1), xlim = c(1e-5, 1e-1), xlab = '', ylab = '', 
     log = 'xy', col = 'red', axes = FALSE)
#arrows(POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$MedLoad-POBR_PredTN_POBRWRTDS$Load05, POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$MedLoad+POBR_PredTN_POBRWRTDS$Load95, length=0.05, angle=90, code=3)
lines(c(1e-7,1e7), c(1e-7,1e7))
legend('topleft', title = 'Detection Limits (mg N/L)', legend = c('0.01', '0.05'), pch = 1, col = c('blue', 'red'))
dev.off()

#Make a map of these water quality sampling locations----
setwd(dir_WChem)

res = 30

world = read.csv("C:\\Users\\js4yd\\Documents\\BaismanSA\\RHESSysRuns\\Run0\\worldfiles\\worldfile.csv", stringsAsFactors = FALSE)

#Taking the unique patch IDs because strata can exist in more than one patch.
Area.basin = length(unique(world$patchID))*res^2

#Get hillslope areas and conversion factor for streamflow in hillslopes
uhills = unique(world$hillID)
Area.Hills = matrix(NA, nrow = length(uhills), ncol = 2)
for (h in 1:length(uhills)){
  Area.Hills[h,1] = h
  #some patches have multiple strata, so their area cannot be counted from the count of cells.
  Area.Hills[h,2] = length(which(world[which(duplicated(world$patchID) == FALSE),]$hillID == h))*res^2
}
rm(h)

coordinates(world) = c('patchX', 'patchY')
proj4string(world) = CRS('+init=epsg:26918')
world=spTransform(world, CRSobj = pCRS)

png('BaismanSynopticWaterQualitySites.png', res = 300, units = 'in', width = 6, height = 6)
par(mar= c(2.5,2.5,1,1))
plot(world, col = 'white')
for (h in 1:length(uhills)){
  plot(world[world$hillID == uhills[h],], pch = 22, add = TRUE, lwd=10, col = 'gray')
  plot(world[world$hillID == uhills[h],], col = 'white', pch = 15, add = TRUE)
}
rm(h)
plot(MDstreams, col = 'skyblue', add = T, lwd = 4)
plot(S_Sites, add = T, col = 'green', cex = 2)
plot(K_Sites, add = T, col = 'red', cex = 2)
text(x = K_Sites@coords[,1], y = K_Sites@coords[,2], K_Sites$NAME, col = 'darkred', cex = 0.4)
text(x = S_Sites@coords[,1], y = S_Sites@coords[,2], S_Sites$NAME, col ='darkgreen', cex = 0.4)
legend('bottomright', legend = c('Kenworth: 2001-2', 'Smith: 2006-7'), col = c('red', 'green'), pch = 3)
box(which = 'figure', lwd = 2)
dev.off()
# degAxis(side = 1, at = seq(-77,-76,.01), labels = FALSE)
# degAxis(side = 1, at = seq(-76.7,-76,.02))
# degAxis(side = 3, at = seq(-77,-76,.01), labels = FALSE)
# degAxis(side = 2, at = seq(39.45, 40,.01))
# degAxis(side = 4, at = seq(39.45, 40,.01), labels = FALSE)
# north.arrow(xb = -76.712, yb = 39.469, len = .0005, lab = 'N', tcol = 'black', col='black')
# text(x = -76.712, y = 39.467, 'WGS84')

#Make a map of the land ID in the world file----
#Check that the number of unique patches equals the number of unique land uses
if (length(unique(world$patchID)) != length(unique(world$patchID[world$patchLandID == 2])) + length(unique(world$patchID[world$patchLandID == 1])) + length(unique(world$patchID[world$patchLandID == 3])) + length(unique(world$patchID[world$patchLandID == 4]))){
  print('Number of unique patches by land use ID is not equal to the number of unique patches.')
}
  
png('Baisman_LandIDs.png', res = 300, units = 'in', width = 6, height = 6)
par(mar= c(2.5,2.5,1,1))
plot(world, col = 'white')
for (h in 1:length(uhills)){
  plot(world[world$hillID == uhills[h],], pch = 22, add = TRUE, lwd=8, col = 'gray')
  #plot the 4 land use types
  plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 1)),], col = 'yellow', pch = 15, add = TRUE, cex = 0.65)
  plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 2)),], col = 'yellow', pch = 15, add = TRUE, cex = 0.65)
  plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 3)),], col = 'purple', pch = 15, add = TRUE, cex = 0.65)
  plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 4)),], col = 'purple', pch = 15, add = TRUE, cex = 0.65)
}
rm(h)
plot(MDstreams, col = 'blue', add = T, lwd = 4)
plot(S_Sites, add = T, col = 'green', cex = 2)
plot(K_Sites, add = T, col = 'red', cex = 2)
legend('topright', legend = c('Kenworth: 2001-2', 'Smith: 2006-7'), col = c('red', 'green'), pch = 3)
legend('bottomright', legend = c('Undeveloped', 'Developed'), col = c('yellow', 'purple'), pch = 15)
box(which = 'figure', lwd = 2)
dev.off()

png('Baisman_LandIDs_grid.png', res = 300, units = 'in', width = 6, height = 6)
par(mar= c(2.5,2.5,1,1))
plot(world, col = 'white')
for (h in 1:length(uhills)){
  plot(world[world$hillID == uhills[h],], pch = 22, add = TRUE, lwd=8, col = 'gray')
  #plot the 4 land use types
  plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 1)),], col = 'yellow', pch = 15, add = TRUE, cex = 0.5)
  plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 2)),], col = 'yellow', pch = 15, add = TRUE, cex = 0.5)
  plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 3)),], col = 'purple', pch = 15, add = TRUE, cex = 0.5)
  plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 4)),], col = 'purple', pch = 15, add = TRUE, cex = 0.5)
}
rm(h)
plot(MDstreams, col = 'blue', add = T, lwd = 4)
plot(S_Sites, add = T, col = 'green', cex = 2)
plot(K_Sites, add = T, col = 'red', cex = 2)
legend('topright', legend = c('Kenworth: 2001-2', 'Smith: 2006-7'), col = c('red', 'green'), pch = 3)
legend('bottomright', legend = c('Undeveloped', 'Developed'), col = c('yellow', 'purple'), pch = 15)
box(which = 'figure', lwd = 2)
dev.off()

png('Baisman_LandIDs_grid_UrbanOnly.png', res = 300, units = 'in', width = 6, height = 6)
par(mar= c(2.5,2.5,1,1))
plot(world, col = 'white')
for (h in 1:length(uhills)){
  plot(world[world$hillID == uhills[h],], pch = 22, add = TRUE, lwd=8, col = 'gray')
  #plot the 4 land use types
  #plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 1)),], col = 'yellow', pch = 15, add = TRUE, cex = 0.5)
  #plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 2)),], col = 'yellow', pch = 15, add = TRUE, cex = 0.5)
  plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 3)),], col = 'purple', pch = 15, add = TRUE, cex = 0.5)
  #plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 4)),], col = 'purple', pch = 15, add = TRUE, cex = 0.5)
}
rm(h)
plot(MDstreams, col = 'blue', add = T, lwd = 4)
plot(S_Sites, add = T, col = 'green', cex = 2)
plot(K_Sites, add = T, col = 'red', cex = 2)
legend('topright', legend = c('Kenworth: 2001-2', 'Smith: 2006-7'), col = c('red', 'green'), pch = 3)
legend('bottomright', legend = c('Undeveloped', 'Developed'), col = c('yellow', 'purple'), pch = 15)
box(which = 'figure', lwd = 2)
dev.off()

png('Baisman_LandIDs_grid_UrbanWithSeptic.png', res = 300, units = 'in', width = 6, height = 6)
par(mar= c(2.5,2.5,1,1))
plot(world, col = 'white')
for (h in 1:length(uhills)){
  plot(world[world$hillID == uhills[h],], pch = 22, add = TRUE, lwd=8, col = 'gray')
  #plot the 4 land use types
  #plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 1)),], col = 'yellow', pch = 15, add = TRUE, cex = 0.5)
  #plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 2)),], col = 'yellow', pch = 15, add = TRUE, cex = 0.5)
  plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 3)),], col = 'purple', pch = 15, add = TRUE, cex = 0.5)
  plot(world[which((world$hillID == uhills[h]) & (world$patchLandID == 4)),], col = 'purple', pch = 15, add = TRUE, cex = 0.5)
}
rm(h)
plot(MDstreams, col = 'blue', add = T, lwd = 4)
plot(S_Sites, add = T, col = 'green', cex = 2)
plot(K_Sites, add = T, col = 'red', cex = 2)
legend('topright', legend = c('Kenworth: 2001-2', 'Smith: 2006-7'), col = c('red', 'green'), pch = 3)
legend('bottomright', legend = c('Undeveloped', 'Developed'), col = c('yellow', 'purple'), pch = 15)
box(which = 'figure', lwd = 2)
dev.off()

# Make a forest model (Pond Branch WRTDS) and a developed model (based on Baisman outlet WRTDS) for TN----
#Obtain the fraction of developed land in each of the hillslopes, and for the basin----
FracDev = vector('numeric', length = length(unique(world$hillID)))
for (i in 1:length(FracDev)){
  FracDev[i] = length(which((world$hillID == uhills[i]) & ((world$patchLandID == 3) | (world$patchLandID == 4))))/length(which((world$hillID == uhills[i])))
}
#Because this is a 2-component mixture model, assume that the remainder is undeveloped land
FracUnDev = 1 - FracDev

#Use export coefficient/source-contribution model----
#This is a signal + background model at BARN because BARN = POBR + all other catchments
Flows_BARN = BES_TN_d$BARN$Flow[which((as.Date(BES_TN_d$BARN$SortDate) <= '2014-01-01'))]
Dates_BARN = BES_TN_d$BARN$SortDate[which((as.Date(BES_TN_d$BARN$SortDate) <= '2014-01-01'))]
BARN_PredTN = matrix(NA, ncol = 3, nrow = length(Flows_BARN))
for(i in 1:length(Flows_BARN)){
  BARN_PredTN[i,] = predictWRTDS(Date = as.character(as.Date(Dates_BARN))[i], Flow = Flows_BARN[i], rowt = as.numeric(rownames(TabInt)), colt = as.numeric(colnames(TabInt)), TabInt = TabInt, TabYear = TabYear, TabLogQ = TabLogQ, TabSinYear = TabSinYear, TabCosYear = TabCosYear, TabLogErr = TabLogErr)
}
rm(i)
colnames(BARN_PredTN) = c('05', 'Med', '95')

BARN_PredTN = as.data.frame(BARN_PredTN)
BARN_PredTN$True = BES_TN_d$BARN$TN..mg.N.L.[which((as.Date(BES_TN_d$BARN$SortDate) <= '2014-01-01'))]

#Load
BARN_PredTN$LowTNLim = NA

#concentration to loads
BARN_PredTN$TrueLoad = BARN_PredTN$True*Flows_BARN*12^3*2.54^3/100^3
BARN_PredTN$Load05 = BARN_PredTN$`05`*Flows_BARN*12^3*2.54^3/100^3
BARN_PredTN$MedLoad = BARN_PredTN$Med*Flows_BARN*12^3*2.54^3/100^3
BARN_PredTN$Load95 = BARN_PredTN$`95`*Flows_BARN*12^3*2.54^3/100^3

#All dates for Pond Branch
Flows_POBR_AllDates = BES_TN_d$POBR$Flow[which((as.Date(BES_TN_d$POBR$SortDate) <= '2014-01-01') & (as.Date(BES_TN_d$POBR$SortDate) >= '1999-11-15'))]
Dates_POBR_AllDates = BES_TN_d$POBR$SortDate[which((as.Date(BES_TN_d$POBR$SortDate) <= '2014-01-01') & (as.Date(BES_TN_d$POBR$SortDate) >= '1999-11-15'))]
POBR_AllDates_PredTN = matrix(NA, ncol = 3, nrow = length(Flows_POBR_AllDates))
for(i in 1:length(Flows_POBR_AllDates)){
  POBR_AllDates_PredTN[i,] = predictWRTDS(Date = as.character(as.Date(Dates_POBR_AllDates))[i], Flow = Flows_POBR_AllDates[i], rowt = as.numeric(rownames(TabInt_POBR)), colt = as.numeric(colnames(TabInt_POBR)), TabInt = TabInt_POBR, TabYear = TabYear_POBR, TabLogQ = TabLogQ_POBR, TabSinYear = TabSinYear_POBR, TabCosYear = TabCosYear_POBR, TabLogErr = TabLogErr_POBR)
}
rm(i)
colnames(POBR_AllDates_PredTN) = c('05', 'Med', '95')

POBR_AllDates_PredTN = as.data.frame(POBR_AllDates_PredTN)
POBR_AllDates_PredTN$True = BES_TN_d$POBR$TN..mg.N.L.[which((as.Date(BES_TN_d$POBR$SortDate) <= '2014-01-01') & (as.Date(BES_TN_d$POBR$SortDate) >= '1999-11-15'))]

#Load
POBR_AllDates_PredTN$LowTNLim = NA
POBR_AllDates_PredTN$LowTNLim[POBR_AllDates_PredTN$True == min(POBR_AllDates_PredTN$True[1:4865])] = 1
POBR_AllDates_PredTN$LowTNLim[POBR_AllDates_PredTN$True == min(POBR_AllDates_PredTN$True[4866:length(POBR_AllDates_PredTN$True)])] = 2

#concentration to loads
POBR_AllDates_PredTN$TrueLoad = POBR_AllDates_PredTN$True*Flows_POBR_AllDates*12^3*2.54^3/100^3
POBR_AllDates_PredTN$Load05 = POBR_AllDates_PredTN$`05`*Flows_POBR_AllDates*12^3*2.54^3/100^3
POBR_AllDates_PredTN$MedLoad = POBR_AllDates_PredTN$Med*Flows_POBR_AllDates*12^3*2.54^3/100^3
POBR_AllDates_PredTN$Load95 = POBR_AllDates_PredTN$`95`*Flows_POBR_AllDates*12^3*2.54^3/100^3

#ECM based on land use type----
#Export coefficient for Pond Branch = undeveloped export coefficient
EC_Undev = POBR_AllDates_PredTN$MedLoad/sum(Area.Hills[c(3,4),2])
#Export coefficient for developed = (Baisman load - undeveloped load)/(developed land area)
EC_Dev = (BARN_PredTN$MedLoad - EC_Undev*sum(Area.Hills[,2]*FracUnDev))/sum(Area.Hills[,2]*FracDev)

#There are some negative values, so maybe this can be attributed to the in-stream losses. Plot dates that these happen
plot(y = EC_Dev, x = Dates_POBR_AllDates, ylim = c(0, 1.5e-6), xlim = c(10500, 16500))
par(new = TRUE)
plot(y = EC_Dev[EC_Dev < 0], x = Dates_POBR_AllDates[EC_Dev < 0], ylim = c(0, 1.5e-6), xlim = c(10500, 16500), col = 'red')

#They are all during a drought! So that's almost surely a result of not accounting for in-stream losses.
Flows_POBR_AllDates[EC_Dev < 0]

#Fixme: setting negative values to 0 for now.
#EC_Undev[EC_Dev < 0] = BARN_PredTN$MedLoad - EC_Dev[EC_Dev < 0]*sum(Area.Hills[,2]*FracDev)
#EC_Dev[EC_Dev < 0] = 0

#ECM based on impervious fraction of land----
impFrac = raster(x = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\ImperviousFraction.tiff")
impFrac = projectRaster(impFrac, crs = CRS(pCRS))

#Add impervious fraction to worldfile information
world$ImpFrac = raster::extract(x = impFrac, y = world)

#Compute the impervious area of each hillslope and for the basin
Area.basin_imp = 0
for (i in 1:length(unique(world$patchID))){
  Area.basin_imp = Area.basin_imp + world$ImpFrac[which(duplicated(world$patchID) == FALSE)[i]]*res^2
}
Area.basin_imp/Area.basin #5.1% impervious
Area.Hills = cbind(Area.Hills, rep(0, nrow(Area.Hills)))
for (h in 1:length(uhills)){
  #Get the impervious surface area in this hillslope
  for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$hillID == h))){
    Area.Hills[h,3] = Area.Hills[h,3] + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$hillID == h)[i]]*res^2
  }
}

#Add fraction of impervious surface to the dataset
Area.Hills = cbind(Area.Hills, Area.Hills[,3]/Area.Hills[,2])

#Estimate the Pond Branch "Undeveloped" model
#Fixme: neglecting the component of Pond Branch that is developed here.
EC_Undev = POBR_AllDates_PredTN$MedLoad/sum(Area.Hills[c(3,4),2])

#Export coefficient for developed = (Baisman load - undeveloped load)/(developed land area)
#Fixme: developed Pond Branch is being subtracted to undeveloped. Added to developed.
EC_Dev = (BARN_PredTN$MedLoad - EC_Undev*sum(Area.Hills[,2] - Area.Hills[,3], 116))/sum(Area.Hills[,3], -116)

#There are some negative values, so maybe this can be attributed to the in-stream losses. Plot dates that these happen
plot(y = EC_Dev, x = Dates_POBR_AllDates, ylim = c(0, 1e-5), xlim = c(10500, 16500))
par(new = TRUE)
plot(y = EC_Dev[EC_Dev < 0], x = Dates_POBR_AllDates[EC_Dev < 0], ylim = c(0, 1e-5), xlim = c(10500, 16500), col = 'red')

#Also all during the drought of 2002 and 2007.

#Fixme: for now setting all of these dates equal to 0 and not adjusting POBR
EC_Dev[EC_Dev < 0] = 0

#Look at fit of model to Pond Branch and Baisman
Load_ECM_Pond = EC_Undev*sum(Area.Hills[c(3,4),2]-Area.Hills[c(3,4),3], 116) + EC_Dev*sum(Area.Hills[c(3,4),3], -116)

plot(as.Date(Dates_POBR_AllDates), POBR_AllDates_PredTN$MedLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,0.03))
par(new = TRUE)
plot(as.Date(Dates_POBR_AllDates), Load_ECM_Pond, col = 'red', xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,0.03))

#Check on method
if(cor(Load_ECM_Pond, POBR_AllDates_PredTN$MedLoad) != 1){
  print('method producing correlation != 1')
}

plot(as.Date(Dates_POBR_AllDates), POBR_AllDates_PredTN$TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,0.03))
par(new = TRUE)
plot(as.Date(Dates_POBR_AllDates), Load_ECM_Pond, col = 'red', xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,0.03))

#in kg/day
plot(as.Date(Dates_POBR_AllDates), POBR_AllDates_PredTN$TrueLoad/1000*3600*24, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,5))
par(new = TRUE)
plot(as.Date(Dates_POBR_AllDates), Load_ECM_Pond/1000*3600*24, col = 'red', xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,5))

plot(POBR_AllDates_PredTN$TrueLoad, Load_ECM_Pond,
     ylim = c(0.000001,0.05), xlim = c(0.000001,0.05), 
     xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Pond Branch Outlet', log = 'xy')
lines(c(0.000001,1.3), c(0.000001,1.3))

Load_ECM_Baisman = EC_Undev*(Area.basin - Area.basin_imp) + EC_Dev*(Area.basin_imp)

plot(as.Date(Dates_BARN), BARN_PredTN$TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,1.3))
par(new = TRUE)
plot(as.Date(Dates_BARN), Load_ECM_Baisman, col = 'red', xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,1.3))

#Load in kg/d
plot(as.Date(Dates_BARN), BARN_PredTN$TrueLoad/1000*3600*24, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,90))
par(new = TRUE)
plot(as.Date(Dates_BARN), Load_ECM_Baisman/1000*3600*24, col = 'red', xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,90))

plot(BARN_PredTN$TrueLoad/1000*3600*24, Load_ECM_Baisman/1000*3600*24,
     ylim = c(0,90), xlim = c(0,90), xlab = 'True TN Load (kg N/d)', ylab = 'ECM Predicted TN Load (kg N/d)', main = 'Baisman Outlet')
lines(c(0,90), c(0,90), col = 'red')

#cor(BARN_PredTN$TrueLoad, Load_ECM_Baisman)

plot(as.Date(Dates_BARN), BARN_PredTN$MedLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,1.3))
par(new = TRUE)
plot(as.Date(Dates_BARN), Load_ECM_Baisman, col = 'red', xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,1.3))

plot(BARN_PredTN$MedLoad, Load_ECM_Baisman,
     ylim = c(0,1.3), xlim = c(0,1.3), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Baisman Outlet')
lines(c(0,1.3), c(0,1.3))

#Check on method
if(cor(Load_ECM_Baisman, BARN_PredTN$MedLoad) != 1){
  print('method producing correlation != 1')
}
#This is okay - some sites were set to 0

Load_ECM_BR3 = EC_Undev*sum(Area.Hills[c(13,14),2]-Area.Hills[c(13,14),3]) + EC_Dev*sum(Area.Hills[c(13,14),3])

plot(as.Date(Dates_BARN), Load_ECM_BR3, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,0.4))
par(new = TRUE)
plot(as.Date(Dates_BR3), BR3_PredTN$TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,0.4), col = 'red')

plot(BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)]*1000, Load_ECM_BR3[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)], 
     xlim = c(0,0.05), ylim = c(0,0.05), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Hillslope 3')
lines(c(0,0.05), c(0,0.05))

Load_ECM_BR3_AllUnDev = EC_Undev*sum(Area.Hills[c(13,14),2])

plot(as.Date(Dates_BARN), Load_ECM_BR3_AllUnDev, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,0.06))
par(new = TRUE)
plot(as.Date(Dates_BR3), BR3_PredTN$TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,0.06), col = 'red')

plot(BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)], Load_ECM_BR3_AllUnDev[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)], 
     xlim = c(0,0.006), ylim = c(0,0.006), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Hillslope 3')
lines(c(0,0.006), c(0,0.006))

# Compare models to where there are datasets----
#Obtain upstream contributing patches for each sampling location
#Add those patch identifiers to the world dataframe, one column for each site

#Evaluate adding flow information and other normalizers to the ECM models----
# Convert back to concentration by subtracting flow that resulted from covered catchments. Then make a model for urban and a model for forest at basin outlet and make sure that forest is same as POBR and that devekoped matches the sampling site data well.
