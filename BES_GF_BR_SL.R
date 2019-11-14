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
#Precipitation
dir_precip = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\Precipitation"
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

# BES Water Quality Gauge Data----
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

#Fixme: Do some diagnostics on the sites to come up with a single timeseries for each site----


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

#Fixme: process nexrad radar info----
setwd("C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\Precipitation\\BES_Nexrad\\Balto2000")
Nexrad = read.table('200004040830.txt', header = FALSE)
colnames(Nexrad) = c('lat', 'long', 'precip')
test = rasterFromXYZ(Nexrad[,c(2,1,3)], crs = '+init=epsg:4326')
coordinates(Nexrad) = c('long', 'lat')
proj4string(Nexrad) = CRS('+init=epsg:4326')
gridded(Nexrad) = TRUE
NexradRas = raster(Nexrad)
NexradRas = projectRaster(NexradRas, crs = CRS(pCRS))

plot(NexradRas, col = terrain.colors(30))
plot(ROI, add = T, border = 'black')

#Fixme: Download climate indices (ENSO, PDO, etc.) and evaluate timeseries relative to those----

#Fixme: add weather station download info to the main repository----

#Create Observation timeseries for RHESSys calibration----
#Baisman USGS Streamflow and TN----
#Streamflow and Nitrogen are both available in the BES_TN_d list
BES_TN_d_load = list.load('C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\WaterChemistry\\Nitrogen\\BES_TN_d.yaml', type = "YAML")

#Earliest possible start date for calibration
min(as.Date(BES_TN_d_load$BARN$SortDate, origin="1970-01-01"))

#Latest possible end date for validation
max(as.Date(BES_TN_d_load$BARN$SortDate, origin="1970-01-01"))

#Total number of samples
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