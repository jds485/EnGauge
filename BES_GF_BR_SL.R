#Loading data from the output of Method1Example to use with additional datasets acquired for a specific site

#Fixme: add satellite background images to maps

#Set directory names----
#EnGauge repository
dir_EnGauge = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\EnGauge\\EnGauge"
#Color functions - from JDS github repo: Geothermal_ESDA
dir_ColFuns = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\USGSGauges"

#Repositories that should match those used for or created by EnGauge 
#Region of interest shapefile
dir_ROI = paste0(dir_EnGauge, "\\DataForExamples")
#dir_ROI = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\BES-Watersheds-Land-Cover-Analysis"  
#Watersheds shapefiles
dir_Sheds = paste0(dir_EnGauge, "\\DataForExamples\\BES_GF_BR_SL_Data")
#dir_Sheds = "C:\\Users\\js4yd\\Documents\\DEMtest"
#DEM
dir_DEM_out = "C:\\Users\\js4yd\\Desktop\\TestEnGauge\\BaismanAndGwynnsFalls"
#dir_DEM_out = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\EnGauge\\EnGauge\\DataForExamples"
#dir_DEM_out = 'C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\DEM'
#Streamflow gauges
wd_sf = "C:\\Users\\js4yd\\Desktop\\TestEnGauge\\BaismanAndGwynnsFalls\\Streamflow"
#wd_sf = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\EnGauge\\EnGauge\\DataForExamples\\Streamflow"
#wd_sf = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\USGSGauges\\Streamflow"
#Nitrogen
wd_N = "C:\\Users\\js4yd\\Desktop\\TestEnGauge\\BaismanAndGwynnsFalls\\Nitrogen"
#wd_N = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\EnGauge\\EnGauge\\DataForExamples\\Nitrogen"
#wd_N = 'C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\USGSGauges\\Nitrogen'
#Phosphorus
wd_P = "C:\\Users\\js4yd\\Desktop\\TestEnGauge\\BaismanAndGwynnsFalls\\Phosphorus"
#wd_P = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\EnGauge\\EnGauge\\DataForExamples\\Phosphorus"
#wd_P = 'C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\USGSGauges\\Phosphorus'

#Water Chemistry from non-EnGauge Datasets
dir_WChem = "C:\\Users\\js4yd\\Desktop\\TestEnGauge\\BaismanAndGwynnsFalls"
#dir_WChem = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\WaterChemistry"
dir_WChemData = paste0(dir_EnGauge, "\\DataForExamples\\BES_GF_BR_SL_Data\\BES_WaterChemistry")
#dir_WChemData = dirWChem

#Synoptic Water Chemistry
dir_SynWChem_Kenworth = paste0(dir_WChemData, "\\BARN_Synoptic\\WaterChemical_Kenworth_01-02")
dir_SynWChem_Smith = paste0(dir_WChemData, "\\BARN_Synoptic\\WaterChemical_Smith_06-07")
dir_ContributingAreas = paste0(dir_WChemData, "\\BARN_Synoptic\\ContributingAreaSamplingSites")

#Baisman Run and Pond Branch Main Directory
wd_BRPOBR = "C:\\Users\\js4yd\\Desktop\\TestEnGauge\\BaismanAndGwynnsFalls\\BR&POBR"
#wd_BRPOBR = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR"

#WRTDS Interpolation Table Directory
wd_WRTDS_BARN = paste0(wd_BRPOBR, '\\WRTDS')
wd_WRTDS_BARN_INFO = paste0(dir_EnGauge, "\\DataForExamples\\BES_GF_BR_SL_Data\\BES_WaterChemistry")
#wd_WRTDS_BARN_INFO = paste0(wd_BRPOBR, '\\WRTDS')
wd_WRTDS_POBR = paste0(wd_WRTDS_BARN, '\\POBR')
wd_WRTDS_POBR_INFO = paste0(dir_EnGauge, "\\DataForExamples\\BES_GF_BR_SL_Data\\BES_WaterChemistry")
#wd_WRTDS_POBR_INFO = paste0(wd_WRTDS_BARN, '\\POBR')

#RHESSys Observation File Directory
wd_RHESSysObs_BARN = paste0(wd_BRPOBR, "\\RHESSysFilePreparation\\obs")
wd_RHESSysObs_POBR = paste0(wd_RHESSysObs_BARN, "\\POBR")

#Precipitation
dir_precip = "C:\\Users\\js4yd\\Desktop\\TestEnGauge\\BaismanAndGwynnsFalls\\Precipitation"
#dir_precip = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\Precipitation"
dir_precipData = paste0(dir_EnGauge, "\\DataForExamples\\BES_GF_BR_SL_Data\\BES_Precipitation")
#dir_precipData = dir_precip

#Streams Shapefile
dir_streams = paste0(dir_EnGauge, "\\DataForExamples")
#dir_streams = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\NHD_H_Maryland_State_Shape\\Shape"

#RHESSys Climate File Directory
wd_clim = paste0(wd_BRPOBR, "\\RHESSysFilePreparation\\clim")
#Laurence Lin's climate files
#wd_LinData = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\JDS_DesktopFiles_13Feb2020\\rhessys30m_Pond\\clim"

#Worldfile directory
dir_worldfile = paste0(dir_EnGauge, "\\DataForExamples\\BES_GF_BR_SL_Data\\BES_Worldfile")
#dir_worldfile = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\RHESSysFiles\\BR&POBR\\SARunReferenceData\\Run0\\worldfiles"

#Check for and make directories that are required to run the script----
dir.create(path = dir_precip, showWarnings = FALSE, recursive = FALSE)
dir.create(path = wd_BRPOBR, showWarnings = FALSE, recursive = FALSE)
dir.create(path = wd_WRTDS_BARN, showWarnings = FALSE, recursive = FALSE)
dir.create(path = wd_WRTDS_POBR, showWarnings = FALSE, recursive = FALSE)
dir.create(path = wd_clim, showWarnings = FALSE, recursive = TRUE)
dir.create(path = wd_RHESSysObs_BARN, showWarnings = FALSE, recursive = FALSE)
dir.create(path = wd_RHESSysObs_POBR, showWarnings = FALSE, recursive = FALSE)

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
#Kenworth Synoptic Water Quality
#Locations
f_KenworthSites = "Kenworth_Chem_Sites"
#Flow
f_KenworthWQ_Q = "Kenworth_Q.csv"
#TN
f_KenworthWQ_TN = "Kenworth_TN.csv"
#Smith Synoptic Water Quality
#Locations
f_SmithSites = "Smith_Chem_Sites"
#Flow
f_SmithWQ_Q = "Smith_discharge.csv"
#TN
f_SmithWQ_TN = "Smith_TN.csv"
#USGS stream gauge matches
f_USGS_GaugeMatch = 'Abbreviations_SampleRecordLengths.csv'
#Precipitation Measurements
f_Precip = 'BES_precipitation.csv'
f_Precip_locs = 'BES_GaugeLocs.csv'
#Stream shapefile
f_streams = 'MDstreams'
#f_streams = 'NHDFlowline'
#Baisman Impervious Fraction
f_BaismanImperviousFrac = "ImperviousFraction.tiff"
#Baisman Slopes
f_BaismanSlope = "slopedem30m.tif"
#Laurence Lin Datasets
#f_LinRain = "Oregon.rain"
#f_LinTmax = "Oregon.tmax"
#f_LinTmin = "Oregon.tmin"
#Baisman Worldfile
f_worldfile = "worldfile.csv"

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
f_NexradDataList_Gauges = 'NexradYearList_Gauges.yaml'
f_BaismanNexradPixels_d = 'BaismanNexradPixels.csv'
f_BaismanNexradPixels_15min = 'BaismanNexradPixels_15minPrecip.csv'
f_BaismanNexradPixelsMap = 'BaismanNexradPixels'
f_BaismanNexradPixels_Gauges_d = 'GaugeNexradPixels.csv'
f_BaismanNexradPixels_Gauges_15min = 'GaugeNexradPixels_15minPrecip.csv'
f_BaismanNexradPixelsMap_Gauges = 'GaugeNexradPixels'
#RHESSys calibration and validation data
f_BaismanStreamflowCal = 'BaismanStreamflow_Feb2020Revised_Cal.txt'
f_BaismanStreamflowVal = 'BaismanStreamflow_Feb2020Revised_Val.txt'
f_POBRStreamflowCal = 'BaismanStreamflow_POBR_Feb2020Revised_Cal.txt'
f_POBRStreamflowVal = 'BaismanStreamflow_POBR_Feb2020Revised_Val.txt'
f_BaismanTNCal = 'TN_Feb2020Revised_Cal.txt'
f_BaismanTNCal_WRTDS = 'TN_Feb2020Revised_Cal_WRTDS.txt'
f_BaismanTNVal = 'TN_Feb2020Revised_Val.txt'
f_BaismanTNVal_WRTDS = 'TN_Feb2020Revised_Val_WRTDS.txt'
f_POBRTNCal = 'TN_POBR_Feb2020Revised_Cal.txt'
f_POBRTNCal_WRTDS = 'TN_POBR_Feb2020Revised_Cal_WRTDS.txt'
f_POBRTNVal = 'TN_POBR_Feb2020Revised_Val.txt'
f_POBRTNVal_WRTDS = 'TN_POBR_Feb2020Revised_Val_WRTDS.txt'
#RHESSys Climate Files
#Precipitation
f_Precip_All = 'AllTN_Feb2020Revised.rain'
f_Precip_SA = 'SA_Feb2020Revised.rain'
f_Precip_Cal = 'Cal_Feb2020Revised.rain'
f_Precip_Val = 'Val_Feb2020Revised.rain'
# Temperature
f_Tmin_All = 'AllTN_Feb2020Revised.tmin'
f_Tmax_All = 'AllTN_Feb2020Revised.tmax'
f_Tmin_SA = 'SA_Feb2020Revised.tmin'
f_Tmax_SA = 'SA_Feb2020Revised.tmax'
f_Tmin_Cal = 'Cal_Feb2020Revised.tmin'
f_Tmax_Cal = 'Cal_Feb2020Revised.tmax'
f_Tmin_Val = 'Val_Feb2020Revised.tmin'
f_Tmax_Val = 'Val_Feb2020Revised.tmax'
#WRTDS datasets
f_BaismanWRTDSInfo = 'WRTDS_INFO.csv'
f_POBRWRTDSInfo = 'WRTDS_INFO_POBR.csv'
f_BaismanWRTDS_TabInt = 'TabIntMod4_p5.txt'
f_BaismanWRTDS_TabYear = 'TabYearMod4_p4.txt'
f_BaismanWRTDS_TabLogQ = 'TabLogQMod4_p4.txt'
f_BaismanWRTDS_TabSinYear = 'TabSinYearMod4_p4.txt'
f_BaismanWRTDS_TabCosYear = 'TabCosYearMod4_p4.txt'
f_BaismanWRTDS_TabLogErr = 'TabLogErrMod4_p5.txt'
f_BaismanWRTDS_TabInt_QLQ = 'TabIntMod5QLQ_p5.txt'
f_BaismanWRTDS_TabYear_QLQ = 'TabYearMod5QLQ_p4.txt'
f_BaismanWRTDS_TabLogQ_QLQ = 'TabLogQMod5QLQ_p4.txt'
f_BaismanWRTDS_TabSinYear_QLQ = 'TabSinYearMod5QLQ_p4.txt'
f_BaismanWRTDS_TabCosYear_QLQ = 'TabCosYearMod5QLQ_p4.txt'
f_BaismanWRTDS_TabLogErr_QLQ = 'TabLogErrMod5QLQ_p5.txt'
f_BaismanWRTDS_TabLogQ2_QLQ = 'TabLogQ2Mod5QLQ_p4.txt'

f_POBRWRTDS_TabInt = 'TabInt_POBRMod5_p5.txt'
f_POBRWRTDS_TabYear = 'TabYear_POBRMod5_p4.txt'
f_POBRWRTDS_TabLogQ = 'TabLogQ_POBRMod5_p4.txt'
f_POBRWRTDS_TabSinYear = 'TabSinYear_POBRMod5_p4.txt'
f_POBRWRTDS_TabCosYear = 'TabCosYear_POBRMod5_p4.txt'
f_POBRWRTDS_TabLogErr = 'TabLogErr_POBRMod5_p5.txt'

#Export coefficient timeseries matrices
f_EC_Undev = "EC_Undev_Baisman.csv"
f_EC_Dev = "EC_Dev_Baisman.csv"

#Processed datasets
f_worldfile_p = 'worldfile_processed'
f_SmithWQ_TN_p = 'S_TN_Processed.csv'
f_SmithWQ_Q_p = 'S_Q_Processed.csv'
f_KenworthWQ_TN_p = 'K_TN_Processed.csv'
f_KenworthWQ_Q_p = 'K_Q_Processed.csv'
f_SmithSites_p = "SmithSites_Processed"
f_KenworthSites_p = "KenworthSites_Processed"

#Set project coordinate system----
#This is the coordinate system that all data will be plotted and written in
# It is not the coordinate system of your data (although it could be)
# EPSG codes from: https://spatialreference.org/ref/?page=2
pCRS = '+init=epsg:26918'

#Define a buffer to use for the ROI to download weather gauges, in map units----
ROIbuffPrecip = 20000
#For weather stations
ROIbuffWeather = 5000

#Set pixel resolution for raster data (m)----
res = 30
#Set plot limits - set to NULL to ignore use----
#BES TN Dates
xlim_BESTN = c(as.POSIXct('2000-01-01'), as.POSIXct('2020-01-01'))
xlim_BESTN_Full = c(as.Date("1995-01-01"), as.Date("2010-01-01"))

#Streamflow y axis limits
ylim_SF = c(0, 7000)

#Set number of Monte Carlo replicates for export coefficient model----
MCreps = 10000

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

#Load functions from repositories----
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
#Load modified WRTDS functions
#Fixme: move this file to the EnGauge repository. It belongs there or in the original EGRET repository.
setwd('C:\\Users\\js4yd\\OneDrive - University of Virginia\\RHESSys_ParameterSA')
source('WRTDS_modifiedFunctions.R')

#Load information from EnGauge downloads and place into the project coordinate system----
DEM = raster(x = paste0(dir_DEM_out, '\\', f_DEM_mosiac))
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
#Save new stream data
#writeOGR(obj = MDstreams, dsn = dir_streams, layer = 'MDstreams', driver = 'ESRI Shapefile')

#Process BES Water Quality Gauge Data----
setwd(dir_WChemData)
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

#   TN----
#Looking that the timeseries and histograms for each gauge, it is possible to see likely detection limits.
#Only interested in correcting sites that have a large proportion of the data as possibly censored. 
#Not able to tell if there are limits for small amounts of possibly censored data.
#Timeseries useful if limits changed over time.
i = 5
hist(log10(BES_WQ_Sites_TN[[i]]$NO3..mg.N.L.), breaks = 10000)
plot(BES_WQ_Sites_TN[[i]]$SortDateTime, BES_WQ_Sites_TN[[i]]$NO3..mg.N.L., log = 'y')
hist(log10(BES_WQ_Sites_TN[[i]]$TN..mg.N.L.), breaks = 10000)
plot(BES_WQ_Sites_TN[[i]]$SortDateTime, BES_WQ_Sites_TN[[i]]$TN..mg.N.L., log = 'y', ylim = c(0.01,5), xlim = xlim_BESTN)
rm(i)
#NoTNDL = c(1, 2, 3, 4, 6, 7, 8)
#TNDL = c(5)

#Add detection limit info to Pond Branch
BES_WQ_Sites_TN$POBR$DetectionLimit = NA
#All that are 0.01 are detection limits
BES_WQ_Sites_TN$POBR$DetectionLimit[BES_WQ_Sites_TN$POBR$TN..mg.N.L. == 0.01] = 1
#All that are 0.05 after 2012 are detection limits
BES_WQ_Sites_TN$POBR$DetectionLimit[which((BES_WQ_Sites_TN$POBR$TN..mg.N.L. == 0.05) & (BES_WQ_Sites_TN$POBR$SortDate > as.Date('2012-01-01')))] = 2

i = 5
plot(BES_WQ_Sites_TN[[i]]$SortDateTime, BES_WQ_Sites_TN[[i]]$TN..mg.N.L., log = 'y', ylim = c(0.01,5), xlim = xlim_BESTN)
par(new = T)
plot(BES_WQ_Sites_TN[[i]]$SortDateTime[is.na(BES_WQ_Sites_TN[[i]]$DetectionLimit) == FALSE], BES_WQ_Sites_TN[[i]]$TN..mg.N.L.[is.na(BES_WQ_Sites_TN[[i]]$DetectionLimit) == FALSE], log = 'y', ylim = c(0.01,5), xlim = xlim_BESTN, col = 'red')
rm(i)

#   TP - quite variable over time----
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
setwd(dir_WChem)
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
north.arrow(xb = 364000, yb = 4348000, len = 700, col = 'black', lab = 'N')
legend('bottomleft', title = 'Water Quality Types', legend = c('T Nitrogen Only', 'T Phosphorus Only', 'Both'), col = c('red', 'blue', 'purple'), pch = c(16,16,16))
legend(x = 334000, y = 4358000, title = 'Water Quality Sites', legend = c('At Stream Gauge', 'Other Location'), col = 'black', pch = c(4, 16))
dev.off()

#  Plot time series for each of the sites----
wd_BESN = paste0(dir_WChem, '/Nitrogen_BES')
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
wd_BESP = paste0(dir_WChem, '/Phosphorus_BES')
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
BES_TN = FillMissingDates_par(Dataset = BES_WQ_Sites_locs_TN, StationList = BES_WQ_Sites_TN_Gauged, Var = 'TN..mg.N.L.', Date = 'SortDate', gapType = 'd', site_no_D = 'Site', site_no_SL = 'Site', NoNAcols = c('USGSgauge', 'site_no'), NumCores = detectCores()-1)
BES_TP = FillMissingDates_par(Dataset = BES_WQ_Sites_locs_TP, StationList = BES_WQ_Sites_TP_Gauged, Var = 'TP..ugP.L.', Date = 'SortDate', gapType = 'd', site_no_D = 'Site', site_no_SL = 'Site', NoNAcols = c('USGSgauge', 'site_no'), NumCores = detectCores()-1)
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

rm(scaleRange, scaleBy, Pal)

#colorbar information
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

rm(scaleRange, scaleBy, Pal)

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
                         site_no_SL = 'Site', NoNAcols = c('USGSgauge', 'site_no'), NumCores = detectCores()-1)
BES_WQ_Sites_locs_TN = BES_TN_d2$Dataset
BES_TN_d = BES_TN_d2$StationList
rm(BES_TN_d2)

BES_TN_m2 = FillMissingDates_par(Dataset = BES_WQ_Sites_locs_TN, StationList = BES_TN_m, Var = 'TN..mg.N.L.', 
                             Date = 'YrMthDy', gapType = 'm', site_no_D = 'Site', 
                             site_no_SL = 'Site', NoNAcols = c('USGSgauge', 'site_no'), NumCores = detectCores()-1)
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
                             site_no_SL = 'Site', NoNAcols = c('USGSgauge', 'site_no'), NumCores = detectCores()-1)
BES_WQ_Sites_locs_TP = BES_TP_d2$Dataset
BES_TP_d = BES_TP_d2$StationList
rm(BES_TP_d2)

BES_TP_m2 = FillMissingDates_par(Dataset = BES_WQ_Sites_locs_TP, StationList = BES_TP_m, Var = 'TP..ugP.L.', 
                             Date = 'YrMthDy', gapType = 'm', site_no_D = 'Site', 
                             site_no_SL = 'Site', NoNAcols = c('USGSgauge', 'site_no'), NumCores = detectCores()-1)
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
       ylim = c(0, 10), xlim = xlim_BESTN_Full)
  if(length(which(is.na(BES_TN[[i]]$TN..mg.N.L.))) > 0){
    #Plot where there are NA days that do not correspond to detection limits
    par(new = TRUE)
    plot(y = rep(0, length(which(is.na(BES_TN[[i]]$TN..mg.N.L.)))), x = as.Date(BES_TN[[i]]$SortDate[which(is.na(BES_TN[[i]]$TN..mg.N.L.))]), pch = 4, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = c(0, 0.5), xlim = xlim_BESTN_Full, axes = FALSE,
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
       ylim = c(0, 10), xlim = xlim_BESTN_Full)
  if(length(which(is.na(BES_TN[[i]]$TN..mg.N.L.))) > 0){
    #Plot where there are NA days that do not correspond to detection limits
    par(new = TRUE)
    plot(y = rep(0, length(which(is.na(BES_TN[[i]]$TN..mg.N.L.)))), x = as.Date(BES_TN[[i]]$SortDate[which(is.na(BES_TN[[i]]$TN..mg.N.L.))]), pch = 4, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = c(0, 0.5), xlim = xlim_BESTN_Full, axes = FALSE,
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
  scatterHist_mod(x = log10(BES_TN_d[[i]]$Flow[BES_TN_d[[i]]$Flow > 0]), y = BES_TN_d[[i]]$TN..mg.N.L.[BES_TN_d[[i]]$Flow > 0], 
                  density = FALSE, freq = FALSE, x.breaks = 100, y.breaks = 100, correl = FALSE, 
                  xlab = expression(paste('log'[10] ~ '(Flow [cfs])')), ylab = 'Total Nitrogen (mg/L as N)', 
                  cex.lab = 1.5, cex.axis = 1.5, title = '', pch = 16, ellipse = FALSE, smooth = FALSE, 
                  col = BES_TN_d[[i]]$cols[BES_TN_d[[i]]$Flow > 0])
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
       ylim = c(0, 10), xlim = xlim_BESTN_Full)
  if(length(which(is.na(BES_TP[[i]]$TP..ugP.L.))) > 0){
    #Plot where there are NA days that do not correspond to detection limits
    par(new = TRUE)
    plot(y = rep(0, length(which(is.na(BES_TP[[i]]$TP..ugP.L.)))), x = as.Date(BES_TP[[i]]$SortDate[which(is.na(BES_TP[[i]]$TP..ugP.L.))]), pch = 4, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = c(0, 0.5), xlim = xlim_BESTN_Full, axes = FALSE,
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
       ylim = c(0, 10), xlim = xlim_BESTN_Full)
  if(length(which(is.na(BES_TP[[i]]$TP..ugP.L.))) > 0){
    #Plot where there are NA days that do not correspond to detection limits
    par(new = TRUE)
    plot(y = rep(0, length(which(is.na(BES_TP[[i]]$TP..ugP.L.)))), x = as.Date(BES_TP[[i]]$SortDate[which(is.na(BES_TP[[i]]$TP..ugP.L.))]), pch = 4, cex = 0.3,
         xlab = '', 
         ylab = "",
         ylim = c(0, 0.5), xlim = xlim_BESTN_Full, axes = FALSE,
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
  scatterHist_mod(x = log10(BES_TP_d[[i]]$Flow[BES_TP_d[[i]]$Flow > 0]), y = BES_TP_d[[i]]$TP..ugP.L.[BES_TP_d[[i]]$Flow > 0], 
                  density = FALSE, freq = FALSE, x.breaks = 100, y.breaks = 100, correl = FALSE, 
                  xlab = expression(paste('log'[10] ~ '(Flow [cfs])')), ylab = expression(paste('Total Phosphorus ( ', mu, 'g/L as P)')), 
                  cex.lab = 1.5, cex.axis = 1.5, title = '', pch = 16, ellipse = FALSE, smooth = FALSE, 
                  col = BES_TP_d[[i]]$cols[BES_TP_d[[i]]$Flow > 0])
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

# Hypsometric plot for DEM in the subwatersheds----
#Fixme: Can shade in the areas above each gauge or select outlet points as horizontal steps
# Allocation of BMPs could be in certain elevation zones
# Couple with map next to it of the elevations above certain thresholds

#Clip DEM to the ROI and make a spatial grid dataframe
DEM_DR = as(mask(DEM, DeadRun_Outlet), 'SpatialGridDataFrame')
DEM_BR = as(mask(DEM, BaisRun_Outlet), 'SpatialGridDataFrame')
DEM_SL = as(mask(DEM, ScottsLevel), 'SpatialGridDataFrame')
#Plot curve
setwd(dir_DEM_out)
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

# Fixme: Search for flow-normalized water chemistry outliers----
# #Select only the non-NA values
# t = BES_TN_d[[i]][which((is.na(BES_TN_d[[i]]$Flow) == FALSE) & (is.na(BES_TN_d[[i]]$TN..mg.N.L.) == FALSE)),c('Flow', 'TN..mg.N.L.')]
# mu = apply(X = t, MARGIN = 2, FUN = mean)
# #Compute the covariance matrix for only the non-NA values
# s = cov(t)
# mahalanobis(x = t, center = FALSE, cov = s)

# Fixme: Check correlation between N, P at all sites to see if it's worthwhile to use both for calibration of model----
# FOR NOW, NOT NEEDED BECAUSE RHESSYS DOESN'T PROVIDE P OUTPUT
# Color by time of year
# for (i in 1:length(BES_TN_d)){
#   #Gather all of the dates that match in both timeseries
#   
# }
#plot(BES_TN_d$POBR$TN..mg.N.L.[362:nrow(BES_TN_d$POBR)][which((is.na(BES_TN_d$POBR$TN..mg.N.L.[362:nrow(BES_TN_d$POBR)]) == FALSE) & (is.na(BES_TP_d$POBR$TP..ugP.L.) == FALSE))], 
#     BES_TP_d$POBR$TP..ugP.L.[which((is.na(BES_TN_d$POBR$TN..mg.N.L.[362:nrow(BES_TN_d$POBR)]) == FALSE) & (is.na(BES_TP_d$POBR$TP..ugP.L.) == FALSE))], log = 'xy')
#plot(BES_TN_d$BARN$TN..mg.N.L.[which((is.na(BES_TN_d$BARN$TN..mg.N.L.) == FALSE) & (is.na(BES_TP_d$BARN$TP..ugP.L.) == FALSE))], 
#     BES_TP_d$BARN$TP..ugP.L.[which((is.na(BES_TN_d$BARN$TN..mg.N.L.) == FALSE) & (is.na(BES_TP_d$BARN$TP..ugP.L.) == FALSE))], log = 'xy')

#Fixme: Check for gaps in the TN/NO3 series that can be filled in with TP correlation with TN

# Fixme: correlation plots of streamflow and N across watershed for different sites----
#       need to have matching timeseries to do this. Hydropairs function may work.
#       can help the spatial prediction
# Write water quality data to files----
setwd(dir_WChem)
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
setwd(dir_precipData)
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
        tvec[t] = paste(as.POSIXct(paste(strsplit(f$Date_Time_.EST.[which(f$DateTime == 'NA')][t], split = ' ', fixed = TRUE)[[1]][1], 
                                         strsplit(as.character(as.POSIXct(strsplit(f$Date_Time_.EST.[which(f$DateTime == 'NA')][t], split = ' ', fixed = TRUE)[[1]][2], format = '%H:%M') + 1*3600), 
                                                  split = ' ', fixed = TRUE)[[1]][2]), format = '%m/%d/%y %H:%M'))
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
        tvec[t] = paste(as.POSIXct(paste(strsplit(f$Date_Time_.EST.[which(f$DateTime == 'NA')][t], split = ' ', fixed = TRUE)[[1]][1], 
                                         strsplit(as.character(as.POSIXct(strsplit(f$Date_Time_.EST.[which(f$DateTime == 'NA')][t], split = ' ', fixed = TRUE)[[1]][2], format = '%H:%M') + 1*3600), 
                                                  split = ' ', fixed = TRUE)[[1]][2]), format = '%m/%d/%Y %H:%M'))
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
setwd(dir_precip)
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

#  Fixme: Compute rainfall rates----
#  compute rates and compare for each gauge. Use for error identification

#  Aggregate into total precip per day, month, year---- 
BES_Precip_agg = aggregateTimesteps(StationList = BES_PrecipList, aggVal = c('d', 'm', 'a'), aggVar = 'Precipitation_.mm.', date = 'SortDate', site = 'Rain_Gauge_ID', fun = 'sum')
BES_Precip_d = BES_Precip_agg$daily
BES_Precip_m = BES_Precip_agg$mthyr
BES_Precip_a = BES_Precip_agg$ann
rm(BES_Precip_agg)

#  Handle missing data in the daily, monthly, and annual aggregated timeseries----
BES_Precip_d2 = FillMissingDates_par(Dataset = BES_Precip_locs, StationList = BES_Precip_d, Var = 'Precipitation_.mm.', 
                             Date = 'SortDate', gapType = 'd', site_no_D = 'Gauge', 
                             site_no_SL = 'Rain_Gauge_ID', NoNAcols = c('Site','Rain_Gauge_ID'), NumCores = detectCores()-1)
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

#  Aggregate the gauges by averaging the 2 gauges per site----
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

#   Compare site daily measurements----
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

#  Fixme: Do some diagnostics on the sites to come up with a single timeseries for each site - completed below for Oregon Ridge----
#  Fixme: check if the high bias days are also biased on the previous and following day----
# Fixme: color the bias by the rainfall rate (> 1.5"/hr is prone to errors)
# NOAA weather station data----
#Get the FIPS code for Maryland
#MDfips = rnoaa::fipscodes[which(rnoaa::fipscodes$state == 'Maryland'),]
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
wd_NOAA = paste0(dir_precip, '\\NOAA')
dir.create(path = wd_NOAA)
file.copy(from = user_cache_dir(appname = 'rnoaa', version = NULL, opinion = TRUE, expand = TRUE), to = wd_NOAA, recursive = TRUE)
setwd(wd_NOAA)

#Select only the met stations in a 5 km buffer for now
ROI_buffer_5km = buffer(ROI, width = ROIbuffWeather)
ROI_buffer_5km_WGS = spTransform(ROI_buffer_5km, CRS('+init=epsg:4326'))
NOAAstations_5km_locs = AllNOAAstations[ROI_buffer_5km_WGS,]
NOAAstations_5km_locs = spTransform(NOAAstations_5km_locs, CRS(pCRS))

MetStations_5km = MetStations[(names(MetStations) %in% NOAAstations_5km_locs$id)]

rm(AllNOAAstations)

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
MetStations_5km_Precip = FillMissingDates_par(Dataset = NOAAstations_5km_locs[NOAAstations_5km_locs$element == 'PRCP',], 
                                              StationList = MetStations_5km[names(MetStations_5km) %in% NOAAstations_5km_locs[NOAAstations_5km_locs$element == 'PRCP',]$id], 
                                              Var = 'prcp', Date = 'date', gapType = 'd', site_no_D = 'id', site_no_SL = 'id', NoNAcols = 'id', NumCores = detectCores()-1)
#Extract data from the function return
NOAAstations_5km_locs_Precip = MetStations_5km_Precip$Dataset
MetStations_5km_Precip = MetStations_5km_Precip$StationList

MetStations_5km_TMAX = FillMissingDates_par(Dataset = NOAAstations_5km_locs[NOAAstations_5km_locs$element == 'TMAX',], 
                                            StationList = MetStations_5km[names(MetStations_5km) %in% NOAAstations_5km_locs[NOAAstations_5km_locs$element == 'TMAX',]$id], 
                                            Var = 'tmax', Date = 'date', gapType = 'd', site_no_D = 'id', site_no_SL = 'id', NoNAcols = 'id', NumCores = detectCores()-1)
#Extract data from the function return
NOAAstations_5km_locs_TMAX = MetStations_5km_TMAX$Dataset
MetStations_5km_TMAX = MetStations_5km_TMAX$StationList

MetStations_5km_TMIN = FillMissingDates_par(Dataset = NOAAstations_5km_locs[NOAAstations_5km_locs$element == 'TMIN',], 
                                            StationList = MetStations_5km[names(MetStations_5km) %in% NOAAstations_5km_locs[NOAAstations_5km_locs$element == 'TMIN',]$id], 
                                            Var = 'tmin', Date = 'date', gapType = 'd', site_no_D = 'id', site_no_SL = 'id', NoNAcols = 'id', NumCores = detectCores()-1)
#Extract data from the function return
NOAAstations_5km_locs_TMIN = MetStations_5km_TMIN$Dataset
MetStations_5km_TMIN = MetStations_5km_TMIN$StationList

#  Plot the met station timeseries----
#   Precip----
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

#   Max and Min Temperature----
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

# Map of record lengths for BES precip stations----
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

rm(scaleRange, scaleBy, Pal)

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

rm(scaleRange, scaleBy, Pal)

# Fixme: ACF for precip and temp datasets----
# Fixme: Spatial predicton of precip and temperature----
# Fixme: evaluate outliers for precip and temp (possible multivariate outliers) - make scatterplots of precip and temp----
# Write met station data to files----
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

#Fixme: Download climate indices (ENSO, PDO, etc.) and evaluate timeseries relative to those----
#Create Observation timeseries for RHESSys runs----
# Baisman USGS Streamflow and TN----
setwd(wd_RHESSysObs_BARN)

#Streamflow and Nitrogen are both available in the BES_TN_d list
#Pond branch had several streamflow revisions in Sept 2019 as a result of a new stage-discharge relation. Updates here are for Feb. 2020.

#Earliest possible start date for calibration
CalEarliestDate = max(min(as.Date(BES_TN_d$BARN$SortDate, origin="1970-01-01")), min(as.Date(BES_TN_d$POBR$SortDate, origin="1970-01-01")))

#Latest possible end date for validation
ValLatestDate = min(max(as.Date(BES_TN_d$BARN$SortDate, origin="1970-01-01")), max(as.Date(BES_TN_d$POBR$SortDate, origin="1970-01-01")))

#Total number of samples
BARN_NumTNSamps = length(which(!is.na(BES_TN_d$BARN$TN..mg.N.L.)))
POBR_NumTNSamps = length(which(!is.na(BES_TN_d$POBR$TN..mg.N.L.)))
#Want to have at least 15-20% of record for validation. BARN has least number of data points, so use that as date reference.
ValDate_15pct = as.Date(BES_TN_d$BARN$SortDate[which(!is.na(BES_TN_d$BARN$TN..mg.N.L.))[floor(BARN_NumTNSamps*.85)]], origin="1970-01-01")
ValDate_20pct = as.Date(BES_TN_d$BARN$SortDate[which(!is.na(BES_TN_d$BARN$TN..mg.N.L.))[floor(BARN_NumTNSamps*.80)]], origin="1970-01-01")

#Select the water year start between these values on Oct. 1, 2013 
#TN about 20% of observations for validation at BARN
1 - length(which((!is.na(BES_TN_d$BARN$TN..mg.N.L.)) & (as.Date(BES_TN_d$BARN$SortDate, origin="1970-01-01") < '2013-10-01')))/length(which((!is.na(BES_TN_d$BARN$TN..mg.N.L.))))

#Flow about 20% of observations for validation at BARN
1 - length(which(as.Date(BES_TN_d$BARN$SortDate, origin="1970-01-01") < '2013-10-01'))/length(BES_TN_d$BARN$SortDate)

#  Save text files for RHESSys calibration and validation----
#Streamflow - Calibration
#Get indices for calibration and validation
IndDatesCalBARN = which((as.Date(BES_TN_d$BARN$SortDate, origin='1970-01-01') < "2013-10-01") & (as.Date(BES_TN_d$BARN$SortDate, origin='1970-01-01') >= CalEarliestDate))
IndDatesValBARN = which((as.Date(BES_TN_d$BARN$SortDate, origin='1970-01-01') >= "2008-11-15") & (as.Date(BES_TN_d$BARN$SortDate, origin='1970-01-01') <= ValLatestDate))
#File for saving
StreamCal = data.frame(Date = as.character(as.Date(BES_TN_d$BARN$SortDate, origin='1970-01-01')[IndDatesCalBARN]), 
                       Flow = BES_TN_d$BARN$Flow[IndDatesCalBARN], 
                       stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = StreamCal, file = f_BaismanStreamflowCal, sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)
#Streamflow - Validation - needs 5 years of spinup before validation starts.
StreamVal = data.frame(Date = as.character(as.Date(BES_TN_d$BARN$SortDate, origin='1970-01-01')[IndDatesValBARN]), 
                       Flow = BES_TN_d$BARN$Flow[IndDatesValBARN], 
                       stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = StreamVal, file = f_BaismanStreamflowVal, sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)

#Streamflow - Calibration - Pond Branch
setwd(wd_RHESSysObs_POBR)
IndDatesCalPOBR = which((as.Date(BES_TN_d$POBR$SortDate, origin='1970-01-01') < "2013-10-01") & (as.Date(BES_TN_d$POBR$SortDate, origin='1970-01-01') >= CalEarliestDate))
IndDatesValPOBR = which((as.Date(BES_TN_d$POBR$SortDate, origin='1970-01-01') >= "2008-11-15") & (as.Date(BES_TN_d$POBR$SortDate, origin='1970-01-01') <= ValLatestDate))
StreamCal_POBR = data.frame(Date = as.character(as.Date(BES_TN_d$POBR$SortDate, origin='1970-01-01')[IndDatesCalPOBR]), 
                            Flow = BES_TN_d$POBR$Flow[IndDatesCalPOBR], 
                            stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = StreamCal_POBR, file = f_POBRStreamflowCal, sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)
#Streamflow - Validation - Pond Branch
StreamVal_POBR = data.frame(Date = as.character(as.Date(BES_TN_d$POBR$SortDate, origin='1970-01-01')[IndDatesValPOBR]), 
                            Flow = BES_TN_d$POBR$Flow[IndDatesValPOBR], 
                            stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = StreamVal_POBR, file = f_POBRStreamflowVal, sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)

#TN - Calibration
setwd(wd_RHESSysObs_BARN)
TNCal = data.frame(Date = as.character(as.Date(BES_TN_d$BARN$SortDate, origin='1970-01-01')[IndDatesCalBARN]), 
                   TN = BES_TN_d$BARN$TN..mg.N.L.[IndDatesCalBARN], 
                   stringsAsFactors = FALSE)
#WRTDS needs another column to run their functions.
TNCal_WRTDS = data.frame(Date = as.character(as.Date(BES_TN_d$BARN$SortDate, origin='1970-01-01')[IndDatesCalBARN]), 
                         remark = rep('', length(BES_TN_d$BARN$TN..mg.N.L.[IndDatesCalBARN])), 
                         TN = BES_TN_d$BARN$TN..mg.N.L.[IndDatesCalBARN], 
                         stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = TNCal, file = f_BaismanTNCal, sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
write.table(x = TNCal_WRTDS, file = f_BaismanTNCal_WRTDS, sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)
#TN - Validation
TNVal = data.frame(Date = as.character(as.Date(BES_TN_d$BARN$SortDate, origin='1970-01-01')[IndDatesValBARN]), 
                   TN = BES_TN_d$BARN$TN..mg.N.L.[IndDatesValBARN], 
                   stringsAsFactors = FALSE)
TNVal_WRTDS = data.frame(Date = as.character(as.Date(BES_TN_d$BARN$SortDate, origin='1970-01-01')[IndDatesValBARN]), 
                         remark = rep('', length(BES_TN_d$BARN$TN..mg.N.L.[IndDatesValBARN])), 
                         TN = BES_TN_d$BARN$TN..mg.N.L.[IndDatesValBARN], 
                         stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = TNVal, file = f_BaismanTNVal, sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
write.table(x = TNVal_WRTDS, file = f_BaismanTNVal_WRTDS, sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)

#TN - Calibration - Pond Branch
setwd(wd_RHESSysObs_POBR)
TNCal_POBR = data.frame(Date = as.character(as.Date(BES_TN_d$POBR$SortDate, origin='1970-01-01')[IndDatesCalPOBR]), 
                        TN = BES_TN_d$POBR$TN..mg.N.L.[IndDatesCalPOBR], 
                        stringsAsFactors = FALSE)
TNCal_POBR_WRTDS = data.frame(Date = as.character(as.Date(BES_TN_d$POBR$SortDate, origin='1970-01-01')[IndDatesCalPOBR]), 
                              remark = rep('', length(BES_TN_d$POBR$TN..mg.N.L.[IndDatesCalPOBR])), 
                              TN = BES_TN_d$POBR$TN..mg.N.L.[IndDatesCalPOBR], 
                              stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = TNCal_POBR, file = f_POBRTNCal, sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
write.table(x = TNCal_POBR_WRTDS, file = f_POBRTNCal_WRTDS, sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)
#TN - Validation - Pond Branch
TNVal_POBR = data.frame(Date = as.character(as.Date(BES_TN_d$POBR$SortDate, origin='1970-01-01')[IndDatesValPOBR]), 
                        TN = BES_TN_d$POBR$TN..mg.N.L.[IndDatesValPOBR], 
                        stringsAsFactors = FALSE)
TNVal_POBR_WRTDS = data.frame(Date = as.character(as.Date(BES_TN_d$POBR$SortDate, origin='1970-01-01')[IndDatesValPOBR]), 
                              remark = rep('', length(BES_TN_d$POBR$TN..mg.N.L.[IndDatesValPOBR])), 
                              TN = BES_TN_d$POBR$TN..mg.N.L.[IndDatesValPOBR], 
                              stringsAsFactors = FALSE)
options(scipen = 999)
write.table(x = TNVal_POBR, file = f_POBRTNVal, sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
write.table(x = TNVal_POBR_WRTDS, file = f_POBRTNVal_WRTDS, sep = '\t', row.names = FALSE, col.names = TRUE, fileEncoding = 'UTF-8')
options(scipen = 0)

#Create climate station timeseries for RHESSys runs----
# Rainfall for Oregon Ridge gauge at Baisman Run----
#Laurence Lin's Oregon Ridge gauge dataset
#OR_Lin = read.table(file = paste0(wd_LinData, "\\", f_LinRain), sep = '\t', stringsAsFactors = FALSE)

#x = as.Date('2006-01-01')
#x = x + seq(0, nrow(OR_Lin)-2, 1)
#plot(x, OR_Lin[-1,], type ='l')

#Compare Lin to BES OR precip
#plot(BES_Precip_Avg_d$`Oregon Ridge Park`$mean/1000, as.numeric(OR_Lin[-1,][1197:4377]), log = 'xy')
#Laurence used OR rain gauge 1, did not average the two gauges
#plot(BES_Precip_d$WXORDG_RG1$Precipitation_.mm./1000, as.numeric(OR_Lin[-1,][1197:4377]), log = 'xy')

#Compare Lin to MD Science Center USW00093784
#plot(MetStations_5km$USW00093784[2742:3937,]$prcp/10000, as.numeric(OR_Lin[-1,][1:1196]), log = 'xy')

#Compare Lin to BALTIMORE WASH INTL AP USW00093721
BWI = meteo_tidy_ghcnd(stationid = 'USW00093721', var = 'all', keep_flags = TRUE)
#BWI_LinDateStart = which(BWI$date == '2006-01-01')
#BWI_LinDateEnd = which(BWI$date == max(x))
#plot(BWI[BWI_LinDateStart:BWI_LinDateEnd,]$prcp/10000, as.numeric(OR_Lin[-1,]), log = 'xy')
#rm(x)

#Compare BWI and MD Sci Center
BWI_CalDateStart = which(BWI$date == CalEarliestDate)
MDSci_CalDateStart = which(MetStations_5km$USW00093784$date == CalEarliestDate)
#Not sure why these specific dates were selected.
plot(BWI[BWI_CalDateStart:(BWI_CalDateStart+6697),]$prcp/10, MetStations_5km$USW00093784[MDSci_CalDateStart:(MDSci_CalDateStart+6697),]$prcp/10, pch = 16, cex = 0.2, log='xy')
lines(c(1,2000), c(1,2000), col = 'red')

#Baisman Oregon Ridge Average vs. MD Science Center precip - OR ridge starts in 2009
MDSci_ORStartDate = which(MetStations_5km$USW00093784$date == BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[1])
MDSci_OREndDate = which(MetStations_5km$USW00093784$date == max(BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate))
plot(MetStations_5km$USW00093784[MDSci_ORStartDate:MDSci_OREndDate,]$prcp/10, BES_Precip_Avg_d$`Oregon Ridge Park`$mean, log = 'xy', xlab = 'MD Science Center', ylab = 'Oregon Ridge Gauge - Average')
lines(c(1,2000), c(1,2000), col = 'red')

#Baisman Oregon Ridge Average vs. BWI
BWI_ORStartDate = which(BWI$date == BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[1])
BWI_OREndDate = which(BWI$date == max(BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate))
plot(y = BES_Precip_Avg_d$`Oregon Ridge Park`$mean, x = BWI[BWI_ORStartDate:BWI_OREndDate,]$prcp/10, log = 'xy', xlab = 'BWI Airport', ylab = 'Oregon Ridge Gauge - Average')
lines(c(1,2000), c(1,2000), col = 'red')

#Oregon Ridge is More correlated with the MD Science Center - Use that station to fill in before OR Ridge data are available.
#From calibration start date to the start of OR ridge data
RainTimeseries = MetStations_5km$USW00093784[MDSci_CalDateStart:(MDSci_ORStartDate-1),]$prcp/10

#Fill in NAs in the MD Science Center data with the BWI data. Most are zeros or a small amount of rain
RainTimeseries[which(is.na(MetStations_5km$USW00093784[MDSci_CalDateStart:(MDSci_ORStartDate-1),]$prcp))] = BWI$prcp[which(BWI$date %in% MetStations_5km$USW00093784[MDSci_CalDateStart:(MDSci_ORStartDate-1),]$date[which(is.na(MetStations_5km$USW00093784[MDSci_CalDateStart:(MDSci_ORStartDate-1),]$prcp))])]/10

#Add the corrected/fixed Oregon Ridge data to this timeseries
#Correct Oregon Ridge dataset
#Rain gauge 2 vs 1
plot(BES_Precip_d$WXORDG_RG1$Precipitation_.mm., BES_Precip_d$WXORDG_RG2$Precipitation_.mm., ylim = c(0,150), xlim = c(0,150), xlab = 'Gauge 1', ylab = 'Gauge 2')
par(new = TRUE)
#More than 15 mm different
plot(BES_Precip_Avg_d$`Oregon Ridge Park`$P_WXORDG_RG1[which(abs(BES_Precip_Avg_d$`Oregon Ridge Park`$GaugeDiff) > 15)], BES_Precip_Avg_d$`Oregon Ridge Park`$P_WXORDG_RG2[which(abs(BES_Precip_Avg_d$`Oregon Ridge Park`$GaugeDiff) > 15)], col = 'red', ylim = c(0,150), xlim = c(0,150), xlab = '', ylab = '', axes = FALSE)
par(new = TRUE)
#More than 20 mm different
plot(BES_Precip_Avg_d$`Oregon Ridge Park`$P_WXORDG_RG1[which(abs(BES_Precip_Avg_d$`Oregon Ridge Park`$GaugeDiff) > 20)], BES_Precip_Avg_d$`Oregon Ridge Park`$P_WXORDG_RG2[which(abs(BES_Precip_Avg_d$`Oregon Ridge Park`$GaugeDiff) > 20)], col = 'blue', ylim = c(0,150), xlim = c(0,150), xlab = '', ylab = '', axes = FALSE)
#Full cloud cover lines
lines(c(.0254*1000,.0254*1000), c(0,150))
lines(y = c(.0254*1000,.0254*1000), x = c(0,150))

#There are two instances of very high precip at one gauge and no precip at the other. 
# Several others differ by more than 20 mm. Use a buffer of "If greater than 20 mm different, take the larger instead of the mean"
BES_Precip_Avg_d$`Oregon Ridge Park`$Selected = BES_Precip_Avg_d$`Oregon Ridge Park`$mean
BES_Precip_Avg_d$`Oregon Ridge Park`$Selected[which(abs(BES_Precip_Avg_d$`Oregon Ridge Park`$GaugeDiff) > 20)] = apply(X = rbind(BES_Precip_Avg_d$`Oregon Ridge Park`$P_WXORDG_RG1[which(abs(BES_Precip_Avg_d$`Oregon Ridge Park`$GaugeDiff) > 20)], 
                                                                                                                                 BES_Precip_Avg_d$`Oregon Ridge Park`$P_WXORDG_RG2[which(abs(BES_Precip_Avg_d$`Oregon Ridge Park`$GaugeDiff) > 20)]), MARGIN = 2, FUN = max)
#Check the NA rain dates with the MD Science Center precip. It's possible that some of these NA dates were days that the gauges did not record.
ORGauge_NADates = BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[which(is.na(BES_Precip_Avg_d$`Oregon Ridge Park`$Selected))]
plot(MetStations_5km$USW00093784[MDSci_ORStartDate:MDSci_OREndDate,]$date[MetStations_5km$USW00093784[MDSci_ORStartDate:MDSci_OREndDate,]$date %in% ORGauge_NADates],
     MetStations_5km$USW00093784[MDSci_ORStartDate:MDSci_OREndDate,]$prcp[MetStations_5km$USW00093784[MDSci_ORStartDate:MDSci_OREndDate,]$date %in% ORGauge_NADates]/10)
#Any storm larger than 1 in of rain will be considered falling in OR Ridge, too.
FillInDates = MetStations_5km$USW00093784[MDSci_ORStartDate:MDSci_OREndDate,]$date[MetStations_5km$USW00093784[MDSci_ORStartDate:MDSci_OREndDate,]$date %in% ORGauge_NADates][which(MetStations_5km$USW00093784[MDSci_ORStartDate:MDSci_OREndDate,]$prcp[MetStations_5km$USW00093784[MDSci_ORStartDate:MDSci_OREndDate,]$date %in% ORGauge_NADates]/10 >= 25.4)]
FillInPrecips = MetStations_5km$USW00093784[MDSci_ORStartDate:MDSci_OREndDate,]$prcp[MetStations_5km$USW00093784[MDSci_ORStartDate:MDSci_OREndDate,]$date %in% ORGauge_NADates][which(MetStations_5km$USW00093784[MDSci_ORStartDate:MDSci_OREndDate,]$prcp[MetStations_5km$USW00093784[MDSci_ORStartDate:MDSci_OREndDate,]$date %in% ORGauge_NADates]/10 >= 25.4)]
BES_Precip_Avg_d$`Oregon Ridge Park`$Selected[which(BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate %in% FillInDates)] = FillInPrecips/10

#Make all other NA values equal to 0
BES_Precip_Avg_d$`Oregon Ridge Park`$Selected[is.na(BES_Precip_Avg_d$`Oregon Ridge Park`$Selected)] = 0

RainTimeseries = c(RainTimeseries, BES_Precip_Avg_d$`Oregon Ridge Park`$Selected)
#Convert to m
RainTimeseries = RainTimeseries/1000
RainDates = seq(as.Date(MetStations_5km$USW00093784[MDSci_CalDateStart,]$date), max(as.Date(BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate)), 1) 

#  Save text file for RHESSys input (m units)----
setwd(wd_clim)
#FullTimeseries
options(scipen = 999)
write.table(x = c("1999 11 15 1", RainTimeseries), file = f_Precip_All, sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote = FALSE)
options(scipen = 0)

#SA
#1999-11-15 to 2010-09-30
options(scipen = 999)
write.table(x = c("1999 11 15 1", RainTimeseries[1:which(RainDates == '2010-09-30')]), file = f_Precip_SA, sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote = FALSE)
options(scipen = 0)

#Calibration
#1999-11-15 to 2013-09-30
options(scipen = 999)
write.table(x = c("1999 11 15 1", RainTimeseries[1:which(RainDates == '2013-09-30')]), file = f_Precip_Cal, sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote = FALSE)
options(scipen = 0)

#Validation
#2008-11-15 to 2017-04-01
options(scipen = 999)
write.table(x = c("2008 11 15 1", RainTimeseries[which(RainDates == '2008-11-15'):which(RainDates == '2017-04-01')]), file = f_Precip_Val, sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote = FALSE)
options(scipen = 0)

# Temperature min and max from MD Science Center and BWI Airport ----
#Laurence's dataset
#OR_Lin_tmin = read.table(file = paste0(wd_LinData, '\\', f_LinTmin), sep = '\t', stringsAsFactors = FALSE)
#OR_Lin_tmax = read.table(file = paste0(wd_LinData, '\\', f_LinTmax), sep = '\t', stringsAsFactors = FALSE)

#x = as.Date('1996-04-01')
#x = x + seq(0, nrow(OR_Lin_tmin)-2, 1)
#plot(x, as.numeric(OR_Lin_tmin[-1,]), type ='l', col = 'blue', ylim = c(-20,50), xlab = 'Years', ylab = expression(paste('Temperature (',degree,'C)')))
#par(new=TRUE)
#plot(x, as.numeric(OR_Lin_tmax[-1,]), type ='l', col = 'red', ylim = c(-20,50), axes=FALSE, xlab = '', ylab = '')

#Compare Lin to MD Science Center USW00093784
#plot(MetStations_5km$USW00093784[MDSci_CalDateStart:(MDSci_CalDateStart+6697),]$tmin/10, as.numeric(OR_Lin_tmin[-1,][1324:(nrow(OR_Lin_tmin)-1)]))
#plot(MetStations_5km$USW00093784[MDSci_CalDateStart:(MDSci_CalDateStart+6697),]$tmax/10, as.numeric(OR_Lin_tmax[-1,][1324:(nrow(OR_Lin_tmin)-1)]))

#Compare Lin to BALTIMORE WASH INTL AP USW00093721
#Laurence used BWI
#plot(BWI[20730:(20730+8020),]$tmin/10, as.numeric(OR_Lin_tmin[-1,]))
#plot(BWI[BWI_CalDateStart:(BWI_CalDateStart+6697),]$tmin/10, as.numeric(OR_Lin_tmin[-1,][1324:(nrow(OR_Lin_tmin)-1)]))
#plot(BWI[20730:(20730+8020),]$tmax/10, as.numeric(OR_Lin_tmax[-1,]))
#plot(BWI[(BWI_CalDateStart):(BWI_CalDateStart+6697),]$tmax/10, as.numeric(OR_Lin_tmax[-1,][1324:(nrow(OR_Lin_tmax)-1)]))


#Compare BWI and MD Sci Center
plot(BWI[(BWI_CalDateStart):(BWI_CalDateStart+6697),]$tmin/10, MetStations_5km$USW00093784[MDSci_CalDateStart:(MDSci_CalDateStart+6697),]$tmin/10, pch = 16, cex = 0.2)
#max better correlated
plot(BWI[(BWI_CalDateStart):(BWI_CalDateStart+6697),]$tmax/10, MetStations_5km$USW00093784[MDSci_CalDateStart:(MDSci_CalDateStart+6697),]$tmax/10, pch = 16, cex = 0.2)

#Make temperature timeseries for MD Sci Center
TminTimeseries = MetStations_5km$USW00093784[MDSci_CalDateStart:MDSci_OREndDate,]$tmin/10
TmaxTimeseries = MetStations_5km$USW00093784[MDSci_CalDateStart:MDSci_OREndDate,]$tmax/10
TempDates = MetStations_5km$USW00093784[MDSci_CalDateStart:MDSci_OREndDate,]$date

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

#  Save text file for RHESSys input (deg. C units)----
setwd(wd_clim)
#Full Timeseries
options(scipen = 999)
write.table(x = c("1999 11 15 1", TminTimeseries), file = f_Tmin_All, sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
write.table(x = c("1999 11 15 1", TmaxTimeseries), file = f_Tmax_All, sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
options(scipen = 0)

#SA
#1999-11-15 to 2010-09-30
options(scipen = 999)
write.table(x = c("1999 11 15 1", TminTimeseries[1:which(TempDates == '2010-09-30')]), file = f_Tmin_SA, sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
write.table(x = c("1999 11 15 1", TmaxTimeseries[1:which(TempDates == '2010-09-30')]), file = f_Tmax_SA, sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
options(scipen = 0)

#Calibration
#1999-11-15 to 2013-09-30
options(scipen = 999)
write.table(x = c("1999 11 15 1", TminTimeseries[1:which(TempDates == '2013-09-30')]), file = f_Tmin_Cal, sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
write.table(x = c("1999 11 15 1", TmaxTimeseries[1:which(TempDates == '2013-09-30')]), file = f_Tmax_Cal, sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
options(scipen = 0)

#Validation
#2008-11-15 to 2017-04-01
options(scipen = 999)
write.table(x = c("2008 11 15 1", TminTimeseries[which(TempDates == '2008-11-15'):which(TempDates == '2017-04-01')]), file = f_Tmin_Val, sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
write.table(x = c("2008 11 15 1", TmaxTimeseries[which(TempDates == '2008-11-15'):which(TempDates == '2017-04-01')]), file = f_Tmax_Val, sep = '\t', row.names = FALSE, col.names = FALSE, fileEncoding = 'UTF-8', quote=FALSE)
options(scipen = 0)

# Elevation of BWI and MD Sci Center for RHESSys----
#No lapse rate correction made for these stations. 
#MD Sci Center - 6.1 m
NOAAstations_locs@data[which(NOAAstations_locs@data$id == "USW00093784")[1],]
#BWI - 47.5 m
NOAAstations_locs@data[which(NOAAstations_locs@data$id == "USW00093721")[1],]
#BES Oregon Ridge: 178.22 according to Laurence. 60 m screen height, which has been changed to 6 m.
#Using these elevations in calculations.

#Baisman rainfall-runoff plots using filled-in precip timeseries----
setwd(wd_sf)
#Find streamflow date that matches the precip date
#Pond Branch
Ind17 = which(StreamStationList$`01583570`$Date == BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[1])
Ind27 = which(StreamStationList$`01583570`$Date == BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[nrow(BES_Precip_Avg_d$`Oregon Ridge Park`)])

# Precip
png('RainfallRunoff_POBR.png', res = 300, units = 'in', width = 6, height = 6)
par(xaxs="i", yaxs="i", mar=c(5,5,5,5))
plot(as.Date(RainDates), RainTimeseries, type="h", ylim=c(max(RainTimeseries)*1.5,0),
     axes=FALSE, xlab=NA, ylab=NA, col="blue",
     lwd=2, lend="square")
axis(4)
mtext("Precipitation (mm)", side=4, line=3)

# Streamflow
par(new=TRUE)
plot(as.Date(StreamStationList$`01583570`$Date[Ind17:Ind27]), StreamStationList$`01583570`$X_00060_00003[Ind17:Ind27], type="l", lwd=1, ylim=c(0, max(StreamStationList$`01583570`$X_00060_00003[Ind17:Ind27])*1.2),
     ylab = 'Streamflow (cfs)', xlab = 'Year', main = 'Pond Branch')
dev.off()

Ind1 = which(StreamStationList$`01583580`$Date == BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[1])
Ind2 = which(StreamStationList$`01583580`$Date == BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate[nrow(BES_Precip_Avg_d$`Oregon Ridge Park`)])

# Baisman Outlet Precip
png('RainfallRunoff_BARN.png', res = 300, units = 'in', width = 6, height = 6)
par(xaxs="i", yaxs="i", mar=c(5,5,5,5))
plot(as.Date(RainDates), RainTimeseries, type="h", ylim=c(max(RainTimeseries)*1.5,0),
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

#Fixme: CO2 from Mauna Loa - just for the general trend ----

#Fixme: Nitrogen deposition from CAST sources----

#Estimate WRTDS interpolation tables----
setwd(wd_WRTDS_BARN)
#Load the streamflow data into WRTDS format
Daily = readUserDaily(filePath = wd_RHESSysObs_BARN, fileName = f_BaismanStreamflowCal, hasHeader = TRUE, separator = '\t', qUnit = 1, verbose = FALSE)
#Read the TN data
Sample = readUserSample(filePath = wd_RHESSysObs_BARN, fileName = f_BaismanTNCal_WRTDS, hasHeader = TRUE, separator = '\t', verbose = FALSE)
#Set the required information
INFO = readUserInfo(filePath = wd_WRTDS_BARN_INFO, fileName = f_BaismanWRTDSInfo, interactive = FALSE)
eList = mergeReport(INFO = INFO, Daily = Daily, Sample = Sample)
saveResults(paste0(wd_WRTDS_BARN, '\\'), eList)

# Default WRTDS parameters----
WRTDSmod = modelEstimation(eList = eList, windowY = 7, windowQ = 2, windowS = .5, minNumObs = 100, minNumUncen = 50, edgeAdjust = TRUE)

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

#Contour plots of the parameters
png('TabInt_mod1.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-300,300,50),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod1.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.2,0.2,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod1.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,.7,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod1.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,.4,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod1.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,.4,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod1.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.3,0.01),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
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

#Contour plots of the parameters
png('TabInt_mod2.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-300,300,50),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod2.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.2,0.2,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod2.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,.7,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod2.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod2.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod2.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.3,0.01),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
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

#Contour plots of the parameters
png('TabInt_mod3.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-300,500,50),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod3.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.2,0.2,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod3.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,.7,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod3.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod3.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod3.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.4,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
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

#Contour plots of the parameters
png('TabInt_mod4.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-300,800,100),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod4.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.5,0.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod4.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,.7,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod4.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod4.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod4.png', units = 'in', res = 300, height = 7, width = 7)
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

#Contour plots of the parameters
png('TabInt_mod5.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-300,900,100),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod5.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.5,0.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod5.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,.7,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod5.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod5.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod5.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

# Default WRTDS parameters with Quadratic Flow Term----
WRTDSmod_QLQ = modelEstimation(eList = eList, windowY = 7, windowQ = 2, windowS = .5, minNumObs = 100, minNumUncen = 50, edgeAdjust = TRUE, QuadLogQ = TRUE)

setwd(wd_WRTDS_BARN)
png('ConcFluxTime_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod_QLQ)
plotFluxTimeDaily(WRTDSmod_QLQ)
dev.off()

png('ConcFluxTimeErrs_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod_QLQ)
plotFluxPred(WRTDSmod_QLQ)
dev.off()

plotResidPred(WRTDSmod_QLQ)
plotResidQ(WRTDSmod_QLQ)
plotResidTime(WRTDSmod_QLQ)
boxResidMonth(WRTDSmod_QLQ)
boxConcThree(WRTDSmod_QLQ)
boxQTwice(WRTDSmod_QLQ)
plotFluxHist(WRTDSmod_QLQ)
plotConcHist(WRTDSmod_QLQ)

png('BiasPlot_QLQ.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod_QLQ, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,2.6,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

#Contour plots of the parameters
png('TabInt_mod1_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-300,100,50),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod1_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.1,0.1,0.02),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod1_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2,5,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod1_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,.5,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod1_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,.5,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod1_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.3,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()
png('TabLogFlow2_mod1_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.3,1.3,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 9)
dev.off()

# Model#2-5 WRTDS with Quadratic Flow Term----
WRTDSmod2_QLQ = modelEstimation(eList = eList, windowY = 7, windowQ = 2, windowS = .25, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE, QuadLogQ = TRUE)
WRTDSmod3_QLQ = modelEstimation(eList = eList, windowY = 4, windowQ = 2, windowS = .25, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE, QuadLogQ = TRUE)
WRTDSmod4_QLQ = modelEstimation(eList = eList, windowY = 2, windowQ = 2, windowS = .25, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE, QuadLogQ = TRUE)
WRTDSmod5_QLQ = modelEstimation(eList = eList, windowY = 1.5, windowQ = 2, windowS = .25, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE, QuadLogQ = TRUE)

png('ConcFluxTime_mod2_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod2_QLQ)
plotFluxTimeDaily(WRTDSmod2_QLQ)
dev.off()

png('ConcFluxTimeErrs_mod2_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod2_QLQ)
plotFluxPred(WRTDSmod2_QLQ)
dev.off()

plotResidPred(WRTDSmod2_QLQ)
plotResidQ(WRTDSmod2_QLQ)
plotResidTime(WRTDSmod2_QLQ)
boxResidMonth(WRTDSmod2_QLQ)
boxConcThree(WRTDSmod2_QLQ)
boxQTwice(WRTDSmod2_QLQ)
plotFluxHist(WRTDSmod2_QLQ)
plotConcHist(WRTDSmod2_QLQ)

png('BiasPlot_mod2_QLQ.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod2_QLQ, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_mod2_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod2_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,5,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr_mod2_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod2_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

#Contour plots of the parameters
png('TabInt_mod2_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-300,100,50),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod2_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.1,0.1,0.02),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod2_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2,7,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod2_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod2_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod2_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.3,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()
png('TabLogFlow2_mod2_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.3,1.3,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 9)
dev.off()

png('ConcFluxTime_mod3_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod3_QLQ)
plotFluxTimeDaily(WRTDSmod3_QLQ)
dev.off()

png('ConcFluxTimeErrs_mod3_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod3_QLQ)
plotFluxPred(WRTDSmod3_QLQ)
dev.off()

plotResidPred(WRTDSmod3_QLQ)
plotResidQ(WRTDSmod3_QLQ)
plotResidTime(WRTDSmod3_QLQ)
boxResidMonth(WRTDSmod3_QLQ)
boxConcThree(WRTDSmod3_QLQ)
boxQTwice(WRTDSmod3_QLQ)
plotFluxHist(WRTDSmod3_QLQ)
plotConcHist(WRTDSmod3_QLQ)

png('BiasPlot_mod3_QLQ.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod3_QLQ, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_mod3_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod3_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,10,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr_mod3_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod3_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

#Contour plots of the parameters
png('TabInt_mod3_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-300,600,50),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod3_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.1,0.15,0.02),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod3_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2,15,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod3_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod3_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod3_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.4,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()
png('TabLogFlow2_mod3_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,3,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 9)
dev.off()

png('ConcFluxTime_mod4_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod4_QLQ)
plotFluxTimeDaily(WRTDSmod4_QLQ)
dev.off()

png('ConcFluxTimeErrs_mod4_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod4_QLQ)
plotFluxPred(WRTDSmod4_QLQ)
dev.off()

plotResidPred(WRTDSmod4_QLQ)
plotResidQ(WRTDSmod4_QLQ)
plotResidTime(WRTDSmod4_QLQ)
boxResidMonth(WRTDSmod4_QLQ)
boxConcThree(WRTDSmod4_QLQ)
boxQTwice(WRTDSmod4_QLQ)
plotFluxHist(WRTDSmod4_QLQ)
plotConcHist(WRTDSmod4_QLQ)

png('BiasPlot_mod4_QLQ.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod4_QLQ, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_mod4_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod4_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,10,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr_mod4_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod4_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

#Contour plots of the parameters
png('TabInt_mod4_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-500,1500,50),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod4_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.1,0.5,0.02),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod4_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2,6,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod4_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod4_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod4_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()
png('TabLogFlow2_mod4_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,1,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 9)
dev.off()


png('ConcFluxTime_mod5_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod5_QLQ)
plotFluxTimeDaily(WRTDSmod5_QLQ)
dev.off()

png('ConcFluxTimeErrs_mod5_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod5_QLQ)
plotFluxPred(WRTDSmod5_QLQ)
dev.off()

plotResidPred(WRTDSmod5_QLQ)
plotResidQ(WRTDSmod5_QLQ)
plotResidTime(WRTDSmod5_QLQ)
boxResidMonth(WRTDSmod5_QLQ)
boxConcThree(WRTDSmod5_QLQ)
boxQTwice(WRTDSmod5_QLQ)
plotFluxHist(WRTDSmod5_QLQ)
plotConcHist(WRTDSmod5_QLQ)

png('BiasPlot_mod5_QLQ.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod5_QLQ, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_mod5_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod5_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,10,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr_mod5_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod5_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

#Contour plots of the parameters
png('TabInt_mod5_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-500,1500,50),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod5_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.1,0.5,0.02),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod5_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2,6,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod5_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod5_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod5_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()
png('TabLogFlow2_mod5_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,1,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 9)
dev.off()

# LOOCV SE for sample and discharge----
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

sum(WRTDSmod_QLQ$Sample$SE^2)
sum(WRTDSmod2_QLQ$Sample$SE^2)
sum(WRTDSmod3_QLQ$Sample$SE^2)
sum(WRTDSmod4_QLQ$Sample$SE^2)
sum(WRTDSmod5_QLQ$Sample$SE^2)

sum(WRTDSmod_QLQ$Daily$SE^2)
sum(WRTDSmod2_QLQ$Daily$SE^2)
sum(WRTDSmod3_QLQ$Daily$SE^2)
sum(WRTDSmod4_QLQ$Daily$SE^2)
sum(WRTDSmod5_QLQ$Daily$SE^2)

# MSE for sample----
mean((WRTDSmod$Sample$ConcHat - WRTDSmod$Sample$ConcAve))^2 + sum((WRTDSmod$Sample$ConcHat - WRTDSmod$Sample$ConcAve)^2)/(nrow(Sample)-1)
mean((WRTDSmod2$Sample$ConcHat - WRTDSmod2$Sample$ConcAve))^2 + sum((WRTDSmod2$Sample$ConcHat - WRTDSmod2$Sample$ConcAve)^2)/(nrow(Sample)-1)
mean((WRTDSmod3$Sample$ConcHat - WRTDSmod3$Sample$ConcAve))^2 + sum((WRTDSmod3$Sample$ConcHat - WRTDSmod3$Sample$ConcAve)^2)/(nrow(Sample)-1)
mean((WRTDSmod4$Sample$ConcHat - WRTDSmod4$Sample$ConcAve))^2 + sum((WRTDSmod4$Sample$ConcHat - WRTDSmod4$Sample$ConcAve)^2)/(nrow(Sample)-1)
mean((WRTDSmod5$Sample$ConcHat - WRTDSmod5$Sample$ConcAve))^2 + sum((WRTDSmod5$Sample$ConcHat - WRTDSmod5$Sample$ConcAve)^2)/(nrow(Sample)-1)

mean((WRTDSmod_QLQ$Sample$ConcHat - WRTDSmod_QLQ$Sample$ConcAve))^2 + sum((WRTDSmod_QLQ$Sample$ConcHat - WRTDSmod_QLQ$Sample$ConcAve)^2)/(nrow(Sample)-1)
mean((WRTDSmod2_QLQ$Sample$ConcHat - WRTDSmod2_QLQ$Sample$ConcAve))^2 + sum((WRTDSmod2_QLQ$Sample$ConcHat - WRTDSmod2_QLQ$Sample$ConcAve)^2)/(nrow(Sample)-1)
mean((WRTDSmod3_QLQ$Sample$ConcHat - WRTDSmod3_QLQ$Sample$ConcAve))^2 + sum((WRTDSmod3_QLQ$Sample$ConcHat - WRTDSmod3_QLQ$Sample$ConcAve)^2)/(nrow(Sample)-1)
mean((WRTDSmod4_QLQ$Sample$ConcHat - WRTDSmod4_QLQ$Sample$ConcAve))^2 + sum((WRTDSmod4_QLQ$Sample$ConcHat - WRTDSmod4_QLQ$Sample$ConcAve)^2)/(nrow(Sample)-1)
mean((WRTDSmod5_QLQ$Sample$ConcHat - WRTDSmod5_QLQ$Sample$ConcAve))^2 + sum((WRTDSmod5_QLQ$Sample$ConcHat - WRTDSmod5_QLQ$Sample$ConcAve)^2)/(nrow(Sample)-1)

# Concentration and Discharge - shift after 2005?----
png('ConcDischarge.png', res = 300, units ='in', width = 5, height = 5)
par(mar=c(4,4,.5,.5))
layout(c(1,2))
plot(Daily$Date, log(Daily$Q), type = 'l', ylab = 'log(Flow [cms])', xlab = 'Year')
plot(Sample$Date, Sample$ConcAve, type = 'l', ylab = 'Total Nitrogen (mg/L)', xlab = 'Year')
dev.off()

# Running selected model with the modified functions that report the parameters of the surfaces----
WRTDSmod4m = modelEstimation(eList = eList, windowY = 2, windowQ = 2, windowS = .25, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE, numTsteps = 50, numQsteps = 100)
WRTDSmod5m_QLQ = modelEstimation(eList = eList, windowY = 1.5, windowQ = 2, windowS = .25, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE, numTsteps = 50, numQsteps = 100, QuadLogQ = TRUE)

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

#  Fill NA values in tables----
#There are 3 cells that have NA values that need to be adjusted using parameter interpolation from the tables.
IndReplace = which(is.na(TabInt))
#6 is one for each of the parameters from WRTDS
RepVal = matrix(0, nrow=length(IndReplace), ncol = 6)
for (ir in 1:length(IndReplace)){
  #Make a new temporary interpolation table that does not include the missing cell in it
  TempColInd = ceiling(IndReplace[ir]/nrow(TabInt))
  TempRowInd = IndReplace[ir] - nrow(TabInt)*(TempColInd-1)
  
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


TabInt_QLQ = WRTDSmod5m_QLQ$surfaces[,,4]
TabYear_QLQ = WRTDSmod5m_QLQ$surfaces[,,5]
TabLogQ_QLQ = WRTDSmod5m_QLQ$surfaces[,,6]
TabSinYear_QLQ = WRTDSmod5m_QLQ$surfaces[,,7]
TabCosYear_QLQ = WRTDSmod5m_QLQ$surfaces[,,8]
TabLogErr_QLQ = WRTDSmod5m_QLQ$surfaces[,,2]
TabLogQ2_QLQ = WRTDSmod5m_QLQ$surfaces[,,9]

rownames(TabInt_QLQ) = rownames(TabYear_QLQ) = rownames(TabLogQ_QLQ) = rownames(TabSinYear_QLQ) = rownames(TabCosYear_QLQ) = rownames(TabLogErr_QLQ) = rownames(TabLogQ2_QLQ) = attr(WRTDSmod5m_QLQ$surfaces, which = 'LogQ')
colnames(TabInt_QLQ) = colnames(TabYear_QLQ) = colnames(TabLogQ_QLQ) = colnames(TabSinYear_QLQ) = colnames(TabCosYear_QLQ) = colnames(TabLogErr_QLQ) = colnames(TabLogQ2_QLQ) = attr(WRTDSmod5m_QLQ$surfaces, which = 'Year')

#  Fill NA values in tables - Quadratic term----
#There are 3 cells that have NA values that need to be adjusted using parameter interpolation from the tables.
IndReplace = which(is.na(TabInt_QLQ))
#7 is one for each of the parameters from WRTDS
RepVal = matrix(0, nrow=length(IndReplace), ncol = 7)
for (ir in 1:length(IndReplace)){
  #Make a new temporary interpolation table that does not include the missing cell in it
  TempColInd = ceiling(IndReplace[ir]/nrow(TabInt_QLQ))
  TempRowInd = IndReplace[ir] - nrow(TabInt_QLQ)*(TempColInd-1)
  
  #Check if this is the first column
  if (TempColInd == 1){
    #Use the parameters from the next column without NA values for both interpolation columns
    count = 1
    while (is.na(TabInt_QLQ[TempRowInd,TempColInd+count])){
      count = count + 1
      if ((TempColInd+count) == (ncol(TabInt_QLQ)+1)){
        stop(paste('No columns in row', TempRowInd, 'of the interpolation table have non-NA values.'))
      }
    }
    Col1 = Col2 = TempColInd+count
  }else if (TempColInd == ncol(TabInt_QLQ)){
    #Used the second to last column's parameters for this column as well
    count = 1
    while (is.na(TabInt_QLQ[TempRowInd,TempColInd-count])){
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
    while (is.na(TabInt_QLQ[TempRowInd,TempColInd-count])){
      count = count - 1
      if ((TempColInd-count) == 0){
        print(paste('Columns in row', TempRowInd, 'of the interpolation table are all NA values for columns less than the target column to be filled in. Checking the greater columns.'))
        count2 = 1
        while (is.na(TabInt_QLQ[TempRowInd,TempColInd+count2])){
          count2 = count2 + 1
          if ((TempColInd+count2) == (ncol(TabInt_QLQ)+1)){
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
      while (is.na(TabInt_QLQ[TempRowInd,TempColInd+count])){
        count = count + 1
        if ((TempColInd+count) == (ncol(TabInt_QLQ)+1)){
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
    while (is.na(TabInt_QLQ[TempRowInd+count, TempColInd])){
      count = count + 1
      if ((TempRowInd+count) == (nrow(TabInt_QLQ)+1)){
        stop(paste('No rows in column', TempColInd, 'of the interpolation table have non-NA values.'))
      }
    }
    Row1 = Row2 = TempRowInd+count
  }else if (TempRowInd == nrow(TabInt_QLQ)){
    #Used the second to last row's parameters for this row as well
    count = 1
    while (is.na(TabInt_QLQ[TempRowInd-count,TempColInd])){
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
    while (is.na(TabInt_QLQ[TempRowInd-count,TempColInd])){
      count = count - 1
      if ((TempRowInd-count) == 0){
        print(paste('Rows in column', TempColInd, 'of the interpolation table are all NA values for rows less than the target row to be filled in. Checking the greater rows.'))
        count2 = 1
        while (is.na(TabInt_QLQ[TempRowInd+count2,TempColInd])){
          count2 = count2 + 1
          if ((TempRowInd+count2) == (nrow(TabInt_QLQ)+1)){
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
      while (is.na(TabInt_QLQ[TempRowInd+count,TempColInd])){
        count = count + 1
        if ((TempRowInd+count) == (nrow(TabInt_QLQ)+1)){
          print(paste('Rows in column', TempColInd, 'of the interpolation table are all NA values for rows greater than the target row to be filled in. Using the row less than the target.'))
          Row2 = Row1
          break
        }
      }
      Row2 = TempRowInd+count
    }
  }
  
  #Get new table using the corners of the tables that bound the NA values
  TempTabInt = TabInt_QLQ[unique(c(Row1:TempRowInd, TempRowInd:Row2)), unique(c(Col1:TempColInd, TempColInd:Col2))]
  TempTabYear = TabYear_QLQ[unique(c(Row1:TempRowInd, TempRowInd:Row2)), unique(c(Col1:TempColInd, TempColInd:Col2))]
  TempTabLogQ = TabLogQ_QLQ[unique(c(Row1:TempRowInd, TempRowInd:Row2)), unique(c(Col1:TempColInd, TempColInd:Col2))]
  TempTabSinYear = TabSinYear_QLQ[unique(c(Row1:TempRowInd, TempRowInd:Row2)), unique(c(Col1:TempColInd, TempColInd:Col2))]
  TempTabCosYear = TabCosYear_QLQ[unique(c(Row1:TempRowInd, TempRowInd:Row2)), unique(c(Col1:TempColInd, TempColInd:Col2))]
  TempTabLogErr = TabLogErr_QLQ[unique(c(Row1:TempRowInd, TempRowInd:Row2)), unique(c(Col1:TempColInd, TempColInd:Col2))]
  TempTabLogQ2 = TabLogQ2_QLQ[unique(c(Row1:TempRowInd, TempRowInd:Row2)), unique(c(Col1:TempColInd, TempColInd:Col2))]
  
  #Extract the new table values for the NA location using the available data in the tables that are not NA
  #Find the column and row index corresponding to the point
  
  RepVal[ir,] = FillTableNAs(DateInd = which(colnames(TempTabInt) == colnames(TabInt_QLQ)[TempColInd]), FlowInd = which(rownames(TempTabInt) == rownames(TabInt_QLQ)[TempRowInd]), QuadLogQ = TRUE)
}
#Replace the table values
TabInt_QLQ[IndReplace] = RepVal[,1]
TabYear_QLQ[IndReplace] = RepVal[,2]
TabLogQ_QLQ[IndReplace] = RepVal[,3]
TabSinYear_QLQ[IndReplace] = RepVal[,4]
TabCosYear_QLQ[IndReplace] = RepVal[,5]
TabLogErr_QLQ[IndReplace] = RepVal[,6]
TabLogQ2_QLQ[IndReplace] = RepVal[,7]

rm(ir, TempTabInt, TempTabCosYear, TempColInd, TempRowInd, TempTabLogErr, TempTabLogQ, TempTabSinYear, TempTabYear, TempTabLogQ2)
rm(Row2, Row1, Col1, Col2, count, Col2Exists, Row2Exists, RepVal, IndReplace)

#  Write the tables----
setwd(wd_WRTDS_BARN)
options(scipen = 999)
write.table(round(TabInt,5), file = f_BaismanWRTDS_TabInt, sep = '\t', row.names = attr(WRTDSmod4m$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod4m$surfaces, which = 'Year'))
write.table(signif(TabYear,4), file = f_BaismanWRTDS_TabYear, sep = '\t', row.names = attr(WRTDSmod4m$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod4m$surfaces, which = 'Year'))
write.table(signif(TabLogQ,4), file = f_BaismanWRTDS_TabLogQ, sep = '\t', row.names = attr(WRTDSmod4m$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod4m$surfaces, which = 'Year'))
write.table(signif(TabSinYear,4), file = f_BaismanWRTDS_TabSinYear, sep = '\t', row.names = attr(WRTDSmod4m$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod4m$surfaces, which = 'Year'))
write.table(signif(TabCosYear,4), file = f_BaismanWRTDS_TabCosYear, sep = '\t', row.names = attr(WRTDSmod4m$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod4m$surfaces, which = 'Year'))
write.table(round(TabLogErr,5), file = f_BaismanWRTDS_TabLogErr, sep = '\t', row.names = attr(WRTDSmod4m$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod4m$surfaces, which = 'Year'))

write.table(round(TabInt_QLQ,5), file = f_BaismanWRTDS_TabInt_QLQ, sep = '\t', row.names = attr(WRTDSmod5m_QLQ$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_QLQ$surfaces, which = 'Year'))
write.table(signif(TabYear_QLQ,4), file = f_BaismanWRTDS_TabYear_QLQ, sep = '\t', row.names = attr(WRTDSmod5m_QLQ$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_QLQ$surfaces, which = 'Year'))
write.table(signif(TabLogQ_QLQ,4), file = f_BaismanWRTDS_TabLogQ_QLQ, sep = '\t', row.names = attr(WRTDSmod5m_QLQ$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_QLQ$surfaces, which = 'Year'))
write.table(signif(TabSinYear_QLQ,4), file = f_BaismanWRTDS_TabSinYear_QLQ, sep = '\t', row.names = attr(WRTDSmod5m_QLQ$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_QLQ$surfaces, which = 'Year'))
write.table(signif(TabCosYear_QLQ,4), file = f_BaismanWRTDS_TabCosYear_QLQ, sep = '\t', row.names = attr(WRTDSmod5m_QLQ$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_QLQ$surfaces, which = 'Year'))
write.table(round(TabLogErr_QLQ,5), file = f_BaismanWRTDS_TabLogErr_QLQ, sep = '\t', row.names = attr(WRTDSmod5m_QLQ$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_QLQ$surfaces, which = 'Year'))
write.table(signif(TabLogQ2_QLQ,4), file = f_BaismanWRTDS_TabLogQ2_QLQ, sep = '\t', row.names = attr(WRTDSmod5m_QLQ$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_QLQ$surfaces, which = 'Year'))
options(scipen = 0)

#Contour plots of the parameters
png('TabInt_SelMod4.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4m, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2500,2100,100),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_SelMod4.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4m, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1.2,1.4,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_SelMod4.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4m, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-8,8,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_SelMod4.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4m, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1.6,1.2,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_SelMod4.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4m, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1.4,1.8,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_SelMod4.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4m, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

png('TabInt_SelMod5QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5m_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2500,2100,100),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_SelMod5QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5m_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1.2,1.4,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_SelMod5QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5m_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-8,8,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_SelMod5QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5m_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1.6,1.2,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_SelMod5QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5m_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1.4,1.8,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_SelMod5QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5m_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()
png('TabLogFlow2_SelMod5QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5m_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 9)
dev.off()

# Load WRTDS Interpolation Tables - rounded to X digits----
setwd(wd_WRTDS_BARN)
#Baisman
TabInt = as.matrix(read.table(file = f_BaismanWRTDS_TabInt, sep = '\t', header = TRUE, check.names = FALSE))
TabYear = as.matrix(read.table(file = f_BaismanWRTDS_TabYear, sep = '\t', header = TRUE, check.names = FALSE))
TabLogQ = as.matrix(read.table(file = f_BaismanWRTDS_TabLogQ, sep = '\t', header = TRUE, check.names = FALSE))
TabSinYear = as.matrix(read.table(file = f_BaismanWRTDS_TabSinYear, sep = '\t', header = TRUE, check.names = FALSE))
TabCosYear = as.matrix(read.table(file = f_BaismanWRTDS_TabCosYear, sep = '\t', header = TRUE, check.names = FALSE))
TabLogErr = as.matrix(read.table(file = f_BaismanWRTDS_TabLogErr, sep = '\t', header = TRUE, check.names = FALSE))

TabInt_QLQ = as.matrix(read.table(file = f_BaismanWRTDS_TabInt_QLQ, sep = '\t', header = TRUE, check.names = FALSE))
TabYear_QLQ = as.matrix(read.table(file = f_BaismanWRTDS_TabYear_QLQ, sep = '\t', header = TRUE, check.names = FALSE))
TabLogQ_QLQ = as.matrix(read.table(file = f_BaismanWRTDS_TabLogQ_QLQ, sep = '\t', header = TRUE, check.names = FALSE))
TabSinYear_QLQ = as.matrix(read.table(file = f_BaismanWRTDS_TabSinYear_QLQ, sep = '\t', header = TRUE, check.names = FALSE))
TabCosYear_QLQ = as.matrix(read.table(file = f_BaismanWRTDS_TabCosYear_QLQ, sep = '\t', header = TRUE, check.names = FALSE))
TabLogErr_QLQ = as.matrix(read.table(file = f_BaismanWRTDS_TabLogErr_QLQ, sep = '\t', header = TRUE, check.names = FALSE))
TabLogQ2_QLQ = as.matrix(read.table(file = f_BaismanWRTDS_TabLogQ2_QLQ, sep = '\t', header = TRUE, check.names = FALSE))

#WRTDS Load Allocation for hillslopes----
# Data: Kenworth 2001-2002----
setwd(dir_SynWChem_Kenworth)
K_Q = read.csv(file = f_KenworthWQ_Q, stringsAsFactors = FALSE)
#Convert to m^3/s
K_Q[,-1] = K_Q[,-1]/1000
K_TN = read.csv(file = f_KenworthWQ_TN, stringsAsFactors = FALSE)
K_Sites = readOGR(dsn = getwd(), layer = f_KenworthSites, stringsAsFactors = FALSE)
K_Sites$NAME[7] = "POBR"
K_Sites$NAME[9] = "PB4"

# Data: Smith 2006-2007----
setwd(dir_SynWChem_Smith)
S_Q = read.csv(file = f_SmithWQ_Q, stringsAsFactors = FALSE)
S_TN = read.csv(file = f_SmithWQ_TN, stringsAsFactors = FALSE)
S_Sites = readOGR(dsn = getwd(), layer = f_SmithSites, stringsAsFactors = FALSE)

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
#cms unit
Flows_BR3 = as.numeric(BR3_Q[which(!is.na(BR3_TN))][-1])

BR3_PredTN = matrix(NA, ncol = 6, nrow = length(Flows_BR3))
for(i in 1:length(Flows_BR3)){
  BR3_PredTN[i,c(1,2,3)] = predictWRTDS(Date = as.character(as.Date(Dates_BR3))[i], Flow = Flows_BR3[i]/.3048^3, 
                                rowt = as.numeric(rownames(TabInt)), colt = as.numeric(colnames(TabInt)), 
                                TabInt = TabInt, TabYear = TabYear, TabLogQ = TabLogQ, 
                                TabSinYear = TabSinYear, TabCosYear = TabCosYear, TabLogErr = TabLogErr)
  BR3_PredTN[i,c(4,5,6)] = predictWRTDS(Date = as.character(as.Date(Dates_BR3))[i], Flow = Flows_BR3[i]/.3048^3, 
                                rowt = as.numeric(rownames(TabInt_QLQ)), colt = as.numeric(colnames(TabInt_QLQ)), 
                                TabInt = TabInt_QLQ, TabYear = TabYear_QLQ, TabLogQ = TabLogQ_QLQ, 
                                TabSinYear = TabSinYear_QLQ, TabCosYear = TabCosYear_QLQ, TabLogErr = TabLogErr_QLQ, TabLogQ2 = TabLogQ2_QLQ)
}
rm(i)
colnames(BR3_PredTN) = c('05', 'Med', '95', '05QLQ', 'MedQLQ', '95QLQ')

BR3_PredTN = as.data.frame(BR3_PredTN)
BR3_PredTN$True = as.numeric(BR3_TN[which(!is.na(BR3_TN))][-1])
BR3_PredTN$DiffMed = BR3_PredTN$True - BR3_PredTN$Med
BR3_PredTN$DiffMedQLQ = BR3_PredTN$True - BR3_PredTN$MedQLQ

setwd(wd_BESN)
png('BR3_TrueTNvsWRTDS.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR3_PredTN$True, BR3_PredTN$Med, ylim = c(0,9), xlim = c(0,9), xlab = 'True BR3 TN (mg N/L)', ylab = 'WRTDS Predicted BR3 TN (mg N/L)')
lines(c(-10,10), c(-10,10))
dev.off()

png('BR3_TrueTNvsWRTDS_QLQ.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR3_PredTN$True, BR3_PredTN$MedQLQ, ylim = c(0,9), xlim = c(0,9), xlab = 'True BR3 TN (mg N/L)', ylab = 'WRTDS Predicted BR3 TN (mg N/L)')
lines(c(-10,10), c(-10,10))
dev.off()

#plot(BR3_PredTN$True, BR3_PredTN$`05`)
#plot(BR3_PredTN$True, BR3_PredTN$`95`)

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
#plot(lm_BR3)

summary(lm_BR3_Dates)
#plot(lm_BR3_Dates)

#plot of regression y vs. x
plot(y = BR3_PredTN$DiffMed, x = BR3_PredTN$Med)
plot(y = BR3_PredTN$DiffMed, x = Flows_BR3)
plot(y = BR3_PredTN$DiffMed, x = sin(2*pi*DatesNums_BR3/366))
plot(y = BR3_PredTN$DiffMed, x = cos(2*pi*DatesNums_BR3/366))

#   Relation dropping suspeccted outlier----
lm_BR3 = lm(BR3_PredTN$DiffMed[BR3_PredTN$True > 1] ~ BR3_PredTN$Med[BR3_PredTN$True > 1])
lm_BR3_Flows = lm(BR3_PredTN$DiffMed[BR3_PredTN$True > 1] ~ BR3_PredTN$Med[BR3_PredTN$True > 1]*I(Flows_BR3[BR3_PredTN$True > 1]*1000))
lm_BR3_Dates = lm(BR3_PredTN$DiffMed[BR3_PredTN$True > 1] ~ BR3_PredTN$Med[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))

summary(lm_BR3)
#plot(lm_BR3)

summary(lm_BR3_Flows)
#plot(lm_BR3_Flows)

summary(lm_BR3_Dates)
#plot(lm_BR3_Dates)

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
#cms unit
Flows_BR5 = as.numeric(BR5_Q[which(!is.na(BR5_TN))][-1])
#Drop the NA flow date
Dates_BR5 = Dates_BR5[-length(Dates_BR5)]
Flows_BR5 = Flows_BR5[-length(Flows_BR5)]

BR5_PredTN = matrix(NA, ncol = 6, nrow = length(Flows_BR5))
for(i in 1:length(Flows_BR5)){
  BR5_PredTN[i,c(1,2,3)] = predictWRTDS(Date = as.character(as.Date(Dates_BR5))[i], Flow = Flows_BR5[i]/.3048^3, 
                                rowt = as.numeric(rownames(TabInt)), colt = as.numeric(colnames(TabInt)), 
                                TabInt = TabInt, TabYear = TabYear, TabLogQ = TabLogQ, 
                                TabSinYear = TabSinYear, TabCosYear = TabCosYear, TabLogErr = TabLogErr)
  BR5_PredTN[i,c(4,5,6)] = predictWRTDS(Date = as.character(as.Date(Dates_BR5))[i], Flow = Flows_BR5[i]/.3048^3, 
                                        rowt = as.numeric(rownames(TabInt_QLQ)), colt = as.numeric(colnames(TabInt_QLQ)), 
                                        TabInt = TabInt_QLQ, TabYear = TabYear_QLQ, TabLogQ = TabLogQ_QLQ, 
                                        TabSinYear = TabSinYear_QLQ, TabCosYear = TabCosYear_QLQ, TabLogErr = TabLogErr_QLQ, TabLogQ2 = TabLogQ2_QLQ)
}
rm(i)
colnames(BR5_PredTN) = c('05', 'Med', '95', '05QLQ', 'MedQLQ', '95QLQ')

BR5_PredTN = as.data.frame(BR5_PredTN)
BR5_PredTN$True = as.numeric(BR5_TN[which(!is.na(BR5_TN))][-1])[-ncol(BR5_TN[which(!is.na(BR5_TN))][-1])]
BR5_PredTN$DiffMed = BR5_PredTN$True - BR5_PredTN$Med
BR5_PredTN$DiffMedQLQ = BR5_PredTN$True - BR5_PredTN$MedQLQ

png('BR5_TrueTNvsWRTDS.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR5_PredTN$True, BR5_PredTN$Med, ylim = c(0,9), xlim = c(0,9), xlab = 'True BR5 TN (mg N/L)', ylab = 'WRTDS Predicted BR5 TN (mg N/L)')
lines(c(-10,10), c(-10,10))
dev.off()

png('BR5_TrueTNvsWRTDS_QLQ.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR5_PredTN$True, BR5_PredTN$MedQLQ, ylim = c(0,9), xlim = c(0,9), xlab = 'True BR5 TN (mg N/L)', ylab = 'WRTDS Predicted BR5 TN (mg N/L)')
lines(c(-10,10), c(-10,10))
dev.off()

#plot(BR5_PredTN$True, BR5_PredTN$`05`)
#plot(BR5_PredTN$True, BR5_PredTN$`95`)

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
#plot(lm_BR5)

summary(lm_BR5_Dates)
#plot(lm_BR5_Dates)

#plot of regression y vs. x
plot(y = BR5_PredTN$DiffMed, x = BR5_PredTN$Med)
plot(y = BR5_PredTN$DiffMed, x = sin(2*pi*DatesNums_BR5/366))
plot(y = BR5_PredTN$DiffMed, x = cos(2*pi*DatesNums_BR5/366))

#   Relation dropping suspeccted outlier----
lm_BR5 = lm(BR5_PredTN$DiffMed[BR5_PredTN$True > 1] ~ BR5_PredTN$Med[BR5_PredTN$True > 1])
lm_BR5_Flows = lm(BR5_PredTN$DiffMed[BR5_PredTN$True > 1] ~ BR5_PredTN$Med[BR5_PredTN$True > 1] + I(Flows_BR5[BR5_PredTN$True > 1]))
lm_BR5_Dates = lm(BR5_PredTN$DiffMed[BR5_PredTN$True > 1] ~ BR5_PredTN$Med[BR5_PredTN$True > 1] + sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + I(Flows_BR5[BR5_PredTN$True > 1]))

summary(lm_BR5)
#plot(lm_BR5)

summary(lm_BR5_Flows)
#plot(lm_BR5_Flows)

summary(lm_BR5_Dates)
#plot(lm_BR5_Dates)

#plot of regression y vs. x
plot(y = BR5_PredTN$DiffMed[BR5_PredTN$True > 1], x = BR5_PredTN$Med[BR5_PredTN$True > 1], xlab = 'Predicted Mean TN (mg N/L)', ylab = 'Difference from True TN (mg N/L)')
plot(y = BR5_PredTN$DiffMed[BR5_PredTN$True > 1], x = Flows_BR5[BR5_PredTN$True > 1])
plot(y = BR5_PredTN$DiffMed[BR5_PredTN$True > 1], x = sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))
plot(y = BR5_PredTN$DiffMed[BR5_PredTN$True > 1], x = cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))

#  Pond Branch----
setwd(wd_BESN)

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
#cfs unit
Flows_POBR = BES_TN_d$POBR$Flow[!is.na(BES_TN_d$POBR$TN..mg.N.L.)][which((as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]) %in% as.Date(BES_TN_d$BARN$SortDate[!is.na(BES_TN_d$BARN$TN..mg.N.L.)])) & (as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]) <= '2014-01-01'))]
Dates_POBR = BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)][which((as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]) %in% as.Date(BES_TN_d$BARN$SortDate[!is.na(BES_TN_d$BARN$TN..mg.N.L.)])) & (as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]) <= '2014-01-01'))]

POBR_PredTN = matrix(NA, ncol = 6, nrow = length(Flows_POBR))
for(i in 1:length(Flows_POBR)){
  POBR_PredTN[i,c(1,2,3)] = predictWRTDS(Date = as.character(as.Date(Dates_POBR))[i], Flow = Flows_POBR[i], 
                                 rowt = as.numeric(rownames(TabInt)), colt = as.numeric(colnames(TabInt)), 
                                 TabInt = TabInt, TabYear = TabYear, TabLogQ = TabLogQ, 
                                 TabSinYear = TabSinYear, TabCosYear = TabCosYear, TabLogErr = TabLogErr)
  POBR_PredTN[i,c(4,5,6)] = predictWRTDS(Date = as.character(as.Date(Dates_POBR))[i], Flow = Flows_POBR[i], 
                                         rowt = as.numeric(rownames(TabInt_QLQ)), colt = as.numeric(colnames(TabInt_QLQ)), 
                                         TabInt = TabInt_QLQ, TabYear = TabYear_QLQ, TabLogQ = TabLogQ_QLQ, 
                                         TabSinYear = TabSinYear_QLQ, TabCosYear = TabCosYear_QLQ, TabLogErr = TabLogErr_QLQ, TabLogQ2 = TabLogQ2_QLQ)
}
rm(i)
colnames(POBR_PredTN) = c('05', 'Med', '95', '05QLQ', 'MedQLQ', '95QLQ')

POBR_PredTN = as.data.frame(POBR_PredTN)
POBR_PredTN$True = BES_TN_d$POBR$TN..mg.N.L.[!is.na(BES_TN_d$POBR$TN..mg.N.L.)][which((as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]) %in% as.Date(BES_TN_d$BARN$SortDate[!is.na(BES_TN_d$BARN$TN..mg.N.L.)])) & (as.Date(BES_TN_d$POBR$SortDate[!is.na(BES_TN_d$POBR$TN..mg.N.L.)]) <= '2014-01-01'))]
POBR_PredTN$DiffMed = POBR_PredTN$True - POBR_PredTN$Med
POBR_PredTN$DiffMedQLQ = POBR_PredTN$True - POBR_PredTN$MedQLQ

plot(POBR_PredTN$True, POBR_PredTN$Med)
lines(c(-10,10), c(-10,10))
#plot(POBR_PredTN$True, POBR_PredTN$`05`)
#plot(POBR_PredTN$True, POBR_PredTN$`95`)

#Assign number values to the dates to be used in regression. 1 - 366
DateNumMat = matrix(NA, nrow = 366, ncol = 2)
DateNumMat[,1] = seq(as.Date("2000-01-01"), as.Date("2000-12-31"), 1)
DateNumMat[,2] = seq(1, 366, 1)

DatesNums_POBR = vector('numeric', length=length(Dates_POBR))
for(i in 1:length(Dates_POBR)){
  DatesNums_POBR[i] = DateNumMat[which((month(as.Date(DateNumMat[,1])) == month(as.Date(Dates_POBR[i]))) & (day(as.Date(DateNumMat[,1])) == day(as.Date(Dates_POBR[i])))), 2]
}

lm_POBR = lm(POBR_PredTN$DiffMed ~ I(Flows_POBR*.3048^3) + POBR_PredTN$Med + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))

summary(lm_POBR)
#plot(lm_POBR)

#   Relation using all available data for Pond Branch----
#cfs unit
Flows_POBR = BES_TN_d$POBR$Flow[which((!is.na(BES_TN_d$POBR$TN..mg.N.L.)) & (as.Date(BES_TN_d$POBR$SortDate) <= '2014-01-01'))]
Dates_POBR = as.Date(BES_TN_d$POBR$SortDate[which((!is.na(BES_TN_d$POBR$TN..mg.N.L.)) & (as.Date(BES_TN_d$POBR$SortDate) <= '2014-01-01'))])

POBR_PredTN = matrix(NA, ncol = 6, nrow = length(Flows_POBR))
for(i in 1:length(Flows_POBR)){
  POBR_PredTN[i,c(1,2,3)] = predictWRTDS(Date = as.character(as.Date(Dates_POBR))[i], Flow = Flows_POBR[i], 
                                         rowt = as.numeric(rownames(TabInt)), colt = as.numeric(colnames(TabInt)), 
                                         TabInt = TabInt, TabYear = TabYear, TabCosYear = TabCosYear, 
                                         TabSinYear = TabSinYear, TabLogQ = TabLogQ, TabLogErr = TabLogErr)
  POBR_PredTN[i,c(4,5,6)] = predictWRTDS(Date = as.character(as.Date(Dates_POBR))[i], Flow = Flows_POBR[i], 
                                         rowt = as.numeric(rownames(TabInt_QLQ)), colt = as.numeric(colnames(TabInt_QLQ)), 
                                         TabInt = TabInt_QLQ, TabYear = TabYear_QLQ, TabCosYear = TabCosYear_QLQ, 
                                         TabSinYear = TabSinYear_QLQ, TabLogQ = TabLogQ_QLQ, TabLogErr = TabLogErr_QLQ, TabLogQ2 = TabLogQ2_QLQ)
}
rm(i)
colnames(POBR_PredTN) = c('05', 'Med', '95', '05QLQ', 'MedQLQ', '95QLQ')

POBR_PredTN = as.data.frame(POBR_PredTN)
POBR_PredTN$True = BES_TN_d$POBR$TN..mg.N.L.[which((!is.na(BES_TN_d$POBR$TN..mg.N.L.)) & (as.Date(BES_TN_d$POBR$SortDate) <= '2014-01-01'))]
POBR_PredTN$DiffMed = POBR_PredTN$True - POBR_PredTN$Med
POBR_PredTN$DiffMedQLQ = POBR_PredTN$True - POBR_PredTN$MedQLQ

png('POBR_TrueTNvsWRTDS.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN$True, POBR_PredTN$Med, ylim = c(0,9), xlim = c(0,9), xlab = 'True Pond Branch TN (mg N/L)', ylab = 'WRTDS Predicted Pond Branch TN (mg N/L)')
lines(c(-10,10), c(-10,10))
dev.off()

png('POBR_TrueTNvsWRTDS_QLQ.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN$True, POBR_PredTN$MedQLQ, ylim = c(0,9), xlim = c(0,9), xlab = 'True Pond Branch TN (mg N/L)', ylab = 'WRTDS Predicted Pond Branch TN (mg N/L)')
lines(c(-10,10), c(-10,10))
dev.off()

#plot(POBR_PredTN$True, POBR_PredTN$`05`)
#plot(POBR_PredTN$True, POBR_PredTN$`95`)

#Assign number values to the dates to be used in regression. 1 - 366
DateNumMat = matrix(NA, nrow = 366, ncol = 2)
DateNumMat[,1] = seq(as.Date("2000-01-01"), as.Date("2000-12-31"), 1)
DateNumMat[,2] = seq(1, 366, 1)

DatesNums_POBR = vector('numeric', length=length(Dates_POBR))
for(i in 1:length(Dates_POBR)){
  DatesNums_POBR[i] = DateNumMat[which((month(as.Date(DateNumMat[,1])) == month(as.Date(Dates_POBR[i]))) & (day(as.Date(DateNumMat[,1])) == day(as.Date(Dates_POBR[i])))), 2]
}

#Fixme - this needs to be a censored regression because of the y censoring. Thresholds change over time. Great.
lm_POBR = lm(POBR_PredTN$DiffMed ~ I(Flows_POBR*.3048^3) + POBR_PredTN$Med)
lm_POBR_Dates = lm(POBR_PredTN$DiffMed ~ I(Flows_POBR*.3048^3) + POBR_PredTN$Med + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))
lm_POBR_log = lm(log10(POBR_PredTN$DiffMed + abs(min(POBR_PredTN$DiffMed)) + .001) ~ POBR_PredTN$Med + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))

summary(lm_POBR)
#plot(lm_POBR)

summary(lm_POBR_Dates)
#plot(lm_POBR_Dates)

summary(lm_POBR_log)

#plot of regression y vs. x
plot(y = POBR_PredTN$DiffMed, x = POBR_PredTN$Med)
plot(y = POBR_PredTN$DiffMed, x = sin(2*pi*DatesNums_POBR/366))
plot(y = POBR_PredTN$DiffMed, x = cos(2*pi*DatesNums_POBR/366))

# Make regressions based on the TN load instead of the TN concentration----
#  Load at Basin Outlet----
#This load must result from the loadings at all other sites, minus the uptake by plants and organisms along the river. Ignoring uptake for now.
BES_TN_d$BARN$Load = BES_TN_d$BARN$TN..mg.N.L.*1000*BES_TN_d$BARN$Flow*.3048^3

#  Load at Pond Branch----
BES_TN_d$POBR$Load = BES_TN_d$POBR$TN..mg.N.L.*1000*BES_TN_d$POBR$Flow*.3048^3

#Make indicators for likely detection limits
BES_TN_d$POBR$LowTNLim[BES_TN_d$POBR$TN..mg.N.L. == min(BES_TN_d$POBR$TN..mg.N.L.[1:5000], na.rm = TRUE)] = 1
BES_TN_d$POBR$LowTNLim[BES_TN_d$POBR$TN..mg.N.L. == min(BES_TN_d$POBR$TN..mg.N.L.[5001:length(BES_TN_d$POBR$TN..mg.N.L.)], na.rm = TRUE)] = 2
POBR_PredTN$LowTNLim = NA
POBR_PredTN$LowTNLim[POBR_PredTN$True == min(POBR_PredTN$True[1:650])] = 1
POBR_PredTN$LowTNLim[POBR_PredTN$True == min(POBR_PredTN$True[651:length(POBR_PredTN$True)])] = 2

#Use WRTDS regression for predicted concentration and transform to loads. Compare predicted to true load.
#mg/s units
POBR_PredTN$TrueLoad = POBR_PredTN$True*1000*Flows_POBR*.3048^3
POBR_PredTN$Load05 = POBR_PredTN$`05`*1000*Flows_POBR*.3048^3
POBR_PredTN$MedLoad = POBR_PredTN$Med*1000*Flows_POBR*.3048^3
POBR_PredTN$Load95 = POBR_PredTN$`95`*1000*Flows_POBR*.3048^3
POBR_PredTN$Load05QLQ = POBR_PredTN$`05QLQ`*1000*Flows_POBR*.3048^3
POBR_PredTN$MedLoadQLQ = POBR_PredTN$MedQLQ*1000*Flows_POBR*.3048^3
POBR_PredTN$Load95QLQ = POBR_PredTN$`95QLQ`*1000*Flows_POBR*.3048^3

png('POBR_TrueTNLoadvsWRTDSLoad.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN$TrueLoad, POBR_PredTN$MedLoad, ylim=range(c(POBR_PredTN$MedLoad-POBR_PredTN$Load05, POBR_PredTN$MedLoad+POBR_PredTN$Load95)), xlim = c(0, 50), xlab = 'True Pond Branch TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of Pond Branch TN (mg N/s)')
arrows(POBR_PredTN$TrueLoad, POBR_PredTN$Load05, POBR_PredTN$TrueLoad, POBR_PredTN$Load95, length=0.05, angle=90, code=3)
lines(c(-10,50), c(-10,50))
dev.off()

png('POBR_TrueTNLoadvsWRTDSLoad_log.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN$TrueLoad, POBR_PredTN$MedLoad, ylim=range(c(POBR_PredTN$MedLoad-POBR_PredTN$Load05, POBR_PredTN$MedLoad+POBR_PredTN$Load95)), xlab = 'True Pond Branch TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of Pond Branch TN (mg N/s)', log = 'xy')
arrows(POBR_PredTN$TrueLoad, POBR_PredTN$Load05, POBR_PredTN$TrueLoad, POBR_PredTN$Load95, length=0.05, angle=90, code=3)
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
#arrows(POBR_PredTN$TrueLoad, POBR_PredTN$Load05, POBR_PredTN$TrueLoad, POBR_PredTN$Load95, length=0.05, angle=90, code=3)
lines(c(1e-7,1e7), c(1e-7,1e7))
legend('bottomright', title = 'Detection Limits (mg N/L)', legend = c('0.01', '0.05'), pch = 1, col = c('blue', 'red'))
dev.off()

png('POBR_TrueTNLoadvsWRTDSLoad_QLQ.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN$TrueLoad, POBR_PredTN$MedLoadQLQ, ylim=range(c(POBR_PredTN$MedLoadQLQ-POBR_PredTN$Load05QLQ, POBR_PredTN$MedLoadQLQ+POBR_PredTN$Load95QLQ)), xlim = c(0, 50), xlab = 'True Pond Branch TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of Pond Branch TN (mg N/s)')
arrows(POBR_PredTN$TrueLoad, POBR_PredTN$Load05QLQ, POBR_PredTN$TrueLoad, POBR_PredTN$Load95QLQ, length=0.05, angle=90, code=3)
lines(c(-10,50), c(-10,50))
dev.off()

png('POBR_TrueTNLoadvsWRTDSLoad_QLQ_log.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN$TrueLoad, POBR_PredTN$MedLoadQLQ, ylim=range(c(POBR_PredTN$MedLoadQLQ-POBR_PredTN$Load05QLQ, POBR_PredTN$MedLoadQLQ+POBR_PredTN$Load95QLQ)), xlab = 'True Pond Branch TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of Pond Branch TN (mg N/s)', log = 'xy')
arrows(POBR_PredTN$TrueLoad, POBR_PredTN$Load05QLQ, POBR_PredTN$TrueLoad, POBR_PredTN$Load95QLQ, length=0.05, angle=90, code=3)
lines(c(1e-7,1e7), c(1e-7,1e7))
dev.off()

png('POBR_TrueTNLoadvsWRTDSLoad_QLQ_log_colors.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN$TrueLoad, POBR_PredTN$MedLoadQLQ, ylim=range(c(POBR_PredTN$MedLoadQLQ-POBR_PredTN$Load05QLQ, POBR_PredTN$MedLoadQLQ+POBR_PredTN$Load95QLQ)),
     xlim=range(POBR_PredTN$TrueLoad),
     xlab = 'True Pond Branch TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of Pond Branch TN (mg N/s)', log = 'xy')
par(new = TRUE)
plot(POBR_PredTN$TrueLoad[POBR_PredTN$LowTNLim == 1], POBR_PredTN$MedLoadQLQ[POBR_PredTN$LowTNLim == 1], 
     xlim=range(POBR_PredTN$TrueLoad),
     ylim=range(c(POBR_PredTN$MedLoadQLQ-POBR_PredTN$Load05QLQ, POBR_PredTN$MedLoadQLQ+POBR_PredTN$Load95QLQ)), xlab = '', ylab = '', 
     log = 'xy', col = 'blue', axes = FALSE)
par(new = TRUE)
plot(POBR_PredTN$TrueLoad[POBR_PredTN$LowTNLim == 2], POBR_PredTN$MedLoadQLQ[POBR_PredTN$LowTNLim == 2], 
     xlim=range(POBR_PredTN$TrueLoad),
     ylim=range(c(POBR_PredTN$MedLoadQLQ-POBR_PredTN$Load05QLQ, POBR_PredTN$MedLoadQLQ+POBR_PredTN$Load95QLQ)), xlab = '', ylab = '', 
     log = 'xy', col = 'red', axes = FALSE)
lines(c(1e-7,1e7), c(1e-7,1e7))
legend('bottomright', title = 'Detection Limits (mg N/L)', legend = c('0.01', '0.05'), pch = 1, col = c('blue', 'red'))
dev.off()

#Add colors for seasons to this dataset
POBR_PredTN$cols = BES_TN_d$POBR$cols[as.Date(BES_TN_d$POBR$SortDate) %in% as.Date(Dates_POBR)]

png('POBR_TN_ScatterHist.png', res = 300, units = 'in', width = 10, height = 10)
scatterHist_mod(x = log10(POBR_PredTN$TrueLoad), y = log10(POBR_PredTN$MedLoad), density = FALSE, freq = FALSE, x.breaks = 100, y.breaks = 100, correl = FALSE, 
                xlab = expression(paste('log'[10] ~ '(True TN Load [mg N/s])')), ylab = expression(paste('log'[10] ~ 'Predicted TN Load (mg N/s)')), 
                cex.lab = 1.5, cex.axis = 1.5, title = '', pch = 16, ellipse = FALSE, smooth = FALSE, col = POBR_PredTN$cols)
dev.off()

png('POBR_TN_ScatterHist_QLQ.png', res = 300, units = 'in', width = 10, height = 10)
scatterHist_mod(x = log10(POBR_PredTN$TrueLoad), y = log10(POBR_PredTN$MedLoadQLQ), density = FALSE, freq = FALSE, x.breaks = 100, y.breaks = 100, correl = FALSE, 
                xlab = expression(paste('log'[10] ~ '(True TN Load [mg N/s])')), ylab = expression(paste('log'[10] ~ 'Predicted TN Load (mg N/s)')), 
                cex.lab = 1.5, cex.axis = 1.5, title = '', pch = 16, ellipse = FALSE, smooth = FALSE, col = POBR_PredTN$cols)
dev.off()

POBR_PredTN$DiffMedLoad = POBR_PredTN$TrueLoad - POBR_PredTN$MedLoad
POBR_PredTN$DiffMedLoadQLQ = POBR_PredTN$TrueLoad - POBR_PredTN$MedLoadQLQ

#plot of regression y vs. x
plot(y = POBR_PredTN$DiffMedLoad, x = POBR_PredTN$MedLoad)
plot(y = POBR_PredTN$DiffMedLoad, x = Flows_POBR*.3048^3)
plot(y = POBR_PredTN$DiffMedLoad, x = sin(2*pi*DatesNums_POBR/366))
plot(y = POBR_PredTN$DiffMedLoad, x = cos(2*pi*DatesNums_POBR/366))

#Fixme - this needs to be a censored regression because of the y censoring. Thresholds change over time. Great.
lm_POBRLoad = lm(POBR_PredTN$DiffMedLoad ~ I(Flows_POBR*.3048^3) + POBR_PredTN$MedLoad)
lm_POBRLoad_Dates = lm(POBR_PredTN$DiffMedLoad ~ I(Flows_POBR*.3048^3) + POBR_PredTN$MedLoad + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))
lm_POBRLoad_Dates_Interaction = lm(POBR_PredTN$DiffMedLoad ~ I(Flows_POBR*.3048^3)*POBR_PredTN$Med + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))
lm_POBRLoad_log = lm(log10(POBR_PredTN$DiffMedLoad + abs(min(POBR_PredTN$DiffMedLoad)) + .001) ~ POBR_PredTN$MedLoad + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))

summary(lm_POBRLoad)
#plot(lm_POBRLoad)

summary(lm_POBRLoad_Dates)
#plot(lm_POBRLoad_Dates)

summary(lm_POBRLoad_Dates_Interaction)
#plot(lm_POBRLoad_Dates_Interaction)

summary(lm_POBRLoad_log)
#plot(lm_POBRLoad_log)

#From the POBRLoad_Dates, use Box-Cox transform to attempt normality
BoxCoxtest = boxcox(I(POBR_PredTN$DiffMedLoad + 40) ~ I(Flows_POBR*.3048^3)*POBR_PredTN$Med + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366), lambda = seq(-2,2, length = 1000))
lambda = BoxCoxtest$x[which(BoxCoxtest$y == max(BoxCoxtest$y))]
POBR_PredTN$DiffMedLoad_BC = bcPower(U = POBR_PredTN$DiffMedLoad+40, gamma = 0, lambda = lambda, jacobian.adjusted = FALSE)
lm_POBRLoad_Dates_BC = lm((((POBR_PredTN$DiffMedLoad + 40)^lambda - 1)/lambda) ~ I(Flows_POBR*.3048^3)*POBR_PredTN$Med + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))

summary(lm_POBRLoad_Dates_BC)
#plot(lm_POBRLoad_Dates_BC)

#plot of regression y vs. x
plot(y = POBR_PredTN$DiffMedLoad_BC, x = POBR_PredTN$MedLoad)
plot(y = POBR_PredTN$DiffMedLoad_BC, x = I(Flows_POBR*.3048^3))
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
#mg/s unit
BR3_PredTN$TrueLoad = BR3_PredTN$True*1000*Flows_BR3
BR3_PredTN$Load05 = BR3_PredTN$`05`*1000*Flows_BR3
BR3_PredTN$MedLoad = BR3_PredTN$Med*1000*Flows_BR3
BR3_PredTN$Load95 = BR3_PredTN$`95`*1000*Flows_BR3
BR3_PredTN$Load05QLQ = BR3_PredTN$`05QLQ`*1000*Flows_BR3
BR3_PredTN$MedLoadQLQ = BR3_PredTN$MedQLQ*1000*Flows_BR3
BR3_PredTN$Load95QLQ = BR3_PredTN$`95QLQ`*1000*Flows_BR3

png('BR3_TrueTNLoadvsWRTDSLoad.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$MedLoad[BR3_PredTN$True > 1], 
     ylim=range(c(BR3_PredTN$MedLoad[BR3_PredTN$True > 1]-BR3_PredTN$Load05[BR3_PredTN$True > 1], BR3_PredTN$MedLoad[BR3_PredTN$True > 1]+BR3_PredTN$Load95[BR3_PredTN$True > 1])), 
     xlab = 'True BR3 TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of BR3 TN (mg N/s)')
arrows(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$Load05[BR3_PredTN$True > 1], 
       BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$Load95[BR3_PredTN$True > 1], 
       length=0.05, angle=90, code=3)
lines(c(-100,100), c(-100,100))
dev.off()

png('BR3_TrueTNLoadvsWRTDSLoad_log.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$MedLoad[BR3_PredTN$True > 1], 
     ylim=range(c(BR3_PredTN$MedLoad[BR3_PredTN$True > 1]-BR3_PredTN$Load05[BR3_PredTN$True > 1], BR3_PredTN$MedLoad[BR3_PredTN$True > 1]+BR3_PredTN$Load95[BR3_PredTN$True > 1])), 
     xlab = expression(paste('log'[10] ~ 'True BR3 TN Load (mg N/s)')), ylab = expression(paste('log'[10] ~ 'WRTDS Predicted Load of BR3 TN (mg N/s)')), log='xy')
arrows(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$Load05[BR3_PredTN$True > 1], 
       BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$Load95[BR3_PredTN$True > 1], 
       length=0.05, angle=90, code=3)
lines(c(1e-7,100), c(1e-7,100))
dev.off()

png('BR3_TrueTNLoadvsWRTDSLoad_QLQ.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$MedLoadQLQ[BR3_PredTN$True > 1], 
     ylim=range(c(BR3_PredTN$MedLoadQLQ[BR3_PredTN$True > 1]-BR3_PredTN$Load05QLQ[BR3_PredTN$True > 1], BR3_PredTN$MedLoadQLQ[BR3_PredTN$True > 1]+BR3_PredTN$Load95QLQ[BR3_PredTN$True > 1])), 
     xlab = 'True BR3 TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of BR3 TN (mg N/s)')
arrows(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$Load05QLQ[BR3_PredTN$True > 1], 
       BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$Load95QLQ[BR3_PredTN$True > 1], 
       length=0.05, angle=90, code=3)
lines(c(-100,100), c(-100,100))
dev.off()

png('BR3_TrueTNLoadvsWRTDSLoad_QLQ_log.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$MedLoadQLQ[BR3_PredTN$True > 1], 
     ylim=range(c(BR3_PredTN$MedLoadQLQ[BR3_PredTN$True > 1]-BR3_PredTN$Load05QLQ[BR3_PredTN$True > 1], BR3_PredTN$MedLoadQLQ[BR3_PredTN$True > 1]+BR3_PredTN$Load95QLQ[BR3_PredTN$True > 1])), 
     xlab = expression(paste('log'[10] ~ 'True BR3 TN Load (mg N/s)')), ylab = expression(paste('log'[10] ~ 'WRTDS Predicted Load of BR3 TN (mg N/s)')), log='xy')
arrows(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$Load05QLQ[BR3_PredTN$True > 1], 
       BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], BR3_PredTN$Load95QLQ[BR3_PredTN$True > 1], 
       length=0.05, angle=90, code=3)
lines(c(1e-7,100), c(1e-7,100))
dev.off()


BR3_PredTN$DiffMedLoad = BR3_PredTN$TrueLoad - BR3_PredTN$MedLoad
BR3_PredTN$DiffMedLoadQLQ = BR3_PredTN$TrueLoad - BR3_PredTN$MedLoadQLQ

#plot of regression y vs. x
plot(y = BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1], x = BR3_PredTN$MedLoad[BR3_PredTN$True > 1])
plot(y = BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1], x = I(Flows_BR3[BR3_PredTN$True > 1]))
plot(y = BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1], x = sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))
plot(y = BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1], x = cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))

#Regression w/o suspected outlier
lm_BR3Load = lm(BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1])
lm_BR3Load_Dates = lm(BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))
lm_BR3Load_Dates_Interaction = lm(BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1])*BR3_PredTN$Med[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))
lm_BR3Load_Interaction = lm(BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1])*BR3_PredTN$Med[BR3_PredTN$True > 1])
lm_BR3Load_log = lm(log10(BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1] + abs(min(BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1])) + .001) ~ BR3_PredTN$MedLoad[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))

summary(lm_BR3Load)
#plot(lm_BR3Load)

summary(lm_BR3Load_Dates)
#plot(lm_BR3Load_Dates)

summary(lm_BR3Load_Dates_Interaction)
#plot(lm_BR3Load_Dates_Interaction)

summary(lm_BR3Load_Interaction)
#plot(lm_BR3Load_Interaction)

summary(lm_BR3Load_log)
#plot(lm_BR3Load_log)


#From the BR3Load_Dates, use Box-Cox transform to attempt normality
BoxCoxtest_BR3 = boxcox(BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366), lambda = seq(0,1, length = 1000))
lambda_BR3 = BoxCoxtest_BR3$x[which(BoxCoxtest_BR3$y == max(BoxCoxtest_BR3$y))]

BR3_PredTN$DiffMedLoad_BC[BR3_PredTN$True > 1] = bcPower(U = BR3_PredTN$DiffMedLoad[BR3_PredTN$True > 1], gamma = 0, lambda = lambda_BR3, jacobian.adjusted = FALSE)
lm_BR3Load_Dates_BC = lm(BR3_PredTN$DiffMedLoad_BC[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))

summary(lm_BR3Load_Dates_BC)
#plot(lm_BR3Load_Dates_BC)

#plot of regression y vs. x
plot(y = BR3_PredTN$DiffMedLoad_BC, x = BR3_PredTN$MedLoad)
plot(y = BR3_PredTN$DiffMedLoad_BC, x = I(Flows_BR3))
plot(y = BR3_PredTN$DiffMedLoad_BC, x = sin(2*pi*DatesNums_BR3/366))
plot(y = BR3_PredTN$DiffMedLoad_BC, x = cos(2*pi*DatesNums_BR3/366))

#  Load at BR5----
#mg/s
BR5_Load = BR5_Q[-1]*BR5_TN[-1]*1000

#Use WRTDS regression for predicted concentration and transform to loads. Compare predicted to true load.
#mg/s unit
BR5_PredTN$TrueLoad = BR5_PredTN$True*1000*Flows_BR5
BR5_PredTN$Load05 = BR5_PredTN$`05`*1000*Flows_BR5
BR5_PredTN$MedLoad = BR5_PredTN$Med*1000*Flows_BR5
BR5_PredTN$Load95 = BR5_PredTN$`95`*1000*Flows_BR5
BR5_PredTN$Load05QLQ = BR5_PredTN$`05QLQ`*1000*Flows_BR5
BR5_PredTN$MedLoadQLQ = BR5_PredTN$MedQLQ*1000*Flows_BR5
BR5_PredTN$Load95QLQ = BR5_PredTN$`95QLQ`*1000*Flows_BR5

png('BR5_TrueTNLoadvsWRTDSLoad.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR5_PredTN$TrueLoad, BR5_PredTN$MedLoad, ylim=range(c(BR5_PredTN$MedLoad-BR5_PredTN$Load05, BR5_PredTN$MedLoad+BR5_PredTN$Load95)), xlab = 'True BR5 TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of BR5 TN (mg N/s)')
arrows(BR5_PredTN$TrueLoad, BR5_PredTN$Load05, BR5_PredTN$TrueLoad, BR5_PredTN$Load95, length=0.05, angle=90, code=3)
lines(c(-10,100), c(-10,100))
dev.off()

png('BR5_TrueTNLoadvsWRTDSLoad_log.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], BR5_PredTN$MedLoad[BR5_PredTN$True > 1], 
     ylim=range(c(BR5_PredTN$MedLoad[BR5_PredTN$True > 1]-BR5_PredTN$Load05[BR5_PredTN$True > 1], BR5_PredTN$MedLoad[BR5_PredTN$True > 1]+BR5_PredTN$Load95[BR5_PredTN$True > 1])), 
     xlab = expression(paste('log'[10] ~ 'True BR5 TN Load (mg N/s)')), ylab = expression(paste('log'[10] ~ 'WRTDS Predicted Load of BR5 TN (mg N/s)')), log='xy')
arrows(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], BR5_PredTN$Load05[BR5_PredTN$True > 1], 
       BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], BR5_PredTN$Load95[BR5_PredTN$True > 1], 
       length=0.05, angle=90, code=3)
lines(c(1e-7,100), c(1e-7,100))
dev.off()

png('BR5_TrueTNLoadvsWRTDSLoad_QLQ.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR5_PredTN$TrueLoad, BR5_PredTN$MedLoadQLQ, ylim=range(c(BR5_PredTN$MedLoadQLQ-BR5_PredTN$Load05QLQ, BR5_PredTN$MedLoadQLQ+BR5_PredTN$Load95QLQ)), xlab = 'True BR5 TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of BR5 TN (mg N/s)')
arrows(BR5_PredTN$TrueLoad, BR5_PredTN$Load05QLQ, BR5_PredTN$TrueLoad, BR5_PredTN$Load95QLQ, length=0.05, angle=90, code=3)
lines(c(-10,100), c(-10,100))
dev.off()

png('BR5_TrueTNLoadvsWRTDSLoad_QLQ_log.png', res = 300, units = 'in', height = 5, width = 5)
plot(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], BR5_PredTN$MedLoadQLQ[BR5_PredTN$True > 1], 
     ylim=range(c(BR5_PredTN$MedLoadQLQ[BR5_PredTN$True > 1]-BR5_PredTN$Load05QLQ[BR5_PredTN$True > 1], BR5_PredTN$MedLoadQLQ[BR5_PredTN$True > 1]+BR5_PredTN$Load95QLQ[BR5_PredTN$True > 1])), 
     xlab = expression(paste('log'[10] ~ 'True BR5 TN Load (mg N/s)')), ylab = expression(paste('log'[10] ~ 'WRTDS Predicted Load of BR5 TN (mg N/s)')), log='xy')
arrows(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], BR5_PredTN$Load05QLQ[BR5_PredTN$True > 1], 
       BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], BR5_PredTN$Load95QLQ[BR5_PredTN$True > 1], 
       length=0.05, angle=90, code=3)
lines(c(1e-7,100), c(1e-7,100))
dev.off()

BR5_PredTN$DiffMedLoad = BR5_PredTN$TrueLoad - BR5_PredTN$MedLoad
BR5_PredTN$DiffMedLoadQLQ = BR5_PredTN$TrueLoad - BR5_PredTN$MedLoadQLQ

#plot of regression y vs. x
plot(y = BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1], x = BR5_PredTN$MedLoad[BR5_PredTN$True > 1])
plot(y = BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1], x = I(Flows_BR5[BR5_PredTN$True > 1]))
plot(y = BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1], x = sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))
plot(y = BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1], x = cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))

#Regression w/o suspected outlier
lm_BR5Load = lm(BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1]) + BR5_PredTN$MedLoad[BR5_PredTN$True > 1])
lm_BR5Load_Interaction = lm(BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1])*BR5_PredTN$MedLoad[BR5_PredTN$True > 1])
lm_BR5Load_Dates = lm(BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1]) + BR5_PredTN$MedLoad[BR5_PredTN$True > 1] + sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))
lm_BR5Load_Dates_Interaction = lm(BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1])*BR5_PredTN$Med[BR5_PredTN$True > 1] + sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))
lm_BR5Load_log = lm(log10(BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1] + abs(min(BR5_PredTN$DiffMedLoad[BR5_PredTN$True > 1])) + .001) ~ BR5_PredTN$MedLoad[BR5_PredTN$True > 1] + sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))

summary(lm_BR5Load)
#plot(lm_BR5Load)

summary(lm_BR5Load_Interaction)
#plot(lm_BR5Load_Interaction)

summary(lm_BR5Load_Dates)
#plot(lm_BR5Load_Dates)

summary(lm_BR5Load_Dates_Interaction)
#plot(lm_BR5Load_Dates_Interaction)

summary(lm_BR5Load_log)
#plot(lm_BR5Load_log)

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
plot(y = POBR_PredTN$TrueLoad, x = POBR_PredTN$MedLoadQLQ)
plot(y = POBR_PredTN$TrueLoad, x = Flows_POBR*.3048^3)
plot(y = POBR_PredTN$TrueLoad, x = sin(2*pi*DatesNums_POBR/366))
plot(y = POBR_PredTN$TrueLoad, x = cos(2*pi*DatesNums_POBR/366))

#Fixme - this needs to be a censored regression because of the y censoring. Thresholds change over time. Great.
lm_POBRLoad = lm(POBR_PredTN$TrueLoad ~ I(Flows_POBR*.3048^3) + POBR_PredTN$MedLoad)
lm_POBRLoad_Dates = lm(POBR_PredTN$TrueLoad ~ I(Flows_POBR*.3048^3) + POBR_PredTN$MedLoad + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))
lm_POBRLoad_Dates_Interaction = lm(POBR_PredTN$TrueLoad ~ I(Flows_POBR*.3048^3)*POBR_PredTN$Med + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))
lm_POBRLoad_log = lm(log10(POBR_PredTN$TrueLoad) ~ POBR_PredTN$MedLoad + sin(2*pi*DatesNums_POBR/366) + cos(2*pi*DatesNums_POBR/366))

summary(lm_POBRLoad)
#plot(lm_POBRLoad)

summary(lm_POBRLoad_Dates)
#plot(lm_POBRLoad_Dates)

summary(lm_POBRLoad_Dates_Interaction)
#plot(lm_POBRLoad_Dates_Interaction)

summary(lm_POBRLoad_log)
#plot(lm_POBRLoad_log)

#  BR3 - a little better than above----
#plot of regression y vs. x
plot(y = BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], x = BR3_PredTN$MedLoad[BR3_PredTN$True > 1])
plot(y = BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], x = Flows_BR3[BR3_PredTN$True > 1])
plot(y = BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], x = sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))
plot(y = BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], x = cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))

#Regression w/o suspected outlier
lm_BR3Load = lm(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1])
lm_BR3Load_Interaction = lm(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1])*BR3_PredTN$MedLoad[BR3_PredTN$True > 1])
lm_BR3Load_Dates = lm(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))
lm_BR3Load_Dates_Interaction = lm(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1])*BR3_PredTN$Med[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))
lm_BR3Load_log = lm(log10(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1]) ~ BR3_PredTN$MedLoad[BR3_PredTN$True > 1] + sin(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR3[BR3_PredTN$True > 1]/366))

summary(lm_BR3Load)
#plot(lm_BR3Load)

summary(lm_BR3Load_Interaction)
#plot(lm_BR3Load_Interaction)

summary(lm_BR3Load_Dates)
#plot(lm_BR3Load_Dates)

summary(lm_BR3Load_Dates_Interaction)
#plot(lm_BR3Load_Dates_Interaction)

summary(lm_BR3Load_log)
#plot(lm_BR3Load_log)

#From the BR3Load, use Box-Cox transform to attempt normality
BoxCoxtest_BR3_True = boxcox(BR3_PredTN$TrueLoad[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1], lambda = seq(0,1, length = 1000))
lambda_BR3_True = BoxCoxtest_BR3_True$x[which(BoxCoxtest_BR3_True$y == max(BoxCoxtest_BR3_True$y))]

BR3_PredTN$DiffMedLoad_BC_True[BR3_PredTN$True > 1] = bcPower(U = BR3_PredTN$TrueLoad[BR3_PredTN$True > 1], gamma = 0, lambda = lambda_BR3_True, jacobian.adjusted = FALSE)
lm_BR3Load_Dates_BC_True = lm(BR3_PredTN$DiffMedLoad_BC_True[BR3_PredTN$True > 1] ~ I(Flows_BR3[BR3_PredTN$True > 1]) + BR3_PredTN$MedLoad[BR3_PredTN$True > 1])

summary(lm_BR3Load_Dates_BC_True)
#plot(lm_BR3Load_Dates_BC)

#plot of regression y vs. x
plot(y = BR3_PredTN$DiffMedLoad_BC_True, x = BR3_PredTN$MedLoad)
plot(y = BR3_PredTN$DiffMedLoad_BC_True, x = I(Flows_BR3))
plot(y = BR3_PredTN$DiffMedLoad_BC_True, x = sin(2*pi*DatesNums_BR3/366))
plot(y = BR3_PredTN$DiffMedLoad_BC_True, x = cos(2*pi*DatesNums_BR3/366))

#  BR5 - essentially same as above----
#plot of regression y vs. x
plot(y = BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], x = BR5_PredTN$MedLoad[BR5_PredTN$True > 1])
plot(y = BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], x = Flows_BR5[BR5_PredTN$True > 1])
plot(y = BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], x = sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))
plot(y = BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], x = cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))

#Regression w/o suspected outlier
lm_BR5Load = lm(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1]) + BR5_PredTN$MedLoad[BR5_PredTN$True > 1])
lm_BR5Load_Interaction = lm(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1])*BR5_PredTN$MedLoad[BR5_PredTN$True > 1])
lm_BR5Load_Dates = lm(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1]) + BR5_PredTN$MedLoad[BR5_PredTN$True > 1] + sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))
lm_BR5Load_Dates_Interaction = lm(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1])*BR5_PredTN$Med[BR5_PredTN$True > 1] + sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))
lm_BR5Load_log = lm(log10(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1]) ~ BR5_PredTN$MedLoad[BR5_PredTN$True > 1] + sin(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366) + cos(2*pi*DatesNums_BR5[BR5_PredTN$True > 1]/366))

summary(lm_BR5Load)
#plot(lm_BR5Load)

summary(lm_BR5Load_Interaction)
#plot(lm_BR5Load_Interaction)

summary(lm_BR5Load_Dates)
#plot(lm_BR5Load_Dates)

summary(lm_BR5Load_Dates_Interaction)
#plot(lm_BR5Load_Dates_Interaction)

summary(lm_BR5Load_log)
#plot(lm_BR5Load_log)

#From the BR5Load, use Box-Cox transform to attempt normality
BoxCoxtest_BR5_True = boxcox(BR5_PredTN$TrueLoad[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1])*BR5_PredTN$MedLoad[BR5_PredTN$True > 1], lambda = seq(0,2, length = 1000))
lambda_BR5_True = BoxCoxtest_BR5_True$x[which(BoxCoxtest_BR5_True$y == max(BoxCoxtest_BR5_True$y))]

BR5_PredTN$DiffMedLoad_BC_True[BR5_PredTN$True > 1] = bcPower(U = BR5_PredTN$TrueLoad[BR5_PredTN$True > 1], gamma = 0, lambda = lambda_BR5_True, jacobian.adjusted = FALSE)
lm_BR5Load_Dates_BC_True = lm(BR5_PredTN$DiffMedLoad_BC_True[BR5_PredTN$True > 1] ~ I(Flows_BR5[BR5_PredTN$True > 1])*BR5_PredTN$MedLoad[BR5_PredTN$True > 1])

summary(lm_BR5Load_Dates_BC_True)
#plot(lm_BR5Load_Dates_BC_True)

#plot of regression y vs. x
plot(y = BR5_PredTN$DiffMedLoad_BC_True, x = BR5_PredTN$MedLoad)
plot(y = BR5_PredTN$DiffMedLoad_BC_True, x = I(Flows_BR5))
plot(y = BR5_PredTN$DiffMedLoad_BC_True, x = sin(2*pi*DatesNums_BR5/366))
plot(y = BR5_PredTN$DiffMedLoad_BC_True, x = cos(2*pi*DatesNums_BR5/366))

plot(y = BR5_PredTN$MedLoad, x = I(Flows_BR5))
cor(y = BR5_PredTN$MedLoad, x = I(Flows_BR5))

# WRTDS Interpolation Tables for Pond Branch----
#Load the streamflow data into WRTDS format
Daily_POBR = readUserDaily(filePath = wd_RHESSysObs_POBR, fileName = f_POBRStreamflowCal, hasHeader = TRUE, separator = '\t', qUnit = 1, verbose = FALSE)
#Read the TN data
Sample_POBR = readUserSample(filePath = wd_RHESSysObs_POBR, fileName = f_POBRTNCal_WRTDS, hasHeader = TRUE, separator = '\t', verbose = FALSE)
#Set the required information
INFO_POBR = readUserInfo(filePath = wd_WRTDS_POBR_INFO, fileName = f_POBRWRTDSInfo, interactive = FALSE)
eList_POBR = mergeReport(INFO = INFO_POBR, Daily = Daily_POBR, Sample = Sample_POBR)

#Make edits to include left-censored data. Nore two different thresholds.
Sample_POBR$ConcLow = ifelse((as.Date(Sample_POBR$Date) >= '2012-03-15') & (Sample_POBR$ConcLow == 0.05), NA, Sample_POBR$ConcLow)
Sample_POBR$ConcLow = ifelse((as.Date(Sample_POBR$Date) < '2012-03-15') & (Sample_POBR$ConcLow == 0.01), NA, Sample_POBR$ConcLow)
Sample_POBR$ConcHigh = ifelse((as.Date(Sample_POBR$Date) >= '2012-03-15') & (Sample_POBR$ConcHigh == 0.05), 0.05, Sample_POBR$ConcHigh)
Sample_POBR$ConcHigh = ifelse((as.Date(Sample_POBR$Date) < '2012-03-15') & (Sample_POBR$ConcHigh == 0.01), 0.01, Sample_POBR$ConcHigh)

#Fix the eList
eList_POBR$Sample = Sample_POBR
eList_POBR = fixSampleFrame(eList_POBR)

saveResults(savePath = paste0(wd_WRTDS_POBR, '\\'), eList = eList_POBR)

#  Default WRTDS parameters----
WRTDSmod_POBR = modelEstimation(eList = eList_POBR, windowY = 7, windowQ = 2, windowS = .5, minNumObs = 100, minNumUncen = 50, edgeAdjust = TRUE)

setwd(wd_WRTDS_POBR)
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


# Default WRTDS parameters with Quadratic Flow Term----
WRTDSmod_POBR_QLQ = modelEstimation(eList = eList_POBR, windowY = 7, windowQ = 2, windowS = .5, minNumObs = 100, minNumUncen = 50, edgeAdjust = TRUE, QuadLogQ = TRUE)

png('ConcFluxTime_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod_POBR_QLQ)
plotFluxTimeDaily(WRTDSmod_POBR_QLQ)
dev.off()

png('ConcFluxTimeErrs_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod_POBR_QLQ)
plotFluxPred(WRTDSmod_POBR_QLQ)
dev.off()

plotResidPred(WRTDSmod_POBR_QLQ)
plotResidQ(WRTDSmod_POBR_QLQ)
plotResidTime(WRTDSmod_POBR_QLQ)
boxResidMonth(WRTDSmod_POBR_QLQ)
boxConcThree(WRTDSmod_POBR_QLQ)
boxQTwice(WRTDSmod_POBR_QLQ)
plotFluxHist(WRTDSmod_POBR_QLQ)
plotConcHist(WRTDSmod_POBR_QLQ)

png('BiasPlot_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod_POBR_QLQ, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_POBR_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,2.6,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr_POBR_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod_POBR, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

#Contour plots of the parameters
png('TabInt_mod1_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-300,100,50),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod1_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.1,0.1,0.02),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod1_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2,5,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod1_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,.5,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod1_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,.5,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod1_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.3,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()
png('TabLogFlow2_mod1_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.3,1.3,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 9)
dev.off()

# Model#2-5 WRTDS with Quadratic Flow Term----
WRTDSmod2_POBR_QLQ = modelEstimation(eList = eList_POBR, windowY = 7, windowQ = 2, windowS = .25, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE, QuadLogQ = TRUE)
WRTDSmod3_POBR_QLQ = modelEstimation(eList = eList_POBR, windowY = 4, windowQ = 3, windowS = .5, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE, QuadLogQ = TRUE)
WRTDSmod4_POBR_QLQ = modelEstimation(eList = eList_POBR, windowY = 2, windowQ = 3, windowS = .5, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE, QuadLogQ = TRUE)
WRTDSmod5_POBR_QLQ = modelEstimation(eList = eList_POBR, windowY = 4, windowQ = 5, windowS = .5, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE, QuadLogQ = TRUE)
WRTDSmod6_POBR_QLQ = modelEstimation(eList = eList_POBR, windowY = 4, windowQ = 5, windowS = .25, minNumObs = 50, minNumUncen = 50, edgeAdjust = TRUE, QuadLogQ = TRUE)

png('ConcFluxTime_mod2_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod2_POBR_QLQ)
plotFluxTimeDaily(WRTDSmod2_POBR_QLQ)
dev.off()

png('ConcFluxTimeErrs_mod2_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod2_POBR_QLQ)
plotFluxPred(WRTDSmod2_POBR_QLQ)
dev.off()

plotResidPred(WRTDSmod2_POBR_QLQ)
plotResidQ(WRTDSmod2_POBR_QLQ)
plotResidTime(WRTDSmod2_POBR_QLQ)
boxResidMonth(WRTDSmod2_POBR_QLQ)
boxConcThree(WRTDSmod2_POBR_QLQ)
boxQTwice(WRTDSmod2_POBR_QLQ)
plotFluxHist(WRTDSmod2_POBR_QLQ)
plotConcHist(WRTDSmod2_POBR_QLQ)

png('BiasPlot_mod2_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod2_POBR_QLQ, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_mod2_POBR_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod2_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,5,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr_mod2_POBR_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod2_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

#Contour plots of the parameters
png('TabInt_mod2_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-300,100,50),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod2_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.1,0.1,0.02),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod2_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2,7,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod2_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod2_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod2_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.3,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()
png('TabLogFlow2_mod2_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod2_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.3,1.3,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 9)
dev.off()

png('ConcFluxTime_mod3_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod3_POBR_QLQ)
plotFluxTimeDaily(WRTDSmod3_POBR_QLQ)
dev.off()

png('ConcFluxTimeErrs_mod3_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod3_POBR_QLQ)
plotFluxPred(WRTDSmod3_POBR_QLQ)
dev.off()

plotResidPred(WRTDSmod3_POBR_QLQ)
plotResidQ(WRTDSmod3_POBR_QLQ)
plotResidTime(WRTDSmod3_POBR_QLQ)
boxResidMonth(WRTDSmod3_POBR_QLQ)
boxConcThree(WRTDSmod3_POBR_QLQ)
boxQTwice(WRTDSmod3_POBR_QLQ)
plotFluxHist(WRTDSmod3_POBR_QLQ)
plotConcHist(WRTDSmod3_POBR_QLQ)

png('BiasPlot_mod3_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod3_POBR_QLQ, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_mod3_POBR_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod3_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,10,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr_mod3_POBR_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod3_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

#Contour plots of the parameters
png('TabInt_mod3_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-300,600,50),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod3_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.1,0.15,0.02),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod3_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2,15,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod3_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod3_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod3_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.4,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()
png('TabLogFlow2_mod3_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod3_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,3,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 9)
dev.off()

png('ConcFluxTime_mod4_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod4_POBR_QLQ)
plotFluxTimeDaily(WRTDSmod4_POBR_QLQ)
dev.off()

png('ConcFluxTimeErrs_mod4_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod4_POBR_QLQ)
plotFluxPred(WRTDSmod4_POBR_QLQ)
dev.off()

plotResidPred(WRTDSmod4_POBR_QLQ)
plotResidQ(WRTDSmod4_POBR_QLQ)
plotResidTime(WRTDSmod4_POBR_QLQ)
boxResidMonth(WRTDSmod4_POBR_QLQ)
boxConcThree(WRTDSmod4_POBR_QLQ)
boxQTwice(WRTDSmod4_POBR_QLQ)
plotFluxHist(WRTDSmod4_POBR_QLQ)
plotConcHist(WRTDSmod4_POBR_QLQ)

png('BiasPlot_mod4_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod4_POBR_QLQ, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_mod4_POBR_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod4_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,10,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr_mod4_POBR_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod4_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

#Contour plots of the parameters
png('TabInt_mod4_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-500,1500,50),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod4_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.1,0.5,0.02),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod4_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2,6,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod4_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod4_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod4_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()
png('TabLogFlow2_mod4_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod4_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,1,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 9)
dev.off()


png('ConcFluxTime_mod5_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod5_POBR_QLQ)
plotFluxTimeDaily(WRTDSmod5_POBR_QLQ)
dev.off()

png('ConcFluxTimeErrs_mod5_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod5_POBR_QLQ)
plotFluxPred(WRTDSmod5_POBR_QLQ)
dev.off()

plotResidPred(WRTDSmod5_POBR_QLQ)
plotResidQ(WRTDSmod5_POBR_QLQ)
plotResidTime(WRTDSmod5_POBR_QLQ)
boxResidMonth(WRTDSmod5_POBR_QLQ)
boxConcThree(WRTDSmod5_POBR_QLQ)
boxQTwice(WRTDSmod5_POBR_QLQ)
plotFluxHist(WRTDSmod5_POBR_QLQ)
plotConcHist(WRTDSmod5_POBR_QLQ)

png('BiasPlot_mod5_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod5_POBR_QLQ, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_mod5_POBR_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod5_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,10,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr_mod5_POBR_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod5_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

#Contour plots of the parameters
png('TabInt_mod5_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-500,1500,50),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod5_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.1,0.5,0.02),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod5_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2,6,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod5_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod5_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod5_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()
png('TabLogFlow2_mod5_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod5_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,1,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 9)
dev.off()


png('ConcFluxTime_mod6_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcTimeDaily(WRTDSmod6_POBR_QLQ)
plotFluxTimeDaily(WRTDSmod6_POBR_QLQ)
dev.off()

png('ConcFluxTimeErrs_mod6_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 6)
layout(rbind(c(1,2)))
plotConcPred(WRTDSmod6_POBR_QLQ)
plotFluxPred(WRTDSmod6_POBR_QLQ)
dev.off()

plotResidPred(WRTDSmod6_POBR_QLQ)
plotResidQ(WRTDSmod6_POBR_QLQ)
plotResidTime(WRTDSmod6_POBR_QLQ)
boxResidMonth(WRTDSmod6_POBR_QLQ)
boxConcThree(WRTDSmod6_POBR_QLQ)
boxQTwice(WRTDSmod6_POBR_QLQ)
plotFluxHist(WRTDSmod6_POBR_QLQ)
plotConcHist(WRTDSmod6_POBR_QLQ)

png('BiasPlot_mod6_POBR_QLQ.png', res = 300, units ='in', width = 12, height = 12)
fluxBiasMulti(WRTDSmod6_POBR_QLQ, cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5)
dev.off()

png('ContourPlotMean_mod6_POBR_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod6_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,10,0.5),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 1)
dev.off()

png('ContourPlotErr_mod6_POBR_QLQ.png', res = 300, units ='in', width = 6, height = 6)
plotContours(WRTDSmod6_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()

#Contour plots of the parameters
png('TabInt_mod6_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod6_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-500,1500,50),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 4)
dev.off()
png('TabYear_mod6_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod6_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-0.1,0.5,0.02),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 5)
dev.off()
png('TabLogFlow_mod6_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod6_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-2,6,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 6)
dev.off()
png('TabSinYear_mod6_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod6_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 7)
dev.off()
png('TabCosYear_mod6_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod6_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-1,1,0.2),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 8)
dev.off()
png('TabLogErr_mod6_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod6_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(0,.5,0.05),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 2)
dev.off()
png('TabLogFlow2_mod6_POBR_QLQ.png', units = 'in', res = 300, height = 7, width = 7)
plotContours(WRTDSmod6_POBR_QLQ, yearStart = 1999, yearEnd = 2011, contourLevels=seq(-.5,1,0.1),qUnit=1, qBottom = 0.001, qTop = 50, whatSurface = 9)
dev.off()

#  LOOCV SE for sample and discharge----
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

sum(WRTDSmod_POBR_QLQ$Sample$SE^2)
sum(WRTDSmod2_POBR_QLQ$Sample$SE^2)
sum(WRTDSmod3_POBR_QLQ$Sample$SE^2)
sum(WRTDSmod4_POBR_QLQ$Sample$SE^2)
sum(WRTDSmod5_POBR_QLQ$Sample$SE^2)
sum(WRTDSmod6_POBR_QLQ$Sample$SE^2)

sum(WRTDSmod_POBR_QLQ$Daily$SE^2)
sum(WRTDSmod2_POBR_QLQ$Daily$SE^2)
sum(WRTDSmod3_POBR_QLQ$Daily$SE^2)
sum(WRTDSmod4_POBR_QLQ$Daily$SE^2)
sum(WRTDSmod5_POBR_QLQ$Daily$SE^2)
sum(WRTDSmod6_POBR_QLQ$Daily$SE^2)

#  MSE for sample----
mean((WRTDSmod_POBR$Sample$ConcHat - WRTDSmod_POBR$Sample$ConcAve))^2 + sum((WRTDSmod_POBR$Sample$ConcHat - WRTDSmod_POBR$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)
mean((WRTDSmod2_POBR$Sample$ConcHat - WRTDSmod2_POBR$Sample$ConcAve))^2 + sum((WRTDSmod2_POBR$Sample$ConcHat - WRTDSmod2_POBR$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)
mean((WRTDSmod3_POBR$Sample$ConcHat - WRTDSmod3_POBR$Sample$ConcAve))^2 + sum((WRTDSmod3_POBR$Sample$ConcHat - WRTDSmod3_POBR$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)
mean((WRTDSmod4_POBR$Sample$ConcHat - WRTDSmod4_POBR$Sample$ConcAve))^2 + sum((WRTDSmod4_POBR$Sample$ConcHat - WRTDSmod4_POBR$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)
mean((WRTDSmod5_POBR$Sample$ConcHat - WRTDSmod5_POBR$Sample$ConcAve))^2 + sum((WRTDSmod5_POBR$Sample$ConcHat - WRTDSmod5_POBR$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)
mean((WRTDSmod6_POBR$Sample$ConcHat - WRTDSmod6_POBR$Sample$ConcAve))^2 + sum((WRTDSmod6_POBR$Sample$ConcHat - WRTDSmod6_POBR$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)

mean((WRTDSmod_POBR_QLQ$Sample$ConcHat - WRTDSmod_POBR_QLQ$Sample$ConcAve))^2 + sum((WRTDSmod_POBR_QLQ$Sample$ConcHat - WRTDSmod_POBR_QLQ$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)
mean((WRTDSmod2_POBR_QLQ$Sample$ConcHat - WRTDSmod2_POBR_QLQ$Sample$ConcAve))^2 + sum((WRTDSmod2_POBR_QLQ$Sample$ConcHat - WRTDSmod2_POBR_QLQ$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)
mean((WRTDSmod3_POBR_QLQ$Sample$ConcHat - WRTDSmod3_POBR_QLQ$Sample$ConcAve))^2 + sum((WRTDSmod3_POBR_QLQ$Sample$ConcHat - WRTDSmod3_POBR_QLQ$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)
mean((WRTDSmod4_POBR_QLQ$Sample$ConcHat - WRTDSmod4_POBR_QLQ$Sample$ConcAve))^2 + sum((WRTDSmod4_POBR_QLQ$Sample$ConcHat - WRTDSmod4_POBR_QLQ$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)
mean((WRTDSmod5_POBR_QLQ$Sample$ConcHat - WRTDSmod5_POBR_QLQ$Sample$ConcAve))^2 + sum((WRTDSmod5_POBR_QLQ$Sample$ConcHat - WRTDSmod5_POBR_QLQ$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)
mean((WRTDSmod6_POBR_QLQ$Sample$ConcHat - WRTDSmod6_POBR_QLQ$Sample$ConcAve))^2 + sum((WRTDSmod6_POBR_QLQ$Sample$ConcHat - WRTDSmod6_POBR_QLQ$Sample$ConcAve)^2)/(nrow(Sample_POBR)-1)

#  Make tables for selected model----
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

#   Write the tables----
setwd(wd_WRTDS_POBR)
options(scipen = 999)
write.table(round(TabInt_POBR,5), file = f_POBRWRTDS_TabInt, sep = '\t', row.names = attr(WRTDSmod5m_POBR$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_POBR$surfaces, which = 'Year'))
write.table(signif(TabYear_POBR,4), file = f_POBRWRTDS_TabYear, sep = '\t', row.names = attr(WRTDSmod5m_POBR$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_POBR$surfaces, which = 'Year'))
write.table(signif(TabLogQ_POBR,4), file = f_POBRWRTDS_TabLogQ, sep = '\t', row.names = attr(WRTDSmod5m_POBR$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_POBR$surfaces, which = 'Year'))
write.table(signif(TabSinYear_POBR,4), file = f_POBRWRTDS_TabSinYear, sep = '\t', row.names = attr(WRTDSmod5m_POBR$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_POBR$surfaces, which = 'Year'))
write.table(signif(TabCosYear_POBR,4), file = f_POBRWRTDS_TabCosYear, sep = '\t', row.names = attr(WRTDSmod5m_POBR$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_POBR$surfaces, which = 'Year'))
write.table(round(TabLogErr_POBR,5), file = f_POBRWRTDS_TabLogErr, sep = '\t', row.names = attr(WRTDSmod5m_POBR$surfaces, which = 'LogQ'), col.names = attr(WRTDSmod5m_POBR$surfaces, which = 'Year'))
options(scipen = 0)

# Load WRTDS Interpolation Tables----
setwd(wd_WRTDS_POBR)
#Pond Branch
TabInt_POBR = as.matrix(read.table(file = f_POBRWRTDS_TabInt, sep = '\t', header = TRUE, check.names = FALSE))
TabYear_POBR = as.matrix(read.table(file = f_POBRWRTDS_TabYear, sep = '\t', header = TRUE, check.names = FALSE))
TabLogQ_POBR = as.matrix(read.table(file = f_POBRWRTDS_TabLogQ, sep = '\t', header = TRUE, check.names = FALSE))
TabSinYear_POBR = as.matrix(read.table(file = f_POBRWRTDS_TabSinYear, sep = '\t', header = TRUE, check.names = FALSE))
TabCosYear_POBR = as.matrix(read.table(file = f_POBRWRTDS_TabCosYear, sep = '\t', header = TRUE, check.names = FALSE))
TabLogErr_POBR = as.matrix(read.table(file = f_POBRWRTDS_TabLogErr, sep = '\t', header = TRUE, check.names = FALSE))

#  Compare predictions of Pond Branch WRTDS and the other regression to true value----
POBR_PredTN_POBRWRTDS = matrix(NA, ncol = 3, nrow = length(Flows_POBR))
for(i in 1:length(Flows_POBR)){
  POBR_PredTN_POBRWRTDS[i,] = predictWRTDS(Date = as.character(as.Date(Dates_POBR))[i], Flow = Flows_POBR[i], 
                                           rowt = as.numeric(rownames(TabInt_POBR)), colt = as.numeric(colnames(TabInt_POBR)), 
                                           TabInt = TabInt_POBR, TabYear = TabYear_POBR, TabLogQ = TabLogQ_POBR, 
                                           TabSinYear = TabSinYear_POBR, TabCosYear = TabCosYear_POBR, TabLogErr = TabLogErr_POBR)
}
rm(i)
colnames(POBR_PredTN_POBRWRTDS) = c('05', 'Med', '95')

POBR_PredTN_POBRWRTDS = as.data.frame(POBR_PredTN_POBRWRTDS)
POBR_PredTN_POBRWRTDS$True = BES_TN_d$POBR$TN..mg.N.L.[which((!is.na(BES_TN_d$POBR$TN..mg.N.L.)) & (as.Date(BES_TN_d$POBR$SortDate) < '2014-01-01'))]
POBR_PredTN_POBRWRTDS$DiffMed = POBR_PredTN_POBRWRTDS$True - POBR_PredTN_POBRWRTDS$Med

png('POBR_TrueTNvsPOBRWRTDS.png', res = 300, units = 'in', height = 5, width = 10)
layout(rbind(c(1,2)))
plot(POBR_PredTN$True, POBR_PredTN$Med, ylim = c(0,9), xlim = c(0,9), xlab = 'True Pond Branch TN (mg N/L)', ylab = 'WRTDS Predicted Pond Branch TN (mg N/L)', main = 'Basin WRTDS Model')
arrows(POBR_PredTN$True, POBR_PredTN$`05`, POBR_PredTN$True, POBR_PredTN$`95`, length=0.05, angle=90, code=3)
lines(c(-10,10), c(-10,10))

plot(POBR_PredTN_POBRWRTDS$True, POBR_PredTN_POBRWRTDS$Med, ylim = c(0,9), xlim = c(0,9), xlab = 'True Pond Branch TN (mg N/L)', ylab = 'WRTDS Predicted Pond Branch TN (mg N/L)', main = 'Pond Branch WRTDS Model')
arrows(POBR_PredTN_POBRWRTDS$True, POBR_PredTN_POBRWRTDS$`05`, POBR_PredTN_POBRWRTDS$True, POBR_PredTN_POBRWRTDS$`95`, length=0.05, angle=90, code=3)
lines(c(-10,10), c(-10,10))
dev.off()

png('POBR_TrueTNvsPOBRWRTDS_noArrows.png', res = 300, units = 'in', height = 5, width = 10)
layout(rbind(c(1,2)))
plot(POBR_PredTN$True, POBR_PredTN$Med, ylim = c(0,9), xlim = c(0,9), xlab = 'True Pond Branch TN (mg N/L)', ylab = 'WRTDS Predicted Pond Branch TN (mg N/L)', main = 'Basin WRTDS Model')
lines(c(-10,10), c(-10,10))

plot(POBR_PredTN_POBRWRTDS$True, POBR_PredTN_POBRWRTDS$Med, ylim = c(0,9), xlim = c(0,9), xlab = 'True Pond Branch TN (mg N/L)', ylab = 'WRTDS Predicted Pond Branch TN (mg N/L)', main = 'Pond Branch WRTDS Model')
lines(c(-10,10), c(-10,10))
dev.off()

#Color by difference from the median
scaleRange = c(-2,4)
scaleBy = 0.5
Pal = rev(rainbow((scaleRange[2] - scaleRange[1])/scaleBy))

#Scatter histogram colored by the difference in prediced concentration
plot(log(Flows_POBR), POBR_PredTN_POBRWRTDS$True, col = colFun(POBR_PredTN_POBRWRTDS$DiffMed))
plot(Dates_POBR, POBR_PredTN_POBRWRTDS$True, col = colFun(POBR_PredTN_POBRWRTDS$DiffMed))
plot(POBR_PredTN_POBRWRTDS$True, POBR_PredTN_POBRWRTDS$Med, col = colFun(POBR_PredTN_POBRWRTDS$DiffMed))
plot(Flows_POBR, POBR_PredTN_POBRWRTDS$DiffMed, col = POBR_PredTN$cols)

#Load
POBR_PredTN_POBRWRTDS$LowTNLim = NA
POBR_PredTN_POBRWRTDS$LowTNLim[POBR_PredTN_POBRWRTDS$True == min(POBR_PredTN_POBRWRTDS$True[1:650])] = 1
POBR_PredTN_POBRWRTDS$LowTNLim[POBR_PredTN_POBRWRTDS$True == min(POBR_PredTN_POBRWRTDS$True[651:length(POBR_PredTN_POBRWRTDS$True)])] = 2

#Use WRTDS regression for predicted concentration and transform to loads. Compare predicted to true load.
#mg/s unit
POBR_PredTN_POBRWRTDS$TrueLoad = POBR_PredTN_POBRWRTDS$True*1000*Flows_POBR*.3048^3
POBR_PredTN_POBRWRTDS$Load05 = POBR_PredTN_POBRWRTDS$`05`*1000*Flows_POBR*.3048^3
POBR_PredTN_POBRWRTDS$MedLoad = POBR_PredTN_POBRWRTDS$Med*1000*Flows_POBR*.3048^3
POBR_PredTN_POBRWRTDS$Load95 = POBR_PredTN_POBRWRTDS$`95`*1000*Flows_POBR*.3048^3

png('POBR_TrueTNLoadvsPOBRWRTDSLoad.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$MedLoad, ylim=range(c(POBR_PredTN_POBRWRTDS$MedLoad-POBR_PredTN_POBRWRTDS$Load05, POBR_PredTN_POBRWRTDS$MedLoad+POBR_PredTN_POBRWRTDS$Load95)), xlab = 'True Pond Branch TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of Pond Branch TN (mg N/s)')
arrows(POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$Load05, POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$Load95, length=0.05, angle=90, code=3)
lines(c(-10,100), c(-10,100))
dev.off()

png('POBR_TrueTNLoadvsPOBRWRTDSLoad_log.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$MedLoad, ylim=range(c(POBR_PredTN_POBRWRTDS$MedLoad-POBR_PredTN_POBRWRTDS$Load05, POBR_PredTN_POBRWRTDS$MedLoad+POBR_PredTN_POBRWRTDS$Load95)), xlab = 'True Pond Branch TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of Pond Branch TN (mg N/s)', log = 'xy')
arrows(POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$Load05, POBR_PredTN_POBRWRTDS$TrueLoad, POBR_PredTN_POBRWRTDS$Load95, length=0.05, angle=90, code=3)
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
scatterHist_mod(x = log10(POBR_PredTN_POBRWRTDS$TrueLoad), y = log10(POBR_PredTN_POBRWRTDS$MedLoad), 
                density = FALSE, freq = FALSE, x.breaks = 100, y.breaks = 100, correl = FALSE, 
                xlab = expression(paste('log'[10] ~ '(True TN Load [mg N/s])')), ylab = expression(paste('log'[10] ~ 'Predicted TN Load (mg N/s)')), 
                cex.lab = 1.5, cex.axis = 1.5, title = '', pch = 16, ellipse = FALSE, smooth = FALSE, col = POBR_PredTN_POBRWRTDS$cols)
dev.off()

POBR_PredTN_POBRWRTDS$DiffMedLoad = POBR_PredTN_POBRWRTDS$TrueLoad - POBR_PredTN_POBRWRTDS$MedLoad

#TN Load as modeled from load correction regression - not recommended model. Used to compare with WRTDS
POBR_PredTN$ModeledTNLoad = POBR_PredTN$MedLoad + lm_POBRLoad_Dates_Interaction$fitted.values

png('POBR_TrueTNLoadvsWRTDSLoad+CorrModel_log_colors.png', res = 300, units = 'in', height = 5, width = 5)
plot(POBR_PredTN$TrueLoad, POBR_PredTN$ModeledTNLoad,
     ylim = c(1e-4, 1e+2), xlim = c(1e-4, 1e+2),
     xlab = 'True Pond Branch TN Load (mg N/s)', ylab = 'WRTDS Predicted Load of Pond Branch TN (mg N/s)', log = 'xy')
par(new = TRUE)
plot(POBR_PredTN$TrueLoad[POBR_PredTN$LowTNLim == 1], POBR_PredTN$ModeledTNLoad[POBR_PredTN$LowTNLim == 1], 
     ylim = c(1e-4, 1e+2), xlim = c(1e-4, 1e+2), xlab = '', ylab = '', 
     log = 'xy', col = 'blue', axes = FALSE)
par(new = TRUE)
plot(POBR_PredTN$TrueLoad[POBR_PredTN$LowTNLim == 2], POBR_PredTN$ModeledTNLoad[POBR_PredTN$LowTNLim == 2], 
     ylim = c(1e-4, 1e+2), xlim = c(1e-4, 1e+2), xlab = '', ylab = '', 
     log = 'xy', col = 'red', axes = FALSE)
lines(c(1e-7,1e7), c(1e-7,1e7))
legend('topleft', title = 'Detection Limits (mg N/L)', legend = c('0.01', '0.05'), pch = 1, col = c('blue', 'red'))
dev.off()

# Make a map of these water quality sampling locations----
setwd(dir_WChem)

world = read.csv(paste0(dir_worldfile, '\\', f_worldfile), stringsAsFactors = FALSE)

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

# Make a map of the land ID in the world file----
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
#  Obtain the fraction of developed land in each of the hillslopes, and for the basin - using RHESSys land use types----
FracDev = vector('numeric', length = length(unique(world$hillID)))
for (i in 1:length(FracDev)){
  FracDev[i] = length(which((world$hillID == uhills[i]) & ((world$patchLandID == 3) | (world$patchLandID == 4))))/length(which((world$hillID == uhills[i])))
}
#Because this is a 2-component mixture model, assume that the remainder is undeveloped land
FracUnDev = 1 - FracDev

#  Use export coefficient/source-contribution model----
#This is a signal + background model at BARN because BARN = POBR + all other catchments

#Get flows and dates for BARN outlet
#cfs unit
Flows_BARN = BES_TN_d$BARN$Flow[which((as.Date(BES_TN_d$BARN$SortDate) <= '2014-01-01'))]
Dates_BARN = BES_TN_d$BARN$SortDate[which((as.Date(BES_TN_d$BARN$SortDate) <= '2014-01-01'))]
BARN_PredTN = matrix(NA, ncol = 6, nrow = length(Flows_BARN))
for(i in 1:length(Flows_BARN)){
  BARN_PredTN[i,c(1,2,3)] = predictWRTDS(Date = as.character(as.Date(Dates_BARN))[i], Flow = Flows_BARN[i], 
                                 rowt = as.numeric(rownames(TabInt)), colt = as.numeric(colnames(TabInt)), 
                                 TabInt = TabInt, TabYear = TabYear, TabLogQ = TabLogQ, 
                                 TabSinYear = TabSinYear, TabCosYear = TabCosYear, TabLogErr = TabLogErr)
  BARN_PredTN[i,c(4,5,6)] = predictWRTDS(Date = as.character(as.Date(Dates_BARN))[i], Flow = Flows_BARN[i], 
                                         rowt = as.numeric(rownames(TabInt_QLQ)), colt = as.numeric(colnames(TabInt_QLQ)), 
                                         TabInt = TabInt_QLQ, TabYear = TabYear_QLQ, TabLogQ = TabLogQ_QLQ, 
                                         TabSinYear = TabSinYear_QLQ, TabCosYear = TabCosYear_QLQ, TabLogErr = TabLogErr_QLQ, TabLogQ2 = TabLogQ2_QLQ)
}
rm(i)
colnames(BARN_PredTN) = c('05', 'Med', '95', '05QLQ', 'MedQLQ', '95QLQ')

BARN_PredTN = as.data.frame(BARN_PredTN)
BARN_PredTN$True = BES_TN_d$BARN$TN..mg.N.L.[which((as.Date(BES_TN_d$BARN$SortDate) <= '2014-01-01'))]

#Load
BARN_PredTN$LowTNLim = NA

#concentration to loads
#mg/s unit
BARN_PredTN$TrueLoad = BARN_PredTN$True*1000*Flows_BARN*.3048^3
BARN_PredTN$Load05 = BARN_PredTN$`05`*1000*Flows_BARN*.3048^3
BARN_PredTN$MedLoad = BARN_PredTN$Med*1000*Flows_BARN*.3048^3
BARN_PredTN$Load95 = BARN_PredTN$`95`*1000*Flows_BARN*.3048^3
BARN_PredTN$Load05QLQ = BARN_PredTN$`05QLQ`*1000*Flows_BARN*.3048^3
BARN_PredTN$MedLoadQLQ = BARN_PredTN$MedQLQ*1000*Flows_BARN*.3048^3
BARN_PredTN$Load95QLQ = BARN_PredTN$`95QLQ`*1000*Flows_BARN*.3048^3

#All dates for Pond Branch
#cfs unit
Flows_POBR_AllDates = BES_TN_d$POBR$Flow[which((as.Date(BES_TN_d$POBR$SortDate) <= '2014-01-01') & (as.Date(BES_TN_d$POBR$SortDate) >= '1999-11-15'))]
Dates_POBR_AllDates = BES_TN_d$POBR$SortDate[which((as.Date(BES_TN_d$POBR$SortDate) <= '2014-01-01') & (as.Date(BES_TN_d$POBR$SortDate) >= '1999-11-15'))]
POBR_AllDates_PredTN = matrix(NA, ncol = 3, nrow = length(Flows_POBR_AllDates))
for(i in 1:length(Flows_POBR_AllDates)){
  POBR_AllDates_PredTN[i,] = predictWRTDS(Date = as.character(as.Date(Dates_POBR_AllDates))[i], Flow = Flows_POBR_AllDates[i], 
                                          rowt = as.numeric(rownames(TabInt_POBR)), colt = as.numeric(colnames(TabInt_POBR)), 
                                          TabInt = TabInt_POBR, TabYear = TabYear_POBR, TabLogQ = TabLogQ_POBR, 
                                          TabSinYear = TabSinYear_POBR, TabCosYear = TabCosYear_POBR, TabLogErr = TabLogErr_POBR)
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
#mg/s unit
POBR_AllDates_PredTN$TrueLoad = POBR_AllDates_PredTN$True*1000*Flows_POBR_AllDates*.3048^3
POBR_AllDates_PredTN$Load05 = POBR_AllDates_PredTN$`05`*1000*Flows_POBR_AllDates*.3048^3
POBR_AllDates_PredTN$MedLoad = POBR_AllDates_PredTN$Med*1000*Flows_POBR_AllDates*.3048^3
POBR_AllDates_PredTN$Load95 = POBR_AllDates_PredTN$`95`*1000*Flows_POBR_AllDates*.3048^3

#   ECM based on RHESSys derived land use type----
#Export coefficient for Pond Branch = undeveloped export coefficient
#All of the area of pond branch is considered undeveloped
EC_Undev = POBR_AllDates_PredTN$MedLoad/sum(Area.Hills[c(3,4),2])
#Export coefficient for developed = (Baisman load - undeveloped load)/(developed land area)
EC_Dev = (BARN_PredTN$MedLoad - EC_Undev*sum(Area.Hills[,2]*FracUnDev))/sum(Area.Hills[,2]*FracDev)
EC_DevQLQ = (BARN_PredTN$MedLoadQLQ - EC_Undev*sum(Area.Hills[,2]*FracUnDev))/sum(Area.Hills[,2]*FracDev)

#There are some negative values, so maybe this can be attributed to the in-stream losses. Plot dates that these happen
plot(y = EC_Dev, x = Dates_POBR_AllDates, ylim = c(-.0005,.002), xlim = c(10500, 16500))
par(new = TRUE)
plot(y = EC_Dev[EC_Dev < 0], x = Dates_POBR_AllDates[EC_Dev < 0], ylim = c(-.0005, .002), xlim = c(10500, 16500), col = 'red')

plot(y = EC_DevQLQ, x = Dates_POBR_AllDates, ylim = c(-.0005,.002), xlim = c(10500, 16500))
par(new = TRUE)
plot(y = EC_DevQLQ[EC_DevQLQ < 0], x = Dates_POBR_AllDates[EC_DevQLQ < 0], ylim = c(-.0005,.002), xlim = c(10500, 16500), col = 'red')

#They are almost all during a drought! Might be a result of not accounting for in-stream losses.
Flows_POBR_AllDates[EC_Dev < 0]
Flows_POBR_AllDates[EC_DevQLQ < 0]

#Fixme: setting negative values to 0 for now.
#EC_Undev[EC_Dev < 0] = BARN_PredTN$MedLoad - EC_Dev[EC_Dev < 0]*sum(Area.Hills[,2]*FracDev)
#EC_Dev[EC_Dev < 0] = 0

#   ECM based on impervious fraction of land (original input to RHESSys land use type)----
impFrac = raster(x = paste0(dir_worldfile, '\\', f_BaismanImperviousFrac))
impFrac = projectRaster(impFrac, crs = CRS(pCRS))

#Add impervious fraction to worldfile information
world$ImpFrac = raster::extract(x = impFrac, y = world)

#    Obtain upstream contributing patches for each sampling location----
#An ArcGIS script was used to get the raster cells corresponding to the upstream areas for each sampling site.
#For that script, the tif file extension did not want to be added. So, need to add that manually to the filename.
#The nice thing about this file.rename function is that it will run and return FALSE if the task has already been completed.
file.rename(from = paste0(dir_ContributingAreas, "\\", grep(grep(list.files(dir_ContributingAreas), pattern = '.', fixed = TRUE, invert = TRUE, value = TRUE), pattern = '_', fixed = TRUE, value = TRUE)),
            to = paste0(dir_ContributingAreas, "\\", grep(grep(list.files(dir_ContributingAreas), pattern = '.', fixed = TRUE, invert = TRUE, value = TRUE), pattern = '_', fixed = TRUE, value = TRUE), '.tif'))

#    Add patch identifiers for each sampling site to the world dataframe, one column for each site----
setwd(wd_BESN)
fs = grep(grep(list.files(dir_ContributingAreas), pattern = '.tif', fixed = TRUE, value = TRUE), pattern = '_', fixed = TRUE, value = TRUE)
for (i in 1:length(fs)){
  temp = raster::extract(x = projectRaster(raster(x = paste0(dir_ContributingAreas, "\\", fs[i])), crs = pCRS), y = world)
  temp[!is.na(temp)] = 1
  temp = as.data.frame(temp)
  colnames(temp) = strsplit(fs[i], split = '.tif', fixed = TRUE)
  world@data = cbind(world@data, temp)
  
  #Create simple maps of the contributing areas for these sampling locations
  png(paste0('ECM_ContributingAreas_', strsplit(fs[i], split = '.', fixed = TRUE)[[1]][1], '.png'), res = 200, height = 5, width = 5, units = 'in')
  plot(BaisRun_Outlet)
  plot(world[which(world@data[colnames(temp)] == 1),], add = T, pch = 15)
  dev.off()
}
rm(i,temp, fs)

#    Compute impervious area for Pond Branch upstream of water quality site----
AreaImp_K2S = 0
Area_K2S = length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_6 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_6 == 1))){
  AreaImp_K2S = AreaImp_K2S + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_6 == 1)[i]]*res^2
}
rm(i)

#    Compute the impervious area of each hillslope and for the basin----
Area.basin_imp = 0
for (i in 1:length(unique(world$patchID))){
  Area.basin_imp = Area.basin_imp + world$ImpFrac[which(duplicated(world$patchID) == FALSE)[i]]*res^2
}
#Correct for Pond Branch impervious correction
Area.basin_imp = Area.basin_imp - AreaImp_K2S
Area.basin_imp/Area.basin #5.1% impervious

Area.Hills = cbind(Area.Hills, rep(0, nrow(Area.Hills)))
for (h in 1:length(uhills)){
  #Get the impervious surface area in this hillslope
  for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$hillID == h))){
    Area.Hills[h,3] = Area.Hills[h,3] + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$hillID == h)[i]]*res^2
  }
}
#Correct for Pond Branch in Hillslope 4
Area.Hills[4,3] = Area.Hills[4,3] - AreaImp_K2S

#Add fraction of impervious surface to the dataset
Area.Hills = cbind(Area.Hills, Area.Hills[,3]/Area.Hills[,2])

#    Estimate the Pond Branch "Undeveloped" model----
#Will treat the Pond Branch upstream of gauge as all undeveloped land because fraction impervious is 7e-5.
#Fixme: neglecting the 28 m^2 of Pond Branch Hillslope 4 that is developed here.
# Also neglecting that the streamflow was measured slightly upstream of water quality.
AreaImp_K2S = 0

#     Method 1: Using quantiles from WRTDS----
EC_Undev05 = POBR_AllDates_PredTN$Load05/Area_K2S
EC_Undev = POBR_AllDates_PredTN$MedLoad/Area_K2S
EC_Undev95 = POBR_AllDates_PredTN$Load95/Area_K2S

#Export coefficient for developed = (Baisman load - undeveloped load)/(developed land area)
EC_Dev05 = (BARN_PredTN$Load05 - EC_Undev05*sum(Area.Hills[,2] - Area.Hills[,3]))/sum(Area.Hills[,3])
EC_Dev = (BARN_PredTN$MedLoad - EC_Undev*sum(Area.Hills[,2] - Area.Hills[,3]))/sum(Area.Hills[,3])
EC_Dev95 = (BARN_PredTN$Load95 - EC_Undev95*sum(Area.Hills[,2] - Area.Hills[,3]))/sum(Area.Hills[,3])
EC_DevQLQ05 = (BARN_PredTN$Load05QLQ - EC_Undev05*sum(Area.Hills[,2] - Area.Hills[,3]))/sum(Area.Hills[,3])
EC_DevQLQ = (BARN_PredTN$MedLoadQLQ - EC_Undev*sum(Area.Hills[,2] - Area.Hills[,3]))/sum(Area.Hills[,3])
EC_DevQLQ95 = (BARN_PredTN$Load95QLQ - EC_Undev95*sum(Area.Hills[,2] - Area.Hills[,3]))/sum(Area.Hills[,3])

#There are some negative values, so maybe this can be attributed to the in-stream losses. Plot dates that these happen
plot(y = EC_Dev05, x = Dates_POBR_AllDates, ylim = c(-0.0005, .007), xlim = c(10500, 16500))
par(new = TRUE)
plot(y = EC_Dev05[EC_Dev05 < 0], x = Dates_POBR_AllDates[EC_Dev05 < 0], ylim = c(-0.0005, .007), xlim = c(10500, 16500), col = 'red')

plot(y = EC_Dev, x = Dates_POBR_AllDates, ylim = c(-0.0005, .007), xlim = c(10500, 16500))
par(new = TRUE)
plot(y = EC_Dev[EC_Dev < 0], x = Dates_POBR_AllDates[EC_Dev < 0], ylim = c(-0.0005, .007), xlim = c(10500, 16500), col = 'red')

plot(y = EC_Dev95, x = Dates_POBR_AllDates, ylim = c(-0.005, .01), xlim = c(10500, 16500))
par(new = TRUE)
plot(y = EC_Dev95[EC_Dev95 < 0], x = Dates_POBR_AllDates[EC_Dev95 < 0], ylim = c(-0.005, .01), xlim = c(10500, 16500), col = 'red')

plot(y = EC_DevQLQ05, x = Dates_POBR_AllDates, ylim = c(-0.0005, .007), xlim = c(10500, 16500))
par(new = TRUE)
plot(y = EC_DevQLQ05[EC_DevQLQ05 < 0], x = Dates_POBR_AllDates[EC_DevQLQ05 < 0], ylim = c(-0.0005, .007), xlim = c(10500, 16500), col = 'red')

plot(y = EC_DevQLQ, x = Dates_POBR_AllDates, ylim = c(-0.0005, .0035), xlim = c(10500, 16500))
par(new = TRUE)
plot(y = EC_DevQLQ[EC_DevQLQ < 0], x = Dates_POBR_AllDates[EC_DevQLQ < 0], ylim = c(-0.0005, .0035), xlim = c(10500, 16500), col = 'red')

#Many negative values for 95th percentile. Probably means that this is not physically realistic.
plot(y = EC_DevQLQ95, x = Dates_POBR_AllDates, ylim = c(-0.005, .01), xlim = c(10500, 16500))
par(new = TRUE)
plot(y = EC_DevQLQ95[EC_DevQLQ95 < 0], x = Dates_POBR_AllDates[EC_DevQLQ95 < 0], ylim = c(-0.005, .01), xlim = c(10500, 16500), col = 'red')

#Also mostly during the drought of 2002 and 2007. 2 are during high flows at Pond Branch
Flows_POBR_AllDates[EC_Dev05 < 0]
Dates_POBR_AllDates[EC_Dev05 < 0]
Flows_POBR_AllDates[EC_Dev < 0]
Dates_POBR_AllDates[EC_Dev < 0]
Flows_POBR_AllDates[EC_DevQLQ05 < 0]
Dates_POBR_AllDates[EC_DevQLQ05 < 0]
Flows_POBR_AllDates[EC_DevQLQ < 0]
Dates_POBR_AllDates[EC_DevQLQ < 0]

#These are basically every low flow day < 0.5 cfs. Some higher flow days, too.
Flows_POBR_AllDates[EC_Dev95 < 0]
Dates_POBR_AllDates[EC_Dev95 < 0]
Flows_POBR_AllDates[EC_DevQLQ95 < 0]
Dates_POBR_AllDates[EC_DevQLQ95 < 0]

#Fixme: for now setting all of these dates equal to 0 and not adjusting POBR
EC_Dev05[EC_Dev05 < 0] = 0
EC_Dev[EC_Dev < 0] = 0
EC_Dev95[EC_Dev95 < 0] = 0
EC_DevQLQ05[EC_DevQLQ05 < 0] = 0
EC_DevQLQ[EC_DevQLQ < 0] = 0
EC_DevQLQ95[EC_DevQLQ95 < 0] = 0

#     Method 2: Using simulation to get quantiles and mean value----
set.seed(3220)
EC_UndevMat = EC_DevMat = EC_UndevMatQLQ = EC_DevMatQLQ = matrix(NA, nrow = length(Dates_POBR_AllDates), ncol = MCreps)
for (i in 1:nrow(EC_UndevMat)){
  #Get pond branch error for this day
  if ((POBR_AllDates_PredTN$Med[i] == 0) & (POBR_AllDates_PredTN$`05`[i] == 0) & (POBR_AllDates_PredTN$`95`[i] == 0)){
    EC_UndevMat[i, ] = 0
  }else{
    e_temp = abs(log(POBR_AllDates_PredTN$Med[i]) - log(POBR_AllDates_PredTN$`05`[i]))*.5
    #Assume this is a standard deviation from a normal distribution of error and simulate possible values, and convert to a load normalized by area
    EC_UndevMat[i, ] = exp(rnorm(n = MCreps, mean = log(POBR_AllDates_PredTN$Med[i]), sd = e_temp))*1000*Flows_POBR_AllDates[i]*.3048^3/Area_K2S
  }
  
  #Now get the developed portion.
  if ((BARN_PredTN$Med[i] == 0) & (BARN_PredTN$`05`[i] == 0) & (BARN_PredTN$`95`[i] == 0)){
    EC_DevMat[i, ] =  -1*EC_UndevMat[i, ]*sum(Area.Hills[,2] - Area.Hills[,3])/sum(Area.Hills[,3])
    
    #Fixme: Not feasible to use rejection sampling here. Should some other method be used?
  }else{
    e_temp_bais = abs(log(BARN_PredTN$Med[i]) - log(BARN_PredTN$`05`[i]))*.5
    EC_DevMat[i, ] = (exp(rnorm(n = MCreps, mean = log(BARN_PredTN$Med[i]), sd = e_temp_bais))*1000*Flows_BARN[i]*.3048^3 - EC_UndevMat[i, ]*sum(Area.Hills[,2] - Area.Hills[,3]))/sum(Area.Hills[,3])
    
    #Rejection sampling for the export coefficients that are less than 0. These are physically impossible.
    while(any(EC_DevMat[i,] < 0)){
      IndL0 = which(EC_DevMat[i,] < 0)
      EC_UndevMat[i, IndL0] = exp(rnorm(n = length(IndL0), mean = log(POBR_AllDates_PredTN$Med[i]), sd = e_temp))*1000*Flows_POBR_AllDates[i]*.3048^3/Area_K2S
      EC_DevMat[i, IndL0] = (exp(rnorm(n = length(IndL0), mean = log(BARN_PredTN$Med[i]), sd = e_temp_bais))*1000*Flows_BARN[i]*.3048^3 - EC_UndevMat[i, IndL0]*sum(Area.Hills[,2] - Area.Hills[,3]))/sum(Area.Hills[,3])
    }
  }
}
rm(e_temp, e_temp_bais, i, IndL0)

for (i in 1:nrow(EC_UndevMatQLQ)){
  #Get pond branch error for this day
  if ((POBR_AllDates_PredTN$Med[i] == 0) & (POBR_AllDates_PredTN$`05`[i] == 0) & (POBR_AllDates_PredTN$`95`[i] == 0)){
    EC_UndevMatQLQ[i, ] = 0
  }else{
    e_temp = abs(log(POBR_AllDates_PredTN$Med[i]) - log(POBR_AllDates_PredTN$`05`[i]))*.5
    #Assume this is a standard deviation from a normal distribution of error and simulate possible values, and convert to a load normalized by area
    EC_UndevMatQLQ[i, ] = exp(rnorm(n = MCreps, mean = log(POBR_AllDates_PredTN$Med[i]), sd = e_temp))*1000*Flows_POBR_AllDates[i]*.3048^3/Area_K2S
  }
  
  #Now get the developed portion.
  if ((BARN_PredTN$MedQLQ[i] == 0) & (BARN_PredTN$`05QLQ`[i] == 0) & (BARN_PredTN$`95QLQ`[i] == 0)){
    EC_DevMatQLQ[i, ] =  -1*EC_UndevMatQLQ[i, ]*sum(Area.Hills[,2] - Area.Hills[,3])/sum(Area.Hills[,3])
    
    #Fixme: Not feasible to use rejection sampling here. Should some other method be used?
  }else{
    e_temp_bais = abs(log(BARN_PredTN$MedQLQ[i]) - log(BARN_PredTN$`05QLQ`[i]))*.5
    EC_DevMatQLQ[i, ] = (exp(rnorm(n = MCreps, mean = log(BARN_PredTN$MedQLQ[i]), sd = e_temp_bais))*1000*Flows_BARN[i]*.3048^3 - EC_UndevMatQLQ[i, ]*sum(Area.Hills[,2] - Area.Hills[,3]))/sum(Area.Hills[,3])
    
    #Rejection sampling for the export coefficients that are less than 0. These are physically impossible.
    while(any(EC_DevMatQLQ[i,] < 0)){
      IndL0 = which(EC_DevMatQLQ[i,] < 0)
      EC_UndevMatQLQ[i, IndL0] = exp(rnorm(n = length(IndL0), mean = log(POBR_AllDates_PredTN$Med[i]), sd = e_temp))*1000*Flows_POBR_AllDates[i]*.3048^3/Area_K2S
      EC_DevMatQLQ[i, IndL0] = (exp(rnorm(n = length(IndL0), mean = log(BARN_PredTN$MedQLQ[i]), sd = e_temp_bais))*1000*Flows_BARN[i]*.3048^3 - EC_UndevMatQLQ[i, IndL0]*sum(Area.Hills[,2] - Area.Hills[,3]))/sum(Area.Hills[,3])
    }
  }
}
rm(e_temp, e_temp_bais, i, IndL0)

#    Pond Branch Load estimation for both methods----
#Method 1
Load_ECM_Pond05 = EC_Undev05*Area_K2S + EC_Dev05*AreaImp_K2S
Load_ECM_Pond = EC_Undev*Area_K2S + EC_Dev*AreaImp_K2S
Load_ECM_Pond95 = EC_Undev95*Area_K2S + EC_Dev95*AreaImp_K2S
Load_ECM_Pond_QLQ05 = EC_Undev05*Area_K2S + EC_DevQLQ05*AreaImp_K2S
Load_ECM_Pond_QLQ = EC_Undev*Area_K2S + EC_DevQLQ*AreaImp_K2S
Load_ECM_Pond_QLQ95 = EC_Undev95*Area_K2S + EC_DevQLQ95*AreaImp_K2S

#Method 2
Load_ECM_Pond_Mat = EC_UndevMat*Area_K2S + EC_DevMat*AreaImp_K2S
#Get mean and 5th, 95th percentiles
Load_ECM_Pond_MatMean = apply(Load_ECM_Pond_Mat, MARGIN = 1, FUN = mean)
Load_ECM_Pond_Mat05 = apply(Load_ECM_Pond_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_Pond_Mat95 = apply(Load_ECM_Pond_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)

Load_ECM_Pond_MatQLQ = EC_UndevMatQLQ*Area_K2S + EC_DevMatQLQ*AreaImp_K2S
#Get mean and 5th, 95th percentiles
Load_ECM_Pond_MatQLQMean = apply(Load_ECM_Pond_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_Pond_MatQLQ05 = apply(Load_ECM_Pond_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_Pond_MatQLQ95 = apply(Load_ECM_Pond_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

#    Look at fit of model to Pond Branch----
plot(as.Date(Dates_POBR_AllDates), POBR_AllDates_PredTN$MedLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,50))
par(new = TRUE)
plot(as.Date(Dates_POBR_AllDates), Load_ECM_Pond, col = 'red', xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,50))

plot(as.Date(Dates_POBR_AllDates), POBR_AllDates_PredTN$MedLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,50))
par(new = TRUE)
plot(as.Date(Dates_POBR_AllDates), Load_ECM_Pond, col = 'red', xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,50))

#Check on method - essentially all 0, but correlation says != 1
if(round(cor(Load_ECM_Pond, POBR_AllDates_PredTN$MedLoad),10) != 1){
  print('method producing correlation != 1')
}
if(round(cor(Load_ECM_Pond_QLQ, POBR_AllDates_PredTN$MedLoad),10) != 1){
  print('method producing correlation != 1')
}

plot(as.Date(Dates_POBR_AllDates), POBR_AllDates_PredTN$TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,50))
par(new = TRUE)
plot(as.Date(Dates_POBR_AllDates), Load_ECM_Pond, col = 'red', xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,50))

#in kg/day
plot(as.Date(Dates_POBR_AllDates), POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,5))
par(new = TRUE)
plot(as.Date(Dates_POBR_AllDates), Load_ECM_Pond/1000000*3600*24, col = 'red', xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,5))

plot(POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, Load_ECM_Pond/1000000*3600*24,
     ylim = c(0.000001,10), xlim = c(0.000001,10), 
     xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Pond Branch Outlet', log = 'xy')
arrows(POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Pond05)/1000000*3600*24, POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Pond95)/1000000*3600*24, length=0.05, angle=90, code=3)
lines(c(0.000001,100), c(0.000001,100), col = 'red')

plot(POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, Load_ECM_Pond_QLQ/1000000*3600*24,
     ylim = c(0.000001,10), xlim = c(0.000001,10), 
     xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Pond Branch Outlet', log = 'xy')
arrows(POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, Load_ECM_Pond_QLQ05/1000000*3600*24, POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, Load_ECM_Pond_QLQ95/1000000*3600*24, length=0.05, angle=90, code=3)
lines(c(0.000001,100), c(0.000001,100), col = 'red')

#With correct error distribution
plot(POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, Load_ECM_Pond_MatMean/1000000*3600*24,
     ylim = c(0.000001,10), xlim = c(0.000001,10), 
     xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Pond Branch Outlet', log = 'xy')
arrows(POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Pond_Mat05)/1000000*3600*24, POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Pond_Mat95)/1000000*3600*24, length=0.05, angle=90, code=3)
lines(c(0.000001,100), c(0.000001,100), col = 'red')

#With correct error distribution
plot(POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, Load_ECM_Pond_MatQLQMean/1000000*3600*24,
     ylim = c(0.000001,10), xlim = c(0.000001,10), 
     xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Pond Branch Outlet', log = 'xy')
arrows(POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Pond_MatQLQ05)/1000000*3600*24, POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Pond_MatQLQ95)/1000000*3600*24, length=0.05, angle=90, code=3)
lines(c(0.000001,100), c(0.000001,100), col = 'red')

#    Baisman Outlet Load estimation for both methods----
Load_ECM_Baisman05 = EC_Undev05*(Area.basin - Area.basin_imp) + EC_Dev05*(Area.basin_imp)
Load_ECM_Baisman = EC_Undev*(Area.basin - Area.basin_imp) + EC_Dev*(Area.basin_imp)
Load_ECM_Baisman95 = EC_Undev95*(Area.basin - Area.basin_imp) + EC_Dev95*(Area.basin_imp)
Load_ECM_Baisman_QLQ05 = EC_Undev05*(Area.basin - Area.basin_imp) + EC_DevQLQ05*(Area.basin_imp)
Load_ECM_Baisman_QLQ = EC_Undev*(Area.basin - Area.basin_imp) + EC_DevQLQ*(Area.basin_imp)
Load_ECM_Baisman_QLQ95 = EC_Undev95*(Area.basin - Area.basin_imp) + EC_DevQLQ95*(Area.basin_imp)

Load_ECM_Baisman_Mat = EC_UndevMat*(Area.basin - Area.basin_imp) + EC_DevMat*Area.basin_imp
Load_ECM_Baisman_MatMean = apply(Load_ECM_Baisman_Mat, MARGIN = 1, FUN = mean)
Load_ECM_Baisman_Mat05 = apply(Load_ECM_Baisman_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_Baisman_Mat95 = apply(Load_ECM_Baisman_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)

Load_ECM_Baisman_MatQLQ = EC_UndevMatQLQ*(Area.basin - Area.basin_imp) + EC_DevMatQLQ*Area.basin_imp
Load_ECM_Baisman_MatQLQMean = apply(Load_ECM_Baisman_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_Baisman_MatQLQ05 = apply(Load_ECM_Baisman_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_Baisman_MatQLQ95 = apply(Load_ECM_Baisman_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

#    Look at fit of model to Baisman----
plot(as.Date(Dates_BARN), BARN_PredTN$TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,1100))
par(new = TRUE)
plot(as.Date(Dates_BARN), Load_ECM_Baisman, col = 'red', xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,1100))

plot(as.Date(Dates_BARN), BARN_PredTN$TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,1100))
par(new = TRUE)
plot(as.Date(Dates_BARN), Load_ECM_Baisman_QLQ, col = 'red', xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,1100))

plot(as.Date(Dates_BARN), BARN_PredTN$MedLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,1100))
par(new = TRUE)
plot(as.Date(Dates_BARN), Load_ECM_Baisman, col = 'red', xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,1100))

plot(BARN_PredTN$MedLoad, Load_ECM_Baisman,
     ylim = c(0,1100), xlim = c(0,1100), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Baisman Outlet')
lines(c(0,1.3), c(0,1.3))

plot(BARN_PredTN$MedLoadQLQ, Load_ECM_Baisman_QLQ,
     ylim = c(0,1100), xlim = c(0,1100), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Baisman Outlet')
lines(c(0,1.3), c(0,1.3))

#Check on method - some sites set to 0, so this will fail.
if(cor(round(Load_ECM_Baisman,7), round(BARN_PredTN$MedLoad,7)) != 1){
  print('method producing correlation != 1')
}
#Load in kg/d
plot(as.Date(Dates_BARN), BARN_PredTN$TrueLoad/1000000*3600*24, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,100))
par(new = TRUE)
plot(as.Date(Dates_BARN), Load_ECM_Baisman/1000000*3600*24, col = 'red', xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,100))

plot(BARN_PredTN$TrueLoad/1000000*3600*24, Load_ECM_Baisman/1000000*3600*24,
     ylim = c(0,100), xlim = c(0,100), xlab = 'True TN Load (kg N/d)', ylab = 'ECM Predicted TN Load (kg N/d)', main = 'Baisman Outlet')
lines(c(0,90), c(0,90), col = 'red')

#With correct error distribution
plot(BARN_PredTN$TrueLoad/1000000*3600*24, Load_ECM_Baisman_MatMean/1000000*3600*24,
     ylim = c(0,100), xlim = c(0,100), 
     xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Baisman Outlet')
arrows(BARN_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Baisman_Mat05)/1000000*3600*24, BARN_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Baisman_Mat95)/1000000*3600*24, length=0.05, angle=90, code=3)
lines(c(0.000001,100), c(0.000001,100), col = 'red')

#Load in kg/d
plot(as.Date(Dates_BARN), BARN_PredTN$TrueLoad/1000000*3600*24, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,100))
par(new = TRUE)
plot(as.Date(Dates_BARN), Load_ECM_Baisman_QLQ/1000000*3600*24, col = 'red', xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,100))

plot(BARN_PredTN$TrueLoad/1000000*3600*24, Load_ECM_Baisman_QLQ/1000000*3600*24, log='xy',
     ylim = c(0.001,100), xlim = c(0.001,100), xlab = 'True TN Load (kg N/d)', ylab = 'ECM Predicted TN Load (kg N/d)', main = 'Baisman Outlet')
arrows(BARN_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Baisman_QLQ05)/1000000*3600*24, BARN_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Baisman_QLQ95)/1000000*3600*24, length=0.05, angle=90, code=3)
lines(c(0.001,90), c(0.001,90), col = 'red')

#With correct error distribution
plot(BARN_PredTN$TrueLoad/1000000*3600*24, Load_ECM_Baisman_MatQLQMean/1000000*3600*24,
     ylim = c(0,100), xlim = c(0,100), 
     xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Baisman Outlet')
arrows(BARN_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Baisman_MatQLQ05)/1000000*3600*24, BARN_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Baisman_MatQLQ95)/1000000*3600*24, length=0.05, angle=90, code=3)
lines(c(0.000001,100), c(0.000001,100), col = 'red')

#    Pond and Baisman on same scale----
setwd(wd_BESN)
png('ECMload_BARN_POBR.png', res = 300, height = 5, width = 10, units = 'in')
layout(rbind(c(1,2)))
plot(POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, Load_ECM_Pond_MatMean/1000000*3600*24,
     ylim = c(0.00001,1000), xlim = c(0.00001,1000), 
     xlab = 'True TN Load (kg/d)', ylab = 'ECM Predicted TN Load (kg/d)', main = 'Pond Branch Outlet', log = 'xy')
arrows(POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Pond_Mat05)/1000000*3600*24, POBR_AllDates_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Pond_Mat95)/1000000*3600*24, length=0.05, angle=90, code=3)
lines(c(0.000001,100), c(0.000001,100), col = 'red')

plot(BARN_PredTN$TrueLoad/1000000*3600*24, Load_ECM_Baisman_MatQLQMean/1000000*3600*24, log = 'xy',
     ylim = c(0.00001,1000), xlim = c(0.00001,1000), 
     xlab = 'True TN Load (kg/d)', ylab = 'ECM Predicted TN Load (kg/d)', main = 'Baisman Outlet')
arrows(BARN_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Baisman_MatQLQ05)/1000000*3600*24, BARN_PredTN$TrueLoad/1000000*3600*24, (Load_ECM_Baisman_MatQLQ95)/1000000*3600*24, length=0.05, angle=90, code=3)
lines(c(0.000001,100), c(0.000001,100), col = 'red')
dev.off()

#cor(BARN_PredTN$TrueLoad, Load_ECM_Baisman)

#  Compare models to where there are datasets----
#   BR3 as collected from above analysis----
Load_ECM_BR305 = EC_Undev05*sum(Area.Hills[c(13,14),2]-Area.Hills[c(13,14),3]) + EC_Dev05*sum(Area.Hills[c(13,14),3])
Load_ECM_BR3 = EC_Undev*sum(Area.Hills[c(13,14),2]-Area.Hills[c(13,14),3]) + EC_Dev*sum(Area.Hills[c(13,14),3])
Load_ECM_BR395 = EC_Undev95*sum(Area.Hills[c(13,14),2]-Area.Hills[c(13,14),3]) + EC_Dev95*sum(Area.Hills[c(13,14),3])
Load_ECM_BR3_QLQ05 = EC_Undev05*sum(Area.Hills[c(13,14),2]-Area.Hills[c(13,14),3]) + EC_DevQLQ05*sum(Area.Hills[c(13,14),3])
Load_ECM_BR3_QLQ = EC_Undev*sum(Area.Hills[c(13,14),2]-Area.Hills[c(13,14),3]) + EC_DevQLQ*sum(Area.Hills[c(13,14),3])
Load_ECM_BR3_QLQ95 = EC_Undev95*sum(Area.Hills[c(13,14),2]-Area.Hills[c(13,14),3]) + EC_DevQLQ95*sum(Area.Hills[c(13,14),3])

Load_ECM_BR3_Mat = EC_UndevMat*sum(Area.Hills[c(13,14),2]-Area.Hills[c(13,14),3]) + EC_DevMat*sum(Area.Hills[c(13,14),3])
Load_ECM_BR3_MatMean = apply(Load_ECM_BR3_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR3_Mat05 = apply(Load_ECM_BR3_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3_Mat95 = apply(Load_ECM_BR3_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)

Load_ECM_BR3_MatQLQ = EC_UndevMatQLQ*sum(Area.Hills[c(13,14),2]-Area.Hills[c(13,14),3]) + EC_DevMatQLQ*sum(Area.Hills[c(13,14),3])
Load_ECM_BR3_MatQLQMean = apply(Load_ECM_BR3_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR3_MatQLQ05 = apply(Load_ECM_BR3_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3_MatQLQ95 = apply(Load_ECM_BR3_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BR3, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR3), BR3_PredTN$TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BR3_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR3), BR3_PredTN$TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

plot(BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)], Load_ECM_BR3[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)], 
     xlim = c(0,200), ylim = c(0,200), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Hillslope 3')
arrows(BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)], (Load_ECM_BR305[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)]), BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)], (Load_ECM_BR395[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')

#With correct error bars
plot(BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)], Load_ECM_BR3_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)], 
     xlim = c(0,200), ylim = c(0,200), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Hillslope 3')
arrows(BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)], (Load_ECM_BR3_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)]), BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)], (Load_ECM_BR3_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')

plot(BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)], Load_ECM_BR3_QLQ[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)], 
     xlim = c(0,200), ylim = c(0,200), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Hillslope 3')
arrows(BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)], (Load_ECM_BR3_QLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)]), BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)], (Load_ECM_BR3_QLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')

#With correct error bars
plot(BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)], Load_ECM_BR3_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)], 
     xlim = c(0,200), ylim = c(0,200), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Hillslope 3')
arrows(BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)], (Load_ECM_BR3_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)]), BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)], (Load_ECM_BR3_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')

plot(BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)], Load_ECM_BR3[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)], 
     xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Hillslope 3', xlim = c(0,100), ylim = c(0,100))
lines(c(0,250), c(0,250))

plot(BR3_PredTN$TrueLoad[as.Date(Dates_BR3) %in% as.Date(Dates_BARN)], Load_ECM_BR3_QLQ[as.Date(Dates_BARN) %in% as.Date(Dates_BR3)], 
     xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'Hillslope 3', xlim = c(0,100), ylim = c(0,100))
lines(c(0,100), c(0,100))

#   Other Smith Sampling Gauges----
#    BR3 S1----
#Get the impervious surface area contributing to this location
AreaImp_S1 = 0
Area_S1 = length(which(world[which(duplicated(world$patchID) == FALSE),]$G_10 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$G_10 == 1))){
  AreaImp_S1 = AreaImp_S1 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$G_10 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
BR3S1_TN = S_TN[S_TN$Site == S_Sites@data$NAME[11],]
BR3S1_Q = S_Q[S_Q$Site == S_Sites@data$NAME[11],]
#Check for NAs and remove
IndRm = which((is.na(BR3S1_TN)) | (is.na(BR3S1_Q)))
BR3S1_Q = BR3S1_Q[-IndRm]
BR3S1_TN = BR3S1_TN[-IndRm]
rm(IndRm)
Dates_BR3S1 = vector('character', ncol(BR3S1_Q)-1)
for (i in 1:(ncol(BR3S1_Q)-1)){
  temp = strsplit(colnames(BR3S1_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_BR3S1[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
rm(temp)
Dates_BR3S1 = as.Date(Dates_BR3S1)

#mg/s unit
BR3S1_TrueLoad = as.numeric(BR3S1_TN[-1]*1000*BR3S1_Q[-1])

Load_ECM_BR3S105 = EC_Undev05*(Area_S1 - AreaImp_S1) + EC_Dev05*AreaImp_S1
Load_ECM_BR3S1 = EC_Undev*(Area_S1 - AreaImp_S1) + EC_Dev*AreaImp_S1
Load_ECM_BR3S195 = EC_Undev95*(Area_S1 - AreaImp_S1) + EC_Dev95*AreaImp_S1
Load_ECM_BR3S1_QLQ05 = EC_Undev05*(Area_S1 - AreaImp_S1) + EC_DevQLQ05*AreaImp_S1
Load_ECM_BR3S1_QLQ = EC_Undev*(Area_S1 - AreaImp_S1) + EC_DevQLQ*AreaImp_S1
Load_ECM_BR3S1_QLQ95 = EC_Undev95*(Area_S1 - AreaImp_S1) + EC_DevQLQ95*AreaImp_S1

Load_ECM_BR3S1_Mat = EC_UndevMat*(Area_S1 - AreaImp_S1) + EC_DevMat*AreaImp_S1
Load_ECM_BR3S1_MatMean = apply(Load_ECM_BR3S1_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR3S1_Mat05 = apply(Load_ECM_BR3S1_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S1_Mat95 = apply(Load_ECM_BR3S1_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR3S1_MatQLQ = EC_UndevMatQLQ*(Area_S1 - AreaImp_S1) + EC_DevMatQLQ*AreaImp_S1
Load_ECM_BR3S1_MatQLQMean = apply(Load_ECM_BR3S1_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR3S1_MatQLQ05 = apply(Load_ECM_BR3S1_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S1_MatQLQ95 = apply(Load_ECM_BR3S1_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BR3S1, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR3S1), BR3S1_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BR3S1_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR3S1), BR3S1_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

png('ECMBR3S1.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR3S1_TrueLoad[as.Date(Dates_BR3S1) %in% as.Date(Dates_BARN)], y = Load_ECM_BR3S1_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S1)], 
     xlim = c(0,30), ylim = c(0,30), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S1')
arrows(BR3S1_TrueLoad[as.Date(Dates_BR3S1) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S1_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S1)]), BR3S1_TrueLoad[as.Date(Dates_BR3S1) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S1_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S1)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR3S1_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR3S1_TrueLoad[as.Date(Dates_BR3S1) %in% as.Date(Dates_BARN)], Load_ECM_BR3S1_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S1)], 
     xlim = c(0,30), ylim = c(0,30), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S1')
arrows(BR3S1_TrueLoad[as.Date(Dates_BR3S1) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S1_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S1)]), BR3S1_TrueLoad[as.Date(Dates_BR3S1) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S1_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S1)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    BR3 S2----
#Get the impervious surface area contributing to this location
AreaImp_S2 = 0
Area_S2 = length(which(world[which(duplicated(world$patchID) == FALSE),]$G_0 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$G_0 == 1))){
  AreaImp_S2 = AreaImp_S2 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$G_0 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
BR3S2_TN = S_TN[S_TN$Site == S_Sites@data$NAME[1],]
BR3S2_Q = S_Q[S_Q$Site == S_Sites@data$NAME[1],]
#Check for NAs and remove
IndRm = which((is.na(BR3S2_TN)) | (is.na(BR3S2_Q)))
BR3S2_Q = BR3S2_Q[-IndRm]
BR3S2_TN = BR3S2_TN[-IndRm]
rm(IndRm)
Dates_BR3S2 = vector('character', ncol(BR3S2_Q)-1)
for (i in 1:(ncol(BR3S2_Q)-1)){
  temp = strsplit(colnames(BR3S2_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_BR3S2[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
rm(temp)
Dates_BR3S2 = as.Date(Dates_BR3S2)

#mg/s unit
BR3S2_TrueLoad = as.numeric(BR3S2_TN[-1]*1000*BR3S2_Q[-1])

Load_ECM_BR3S205 = EC_Undev05*(Area_S2 - AreaImp_S2) + EC_Dev05*AreaImp_S2
Load_ECM_BR3S2 = EC_Undev*(Area_S2 - AreaImp_S2) + EC_Dev*AreaImp_S2
Load_ECM_BR3S295 = EC_Undev95*(Area_S2 - AreaImp_S2) + EC_Dev95*AreaImp_S2
Load_ECM_BR3S2_QLQ05 = EC_Undev05*(Area_S2 - AreaImp_S2) + EC_DevQLQ05*AreaImp_S2
Load_ECM_BR3S2_QLQ = EC_Undev*(Area_S2 - AreaImp_S2) + EC_DevQLQ*AreaImp_S2
Load_ECM_BR3S2_QLQ95 = EC_Undev95*(Area_S2 - AreaImp_S2) + EC_DevQLQ95*AreaImp_S2

Load_ECM_BR3S2_Mat = EC_UndevMat*(Area_S2 - AreaImp_S2) + EC_DevMat*AreaImp_S2
Load_ECM_BR3S2_MatMean = apply(Load_ECM_BR3S2_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR3S2_Mat05 = apply(Load_ECM_BR3S2_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S2_Mat95 = apply(Load_ECM_BR3S2_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR3S2_MatQLQ = EC_UndevMatQLQ*(Area_S2 - AreaImp_S2) + EC_DevMatQLQ*AreaImp_S2
Load_ECM_BR3S2_MatQLQMean = apply(Load_ECM_BR3S2_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR3S2_MatQLQ05 = apply(Load_ECM_BR3S2_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S2_MatQLQ95 = apply(Load_ECM_BR3S2_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BR3S2, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR3S2), BR3S2_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BR3S2_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR3S2), BR3S2_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

png('ECMBR3S2.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR3S2_TrueLoad[as.Date(Dates_BR3S2) %in% as.Date(Dates_BARN)], y = Load_ECM_BR3S2_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S2)], 
     xlim = c(0,40), ylim = c(0,40), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S2')
arrows(BR3S2_TrueLoad[as.Date(Dates_BR3S2) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S2_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S2)]), BR3S2_TrueLoad[as.Date(Dates_BR3S2) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S2_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S2)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR3S2_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR3S2_TrueLoad[as.Date(Dates_BR3S2) %in% as.Date(Dates_BARN)], Load_ECM_BR3S2_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S2)], 
     xlim = c(0,40), ylim = c(0,40), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S2')
arrows(BR3S2_TrueLoad[as.Date(Dates_BR3S2) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S2_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S2)]), BR3S2_TrueLoad[as.Date(Dates_BR3S2) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S2_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S2)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    BR3 S3----
#Get the impervious surface area contributing to this location
AreaImp_S3 = 0
Area_S3 = length(which(world[which(duplicated(world$patchID) == FALSE),]$G_1 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$G_1 == 1))){
  AreaImp_S3 = AreaImp_S3 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$G_1 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
BR3S3_TN = S_TN[S_TN$Site == S_Sites@data$NAME[2],]
BR3S3_Q = S_Q[S_Q$Site == S_Sites@data$NAME[2],]
#Check for NAs and remove
IndRm = which((is.na(BR3S3_TN)) | (is.na(BR3S3_Q)))
BR3S3_Q = BR3S3_Q[-IndRm]
BR3S3_TN = BR3S3_TN[-IndRm]
rm(IndRm)
Dates_BR3S3 = vector('character', ncol(BR3S3_Q)-1)
for (i in 1:(ncol(BR3S3_Q)-1)){
  temp = strsplit(colnames(BR3S3_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_BR3S3[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
rm(temp)
Dates_BR3S3 = as.Date(Dates_BR3S3)

#mg/s unit
BR3S3_TrueLoad = as.numeric(BR3S3_TN[-1]*1000*BR3S3_Q[-1])

Load_ECM_BR3S305 = EC_Undev05*(Area_S3 - AreaImp_S3) + EC_Dev05*AreaImp_S3
Load_ECM_BR3S3 = EC_Undev*(Area_S3 - AreaImp_S3) + EC_Dev*AreaImp_S3
Load_ECM_BR3S395 = EC_Undev95*(Area_S3 - AreaImp_S3) + EC_Dev95*AreaImp_S3
Load_ECM_BR3S3_QLQ05 = EC_Undev05*(Area_S3 - AreaImp_S3) + EC_DevQLQ05*AreaImp_S3
Load_ECM_BR3S3_QLQ = EC_Undev*(Area_S3 - AreaImp_S3) + EC_DevQLQ*AreaImp_S3
Load_ECM_BR3S3_QLQ95 = EC_Undev95*(Area_S3 - AreaImp_S3) + EC_DevQLQ95*AreaImp_S3

Load_ECM_BR3S3_Mat = EC_UndevMat*(Area_S3 - AreaImp_S3) + EC_DevMat*AreaImp_S3
Load_ECM_BR3S3_MatMean = apply(Load_ECM_BR3S3_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR3S3_Mat05 = apply(Load_ECM_BR3S3_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S3_Mat95 = apply(Load_ECM_BR3S3_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR3S3_MatQLQ = EC_UndevMatQLQ*(Area_S3 - AreaImp_S3) + EC_DevMatQLQ*AreaImp_S3
Load_ECM_BR3S3_MatQLQMean = apply(Load_ECM_BR3S3_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR3S3_MatQLQ05 = apply(Load_ECM_BR3S3_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S3_MatQLQ95 = apply(Load_ECM_BR3S3_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BR3S3, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR3S3), BR3S3_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BR3S3_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR3S3), BR3S3_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

png('ECMBR3S3.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR3S3_TrueLoad[as.Date(Dates_BR3S3) %in% as.Date(Dates_BARN)], y = Load_ECM_BR3S3_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S3)], 
     xlim = c(0,30), ylim = c(0,30), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S3')
arrows(BR3S3_TrueLoad[as.Date(Dates_BR3S3) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S3_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S3)]), BR3S3_TrueLoad[as.Date(Dates_BR3S3) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S3_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S3)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR3S3_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR3S3_TrueLoad[as.Date(Dates_BR3S3) %in% as.Date(Dates_BARN)], Load_ECM_BR3S3_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S3)], 
     xlim = c(0,30), ylim = c(0,30), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S3')
arrows(BR3S3_TrueLoad[as.Date(Dates_BR3S3) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S3_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S3)]), BR3S3_TrueLoad[as.Date(Dates_BR3S3) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S3_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S3)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    BR3 S4----
#Get the impervious surface area contributing to this location
AreaImp_S4 = 0
Area_S4 = length(which(world[which(duplicated(world$patchID) == FALSE),]$G_2 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$G_2 == 1))){
  AreaImp_S4 = AreaImp_S4 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$G_2 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
BR3S4_TN = S_TN[S_TN$Site == S_Sites@data$NAME[3],]
BR3S4_Q = S_Q[S_Q$Site == S_Sites@data$NAME[3],]
#Check for NAs and remove
IndRm = which((is.na(BR3S4_TN)) | (is.na(BR3S4_Q)))
BR3S4_Q = BR3S4_Q[-IndRm]
BR3S4_TN = BR3S4_TN[-IndRm]
rm(IndRm)
Dates_BR3S4 = vector('character', ncol(BR3S4_Q)-1)
for (i in 1:(ncol(BR3S4_Q)-1)){
  temp = strsplit(colnames(BR3S4_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_BR3S4[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_BR3S4 = as.Date(Dates_BR3S4)

#mg/s unit
BR3S4_TrueLoad = as.numeric(BR3S4_TN[-1]*1000*BR3S4_Q[-1])

Load_ECM_BR3S405 = EC_Undev05*(Area_S4 - AreaImp_S4) + EC_Dev05*AreaImp_S4
Load_ECM_BR3S4 = EC_Undev*(Area_S4 - AreaImp_S4) + EC_Dev*AreaImp_S4
Load_ECM_BR3S495 = EC_Undev95*(Area_S4 - AreaImp_S4) + EC_Dev95*AreaImp_S4
Load_ECM_BR3S4_QLQ05 = EC_Undev05*(Area_S4 - AreaImp_S4) + EC_DevQLQ05*AreaImp_S4
Load_ECM_BR3S4_QLQ = EC_Undev*(Area_S4 - AreaImp_S4) + EC_DevQLQ*AreaImp_S4
Load_ECM_BR3S4_QLQ95 = EC_Undev95*(Area_S4 - AreaImp_S4) + EC_DevQLQ95*AreaImp_S4

Load_ECM_BR3S4_Mat = EC_UndevMat*(Area_S4 - AreaImp_S4) + EC_DevMat*AreaImp_S4
Load_ECM_BR3S4_MatMean = apply(Load_ECM_BR3S4_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR3S4_Mat05 = apply(Load_ECM_BR3S4_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S4_Mat95 = apply(Load_ECM_BR3S4_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR3S4_MatQLQ = EC_UndevMatQLQ*(Area_S4 - AreaImp_S4) + EC_DevMatQLQ*AreaImp_S4
Load_ECM_BR3S4_MatQLQMean = apply(Load_ECM_BR3S4_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR3S4_MatQLQ05 = apply(Load_ECM_BR3S4_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S4_MatQLQ95 = apply(Load_ECM_BR3S4_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BR3S4, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR3S4), BR3S4_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BR3S4_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR3S4), BR3S4_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

png('ECMBR3S4.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR3S4_TrueLoad[as.Date(Dates_BR3S4) %in% as.Date(Dates_BARN)], y = Load_ECM_BR3S4_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S4)], 
     xlim = c(0,30), ylim = c(0,30), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S4')
arrows(BR3S4_TrueLoad[as.Date(Dates_BR3S4) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S4_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S4)]), BR3S4_TrueLoad[as.Date(Dates_BR3S4) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S4_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S4)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR3S4_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR3S4_TrueLoad[as.Date(Dates_BR3S4) %in% as.Date(Dates_BARN)], Load_ECM_BR3S4_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S4)], 
     xlim = c(0,30), ylim = c(0,30), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S4')
arrows(BR3S4_TrueLoad[as.Date(Dates_BR3S4) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S4_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S4)]), BR3S4_TrueLoad[as.Date(Dates_BR3S4) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S4_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S4)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    BR3 S5 - poor fit - located in headwaters, might be Arc problem----
#Get the impervious surface area contributing to this location
AreaImp_S5 = 0
Area_S5 = length(which(world[which(duplicated(world$patchID) == FALSE),]$G_3 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$G_3 == 1))){
  AreaImp_S5 = AreaImp_S5 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$G_3 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
BR3S5_TN = S_TN[S_TN$Site == S_Sites@data$NAME[4],]
BR3S5_Q = S_Q[S_Q$Site == S_Sites@data$NAME[4],]
#Check for NAs and remove
IndRm = which((is.na(BR3S5_TN)) | (is.na(BR3S5_Q)))
BR3S5_Q = BR3S5_Q[-IndRm]
BR3S5_TN = BR3S5_TN[-IndRm]
rm(IndRm)
Dates_BR3S5 = vector('character', ncol(BR3S5_Q)-1)
for (i in 1:(ncol(BR3S5_Q)-1)){
  temp = strsplit(colnames(BR3S5_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_BR3S5[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_BR3S5 = as.Date(Dates_BR3S5)

#mg/s unit
BR3S5_TrueLoad = as.numeric(BR3S5_TN[-1]*1000*BR3S5_Q[-1])

Load_ECM_BR3S505 = EC_Undev05*(Area_S5 - AreaImp_S5) + EC_Dev05*AreaImp_S5
Load_ECM_BR3S5 = EC_Undev*(Area_S5 - AreaImp_S5) + EC_Dev*AreaImp_S5
Load_ECM_BR3S595 = EC_Undev95*(Area_S5 - AreaImp_S5) + EC_Dev95*AreaImp_S5
Load_ECM_BR3S5_QLQ05 = EC_Undev05*(Area_S5 - AreaImp_S5) + EC_DevQLQ05*AreaImp_S5
Load_ECM_BR3S5_QLQ = EC_Undev*(Area_S5 - AreaImp_S5) + EC_DevQLQ*AreaImp_S5
Load_ECM_BR3S5_QLQ95 = EC_Undev95*(Area_S5 - AreaImp_S5) + EC_DevQLQ95*AreaImp_S5

Load_ECM_BR3S5_Mat = EC_UndevMat*(Area_S5 - AreaImp_S5) + EC_DevMat*AreaImp_S5
Load_ECM_BR3S5_MatMean = apply(Load_ECM_BR3S5_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR3S5_Mat05 = apply(Load_ECM_BR3S5_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S5_Mat95 = apply(Load_ECM_BR3S5_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR3S5_MatQLQ = EC_UndevMatQLQ*(Area_S5 - AreaImp_S5) + EC_DevMatQLQ*AreaImp_S5
Load_ECM_BR3S5_MatQLQMean = apply(Load_ECM_BR3S5_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR3S5_MatQLQ05 = apply(Load_ECM_BR3S5_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S5_MatQLQ95 = apply(Load_ECM_BR3S5_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BR3S5, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_BR3S5), BR3S5_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BR3S5_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_BR3S5), BR3S5_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

png('ECMBR3S5.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR3S5_TrueLoad[as.Date(Dates_BR3S5) %in% as.Date(Dates_BARN)], y = Load_ECM_BR3S5_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S5)], 
     xlim = c(0,30), ylim = c(0,30), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S5')
arrows(BR3S5_TrueLoad[as.Date(Dates_BR3S5) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S5_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S5)]), BR3S5_TrueLoad[as.Date(Dates_BR3S5) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S5_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S5)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR3S5_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR3S5_TrueLoad[as.Date(Dates_BR3S5) %in% as.Date(Dates_BARN)], Load_ECM_BR3S5_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S5)], 
     xlim = c(0,30), ylim = c(0,30), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S5')
arrows(BR3S5_TrueLoad[as.Date(Dates_BR3S5) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S5_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S5)]), BR3S5_TrueLoad[as.Date(Dates_BR3S5) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S5_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S5)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    BR3 S6 - poor fit - located in headwaters, might be Arc problem----
#Get the impervious surface area contributing to this location
AreaImp_S6 = 0
Area_S6 = length(which(world[which(duplicated(world$patchID) == FALSE),]$G_4 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$G_4 == 1))){
  AreaImp_S6 = AreaImp_S6 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$G_4 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
BR3S6_TN = S_TN[S_TN$Site == S_Sites@data$NAME[5],]
BR3S6_Q = S_Q[S_Q$Site == S_Sites@data$NAME[5],]
#Check for NAs and remove
IndRm = which((is.na(BR3S6_TN)) | (is.na(BR3S6_Q)))
BR3S6_Q = BR3S6_Q[-IndRm]
BR3S6_TN = BR3S6_TN[-IndRm]
rm(IndRm)
Dates_BR3S6 = vector('character', ncol(BR3S6_Q)-1)
for (i in 1:(ncol(BR3S6_Q)-1)){
  temp = strsplit(colnames(BR3S6_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_BR3S6[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_BR3S6 = as.Date(Dates_BR3S6)

#mg/s unit
BR3S6_TrueLoad = as.numeric(BR3S6_TN[-1]*1000*BR3S6_Q[-1])

Load_ECM_BR3S605 = EC_Undev05*(Area_S6 - AreaImp_S6) + EC_Dev05*AreaImp_S6
Load_ECM_BR3S6 = EC_Undev*(Area_S6 - AreaImp_S6) + EC_Dev*AreaImp_S6
Load_ECM_BR3S695 = EC_Undev95*(Area_S6 - AreaImp_S6) + EC_Dev95*AreaImp_S6
Load_ECM_BR3S6_QLQ05 = EC_Undev05*(Area_S6 - AreaImp_S6) + EC_DevQLQ05*AreaImp_S6
Load_ECM_BR3S6_QLQ = EC_Undev*(Area_S6 - AreaImp_S6) + EC_DevQLQ*AreaImp_S6
Load_ECM_BR3S6_QLQ95 = EC_Undev95*(Area_S6 - AreaImp_S6) + EC_DevQLQ95*AreaImp_S6

Load_ECM_BR3S6_Mat = EC_UndevMat*(Area_S6 - AreaImp_S6) + EC_DevMat*AreaImp_S6
Load_ECM_BR3S6_MatMean = apply(Load_ECM_BR3S6_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR3S6_Mat05 = apply(Load_ECM_BR3S6_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S6_Mat95 = apply(Load_ECM_BR3S6_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR3S6_MatQLQ = EC_UndevMatQLQ*(Area_S6 - AreaImp_S6) + EC_DevMatQLQ*AreaImp_S6
Load_ECM_BR3S6_MatQLQMean = apply(Load_ECM_BR3S6_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR3S6_MatQLQ05 = apply(Load_ECM_BR3S6_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S6_MatQLQ95 = apply(Load_ECM_BR3S6_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BR3S6, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR3S6), BR3S6_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BR3S6_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR3S6), BR3S6_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

png('ECMBR3S6.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR3S6_TrueLoad[as.Date(Dates_BR3S6) %in% as.Date(Dates_BARN)], y = Load_ECM_BR3S6_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S6)], 
     xlim = c(0,30), ylim = c(0,30), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S6')
arrows(BR3S6_TrueLoad[as.Date(Dates_BR3S6) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S6_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S6)]), BR3S6_TrueLoad[as.Date(Dates_BR3S6) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S6_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S6)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR3S6_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR3S6_TrueLoad[as.Date(Dates_BR3S6) %in% as.Date(Dates_BARN)], Load_ECM_BR3S6_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S6)], 
     xlim = c(0,30), ylim = c(0,30), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S6')
arrows(BR3S6_TrueLoad[as.Date(Dates_BR3S6) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S6_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S6)]), BR3S6_TrueLoad[as.Date(Dates_BR3S6) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S6_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S6)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    BR3 S7 - no data----
#    BR3 S8----
#Get the impervious surface area contributing to this location
AreaImp_S8 = 0
Area_S8 = length(which(world[which(duplicated(world$patchID) == FALSE),]$G_6 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$G_6 == 1))){
  AreaImp_S8 = AreaImp_S8 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$G_6 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
BR3S8_TN = S_TN[S_TN$Site == S_Sites@data$NAME[7],]
BR3S8_Q = S_Q[S_Q$Site == S_Sites@data$NAME[7],]
#Check for NAs and remove
IndRm = which((is.na(BR3S8_TN)) | (is.na(BR3S8_Q)))
BR3S8_Q = BR3S8_Q[-IndRm]
BR3S8_TN = BR3S8_TN[-IndRm]
rm(IndRm)
Dates_BR3S8 = vector('character', ncol(BR3S8_Q)-1)
for (i in 1:(ncol(BR3S8_Q)-1)){
  temp = strsplit(colnames(BR3S8_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_BR3S8[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_BR3S8 = as.Date(Dates_BR3S8)

#mg/s unit
BR3S8_TrueLoad = as.numeric(BR3S8_TN[-1]*1000*BR3S8_Q[-1])

Load_ECM_BR3S805 = EC_Undev05*(Area_S8 - AreaImp_S8) + EC_Dev05*AreaImp_S8
Load_ECM_BR3S8 = EC_Undev*(Area_S8 - AreaImp_S8) + EC_Dev*AreaImp_S8
Load_ECM_BR3S895 = EC_Undev95*(Area_S8 - AreaImp_S8) + EC_Dev95*AreaImp_S8
Load_ECM_BR3S8_QLQ05 = EC_Undev05*(Area_S8 - AreaImp_S8) + EC_DevQLQ05*AreaImp_S8
Load_ECM_BR3S8_QLQ = EC_Undev*(Area_S8 - AreaImp_S8) + EC_DevQLQ*AreaImp_S8
Load_ECM_BR3S8_QLQ95 = EC_Undev95*(Area_S8 - AreaImp_S8) + EC_DevQLQ95*AreaImp_S8

Load_ECM_BR3S8_Mat = EC_UndevMat*(Area_S8 - AreaImp_S8) + EC_DevMat*AreaImp_S8
Load_ECM_BR3S8_MatMean = apply(Load_ECM_BR3S8_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR3S8_Mat05 = apply(Load_ECM_BR3S8_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S8_Mat95 = apply(Load_ECM_BR3S8_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR3S8_MatQLQ = EC_UndevMatQLQ*(Area_S8 - AreaImp_S8) + EC_DevMatQLQ*AreaImp_S8
Load_ECM_BR3S8_MatQLQMean = apply(Load_ECM_BR3S8_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR3S8_MatQLQ05 = apply(Load_ECM_BR3S8_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S8_MatQLQ95 = apply(Load_ECM_BR3S8_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BR3S8, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR3S8), BR3S8_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BR3S8_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR3S8), BR3S8_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

png('ECMBR3S8.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR3S8_TrueLoad[as.Date(Dates_BR3S8) %in% as.Date(Dates_BARN)], y = Load_ECM_BR3S8_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S8)], 
     xlim = c(0,30), ylim = c(0,30), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S8')
arrows(BR3S8_TrueLoad[as.Date(Dates_BR3S8) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S8_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S8)]), BR3S8_TrueLoad[as.Date(Dates_BR3S8) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S8_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S8)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR3S8_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR3S8_TrueLoad[as.Date(Dates_BR3S8) %in% as.Date(Dates_BARN)], Load_ECM_BR3S8_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S8)], 
     xlim = c(0,30), ylim = c(0,30), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S8')
arrows(BR3S8_TrueLoad[as.Date(Dates_BR3S8) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S8_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S8)]), BR3S8_TrueLoad[as.Date(Dates_BR3S8) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S8_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S8)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    BR5a S10 - poor fit----
#Get the impervious surface area contributing to this location
AreaImp_S10 = 0
Area_S10 = length(which(world[which(duplicated(world$patchID) == FALSE),]$G_8 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$G_8 == 1))){
  AreaImp_S10 = AreaImp_S10 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$G_8 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
BR5S10_TN = S_TN[S_TN$Site == toupper(S_Sites@data$NAME[9]),]
BR5S10_Q = S_Q[S_Q$Site == toupper(S_Sites@data$NAME[9]),]
#Check for NAs and remove
IndRm = which((is.na(BR5S10_TN)) | (is.na(BR5S10_Q)))
BR5S10_Q = BR5S10_Q[-IndRm]
BR5S10_TN = BR5S10_TN[-IndRm]
rm(IndRm)
Dates_BR5S10 = vector('character', ncol(BR5S10_Q)-1)
for (i in 1:(ncol(BR5S10_Q)-1)){
  temp = strsplit(colnames(BR5S10_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_BR5S10[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_BR5S10 = as.Date(Dates_BR5S10)

#mg/s unit
BR5S10_TrueLoad = as.numeric(BR5S10_TN[-1]*1000*BR5S10_Q[-1])

Load_ECM_BR5S1005 = EC_Undev05*(Area_S10 - AreaImp_S10) + EC_Dev05*AreaImp_S10
Load_ECM_BR5S10 = EC_Undev*(Area_S10 - AreaImp_S10) + EC_Dev*AreaImp_S10
Load_ECM_BR5S1095 = EC_Undev95*(Area_S10 - AreaImp_S10) + EC_Dev95*AreaImp_S10
Load_ECM_BR5S10_QLQ05 = EC_Undev05*(Area_S10 - AreaImp_S10) + EC_DevQLQ05*AreaImp_S10
Load_ECM_BR5S10_QLQ = EC_Undev*(Area_S10 - AreaImp_S10) + EC_DevQLQ*AreaImp_S10
Load_ECM_BR5S10_QLQ95 = EC_Undev95*(Area_S10 - AreaImp_S10) + EC_DevQLQ95*AreaImp_S10

Load_ECM_BR5S10_Mat = EC_UndevMat*(Area_S10 - AreaImp_S10) + EC_DevMat*AreaImp_S10
Load_ECM_BR5S10_MatMean = apply(Load_ECM_BR5S10_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR5S10_Mat05 = apply(Load_ECM_BR5S10_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR5S10_Mat95 = apply(Load_ECM_BR5S10_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR5S10_MatQLQ = EC_UndevMatQLQ*(Area_S10 - AreaImp_S10) + EC_DevMatQLQ*AreaImp_S10
Load_ECM_BR5S10_MatQLQMean = apply(Load_ECM_BR5S10_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR5S10_MatQLQ05 = apply(Load_ECM_BR5S10_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR5S10_MatQLQ95 = apply(Load_ECM_BR5S10_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BR5S10, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR5S10), BR5S10_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BR5S10_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR5S10), BR5S10_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

png('ECMBR5aS10.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR5S10_TrueLoad[as.Date(Dates_BR5S10) %in% as.Date(Dates_BARN)], y = Load_ECM_BR5S10_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR5S10)], 
     xlim = c(0,35), ylim = c(0,35), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR5 S10')
arrows(BR5S10_TrueLoad[as.Date(Dates_BR5S10) %in% as.Date(Dates_BARN)], (Load_ECM_BR5S10_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR5S10)]), BR5S10_TrueLoad[as.Date(Dates_BR5S10) %in% as.Date(Dates_BARN)], (Load_ECM_BR5S10_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR5S10)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR5aS10_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR5S10_TrueLoad[as.Date(Dates_BR5S10) %in% as.Date(Dates_BARN)], Load_ECM_BR5S10_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR5S10)], 
     xlim = c(0,35), ylim = c(0,35), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR5 S10')
arrows(BR5S10_TrueLoad[as.Date(Dates_BR5S10) %in% as.Date(Dates_BARN)], (Load_ECM_BR5S10_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR5S10)]), BR5S10_TrueLoad[as.Date(Dates_BR5S10) %in% as.Date(Dates_BARN)], (Load_ECM_BR5S10_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR5S10)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    BR5a S11 - poor fit----
#Get the impervious surface area contributing to this location
AreaImp_S11 = 0
Area_S11 = length(which(world[which(duplicated(world$patchID) == FALSE),]$G_9 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$G_9 == 1))){
  AreaImp_S11 = AreaImp_S11 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$G_9 == 1)[i]]*res^2
}
rm(i)

#Observed data for this site
BR5S11_TN = S_TN[S_TN$Site == toupper(S_Sites@data$NAME[10]),]
BR5S11_Q = S_Q[S_Q$Site == toupper(S_Sites@data$NAME[10]),]
#Check for NAs and remove
IndRm = which((is.na(BR5S11_TN)) | (is.na(BR5S11_Q)))
BR5S11_Q = BR5S11_Q[-IndRm]
BR5S11_TN = BR5S11_TN[-IndRm]
rm(IndRm)
Dates_BR5S11 = vector('character', ncol(BR5S11_Q)-1)
for (i in 1:(ncol(BR5S11_Q)-1)){
  temp = strsplit(colnames(BR5S11_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_BR5S11[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_BR5S11 = as.Date(Dates_BR5S11)

#mg/s unit
BR5S11_TrueLoad = as.numeric(BR5S11_TN[-1]*1000*BR5S11_Q[-1])

Load_ECM_BR5S1105 = EC_Undev05*(Area_S11 - AreaImp_S11) + EC_Dev05*AreaImp_S11
Load_ECM_BR5S11 = EC_Undev*(Area_S11 - AreaImp_S11) + EC_Dev*AreaImp_S11
Load_ECM_BR5S1195 = EC_Undev95*(Area_S11 - AreaImp_S11) + EC_Dev95*AreaImp_S11
Load_ECM_BR5S11_QLQ05 = EC_Undev05*(Area_S11 - AreaImp_S11) + EC_DevQLQ05*AreaImp_S11
Load_ECM_BR5S11_QLQ = EC_Undev*(Area_S11 - AreaImp_S11) + EC_DevQLQ*AreaImp_S11
Load_ECM_BR5S11_QLQ95 = EC_Undev95*(Area_S11 - AreaImp_S11) + EC_DevQLQ95*AreaImp_S11

Load_ECM_BR5S11_Mat = EC_UndevMat*(Area_S11 - AreaImp_S11) + EC_DevMat*AreaImp_S11
Load_ECM_BR5S11_MatMean = apply(Load_ECM_BR5S11_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR5S11_Mat05 = apply(Load_ECM_BR5S11_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR5S11_Mat95 = apply(Load_ECM_BR5S11_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR5S11_MatQLQ = EC_UndevMatQLQ*(Area_S11 - AreaImp_S11) + EC_DevMatQLQ*AreaImp_S11
Load_ECM_BR5S11_MatQLQMean = apply(Load_ECM_BR5S11_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR5S11_MatQLQ05 = apply(Load_ECM_BR5S11_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR5S11_MatQLQ95 = apply(Load_ECM_BR5S11_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BR5S11, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR5S11), BR5S11_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BR5S11_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR5S11), BR5S11_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

png('ECMBR5aS11.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR5S11_TrueLoad[as.Date(Dates_BR5S11) %in% as.Date(Dates_BARN)], y = Load_ECM_BR5S11_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR5S11)], 
     xlim = c(0,80), ylim = c(0,80), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR5 S11')
arrows(BR5S11_TrueLoad[as.Date(Dates_BR5S11) %in% as.Date(Dates_BARN)], (Load_ECM_BR5S11_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR5S11)]), BR5S11_TrueLoad[as.Date(Dates_BR5S11) %in% as.Date(Dates_BARN)], (Load_ECM_BR5S11_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR5S11)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR5aS11_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR5S11_TrueLoad[as.Date(Dates_BR5S11) %in% as.Date(Dates_BARN)], Load_ECM_BR5S11_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR5S11)], 
     xlim = c(0,80), ylim = c(0,80), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR5 S11')
arrows(BR5S11_TrueLoad[as.Date(Dates_BR5S11) %in% as.Date(Dates_BARN)], (Load_ECM_BR5S11_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR5S11)]), BR5S11_TrueLoad[as.Date(Dates_BR5S11) %in% as.Date(Dates_BARN)], (Load_ECM_BR5S11_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR5S11)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    BARN K1 Synoptic - poor fit because of an Arc problem - should be whole watershed contributing, but it's only a small portion----
#Fixme: Get the impervious surface area contributing to this location
#AreaImp_K1S = 0
#Area_K1S = length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_7 == 1))*res^2
#for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_7 == 1))){
#  AreaImp_K1S = AreaImp_K1S + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_7 == 1)[i]]*res^2
#}

#Using basin impervious for now
Area_K1S = Area.basin
AreaImp_K1S = Area.basin_imp

#Observed data for this site
BARNK1S_TN = K_TN[K_TN$Site == K_Sites@data$NAME[8],]
BARNK1S_Q = K_Q[K_Q$Site == K_Sites@data$NAME[8],]
#Check for NAs and remove
IndRm = which((is.na(BARNK1S_TN)) | (is.na(BARNK1S_Q)))
BARNK1S_Q = BARNK1S_Q[-IndRm]
BARNK1S_TN = BARNK1S_TN[-IndRm]
rm(IndRm)
Dates_BARNK1S = vector('character', ncol(BARNK1S_Q)-1)
for (i in 1:(ncol(BARNK1S_Q)-1)){
  temp = strsplit(colnames(BARNK1S_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_BARNK1S[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_BARNK1S = as.Date(Dates_BARNK1S)

#mg/s unit
BARNK1S_TrueLoad = as.numeric(BARNK1S_TN[-1]*1000*BARNK1S_Q[-1])

Load_ECM_BARNK1S05 = EC_Undev05*(Area_K1S - AreaImp_K1S) + EC_Dev05*AreaImp_K1S
Load_ECM_BARNK1S = EC_Undev*(Area_K1S - AreaImp_K1S) + EC_Dev*AreaImp_K1S
Load_ECM_BARNK1S95 = EC_Undev95*(Area_K1S - AreaImp_K1S) + EC_Dev95*AreaImp_K1S
Load_ECM_BARNK1S_QLQ05 = EC_Undev05*(Area_K1S - AreaImp_K1S) + EC_DevQLQ05*AreaImp_K1S
Load_ECM_BARNK1S_QLQ = EC_Undev*(Area_K1S - AreaImp_K1S) + EC_DevQLQ*AreaImp_K1S
Load_ECM_BARNK1S_QLQ95 = EC_Undev95*(Area_K1S - AreaImp_K1S) + EC_DevQLQ95*AreaImp_K1S

Load_ECM_BARNK1S_Mat = EC_UndevMat*(Area_K1S - AreaImp_K1S) + EC_DevMat*AreaImp_K1S
Load_ECM_BARNK1S_MatMean = apply(Load_ECM_BARNK1S_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BARNK1S_Mat05 = apply(Load_ECM_BARNK1S_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BARNK1S_Mat95 = apply(Load_ECM_BARNK1S_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BARNK1S_MatQLQ = EC_UndevMatQLQ*(Area_K1S - AreaImp_K1S) + EC_DevMatQLQ*AreaImp_K1S
Load_ECM_BARNK1S_MatQLQMean = apply(Load_ECM_BARNK1S_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BARNK1S_MatQLQ05 = apply(Load_ECM_BARNK1S_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BARNK1S_MatQLQ95 = apply(Load_ECM_BARNK1S_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BARNK1S, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_BARNK1S), BARNK1S_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BARNK1S_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_BARNK1S), BARNK1S_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

png('ECMBARNK1s.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BARNK1S_TrueLoad[as.Date(Dates_BARNK1S) %in% as.Date(Dates_BARN)], y = Load_ECM_BARNK1S_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BARNK1S)], 
     xlim = c(0,50), ylim = c(0,50), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BARN K1S')
arrows(BARNK1S_TrueLoad[as.Date(Dates_BARNK1S) %in% as.Date(Dates_BARN)], (Load_ECM_BARNK1S_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BARNK1S)]), BARNK1S_TrueLoad[as.Date(Dates_BARNK1S) %in% as.Date(Dates_BARN)], (Load_ECM_BARNK1S_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BARNK1S)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBARNK1s_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BARNK1S_TrueLoad[as.Date(Dates_BARNK1S) %in% as.Date(Dates_BARN)], Load_ECM_BARNK1S_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BARNK1S)], 
     xlim = c(0,50), ylim = c(0,50), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BARN K1S')
arrows(BARNK1S_TrueLoad[as.Date(Dates_BARNK1S) %in% as.Date(Dates_BARN)], (Load_ECM_BARNK1S_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BARNK1S)]), BARNK1S_TrueLoad[as.Date(Dates_BARNK1S) %in% as.Date(Dates_BARN)], (Load_ECM_BARNK1S_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BARNK1S)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    POBR K2----
#Get the impervious surface area contributing to this location
AreaImp_K2 = 0
Area_K2 = length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_5 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_5 == 1))){
  AreaImp_K2 = AreaImp_K2 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_5 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
POBRK2_TN = K_TN[K_TN$Site == K_Sites@data$NAME[6],]
POBRK2_Q = K_Q[K_Q$Site == K_Sites@data$NAME[6],]
#Check for NAs and remove
IndRm = which((is.na(POBRK2_TN)) | (is.na(POBRK2_Q)))
POBRK2_Q = POBRK2_Q[-IndRm]
POBRK2_TN = POBRK2_TN[-IndRm]
rm(IndRm)
Dates_POBRK2 = vector('character', ncol(POBRK2_Q)-1)
for (i in 1:(ncol(POBRK2_Q)-1)){
  temp = strsplit(colnames(POBRK2_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_POBRK2[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_POBRK2 = as.Date(Dates_POBRK2)

#mg/s unit
POBRK2_TrueLoad = as.numeric(POBRK2_TN[-1]*1000*POBRK2_Q[-1])

Load_ECM_POBRK205 = EC_Undev05*(Area_K2 - AreaImp_K2) + EC_Dev05*AreaImp_K2
Load_ECM_POBRK2 = EC_Undev*(Area_K2 - AreaImp_K2) + EC_Dev*AreaImp_K2
Load_ECM_POBRK295 = EC_Undev95*(Area_K2 - AreaImp_K2) + EC_Dev95*AreaImp_K2
Load_ECM_POBRK2_QLQ05 = EC_Undev05*(Area_K2 - AreaImp_K2) + EC_DevQLQ05*AreaImp_K2
Load_ECM_POBRK2_QLQ = EC_Undev*(Area_K2 - AreaImp_K2) + EC_DevQLQ*AreaImp_K2
Load_ECM_POBRK2_QLQ95 = EC_Undev95*(Area_K2 - AreaImp_K2) + EC_DevQLQ95*AreaImp_K2

Load_ECM_POBRK2_Mat = EC_UndevMat*(Area_K2 - AreaImp_K2) + EC_DevMat*AreaImp_K2
Load_ECM_POBRK2_MatMean = apply(Load_ECM_POBRK2_Mat, MARGIN = 1, FUN = mean)
Load_ECM_POBRK2_Mat05 = apply(Load_ECM_POBRK2_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_POBRK2_Mat95 = apply(Load_ECM_POBRK2_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_POBRK2_MatQLQ = EC_UndevMatQLQ*(Area_K2 - AreaImp_K2) + EC_DevMatQLQ*AreaImp_K2
Load_ECM_POBRK2_MatQLQMean = apply(Load_ECM_POBRK2_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_POBRK2_MatQLQ05 = apply(Load_ECM_POBRK2_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_POBRK2_MatQLQ95 = apply(Load_ECM_POBRK2_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_POBRK2, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_POBRK2), POBRK2_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_POBRK2_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_POBRK2), POBRK2_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

png('ECMPOBRK2.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = POBRK2_TrueLoad[as.Date(Dates_POBRK2) %in% as.Date(Dates_BARN)], y = Load_ECM_POBRK2_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK2)], 
     xlim = c(0,5), ylim = c(0,5), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'POBR K2')
arrows(POBRK2_TrueLoad[as.Date(Dates_POBRK2) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK2_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK2)]), POBRK2_TrueLoad[as.Date(Dates_POBRK2) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK2_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK2)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMPOBRK2_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(POBRK2_TrueLoad[as.Date(Dates_POBRK2) %in% as.Date(Dates_BARN)], Load_ECM_POBRK2_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK2)], 
     xlim = c(0,5), ylim = c(0,5), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'POBR K2')
arrows(POBRK2_TrueLoad[as.Date(Dates_POBRK2) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK2_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK2)]), POBRK2_TrueLoad[as.Date(Dates_POBRK2) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK2_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK2)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    POBR K2 Synoptic----
#Get the impervious surface area contributing to this location
#This is same as used above for Pond Branch
#AreaImp_K2S = 0
#Area_K2S = length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_6 == 1))*res^2
#for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_6 == 1))){
#  AreaImp_K2S = AreaImp_K2S + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_6 == 1)[i]]*res^2
#}
#rm(i)
#Will treat this as all undeveloped land because fraction impervious is 7e-5.
#AreaImp_K2S = 0

#Observed data for this site
POBRK2S_TN = K_TN[K_TN$Site == K_Sites@data$NAME[7],]
POBRK2S_Q = K_Q[K_Q$Site == K_Sites@data$NAME[7],]
#Check for NAs and remove
IndRm = which((is.na(POBRK2S_TN)) | (is.na(POBRK2S_Q)))
POBRK2S_Q = POBRK2S_Q[-IndRm]
POBRK2S_TN = POBRK2S_TN[-IndRm]
rm(IndRm)
Dates_POBRK2S = vector('character', ncol(POBRK2S_Q)-1)
for (i in 1:(ncol(POBRK2S_Q)-1)){
  temp = strsplit(colnames(POBRK2S_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_POBRK2S[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_POBRK2S = as.Date(Dates_POBRK2S)

#mg/s unit
POBRK2S_TrueLoad = as.numeric(POBRK2S_TN[-1]*1000*POBRK2S_Q[-1])

Load_ECM_POBRK2S05 = EC_Undev05*(Area_K2S - AreaImp_K2S) + EC_Dev05*AreaImp_K2S
Load_ECM_POBRK2S = EC_Undev*(Area_K2S - AreaImp_K2S) + EC_Dev*AreaImp_K2S
Load_ECM_POBRK2S95 = EC_Undev95*(Area_K2S - AreaImp_K2S) + EC_Dev95*AreaImp_K2S
Load_ECM_POBRK2S_QLQ05 = EC_Undev05*(Area_K2S - AreaImp_K2S) + EC_DevQLQ05*AreaImp_K2S
Load_ECM_POBRK2S_QLQ = EC_Undev*(Area_K2S - AreaImp_K2S) + EC_DevQLQ*AreaImp_K2S
Load_ECM_POBRK2S_QLQ95 = EC_Undev95*(Area_K2S - AreaImp_K2S) + EC_DevQLQ95*AreaImp_K2S

Load_ECM_POBRK2S_Mat = EC_UndevMat*(Area_K2S - AreaImp_K2S) + EC_DevMat*AreaImp_K2S
Load_ECM_POBRK2S_MatMean = apply(Load_ECM_POBRK2S_Mat, MARGIN = 1, FUN = mean)
Load_ECM_POBRK2S_Mat05 = apply(Load_ECM_POBRK2S_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_POBRK2S_Mat95 = apply(Load_ECM_POBRK2S_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_POBRK2S_MatQLQ = EC_UndevMatQLQ*(Area_K2S - AreaImp_K2S) + EC_DevMatQLQ*AreaImp_K2S
Load_ECM_POBRK2S_MatQLQMean = apply(Load_ECM_POBRK2S_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_POBRK2S_MatQLQ05 = apply(Load_ECM_POBRK2S_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_POBRK2S_MatQLQ95 = apply(Load_ECM_POBRK2S_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_POBRK2S, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_POBRK2S), POBRK2S_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_POBRK2S_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_POBRK2S), POBRK2S_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

png('ECMPOBRK2s.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = POBRK2S_TrueLoad[as.Date(Dates_POBRK2S) %in% as.Date(Dates_BARN)], y = Load_ECM_POBRK2S_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK2S)], 
     xlim = c(0,5), ylim = c(0,5), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'POBR K2S')
arrows(POBRK2S_TrueLoad[as.Date(Dates_POBRK2S) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK2S_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK2S)]), POBRK2S_TrueLoad[as.Date(Dates_POBRK2S) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK2S_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK2S)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMPOBRK2s_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(POBRK2S_TrueLoad[as.Date(Dates_POBRK2S) %in% as.Date(Dates_BARN)], Load_ECM_POBRK2S_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK2S)], 
     xlim = c(0,5), ylim = c(0,5), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'POBR K2S')
arrows(POBRK2S_TrueLoad[as.Date(Dates_POBRK2S) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK2S_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK2S)]), POBRK2S_TrueLoad[as.Date(Dates_POBRK2S) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK2S_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK2S)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    POBR K3a----
#Get the impervious surface area contributing to this location
AreaImp_K3a = 0
Area_K3a = length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_9 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_9 == 1))){
  AreaImp_K3a = AreaImp_K3a + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_9 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
POBRK3a_TN = K_TN[K_TN$Site == K_Sites@data$NAME[10],]
POBRK3a_Q = K_Q[K_Q$Site == K_Sites@data$NAME[10],]
#Check for NAs and remove
IndRm = which((is.na(POBRK3a_TN)) | (is.na(POBRK3a_Q)))
POBRK3a_Q = POBRK3a_Q[-IndRm]
POBRK3a_TN = POBRK3a_TN[-IndRm]
rm(IndRm)
Dates_POBRK3a = vector('character', ncol(POBRK3a_Q)-1)
for (i in 1:(ncol(POBRK3a_Q)-1)){
  temp = strsplit(colnames(POBRK3a_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_POBRK3a[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_POBRK3a = as.Date(Dates_POBRK3a)

#mg/s unit
POBRK3a_TrueLoad = as.numeric(POBRK3a_TN[-1]*1000*POBRK3a_Q[-1])

Load_ECM_POBRK3a05 = EC_Undev05*(Area_K3a - AreaImp_K3a) + EC_Dev05*AreaImp_K3a
Load_ECM_POBRK3a = EC_Undev*(Area_K3a - AreaImp_K3a) + EC_Dev*AreaImp_K3a
Load_ECM_POBRK3a95 = EC_Undev95*(Area_K3a - AreaImp_K3a) + EC_Dev95*AreaImp_K3a
Load_ECM_POBRK3a_QLQ05 = EC_Undev05*(Area_K3a - AreaImp_K3a) + EC_DevQLQ05*AreaImp_K3a
Load_ECM_POBRK3a_QLQ = EC_Undev*(Area_K3a - AreaImp_K3a) + EC_DevQLQ*AreaImp_K3a
Load_ECM_POBRK3a_QLQ95 = EC_Undev95*(Area_K3a - AreaImp_K3a) + EC_DevQLQ95*AreaImp_K3a

Load_ECM_POBRK3a_Mat = EC_UndevMat*(Area_K3a - AreaImp_K3a) + EC_DevMat*AreaImp_K3a
Load_ECM_POBRK3a_MatMean = apply(Load_ECM_POBRK3a_Mat, MARGIN = 1, FUN = mean)
Load_ECM_POBRK3a_Mat05 = apply(Load_ECM_POBRK3a_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_POBRK3a_Mat95 = apply(Load_ECM_POBRK3a_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_POBRK3a_MatQLQ = EC_UndevMatQLQ*(Area_K3a - AreaImp_K3a) + EC_DevMatQLQ*AreaImp_K3a
Load_ECM_POBRK3a_MatQLQMean = apply(Load_ECM_POBRK3a_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_POBRK3a_MatQLQ05 = apply(Load_ECM_POBRK3a_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_POBRK3a_MatQLQ95 = apply(Load_ECM_POBRK3a_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_POBRK3a, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_POBRK3a), POBRK3a_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_POBRK3a_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_POBRK3a), POBRK3a_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

png('ECMPOBRK3a.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = POBRK3a_TrueLoad[as.Date(Dates_POBRK3a) %in% as.Date(Dates_BARN)], y = Load_ECM_POBRK3a_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK3a)], 
     xlim = c(0,1), ylim = c(0,1), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'POBR K3a')
arrows(POBRK3a_TrueLoad[as.Date(Dates_POBRK3a) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK3a_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK3a)]), POBRK3a_TrueLoad[as.Date(Dates_POBRK3a) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK3a_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK3a)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMPOBRK3a_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(POBRK3a_TrueLoad[as.Date(Dates_POBRK3a) %in% as.Date(Dates_BARN)], Load_ECM_POBRK3a_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK3a)], 
     xlim = c(0,1), ylim = c(0,1), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'POBR K3a')
arrows(POBRK3a_TrueLoad[as.Date(Dates_POBRK3a) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK3a_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK3a)]), POBRK3a_TrueLoad[as.Date(Dates_POBRK3a) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK3a_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK3a)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    POBR K3----
#Get the impervious surface area contributing to this location
AreaImp_K3 = 0
Area_K3 = length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_10 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_10 == 1))){
  AreaImp_K3 = AreaImp_K3 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_10 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
POBRK3_TN = K_TN[K_TN$Site == K_Sites@data$NAME[11],]
POBRK3_Q = K_Q[K_Q$Site == K_Sites@data$NAME[11],]
#Check for NAs and remove
IndRm = which((is.na(POBRK3_TN)) | (is.na(POBRK3_Q)))
POBRK3_Q = POBRK3_Q[-IndRm]
POBRK3_TN = POBRK3_TN[-IndRm]
rm(IndRm)
Dates_POBRK3 = vector('character', ncol(POBRK3_Q)-1)
for (i in 1:(ncol(POBRK3_Q)-1)){
  temp = strsplit(colnames(POBRK3_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_POBRK3[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_POBRK3 = as.Date(Dates_POBRK3)

#mg/s unit
POBRK3_TrueLoad = as.numeric(POBRK3_TN[-1]*1000*POBRK3_Q[-1])

Load_ECM_POBRK305 = EC_Undev05*(Area_K3 - AreaImp_K3) + EC_Dev05*AreaImp_K3
Load_ECM_POBRK3 = EC_Undev*(Area_K3 - AreaImp_K3) + EC_Dev*AreaImp_K3
Load_ECM_POBRK395 = EC_Undev95*(Area_K3 - AreaImp_K3) + EC_Dev95*AreaImp_K3
Load_ECM_POBRK3_QLQ05 = EC_Undev05*(Area_K3 - AreaImp_K3) + EC_DevQLQ05*AreaImp_K3
Load_ECM_POBRK3_QLQ = EC_Undev*(Area_K3 - AreaImp_K3) + EC_DevQLQ*AreaImp_K3
Load_ECM_POBRK3_QLQ95 = EC_Undev95*(Area_K3 - AreaImp_K3) + EC_DevQLQ95*AreaImp_K3

Load_ECM_POBRK3_Mat = EC_UndevMat*(Area_K3 - AreaImp_K3) + EC_DevMat*AreaImp_K3
Load_ECM_POBRK3_MatMean = apply(Load_ECM_POBRK3_Mat, MARGIN = 1, FUN = mean)
Load_ECM_POBRK3_Mat05 = apply(Load_ECM_POBRK3_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_POBRK3_Mat95 = apply(Load_ECM_POBRK3_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_POBRK3_MatQLQ = EC_UndevMatQLQ*(Area_K3 - AreaImp_K3) + EC_DevMatQLQ*AreaImp_K3
Load_ECM_POBRK3_MatQLQMean = apply(Load_ECM_POBRK3_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_POBRK3_MatQLQ05 = apply(Load_ECM_POBRK3_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_POBRK3_MatQLQ95 = apply(Load_ECM_POBRK3_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_POBRK3, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_POBRK3), POBRK3_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_POBRK3_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_POBRK3), POBRK3_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

png('ECMPOBRK3.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = POBRK3_TrueLoad[as.Date(Dates_POBRK3) %in% as.Date(Dates_BARN)], y = Load_ECM_POBRK3_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK3)], 
     xlim = c(0,1), ylim = c(0,1), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'POBR K3')
arrows(POBRK3_TrueLoad[as.Date(Dates_POBRK3) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK3_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK3)]), POBRK3_TrueLoad[as.Date(Dates_POBRK3) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK3_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK3)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMPOBRK3_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(POBRK3_TrueLoad[as.Date(Dates_POBRK3) %in% as.Date(Dates_BARN)], Load_ECM_POBRK3_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK3)], 
     xlim = c(0,1), ylim = c(0,1), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'POBR K3')
arrows(POBRK3_TrueLoad[as.Date(Dates_POBRK3) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK3_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK3)]), POBRK3_TrueLoad[as.Date(Dates_POBRK3) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK3_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK3)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    POBR K11----
#Get the impervious surface area contributing to this location
AreaImp_K11 = 0
Area_K11 = length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_11 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_11 == 1))){
  AreaImp_K11 = AreaImp_K11 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_11 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
POBRK11_TN = K_TN[K_TN$Site == K_Sites@data$NAME[12],]
POBRK11_Q = K_Q[K_Q$Site == K_Sites@data$NAME[12],]
#Check for NAs and remove
IndRm = which((is.na(POBRK11_TN)) | (is.na(POBRK11_Q)))
POBRK11_Q = POBRK11_Q[-IndRm]
POBRK11_TN = POBRK11_TN[-IndRm]
rm(IndRm)
Dates_POBRK11 = vector('character', ncol(POBRK11_Q)-1)
for (i in 1:(ncol(POBRK11_Q)-1)){
  temp = strsplit(colnames(POBRK11_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_POBRK11[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_POBRK11 = as.Date(Dates_POBRK11)

#mg/s unit
POBRK11_TrueLoad = as.numeric(POBRK11_TN[-1]*1000*POBRK11_Q[-1])

Load_ECM_POBRK1105 = EC_Undev05*(Area_K11 - AreaImp_K11) + EC_Dev05*AreaImp_K11
Load_ECM_POBRK11 = EC_Undev*(Area_K11 - AreaImp_K11) + EC_Dev*AreaImp_K11
Load_ECM_POBRK1195 = EC_Undev95*(Area_K11 - AreaImp_K11) + EC_Dev95*AreaImp_K11
Load_ECM_POBRK11_QLQ05 = EC_Undev05*(Area_K11 - AreaImp_K11) + EC_DevQLQ05*AreaImp_K11
Load_ECM_POBRK11_QLQ = EC_Undev*(Area_K11 - AreaImp_K11) + EC_DevQLQ*AreaImp_K11
Load_ECM_POBRK11_QLQ95 = EC_Undev95*(Area_K11 - AreaImp_K11) + EC_DevQLQ95*AreaImp_K11

Load_ECM_POBRK11_Mat = EC_UndevMat*(Area_K11 - AreaImp_K11) + EC_DevMat*AreaImp_K11
Load_ECM_POBRK11_MatMean = apply(Load_ECM_POBRK11_Mat, MARGIN = 1, FUN = mean)
Load_ECM_POBRK11_Mat05 = apply(Load_ECM_POBRK11_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_POBRK11_Mat95 = apply(Load_ECM_POBRK11_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_POBRK11_MatQLQ = EC_UndevMatQLQ*(Area_K11 - AreaImp_K11) + EC_DevMatQLQ*AreaImp_K11
Load_ECM_POBRK11_MatQLQMean = apply(Load_ECM_POBRK11_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_POBRK11_MatQLQ05 = apply(Load_ECM_POBRK11_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_POBRK11_MatQLQ95 = apply(Load_ECM_POBRK11_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_POBRK11, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_POBRK11), POBRK11_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_POBRK11_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_POBRK11), POBRK11_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

png('ECMPOBRK11.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = POBRK11_TrueLoad[as.Date(Dates_POBRK11) %in% as.Date(Dates_BARN)], y = Load_ECM_POBRK11_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK11)], 
     xlim = c(0,1), ylim = c(0,1), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'POBR K11')
arrows(POBRK11_TrueLoad[as.Date(Dates_POBRK11) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK11_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK11)]), POBRK11_TrueLoad[as.Date(Dates_POBRK11) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK11_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK11)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMPOBRK11_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(POBRK11_TrueLoad[as.Date(Dates_POBRK11) %in% as.Date(Dates_BARN)], Load_ECM_POBRK11_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK11)], 
     xlim = c(0,1), ylim = c(0,1), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'POBR K11')
arrows(POBRK11_TrueLoad[as.Date(Dates_POBRK11) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK11_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK11)]), POBRK11_TrueLoad[as.Date(Dates_POBRK11) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK11_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK11)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    POBR K8----
#Get the impervious surface area contributing to this location
AreaImp_K8 = 0
Area_K8 = length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_8 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_8 == 1))){
  AreaImp_K8 = AreaImp_K8 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_8 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
POBRK8_TN = K_TN[K_TN$Site == K_Sites@data$NAME[9],]
POBRK8_Q = K_Q[K_Q$Site == K_Sites@data$NAME[9],]
#Check for NAs and remove
IndRm = which((is.na(POBRK8_TN)) | (is.na(POBRK8_Q)))
POBRK8_Q = POBRK8_Q[-IndRm]
POBRK8_TN = POBRK8_TN[-IndRm]
rm(IndRm)
Dates_POBRK8 = vector('character', ncol(POBRK8_Q)-1)
for (i in 1:(ncol(POBRK8_Q)-1)){
  temp = strsplit(colnames(POBRK8_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_POBRK8[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_POBRK8 = as.Date(Dates_POBRK8)

#mg/s unit
POBRK8_TrueLoad = as.numeric(POBRK8_TN[-1]*1000*POBRK8_Q[-1])

Load_ECM_POBRK805 = EC_Undev05*(Area_K8 - AreaImp_K8) + EC_Dev05*AreaImp_K8
Load_ECM_POBRK8 = EC_Undev*(Area_K8 - AreaImp_K8) + EC_Dev*AreaImp_K8
Load_ECM_POBRK895 = EC_Undev95*(Area_K8 - AreaImp_K8) + EC_Dev95*AreaImp_K8
Load_ECM_POBRK8_QLQ05 = EC_Undev05*(Area_K8 - AreaImp_K8) + EC_DevQLQ05*AreaImp_K8
Load_ECM_POBRK8_QLQ = EC_Undev*(Area_K8 - AreaImp_K8) + EC_DevQLQ*AreaImp_K8
Load_ECM_POBRK8_QLQ95 = EC_Undev95*(Area_K8 - AreaImp_K8) + EC_DevQLQ95*AreaImp_K8

Load_ECM_POBRK8_Mat = EC_UndevMat*(Area_K8 - AreaImp_K8) + EC_DevMat*AreaImp_K8
Load_ECM_POBRK8_MatMean = apply(Load_ECM_POBRK8_Mat, MARGIN = 1, FUN = mean)
Load_ECM_POBRK8_Mat05 = apply(Load_ECM_POBRK8_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_POBRK8_Mat95 = apply(Load_ECM_POBRK8_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_POBRK8_MatQLQ = EC_UndevMatQLQ*(Area_K8 - AreaImp_K8) + EC_DevMatQLQ*AreaImp_K8
Load_ECM_POBRK8_MatQLQMean = apply(Load_ECM_POBRK8_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_POBRK8_MatQLQ05 = apply(Load_ECM_POBRK8_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_POBRK8_MatQLQ95 = apply(Load_ECM_POBRK8_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_POBRK8, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_POBRK8), POBRK8_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_POBRK8_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40))
par(new = TRUE)
plot(as.Date(Dates_POBRK8), POBRK8_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,40), col = 'red')

png('ECMPOBRK8.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = POBRK8_TrueLoad[as.Date(Dates_POBRK8) %in% as.Date(Dates_BARN)], y = Load_ECM_POBRK8_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK8)], 
     xlim = c(0,1), ylim = c(0,1), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'POBR K8')
arrows(POBRK8_TrueLoad[as.Date(Dates_POBRK8) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK8_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK8)]), POBRK8_TrueLoad[as.Date(Dates_POBRK8) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK8_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK8)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMPOBRK8_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(POBRK8_TrueLoad[as.Date(Dates_POBRK8) %in% as.Date(Dates_BARN)], Load_ECM_POBRK8_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK8)], 
     xlim = c(0,1), ylim = c(0,1), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'POBR K8')
arrows(POBRK8_TrueLoad[as.Date(Dates_POBRK8) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK8_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK8)]), POBRK8_TrueLoad[as.Date(Dates_POBRK8) %in% as.Date(Dates_BARN)], (Load_ECM_POBRK8_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_POBRK8)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    BR4 K4----
#Get the impervious surface area contributing to this location
AreaImp_K4 = 0
Area_K4 = length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_2 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_2 == 1))){
  AreaImp_K4 = AreaImp_K4 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_2 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
BR4K4_TN = K_TN[K_TN$Site == K_Sites@data$NAME[3],]
BR4K4_Q = K_Q[K_Q$Site == K_Sites@data$NAME[3],]
#Check for NAs and remove
IndRm = which((is.na(BR4K4_TN)) | (is.na(BR4K4_Q)))
BR4K4_Q = BR4K4_Q[-IndRm]
BR4K4_TN = BR4K4_TN[-IndRm]
rm(IndRm)
Dates_BR4K4 = vector('character', ncol(BR4K4_Q)-1)
for (i in 1:(ncol(BR4K4_Q)-1)){
  temp = strsplit(colnames(BR4K4_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_BR4K4[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_BR4K4 = as.Date(Dates_BR4K4)

#mg/s unit
BR4K4_TrueLoad = as.numeric(BR4K4_TN[-1]*1000*BR4K4_Q[-1])

Load_ECM_BR4K405 = EC_Undev05*(Area_K4 - AreaImp_K4) + EC_Dev05*AreaImp_K4
Load_ECM_BR4K4 = EC_Undev*(Area_K4 - AreaImp_K4) + EC_Dev*AreaImp_K4
Load_ECM_BR4K495 = EC_Undev95*(Area_K4 - AreaImp_K4) + EC_Dev95*AreaImp_K4
Load_ECM_BR4K4_QLQ05 = EC_Undev05*(Area_K4 - AreaImp_K4) + EC_DevQLQ05*AreaImp_K4
Load_ECM_BR4K4_QLQ = EC_Undev*(Area_K4 - AreaImp_K4) + EC_DevQLQ*AreaImp_K4
Load_ECM_BR4K4_QLQ95 = EC_Undev95*(Area_K4 - AreaImp_K4) + EC_DevQLQ95*AreaImp_K4

Load_ECM_BR4K4_Mat = EC_UndevMat*(Area_K4 - AreaImp_K4) + EC_DevMat*AreaImp_K4
Load_ECM_BR4K4_MatMean = apply(Load_ECM_BR4K4_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR4K4_Mat05 = apply(Load_ECM_BR4K4_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR4K4_Mat95 = apply(Load_ECM_BR4K4_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR4K4_MatQLQ = EC_UndevMatQLQ*(Area_K4 - AreaImp_K4) + EC_DevMatQLQ*AreaImp_K4
Load_ECM_BR4K4_MatQLQMean = apply(Load_ECM_BR4K4_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR4K4_MatQLQ05 = apply(Load_ECM_BR4K4_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR4K4_MatQLQ95 = apply(Load_ECM_BR4K4_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BR4K4, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR4K4), BR4K4_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BR4K4_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR4K4), BR4K4_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

png('ECMBR4K4.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR4K4_TrueLoad[as.Date(Dates_BR4K4) %in% as.Date(Dates_BARN)], y = Load_ECM_BR4K4_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR4K4)], 
     xlim = c(0,5), ylim = c(0,5), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR4 K4')
arrows(BR4K4_TrueLoad[as.Date(Dates_BR4K4) %in% as.Date(Dates_BARN)], (Load_ECM_BR4K4_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR4K4)]), BR4K4_TrueLoad[as.Date(Dates_BR4K4) %in% as.Date(Dates_BARN)], (Load_ECM_BR4K4_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR4K4)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR4K4_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR4K4_TrueLoad[as.Date(Dates_BR4K4) %in% as.Date(Dates_BARN)], Load_ECM_BR4K4_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR4K4)], 
     xlim = c(0,5), ylim = c(0,5), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR4 K4')
arrows(BR4K4_TrueLoad[as.Date(Dates_BR4K4) %in% as.Date(Dates_BARN)], (Load_ECM_BR4K4_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR4K4)]), BR4K4_TrueLoad[as.Date(Dates_BR4K4) %in% as.Date(Dates_BARN)], (Load_ECM_BR4K4_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR4K4)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    BR5 K5----
#Get the impervious surface area contributing to this location
AreaImp_K5 = 0
Area_K5 = length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_3 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_3 == 1))){
  AreaImp_K5 = AreaImp_K5 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_3 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
BR5K5_TN = K_TN[K_TN$Site == K_Sites@data$NAME[4],]
BR5K5_Q = K_Q[K_Q$Site == K_Sites@data$NAME[4],]
#Check for NAs and remove
IndRm = which((is.na(BR5K5_TN)) | (is.na(BR5K5_Q)))
BR5K5_Q = BR5K5_Q[-IndRm]
BR5K5_TN = BR5K5_TN[-IndRm]
rm(IndRm)
Dates_BR5K5 = vector('character', ncol(BR5K5_Q)-1)
for (i in 1:(ncol(BR5K5_Q)-1)){
  temp = strsplit(colnames(BR5K5_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_BR5K5[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_BR5K5 = as.Date(Dates_BR5K5)

#mg/s unit
BR5K5_TrueLoad = as.numeric(BR5K5_TN[-1]*1000*BR5K5_Q[-1])

Load_ECM_BR5K505 = EC_Undev05*(Area_K5 - AreaImp_K5) + EC_Dev05*AreaImp_K5
Load_ECM_BR5K5 = EC_Undev*(Area_K5 - AreaImp_K5) + EC_Dev*AreaImp_K5
Load_ECM_BR5K595 = EC_Undev95*(Area_K5 - AreaImp_K5) + EC_Dev95*AreaImp_K5
Load_ECM_BR5K5_QLQ05 = EC_Undev05*(Area_K5 - AreaImp_K5) + EC_DevQLQ05*AreaImp_K5
Load_ECM_BR5K5_QLQ = EC_Undev*(Area_K5 - AreaImp_K5) + EC_DevQLQ*AreaImp_K5
Load_ECM_BR5K5_QLQ95 = EC_Undev95*(Area_K5 - AreaImp_K5) + EC_DevQLQ95*AreaImp_K5

Load_ECM_BR5K5_Mat = EC_UndevMat*(Area_K5 - AreaImp_K5) + EC_DevMat*AreaImp_K5
Load_ECM_BR5K5_MatMean = apply(Load_ECM_BR5K5_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR5K5_Mat05 = apply(Load_ECM_BR5K5_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR5K5_Mat95 = apply(Load_ECM_BR5K5_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR5K5_MatQLQ = EC_UndevMatQLQ*(Area_K5 - AreaImp_K5) + EC_DevMatQLQ*AreaImp_K5
Load_ECM_BR5K5_MatQLQMean = apply(Load_ECM_BR5K5_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR5K5_MatQLQ05 = apply(Load_ECM_BR5K5_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR5K5_MatQLQ95 = apply(Load_ECM_BR5K5_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BR5K5, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR5K5), BR5K5_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BR5K5_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR5K5), BR5K5_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

png('ECMBR5K5.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR5K5_TrueLoad[as.Date(Dates_BR5K5) %in% as.Date(Dates_BARN)], y = Load_ECM_BR5K5_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5)], 
     xlim = c(0,20), ylim = c(0,20), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR5 K5')
arrows(BR5K5_TrueLoad[as.Date(Dates_BR5K5) %in% as.Date(Dates_BARN)], (Load_ECM_BR5K5_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5)]), BR5K5_TrueLoad[as.Date(Dates_BR5K5) %in% as.Date(Dates_BARN)], (Load_ECM_BR5K5_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR5K5_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR5K5_TrueLoad[as.Date(Dates_BR5K5) %in% as.Date(Dates_BARN)], Load_ECM_BR5K5_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5)], 
     xlim = c(0,20), ylim = c(0,20), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR5 K5')
arrows(BR5K5_TrueLoad[as.Date(Dates_BR5K5) %in% as.Date(Dates_BARN)], (Load_ECM_BR5K5_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5)]), BR5K5_TrueLoad[as.Date(Dates_BR5K5) %in% as.Date(Dates_BARN)], (Load_ECM_BR5K5_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    BR5a K5a----
#Get the impervious surface area contributing to this location
AreaImp_K5a = 0
Area_K5a = length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_12 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_12 == 1))){
  AreaImp_K5a = AreaImp_K5a + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_12 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
BR5K5a_TN = K_TN[K_TN$Site == K_Sites@data$NAME[4],]
BR5K5a_Q = K_Q[K_Q$Site == K_Sites@data$NAME[4],]
#Check for NAs and remove
IndRm = which((is.na(BR5K5a_TN)) | (is.na(BR5K5a_Q)))
BR5K5a_Q = BR5K5a_Q[-IndRm]
BR5K5a_TN = BR5K5a_TN[-IndRm]
rm(IndRm)
Dates_BR5K5a = vector('character', ncol(BR5K5a_Q)-1)
for (i in 1:(ncol(BR5K5a_Q)-1)){
  temp = strsplit(colnames(BR5K5a_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_BR5K5a[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_BR5K5a = as.Date(Dates_BR5K5a)

#mg/s unit
BR5K5a_TrueLoad = as.numeric(BR5K5a_TN[-1]*1000*BR5K5a_Q[-1])

Load_ECM_BR5K5a05 = EC_Undev05*(Area_K5a - AreaImp_K5a) + EC_Dev05*AreaImp_K5a
Load_ECM_BR5K5a = EC_Undev*(Area_K5a - AreaImp_K5a) + EC_Dev*AreaImp_K5a
Load_ECM_BR5K5a95 = EC_Undev95*(Area_K5a - AreaImp_K5a) + EC_Dev95*AreaImp_K5a
Load_ECM_BR5K5a_QLQ05 = EC_Undev05*(Area_K5a - AreaImp_K5a) + EC_DevQLQ05*AreaImp_K5a
Load_ECM_BR5K5a_QLQ = EC_Undev*(Area_K5a - AreaImp_K5a) + EC_DevQLQ*AreaImp_K5a
Load_ECM_BR5K5a_QLQ95 = EC_Undev95*(Area_K5a - AreaImp_K5a) + EC_DevQLQ95*AreaImp_K5a

Load_ECM_BR5K5a_Mat = EC_UndevMat*(Area_K5a - AreaImp_K5a) + EC_DevMat*AreaImp_K5a
Load_ECM_BR5K5a_MatMean = apply(Load_ECM_BR5K5a_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR5K5a_Mat05 = apply(Load_ECM_BR5K5a_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR5K5a_Mat95 = apply(Load_ECM_BR5K5a_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR5K5a_MatQLQ = EC_UndevMatQLQ*(Area_K5a - AreaImp_K5a) + EC_DevMatQLQ*AreaImp_K5a
Load_ECM_BR5K5a_MatQLQMean = apply(Load_ECM_BR5K5a_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR5K5a_MatQLQ05 = apply(Load_ECM_BR5K5a_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR5K5a_MatQLQ95 = apply(Load_ECM_BR5K5a_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BR5K5a, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR5K5a), BR5K5a_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BR5K5a_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR5K5a), BR5K5a_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

png('ECMBR5K5a.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR5K5a_TrueLoad[as.Date(Dates_BR5K5a) %in% as.Date(Dates_BARN)], y = Load_ECM_BR5K5a_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5a)], 
     xlim = c(0,20), ylim = c(0,20), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR5 K5a')
arrows(BR5K5a_TrueLoad[as.Date(Dates_BR5K5a) %in% as.Date(Dates_BARN)], (Load_ECM_BR5K5a_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5a)]), BR5K5a_TrueLoad[as.Date(Dates_BR5K5a) %in% as.Date(Dates_BARN)], (Load_ECM_BR5K5a_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5a)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR5K5a_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR5K5a_TrueLoad[as.Date(Dates_BR5K5a) %in% as.Date(Dates_BARN)], Load_ECM_BR5K5a_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5a)], 
     xlim = c(0,20), ylim = c(0,20), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR5 K5a')
arrows(BR5K5a_TrueLoad[as.Date(Dates_BR5K5a) %in% as.Date(Dates_BARN)], (Load_ECM_BR5K5a_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5a)]), BR5K5a_TrueLoad[as.Date(Dates_BR5K5a) %in% as.Date(Dates_BARN)], (Load_ECM_BR5K5a_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5a)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    BR6 K6----
#Get the impervious surface area contributing to this location
AreaImp_K6 = 0
Area_K6 = length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_4 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_4 == 1))){
  AreaImp_K6 = AreaImp_K6 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_4 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
BR6K6_TN = K_TN[K_TN$Site == K_Sites@data$NAME[5],]
BR6K6_Q = K_Q[K_Q$Site == K_Sites@data$NAME[5],]
#Check for NAs and remove
IndRm = which((is.na(BR6K6_TN)) | (is.na(BR6K6_Q)))
BR6K6_Q = BR6K6_Q[-IndRm]
BR6K6_TN = BR6K6_TN[-IndRm]
rm(IndRm)
Dates_BR6K6 = vector('character', ncol(BR6K6_Q)-1)
for (i in 1:(ncol(BR6K6_Q)-1)){
  temp = strsplit(colnames(BR6K6_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_BR6K6[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_BR6K6 = as.Date(Dates_BR6K6)

#mg/s unit
BR6K6_TrueLoad = as.numeric(BR6K6_TN[-1]*1000*BR6K6_Q[-1])

Load_ECM_BR6K605 = EC_Undev05*(Area_K6 - AreaImp_K6) + EC_Dev05*AreaImp_K6
Load_ECM_BR6K6 = EC_Undev*(Area_K6 - AreaImp_K6) + EC_Dev*AreaImp_K6
Load_ECM_BR6K695 = EC_Undev95*(Area_K6 - AreaImp_K6) + EC_Dev95*AreaImp_K6
Load_ECM_BR6K6_QLQ05 = EC_Undev05*(Area_K6 - AreaImp_K6) + EC_DevQLQ05*AreaImp_K6
Load_ECM_BR6K6_QLQ = EC_Undev*(Area_K6 - AreaImp_K6) + EC_DevQLQ*AreaImp_K6
Load_ECM_BR6K6_QLQ95 = EC_Undev95*(Area_K6 - AreaImp_K6) + EC_DevQLQ95*AreaImp_K6

Load_ECM_BR6K6_Mat = EC_UndevMat*(Area_K6 - AreaImp_K6) + EC_DevMat*AreaImp_K6
Load_ECM_BR6K6_MatMean = apply(Load_ECM_BR6K6_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR6K6_Mat05 = apply(Load_ECM_BR6K6_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR6K6_Mat95 = apply(Load_ECM_BR6K6_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR6K6_MatQLQ = EC_UndevMatQLQ*(Area_K6 - AreaImp_K6) + EC_DevMatQLQ*AreaImp_K6
Load_ECM_BR6K6_MatQLQMean = apply(Load_ECM_BR6K6_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR6K6_MatQLQ05 = apply(Load_ECM_BR6K6_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR6K6_MatQLQ95 = apply(Load_ECM_BR6K6_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BR6K6, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR6K6), BR6K6_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BR6K6_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR6K6), BR6K6_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

png('ECMBR6K6.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR6K6_TrueLoad[as.Date(Dates_BR6K6) %in% as.Date(Dates_BARN)], y = Load_ECM_BR6K6_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR6K6)], 
    xlim = c(0,15), ylim = c(0,15), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR6 K6')
arrows(BR6K6_TrueLoad[as.Date(Dates_BR6K6) %in% as.Date(Dates_BARN)], (Load_ECM_BR6K6_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR6K6)]), BR6K6_TrueLoad[as.Date(Dates_BR6K6) %in% as.Date(Dates_BARN)], (Load_ECM_BR6K6_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR6K6)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR6K6_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR6K6_TrueLoad[as.Date(Dates_BR6K6) %in% as.Date(Dates_BARN)], Load_ECM_BR6K6_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR6K6)], 
     xlim = c(0,15), ylim = c(0,15), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR6 K6')
arrows(BR6K6_TrueLoad[as.Date(Dates_BR6K6) %in% as.Date(Dates_BARN)], (Load_ECM_BR6K6_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR6K6)]), BR6K6_TrueLoad[as.Date(Dates_BR6K6) %in% as.Date(Dates_BARN)], (Load_ECM_BR6K6_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR6K6)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    BR7 K7 - poor fit might result from Arc method----
#Get the impervious surface area contributing to this location
AreaImp_K7 = 0
Area_K7 = length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_0 == 1))*res^2
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_0 == 1))){
  AreaImp_K7 = AreaImp_K7 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_0 == 1)[i]]*res^2
}
rm(i)
#Observed data for this site
BR7K7_TN = K_TN[K_TN$Site == K_Sites@data$NAME[1],]
BR7K7_Q = K_Q[K_Q$Site == K_Sites@data$NAME[1],]
#Check for NAs and remove
IndRm = which((is.na(BR7K7_TN)) | (is.na(BR7K7_Q)))
BR7K7_Q = BR7K7_Q[-IndRm]
BR7K7_TN = BR7K7_TN[-IndRm]
rm(IndRm)
Dates_BR7K7 = vector('character', ncol(BR7K7_Q)-1)
for (i in 1:(ncol(BR7K7_Q)-1)){
  temp = strsplit(colnames(BR7K7_Q)[i+1], split = '.', fixed = TRUE)[[1]]
  Dates_BR7K7[i] = paste(temp[2], temp[3], temp[4], sep = '-')
}
Dates_BR7K7 = as.Date(Dates_BR7K7)

#mg/s unit
BR7K7_TrueLoad = as.numeric(BR7K7_TN[-1]*1000*BR7K7_Q[-1])

Load_ECM_BR7K705 = EC_Undev05*(Area_K7 - AreaImp_K7) + EC_Dev05*AreaImp_K7
Load_ECM_BR7K7 = EC_Undev*(Area_K7 - AreaImp_K7) + EC_Dev*AreaImp_K7
Load_ECM_BR7K795 = EC_Undev95*(Area_K7 - AreaImp_K7) + EC_Dev95*AreaImp_K7
Load_ECM_BR7K7_QLQ05 = EC_Undev05*(Area_K7 - AreaImp_K7) + EC_DevQLQ05*AreaImp_K7
Load_ECM_BR7K7_QLQ = EC_Undev*(Area_K7 - AreaImp_K7) + EC_DevQLQ*AreaImp_K7
Load_ECM_BR7K7_QLQ95 = EC_Undev95*(Area_K7 - AreaImp_K7) + EC_DevQLQ95*AreaImp_K7

Load_ECM_BR7K7_Mat = EC_UndevMat*(Area_K7 - AreaImp_K7) + EC_DevMat*AreaImp_K7
Load_ECM_BR7K7_MatMean = apply(Load_ECM_BR7K7_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR7K7_Mat05 = apply(Load_ECM_BR7K7_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR7K7_Mat95 = apply(Load_ECM_BR7K7_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR7K7_MatQLQ = EC_UndevMatQLQ*(Area_K7 - AreaImp_K7) + EC_DevMatQLQ*AreaImp_K7
Load_ECM_BR7K7_MatQLQMean = apply(Load_ECM_BR7K7_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR7K7_MatQLQ05 = apply(Load_ECM_BR7K7_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR7K7_MatQLQ95 = apply(Load_ECM_BR7K7_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

plot(as.Date(Dates_BARN), Load_ECM_BR7K7, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR7K7), BR7K7_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

plot(as.Date(Dates_BARN), Load_ECM_BR7K7_QLQ, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400))
par(new = TRUE)
plot(as.Date(Dates_BR7K7), BR7K7_TrueLoad, xlim = c(as.Date('1999-01-01'), as.Date('2014-01-01')), ylim = c(0,400), col = 'red')

png('ECMBR7K7.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR7K7_TrueLoad[as.Date(Dates_BR7K7) %in% as.Date(Dates_BARN)], y = Load_ECM_BR7K7_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR7K7)], 
     xlim = c(0,50), ylim = c(0,50), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR7 K7')
arrows(BR7K7_TrueLoad[as.Date(Dates_BR7K7) %in% as.Date(Dates_BARN)], (Load_ECM_BR7K7_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR7K7)]), BR7K7_TrueLoad[as.Date(Dates_BR7K7) %in% as.Date(Dates_BARN)], (Load_ECM_BR7K7_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR7K7)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR7K7_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR7K7_TrueLoad[as.Date(Dates_BR7K7) %in% as.Date(Dates_BARN)], Load_ECM_BR7K7_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR7K7)], 
     xlim = c(0,50), ylim = c(0,50), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR7 K7')
arrows(BR7K7_TrueLoad[as.Date(Dates_BR7K7) %in% as.Date(Dates_BARN)], (Load_ECM_BR7K7_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR7K7)]), BR7K7_TrueLoad[as.Date(Dates_BR7K7) %in% as.Date(Dates_BARN)], (Load_ECM_BR7K7_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR7K7)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#  Evaluate adding other normalizers to the ECM models----
#   Fixme: Population----

#   Slope gradient----
Slopes = raster(x = paste0(dir_worldfile, '\\', f_BaismanSlope))
Slopes = projectRaster(Slopes, crs = CRS(pCRS))

#Add impervious fraction to worldfile information
world$Slope = raster::extract(x = Slopes, y = world)

#Get the slope / average slope upstream of the sampling location
world$SlopeG0 = world$SlopeG1 = world$SlopeG2 = world$SlopeG3 = world$SlopeG4 = world$SlopeG5 = world$SlopeG6 = world$SlopeG7 = world$SlopeG8 = world$SlopeG9 = world$SlopeG10 = NA
world$SlopeGK0 = world$SlopeGK1 = world$SlopeGK2 = world$SlopeGK3 = world$SlopeGK4 = world$SlopeGK5 = world$SlopeGK6 = world$SlopeGK7 = world$SlopeGK8 = world$SlopeGK9 = world$SlopeGK10 = world$SlopeGK11 = NA
world$SlopeG0[which(!is.na(world$G_0))] = world$Slope[which(!is.na(world$G_0))]/mean(world$Slope[which(!is.na(world$G_0))])
world$SlopeG1[which(!is.na(world$G_1))] = world$Slope[which(!is.na(world$G_1))]/mean(world$Slope[which(!is.na(world$G_1))])
world$SlopeG2[which(!is.na(world$G_2))] = world$Slope[which(!is.na(world$G_2))]/mean(world$Slope[which(!is.na(world$G_2))])
world$SlopeG3[which(!is.na(world$G_3))] = world$Slope[which(!is.na(world$G_3))]/mean(world$Slope[which(!is.na(world$G_3))])
world$SlopeG4[which(!is.na(world$G_4))] = world$Slope[which(!is.na(world$G_4))]/mean(world$Slope[which(!is.na(world$G_4))])
world$SlopeG5[which(!is.na(world$G_5))] = world$Slope[which(!is.na(world$G_5))]/mean(world$Slope[which(!is.na(world$G_5))])
world$SlopeG6[which(!is.na(world$G_6))] = world$Slope[which(!is.na(world$G_6))]/mean(world$Slope[which(!is.na(world$G_6))])
world$SlopeG7[which(!is.na(world$G_7))] = world$Slope[which(!is.na(world$G_7))]/mean(world$Slope[which(!is.na(world$G_7))])
world$SlopeG8[which(!is.na(world$G_8))] = world$Slope[which(!is.na(world$G_8))]/mean(world$Slope[which(!is.na(world$G_8))])
world$SlopeG9[which(!is.na(world$G_9))] = world$Slope[which(!is.na(world$G_9))]/mean(world$Slope[which(!is.na(world$G_9))])
world$SlopeG10[which(!is.na(world$G_10))] = world$Slope[which(!is.na(world$G_10))]/mean(world$Slope[which(!is.na(world$G_10))])
world$SlopeGK0[which(!is.na(world$GK_0))] = world$Slope[which(!is.na(world$GK_0))]/mean(world$Slope[which(!is.na(world$GK_0))])
world$SlopeGK1[which(!is.na(world$GK_1))] = world$Slope[which(!is.na(world$GK_1))]/mean(world$Slope[which(!is.na(world$GK_1))])
world$SlopeGK2[which(!is.na(world$GK_2))] = world$Slope[which(!is.na(world$GK_2))]/mean(world$Slope[which(!is.na(world$GK_2))])
world$SlopeGK3[which(!is.na(world$GK_3))] = world$Slope[which(!is.na(world$GK_3))]/mean(world$Slope[which(!is.na(world$GK_3))])
world$SlopeGK4[which(!is.na(world$GK_4))] = world$Slope[which(!is.na(world$GK_4))]/mean(world$Slope[which(!is.na(world$GK_4))])
world$SlopeGK5[which(!is.na(world$GK_5))] = world$Slope[which(!is.na(world$GK_5))]/mean(world$Slope[which(!is.na(world$GK_5))])
world$SlopeGK6[which(!is.na(world$GK_6))] = world$Slope[which(!is.na(world$GK_6))]/mean(world$Slope[which(!is.na(world$GK_6))])
world$SlopeGK7[which(!is.na(world$GK_7))] = world$Slope[which(!is.na(world$GK_7))]/mean(world$Slope[which(!is.na(world$GK_7))])
world$SlopeGK8[which(!is.na(world$GK_8))] = world$Slope[which(!is.na(world$GK_8))]/mean(world$Slope[which(!is.na(world$GK_8))])
world$SlopeGK9[which(!is.na(world$GK_9))] = world$Slope[which(!is.na(world$GK_9))]/mean(world$Slope[which(!is.na(world$GK_9))])
world$SlopeGK10[which(!is.na(world$GK_10))] = world$Slope[which(!is.na(world$GK_10))]/mean(world$Slope[which(!is.na(world$GK_10))])
world$SlopeGK11[which(!is.na(world$GK_11))] = world$Slope[which(!is.na(world$GK_11))]/mean(world$Slope[which(!is.na(world$GK_11))])

#    BR3 S1----
#Get the impervious surface area contributing to this location
AreaImpSlp_S1 = 0
AreaSlp_S1 = 0
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$G_10 == 1))){
  AreaImpSlp_S1 = AreaImpSlp_S1 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$G_10 == 1)[i]]*res^2*world$SlopeG10[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$G_10 == 1)[i]]
  AreaSlp_S1 = AreaSlp_S1 + (1-world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$G_10 == 1)[i]])*res^2*world$SlopeG10[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$G_10 == 1)[i]]
}
rm(i)

Load_ECM_BR3S1Slp_Mat = EC_UndevMat*(AreaSlp_S1) + EC_DevMat*AreaImpSlp_S1
Load_ECM_BR3S1Slp_MatMean = apply(Load_ECM_BR3S1Slp_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR3S1Slp_Mat05 = apply(Load_ECM_BR3S1Slp_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S1Slp_Mat95 = apply(Load_ECM_BR3S1Slp_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR3S1Slp_MatQLQ = EC_UndevMatQLQ*(AreaSlp_S1) + EC_DevMatQLQ*AreaImpSlp_S1
Load_ECM_BR3S1Slp_MatQLQMean = apply(Load_ECM_BR3S1Slp_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR3S1Slp_MatQLQ05 = apply(Load_ECM_BR3S1Slp_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR3S1Slp_MatQLQ95 = apply(Load_ECM_BR3S1Slp_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

png('ECMBR3S1Slp.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR3S1_TrueLoad[as.Date(Dates_BR3S1) %in% as.Date(Dates_BARN)], y = Load_ECM_BR3S1Slp_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S1)], 
     xlim = c(0,30), ylim = c(0,30), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S1')
arrows(BR3S1_TrueLoad[as.Date(Dates_BR3S1) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S1Slp_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S1)]), BR3S1_TrueLoad[as.Date(Dates_BR3S1) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S1Slp_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S1)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR3S1Slp_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR3S1_TrueLoad[as.Date(Dates_BR3S1) %in% as.Date(Dates_BARN)], Load_ECM_BR3S1Slp_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S1)], 
     xlim = c(0,30), ylim = c(0,30), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR3 S1')
arrows(BR3S1_TrueLoad[as.Date(Dates_BR3S1) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S1Slp_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S1)]), BR3S1_TrueLoad[as.Date(Dates_BR3S1) %in% as.Date(Dates_BARN)], (Load_ECM_BR3S1Slp_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR3S1)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#    BR5 K5----
#Get the impervious surface area contributing to this location
AreaImpSlp_K5 = 0
AreaSlp_K5 = 0
for (i in 1:length(which(world[which(duplicated(world$patchID) == FALSE),]$GK_3 == 1))){
  AreaImpSlp_K5 = AreaImpSlp_K5 + world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_3 == 1)[i]]*res^2*world$SlopeGK3[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_3 == 1)[i]]
  AreaSlp_K5 = AreaSlp_K5 + (1-world$ImpFrac[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_3 == 1)[i]])*res^2*world$SlopeGK3[which(duplicated(world$patchID) == FALSE)][which(world[which(duplicated(world$patchID) == FALSE),]$GK_3 == 1)[i]]
}
rm(i)

Load_ECM_BR5K5Slp_Mat = EC_UndevMat*(AreaSlp_K5) + EC_DevMat*AreaImpSlp_K5
Load_ECM_BR5K5Slp_MatMean = apply(Load_ECM_BR5K5Slp_Mat, MARGIN = 1, FUN = mean)
Load_ECM_BR5K5Slp_Mat05 = apply(Load_ECM_BR5K5Slp_Mat, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR5K5Slp_Mat95 = apply(Load_ECM_BR5K5Slp_Mat, MARGIN = 1, FUN = quantile, probs = 0.95)
Load_ECM_BR5K5Slp_MatQLQ = EC_UndevMatQLQ*(AreaSlp_K5) + EC_DevMatQLQ*AreaImpSlp_K5
Load_ECM_BR5K5Slp_MatQLQMean = apply(Load_ECM_BR5K5Slp_MatQLQ, MARGIN = 1, FUN = mean)
Load_ECM_BR5K5Slp_MatQLQ05 = apply(Load_ECM_BR5K5Slp_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.05)
Load_ECM_BR5K5Slp_MatQLQ95 = apply(Load_ECM_BR5K5Slp_MatQLQ, MARGIN = 1, FUN = quantile, probs = 0.95)

png('ECMBR5K5Slp.png', res = 300, units = 'in', width = 5, height = 5)
plot(x = BR5K5_TrueLoad[as.Date(Dates_BR5K5) %in% as.Date(Dates_BARN)], y = Load_ECM_BR5K5Slp_MatMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5)], 
     xlim = c(0,20), ylim = c(0,20), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR5 K5')
arrows(BR5K5_TrueLoad[as.Date(Dates_BR5K5) %in% as.Date(Dates_BARN)], (Load_ECM_BR5K5Slp_Mat05[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5)]), BR5K5_TrueLoad[as.Date(Dates_BR5K5) %in% as.Date(Dates_BARN)], (Load_ECM_BR5K5Slp_Mat95[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col = 'red')
dev.off()

png('ECMBR5K5Slp_QLQ.png', res = 300, units = 'in', width = 5, height = 5)
plot(BR5K5_TrueLoad[as.Date(Dates_BR5K5) %in% as.Date(Dates_BARN)], Load_ECM_BR5K5Slp_MatQLQMean[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5)], 
     xlim = c(0,20), ylim = c(0,20), xlab = 'True TN Load', ylab = 'ECM Predicted TN Load', main = 'BR5 K5')
arrows(BR5K5_TrueLoad[as.Date(Dates_BR5K5) %in% as.Date(Dates_BARN)], (Load_ECM_BR5K5Slp_MatQLQ05[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5)]), BR5K5_TrueLoad[as.Date(Dates_BR5K5) %in% as.Date(Dates_BARN)], (Load_ECM_BR5K5Slp_MatQLQ95[as.Date(Dates_BARN) %in% as.Date(Dates_BR5K5)]), length=0.05, angle=90, code=3)
lines(c(0,200), c(0,200), col ='red')
dev.off()

#  Evaluate temporal correlation in ECM TN load - true TN load residuals----
setwd(wd_BESN)
acf(x = Load_ECM_Baisman_MatMean - BARN_PredTN$TrueLoad, na.action = na.pass)

#Variogram method
#Make the variogram cloud matrix
Dists = Vals = matrix(NA, nrow = length(BARN_PredTN$TrueLoad[!is.na(BARN_PredTN$TrueLoad)]), ncol = length(BARN_PredTN$TrueLoad[!is.na(BARN_PredTN$TrueLoad)]))
for (i in 1:length(BARN_PredTN$TrueLoad[!is.na(BARN_PredTN$TrueLoad)])){
  #Get the distance of this point to all other points
  Dists[i,] = abs(as.numeric(Dates_BARN[!is.na(BARN_PredTN$TrueLoad)] - Dates_BARN[!is.na(BARN_PredTN$TrueLoad)][i]))
  #Get the squared difference of this point and all other points, and add to matrix
  Vals[i,] = 0.5*((Load_ECM_Baisman_MatMean[!is.na(BARN_PredTN$TrueLoad)] - BARN_PredTN$TrueLoad[!is.na(BARN_PredTN$TrueLoad)]) - (Load_ECM_Baisman_MatMean[!is.na(BARN_PredTN$TrueLoad)][i] - BARN_PredTN$TrueLoad[!is.na(BARN_PredTN$TrueLoad)][i]))^2
} 

#Estimate bin lags
h = seq(0,2500,100)
Bins = n = vector('numeric', length=length(h))
tol=50

for (j in 1:length(h)){
  #Assign the value of semivariance
  n[j] = length(which(Dists[upper.tri(Dists)] >= (h[j]-tol) & Dists[upper.tri(Dists)] <= (h[j]+tol)))
  Bins[j] = sum(Vals[upper.tri(Vals)][which((Dists[upper.tri(Dists)] >= (h[j]-tol)) & (Dists[upper.tri(Dists)] <= (h[j]+tol)))])/n[j]
}

#Plot variogram cloud and bin means
#Real space
png('Variogram_RealSpace_BARNecm.png', res = 300, units = 'in', width = 5, height = 5)
par(mar=c(5,5,2,1))
plot(x = Dists[upper.tri(Dists)], y = Vals[upper.tri(Vals)], xlab='Time Distance [days]', ylab=expression(Semivariance ~ (mg/s)^2), cex.lab=1.5, cex.axis=1.5,
     ylim=c(0, 1000), xlim=c(0,6000), main = '100 day lag means, y-axis trimmed')
par(new=TRUE)
plot(h, Bins, pch=16, col='blue', ylim=c(0, 1000), xlim=c(0,6000), axes=FALSE, xlab='', ylab='')
legend('bottomright', legend = c('Variogram Cloud', 'Bin Estimates'), pch=c(1,16), col=c('black', 'blue'))
dev.off()

#log space
png('Variogram_LogSpace_BARNecm.png', res = 300, units = 'in', width = 5, height = 5)
par(mar=c(5,5,2,1))
plot(x = Dists[upper.tri(Dists)], y = Vals[upper.tri(Vals)], xlab='Time Distance [days]', ylab=expression(Semivariance ~ (mg/s)^2), cex.lab=1.5, cex.axis=1.5,
     ylim=c(0.001, 100000), xlim=c(0,6000), log = 'y', main = '100 day lag means')
par(new=TRUE)
plot(h, Bins, pch=16, col='blue', ylim=c(0.001, 100000), xlim=c(0,6000), axes=FALSE, xlab='', ylab='', log = 'y')
legend('bottomright', legend = c('Variogram Cloud', 'Bin Estimates'), pch=c(1,16), col=c('black', 'blue'))
dev.off()

#Smaller time distances
h = seq(0,500,20)
Bins = n = vector('numeric', length=length(h))
tol=10

for (j in 1:length(h)){
  #Assign the value of semivariance
  n[j] = length(which(Dists[upper.tri(Dists)] >= (h[j]-tol) & Dists[upper.tri(Dists)] <= (h[j]+tol)))
  Bins[j] = sum(Vals[upper.tri(Vals)][which((Dists[upper.tri(Dists)] >= (h[j]-tol)) & (Dists[upper.tri(Dists)] <= (h[j]+tol)))])/n[j]
}

png('Variogram_LogSpace_BARNecm_500days.png', res = 300, units = 'in', width = 5, height = 5)
par(mar=c(5,5,2,1))
plot(x = Dists[upper.tri(Dists)], y = Vals[upper.tri(Vals)], xlab='Time Distance [days]', ylab=expression(Semivariance ~ (mg/s)^2), cex.lab=1.5, cex.axis=1.5,
     ylim=c(0, 1000), xlim=c(0,500), main = '20 day lags')
par(new=TRUE)
plot(h, Bins, pch=16, col='blue', ylim=c(0, 1000), xlim=c(0,500), axes=FALSE, xlab='', ylab='')
legend('bottomright', legend = c('Variogram Cloud', 'Bin Estimates'), pch=c(1,16), col=c('black', 'blue'))
dev.off()

#  Save ECM Timeseries for undeveloped and developed land----
write.csv(x = EC_UndevMat, file = f_EC_Undev, row.names = FALSE)
write.csv(x = EC_DevMat, file = f_EC_Dev, row.names = FALSE)

#Save updated site water chem data----
writeOGR(obj = K_Sites, dsn = wd_BESN, layer = f_KenworthSites_p, driver = "ESRI Shapefile")
write.csv(x = K_Q, file = f_KenworthWQ_Q_p, row.names = FALSE)
write.csv(x = K_TN, file = f_KenworthWQ_TN_p, row.names = FALSE)
writeOGR(obj = S_Sites, dsn = wd_BESN, layer = f_SmithSites_p, driver = "ESRI Shapefile")
write.csv(x = S_Q, file = f_SmithWQ_Q_p, row.names = FALSE)
write.csv(x = S_TN, file = f_SmithWQ_TN_p, row.names = FALSE)

#Save updated worldfile csv----
writeOGR(obj = world, dsn = wd_BRPOBR, layer = f_worldfile_p, driver = "ESRI Shapefile")

#extra plots----
# plot(POBR_PredTN_POBRWRTDS$MedLoad, POBR_PredTN_POBRWRTDS$TrueLoad, log = 'xy', xlim = c(0.001, 1000), ylim = c(0.001, 1000), xlab = 'Estimated TN Load (kg N/s)', ylab = 'True TN Load (kg N/s)')
# par(new = T)
# plot(POBR_PredTN_POBRWRTDS$MedLoad[POBR_PredTN_POBRWRTDS$LowTNLim == 1], POBR_PredTN_POBRWRTDS$TrueLoad[POBR_PredTN_POBRWRTDS$LowTNLim == 1], log = 'xy', xlim = c(0.001, 1000), ylim = c(0.001, 1000), col ='blue', xlab = '', ylab = '', axes = FALSE)
# par(new = T)
# plot(POBR_PredTN_POBRWRTDS$MedLoad[POBR_PredTN_POBRWRTDS$LowTNLim == 2], POBR_PredTN_POBRWRTDS$TrueLoad[POBR_PredTN_POBRWRTDS$LowTNLim == 2], log = 'xy', xlim = c(0.001, 1000), ylim = c(0.001, 1000), col = 'red', xlab = '', ylab = '', axes = FALSE)
# lines(c(0.00001,10000), c(0.00001,10000))
# legend('bottomright', title = 'Detection Limits', legend = c('0.01', '0.05'), pch = 1, col = c('blue', 'red'))
# 
# plot(BARN_PredTN$MedLoad, BARN_PredTN$TrueLoad, log = 'xy', xlim = c(0.001, 1000), ylim = c(0.001, 1000), xlab = 'Estimated TN Load (kg N/s)', ylab = 'True TN Load (kg N/s)')
# lines(c(0.00001,10000), c(0.00001,10000))
