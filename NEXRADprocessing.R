#Script for processing Hydro NEXRAD data for the Baltimore area

#Set directory names----
#EnGauge repository
dir_EnGauge = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\EnGauge\\EnGauge"
#Color functions - from JDS github repo: Geothermal_ESDA
dir_ColFuns = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\USGSGauges"

#Region of interest shapefile
dir_ROI = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\BES-Watersheds-Land-Cover-Analysis"

#Precipitation
dir_precip = "C:\\Users\\js4yd\\OneDrive - University of Virginia\\BES_Data\\BES_Data\\Hydrology\\Precipitation"
dir_precipData = dir_precip

dir_Nexrad = paste0(dir_precip, "\\BES_Nexrad")
wd_NOAA = paste0(dir_precip, '\\NOAA2')

#Watersheds shapefiles
dir_Sheds = "C:\\Users\\js4yd\\Documents\\DEMtest"

#Set input filenames----
#Region of interest shapefile name
f_ROI = "BaismanGFMerge2"
#NOAA station data - all
f_NOAAstationsROI = 'NOAA_StationLocs_20km'
f_NOAAstationsDataList = 'NOAA_MetStations.yaml'
#Watersheds
f_BR = 'BaisRun5'
# Aggregated to annual, monthly, daily
f_BESPrecipAggDay = 'BES_Precip_d_AllGauges.yaml'
# Daily + site average of 2 gauges
f_BESPrecipAggDay_SiteAverage = 'BES_Precip_d_Sites.yaml'

#Set output filenames----
#Nexrad Data
f_NexradDataList = 'NexradYearList.yaml'
f_NexradDataList_Gauges = 'NexradYearList_Gauges.yaml'
f_BaismanNexradPixels_d = 'BaismanNexradPixels.csv'
f_BaismanNexradPixels_15min = 'BaismanNexradPixels_15minPrecip.csv'
f_BaismanNexradPixelsMap = 'BaismanNexradPixels'
f_BaismanNexradPixels_Gauges_d = 'GaugeNexradPixels.csv'
f_BaismanNexradPixels_Gauges_15min = 'GaugeNexradPixels_15minPrecip.csv'
f_BaismanNexradPixelsMap_Gauges = 'GaugeNexradPixels'

#Set project coordinate system----
#This is the coordinate system that all data will be plotted and written in
# It is not the coordinate system of your data (although it could be)
# EPSG codes from: https://spatialreference.org/ref/?page=2
pCRS = '+init=epsg:26918'

#Load libraries----
library(sp)
library(rgdal)
library(raster)
library(foreach)
library(doParallel)
library(stringr)
library(stringi)
library(rlist)
library(lubridate)
library(zoo)
library(psych)

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

#Load information from EnGauge downloads and place into the project coordinate system----
#Region of interest
ROI = readOGR(dsn = dir_ROI, layer = f_ROI, stringsAsFactors = FALSE)
ROI = spTransform(ROI, CRS(pCRS))

NOAAstations_locs = readOGR(dsn = wd_NOAA, layer = f_NOAAstationsROI, stringsAsFactors = FALSE)
BES_Precip_d = list.load(file = paste0(dir_precip, "\\", f_BESPrecipAggDay), type = 'YAML')
BES_Precip_Avg_d = list.load(file = paste0(dir_precip, "\\", f_BESPrecipAggDay_SiteAverage, '.RData'), type = 'RData')

MetStations = list.load(file = paste0(wd_NOAA, "\\", f_NOAAstationsDataList), type = 'YAML')

#Load watersheds----
setwd(dir_Sheds)
BaisRun_Outlet = readOGR(dsn = getwd(), layer = f_BR, stringsAsFactors = FALSE)
BaisRun_Outlet = spTransform(BaisRun_Outlet, CRS(pCRS))

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
if (any((Nexrad@coords == Nexrad2@coords) == FALSE)){
  print('Some Nexrad coordinates are not the same from one time step to the next.')
}

#Plot the location of the MD Science Center and BWI Airport gauges. 
#These were used to fill in gaps. Want to check agreement with Nexrad data
plot(Nexrad, col = 'red', pch = 15, cex = 0.3)
plot(ROI, add = T)
plot(NOAAstations_locs[NOAAstations_locs$id == "USW00093721",], add = T, col = 'blue')
plot(NOAAstations_locs[NOAAstations_locs$id == "USW00093784",], add = T, col = 'green')

#Get the pixel corresponding to the BWI and MD Science Center gauges
#Find the minimum distance from all pixels to the gauge location
NexradDistsMDSci = sqrt((Nexrad@coords[,1] - NOAAstations_locs[NOAAstations_locs$id == "USW00093784",]@coords[1,1])^2 + (Nexrad@coords[,2] - NOAAstations_locs[NOAAstations_locs$id == "USW00093784",]@coords[1,2])^2) 
NexradDistsBWI = sqrt((Nexrad@coords[,1] - NOAAstations_locs[NOAAstations_locs$id == "USW00093721",]@coords[1,1])^2 + (Nexrad@coords[,2] - NOAAstations_locs[NOAAstations_locs$id == "USW00093721",]@coords[1,2])^2) 
NexradPix_MDSci = Nexrad[which(NexradDistsMDSci == min(NexradDistsMDSci)), ]
NexradPix_BWI = Nexrad[which(NexradDistsBWI == min(NexradDistsBWI)), ]

#Plot the raster point (cell center) locations for ROI
plot(ROI)
plot(Nexrad, col = 'red', pch = 15, cex = 0.3, add = T)
plot(NexradPix_MDSci, col = 'green', pch = 15, cex = 0.3, add = T)
#BWI off screen
plot(NexradPix_BWI, col = 'blue', pch = 15, cex = 0.3, add = T)

BaispCRS = spTransform(BaisRun_Outlet, CRSobj = pCRS)
plot(BaispCRS)
plot(Nexrad, col = 'red', pch = 15, cex = 0.5, add = T)

#Extract the Baisman Run pixels within specified buffer
BaispCRS_buff = buffer(BaispCRS, width = 800)
BaisPix = Nexrad[BaispCRS_buff,]
BaisPix$Pix = seq(1, nrow(BaisPix@data), 1)

# Map of which pixel is which----
scaleRange = c(1,nrow(BaisPix@data))
scaleBy = 1
Pal = colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue'))((scaleRange[2] - scaleRange[1])/scaleBy)

setwd(dir_Nexrad)
png('BaismanNexradPixelsMap.png', res = 300, width = 5, height = 5, units = 'in')
plot(BaispCRS_buff, col = 'white', border = 'white')
plot(BaispCRS, add = TRUE)
plot(Nexrad, pch = 15, add = TRUE)
plot(BaisPix, col = colFun(BaisPix$Pix), pch = 15, add = TRUE)
#add labels
text(x = BaisPix@coords[,1]+100, y = BaisPix@coords[,2], BaisPix$Pix)
dev.off()

rm(scaleRange, scaleBy, Pal)

# Loop over the available NEXRAD txt files and extract the pixels for Baisman Run----
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
    f_BaisPix = f[BaispCRS_buff,]
    rm(f)
    
    #Extract the date and time information
    Date = paste(paste(substr(fs[i], start = 1, stop = 4), substr(fs[i], start = 5, stop = 6), substr(fs[i], start = 7, stop = 8), sep = '-'), paste(substr(fs[i], start = 9, stop = 10), substr(fs[i], start = 11, stop = 12), sep = ':'), sep = ' ') 
    
    #Store the date and time as Posix
    Date = as.character(as.POSIXct(Date))
    
    #Store the precip info (mm/ha) and return as vector
    retvals = c(Date, f_BaisPix$precip)
  }
  
  rownames(BaisNexMat) = NULL
  colnames(BaisNexMat) = c('Date', paste0('Pix', seq(1,length(BaisPix$Pix),1)))
  
  #Save to list
  NexradList = c(NexradList, list(BaisNexMat))
  names(NexradList) = c(names(NexradList)[1:(j-1)], folsNexrad[j])
}
stopCluster(cl)
rm(fs, j, cl)

#Save the Nexrad product list
list.save(x = NexradList, file = f_NexradDataList, type = 'YAML')

#Go from a list of years to a matrix / data.frame to make a continuous timeseries
NexradMat = as.data.frame(NexradList[[1]])
NexradMat$Date = as.character(NexradMat$Date)
for(i in 2:ncol(NexradMat)){
  NexradMat[,i] = as.numeric(matrix(NexradMat[,i]))
}
rm(i)
#bind list elements
if(length(NexradList) > 1){
  for(i in 2:length(NexradList)){
    NexradMat = rbind(NexradMat, NexradList[[i]])
  }
  rm(i)
}
for(i in 2:ncol(NexradMat)){
  NexradMat[,i] = as.numeric(matrix(NexradMat[,i]))
}
rm(i)

#Convert all factor and character to numeric values
for(i in 2:ncol(NexradMat)){
  NexradMat[,i] = as.numeric(as.character(NexradMat[,i]))
}
rm(i)
rm(NexradList)

#   Check for duplicate records----
NexradDuplicates = which(duplicated(NexradMat))
#All of these are for the timestamp 2:00 - 2:45 AM. 
#I'm not sure why only a few dates have this problem. It does not affect further analysis, so leaving as is.

#   Check for negative values----
#Add new columns to store negative-removed information
for (i in 2:ncol(NexradMat)){
  NexradMat = cbind(NexradMat, NexradMat[,i])
}
rm(i)
colnames(NexradMat) = c(colnames(NexradMat)[1:((ncol(NexradMat)-1)/2+1)], paste0(colnames(NexradMat)[2:((ncol(NexradMat)-1)/2+1)], 'NegRm'))

#make all the negative values NAs
for (i in ((ncol(NexradMat)-1)/2 + 2):ncol(NexradMat)){
  NexradMat[NexradMat[,i] < 0, i] = NA
}
rm(i)

#   Aggregate to daily timeseries----
#Make a list of the pixels to use the aggregate Timeseries function
BaisPixList = list()
for (i in ((ncol(NexradMat)-1)/2+2):ncol(NexradMat)){
  BaisPixList = c(BaisPixList, list(NexradMat[,c(1,i)]))
}
#Add columns for the pixel number, date, and date with time
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
                                         site_no_SL = 'Pix', NoNAcols = c('Pix'), NumCores = detectCores()-1)
BaisPix_Processed = BaisNex_Precip_d2$Dataset
BaisNex_Precip_d = BaisNex_Precip_d2$StationList
rm(BaisNex_Precip_d2)

names(BaisNex_Precip_d) = paste0('Pix', seq(1,length(BaisNex_Precip_d),1))

#   Fixme: Flag all dates on record that had gaps in 15-minute recording----
#Making a label for 15 minute gaps
#for (i in 1:length(BaisNex_Precip_d)){
#  BaisNex_Precip_d[[i]]$Gap15 = NA
#  #Checking the matrix of rain data for 15 minute gaps
#  #Check that the first entry starts at "2000-01-01 00:00:00"
#  if (as.POSIXct(NexradMat$Date[1]) != as.POSIXct("2000-01-01 00:00:00")){
#    #Flag 2000-01-01 as having a gap
#    BaisNex_Precip_d[[i]]$Gap15[which(BaisNex_Precip_d[[i]]$Date == '2000-01-01')] = 1
#  }
#}

#   Fixme?: Convert units to mm/km^2 from mm/ha----
# Based on the BES website, the units are simply mm, not mm/ha. But, the readme file says mm/ha.
#1 ha = .01 km^2
#for(i in 1:length(BaisNex_Precip_d)){
#  BaisNex_Precip_d[[i]]$Precip = BaisNex_Precip_d[[i]]$Precip/.01
#}
#rm(i)

# Loop over the available NEXRAD txt files and extract the pixels for BWI and MD Science Center----
#  Parallel over all available Nexrad years----
setwd(dir_Nexrad)
cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
NexradList_Gauges = list()
BaisPix_Gauges = rbind(NexradPix_MDSci, NexradPix_BWI)
BaisPix_Gauges$Pix = c('MDSci', 'BWI')
for (j in 1:length(folsNexrad)){
  fs = list.files(paste0(getwd(), '/', folsNexrad[j]))
  #Get all of the nexrad data for this year.
  BaisNexMat_Gauges = foreach (i = 1:length(fs), .combine = rbind, .inorder = FALSE, .packages = 'sp') %dopar% {
    #Open the file and project into spatial dataframe
    f = read.table(paste0(getwd(), '/', folsNexrad[j], '/', fs[i]), header = FALSE)
    colnames(f) = c('lat', 'long', 'precip')
    coordinates(f) = c('long', 'lat')
    proj4string(f) = CRS('+init=epsg:4326')
    f = spTransform(f, CRSobj = pCRS)
    
    #Extract the BWI and MD Science Center pixels
    f_Pix = f[BaisPix_Gauges,]
    rm(f)
    
    #Extract the date and time information
    Date = paste(paste(substr(fs[i], start = 1, stop = 4), substr(fs[i], start = 5, stop = 6), substr(fs[i], start = 7, stop = 8), sep = '-'), paste(substr(fs[i], start = 9, stop = 10), substr(fs[i], start = 11, stop = 12), sep = ':'), sep = ' ') 
    
    #Store the date and time as Posix
    Date = as.character(as.POSIXct(Date))
    
    #Store the precip info (mm/ha) and return as vector
    retvals = c(Date, f_Pix$precip)
  }
  
  rownames(BaisNexMat_Gauges) = NULL
  colnames(BaisNexMat_Gauges) = c('Date', 'MDSci', 'BWI')
  
  #Save to list
  NexradList_Gauges = c(NexradList_Gauges, list(BaisNexMat_Gauges))
  names(NexradList_Gauges) = c(names(NexradList_Gauges)[1:(j-1)], folsNexrad[j])
}
stopCluster(cl)
rm(fs, j, cl)

#Save the Nexrad product list
list.save(x = NexradList_Gauges, file = f_NexradDataList_Gauges, type = 'YAML')

#Go from a list of years to a matrix / data.frame to make a continuous timeseries
#For the first year: character dates, numeric precipitation
NexradMat_Gauges = as.data.frame(NexradList_Gauges[[1]])
NexradMat_Gauges$Date = as.character(NexradMat_Gauges$Date)
for(i in 2:length(NexradMat_Gauges)){
  NexradMat_Gauges[,i] = as.numeric(matrix(NexradMat_Gauges[,i]))
}
rm(i)
if (length(NexradList_Gauges) > 1){
  #bind list elements to this matrix
  for(i in 2:length(NexradList_Gauges)){
    NexradMat_Gauges = rbind(NexradMat_Gauges, NexradList_Gauges[[i]])
  }
  rm(i)
}
#Convert precip values to numeric
for(i in 2:ncol(NexradMat_Gauges)){
  NexradMat_Gauges[,i] = as.numeric(matrix(NexradMat_Gauges[,i]))
}
rm(i)

rm(NexradList_Gauges)

#   Check for duplicate records----
NexradDuplicates_Gauges = which(duplicated(NexradMat_Gauges))
#All of these are for the timestamp 2:00 - 2:45 AM. 
#I think this is for daylight savings time.
#But, I'm not sure why only a few dates have this problem. It does not affect further analysis, so leaving as is.

#   Check for negative values----
#Add new columns to store negative-removed information
for (i in 2:ncol(NexradMat_Gauges)){
  NexradMat_Gauges = cbind(NexradMat_Gauges, NexradMat_Gauges[,i])
}
rm(i)
colnames(NexradMat_Gauges) = c(colnames(NexradMat_Gauges)[1:((ncol(NexradMat_Gauges)-1)/2+1)], paste0(colnames(NexradMat_Gauges)[2:((ncol(NexradMat_Gauges)-1)/2+1)], 'NegRm'))

#make all the negative values NAs
for (i in ((ncol(NexradMat_Gauges)-1)/2 + 2):ncol(NexradMat_Gauges)){
  NexradMat_Gauges[NexradMat_Gauges[,i] < 0, i] = NA
}
rm(i)

#   Aggregate to daily timeseries----
#Make a list of the pixels to use the aggregate Timeseries function
BaisPixList_Gauges = list()
for (i in ((ncol(NexradMat_Gauges)-1)/2+2):ncol(NexradMat_Gauges)){
  BaisPixList_Gauges = c(BaisPixList_Gauges, list(NexradMat_Gauges[,c(1,i)]))
}
rm(i)
#Add columns for the pixel number, date, and date with time
for(i in 1:length(BaisPixList_Gauges)){
  BaisPixList_Gauges[[i]]$Pix = ifelse(i == 1, "MDSci", "BWI")
  BaisPixList_Gauges[[i]]$DateTime = BaisPixList_Gauges[[i]]$Date
  colnames(BaisPixList_Gauges[[i]]) = c("Date", 'Precip', 'Pix', "DateTime")
  BaisPixList_Gauges[[i]]$Date = as.Date(BaisPixList_Gauges[[i]]$Date)
}
rm(i)
BaisNexGauges_Precip_agg = aggregateTimesteps(StationList = BaisPixList_Gauges, aggVal = 'd', aggVar = 'Precip', date = 'Date', site = 'Pix', fun = 'sum')
BaisNexGauges_Precip_d = BaisNexGauges_Precip_agg$daily
rm(BaisNexGauges_Precip_agg)

#   Fill in missing dates in the aggregated timeseries----
BaisNexGauges_Precip_d2 = FillMissingDates_par(Dataset = BaisPix_Gauges, StationList = BaisNexGauges_Precip_d, Var = 'Precip', 
                                               Date = 'Date', gapType = 'd', site_no_D = 'Pix', 
                                               site_no_SL = 'Pix', NoNAcols = c('Pix'), NumCores = detectCores()-1)
BaisPix_Gauges_Processed = BaisNexGauges_Precip_d2$Dataset
BaisNexGauges_Precip_d = BaisNexGauges_Precip_d2$StationList
rm(BaisNexGauges_Precip_d2)

names(BaisNexGauges_Precip_d) = c('MDSci', 'BWI')
#   Fixme: Flag all dates on record that had gaps in 15-minute recording----
#Making a label for 15 minute gaps
#for (i in 1:length(BaisNexGauges_Precip_d)){
#  BaisNexGauges_Precip_d[[i]]$Gap15 = NA
#  #Checking the matrix of rain data for 15 minute gaps
#  #Check that the first entry starts at "2000-01-01 00:00:00"
#  if (as.POSIXct(NexradMat_Gauges$Date[1]) != as.POSIXct("2000-01-01 00:00:00")){
#    #Flag 2000-01-01 as having a gap
#    BaisNexGauges_Precip_d[[i]]$Gap15[which(BaisNexGauges_Precip_d[[i]]$Date == '2000-01-01')] = 1
#  }
#}

# Compare the pixels with each other----
#make a matrix of the data to get a scatterplot matrix
BaisNexPrecipMat = matrix(NA, nrow = nrow(BaisNex_Precip_d[[1]]), ncol = 1+length(BaisPix))
BaisNexPrecipMat[,1] = BaisNex_Precip_d[[1]]$Date
for (i in 1:length(BaisPix)){
  BaisNexPrecipMat[,1+i] = BaisNex_Precip_d[[i]]$Precip
}
BaisNexPrecipMat = as.data.frame(BaisNexPrecipMat)
colnames(BaisNexPrecipMat) = c('Date', paste0('Pix', seq(1,length(BaisPix),1)))

#For rain gauges
BaisNexPrecipMat_Gauges = matrix(NA, nrow = nrow(BaisNexGauges_Precip_d[[1]]), ncol = 1+length(BaisPix_Gauges))
BaisNexPrecipMat_Gauges[,1] = BaisNexGauges_Precip_d[[1]]$Date
for (i in 1:length(BaisPix_Gauges)){
  BaisNexPrecipMat_Gauges[,1+i] = BaisNexGauges_Precip_d[[i]]$Precip
}
BaisNexPrecipMat_Gauges = as.data.frame(BaisNexPrecipMat_Gauges)
colnames(BaisNexPrecipMat_Gauges) = c('Date', 'MDSci', 'BWI')

#Figures for Baisman Only
png('BaismanNexradPixPrecipTimeseries.png', res = 300, height = 5, width = 5, units = 'in')
par(mar = c(5,5,2,2))
matplotDates(as.Date(BaisNexPrecipMat[,1]), BaisNexPrecipMat[,c(5,6,7,10)], type = 'l', col = c('red', 'orange', 'green', 'blue'), 
             xlab = 'Time', ylab = expression(paste('Precipitation (mm/km'^2,')')))
dev.off()

png('BaismanNexradPixPrecipScatterPlotMatrix.png', res = 300, height = 8, width = 8, units = 'in')
pairs.panels(x = log10(BaisNexPrecipMat[,c(6,7,8,11)]+.01), scale = FALSE, density = FALSE, ellipses = FALSE, smooth = FALSE, 
             digits = 3, lm = TRUE, jiggle = FALSE, rug = FALSE, cex.cor = .7, method = 'spearman')
dev.off()

png('BaismanNexradPixPrecipScatterPlotMatrix_BufferPix.png', res = 300, height = 16, width = 16, units = 'in')
pairs.panels(x = log10(BaisNexPrecipMat[,-1]+.01), scale = FALSE, density = FALSE, ellipses = FALSE, smooth = FALSE, 
             digits = 3, lm = TRUE, jiggle = FALSE, rug = FALSE, cex.cor = .7, method = 'spearman')
dev.off()

#Figures with the other rain gauge Nexrad Data
png('BaismanNexradPixPrecipScatterPlotMatrix_GaugePixels.png', res = 300, height = 8, width = 8, units = 'in')
pairs.panels(x = log10(cbind(BaisNexPrecipMat[,c(6,7,8,11)], BaisNexPrecipMat_Gauges[,-1])+.01), scale = FALSE, density = FALSE, ellipses = FALSE, smooth = FALSE, 
             digits = 3, lm = TRUE, jiggle = FALSE, rug = FALSE, cex.cor = .7, method = 'spearman')
dev.off()

# Compare to the Oregon Ridge rain gauge----
min(as.Date(BES_Precip_d$WXORDG_RG1$SortDate))
max(as.Date(BaisNex_Precip_d[[1]]$Date))
IndStart = which(as.Date(BaisNex_Precip_d[[1]]$Date) == min(as.Date(BES_Precip_d$WXORDG_RG1$SortDate)))

#Plot the rain gauge vs. this gauge
#BaisNex_Precip_d[[1]][IndStart:nrow(BaisNex_Precip_d[[1]]),]

BaisNexPrecipMat$ORGauge = NA 
BaisNexPrecipMat$ORGauge[IndStart:nrow(BaisNexPrecipMat)] = BES_Precip_Avg_d$`Oregon Ridge Park`$mean[as.Date(BES_Precip_Avg_d$`Oregon Ridge Park`$SortDate) <= max(as.Date(BaisNex_Precip_d[[1]]$Date))]

png('BaismanNexradPixPrecipScatterPlotMatrix_WithRainGauge.png', res = 300, height = 8, width = 8, units = 'in')
pairs.panels(x = log10(BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),c(6,7,8,11,14)][which((BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),6] > 0) & !is.nan(BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),14])),]), scale = FALSE, density = FALSE, ellipses = FALSE, smooth = FALSE, 
             digits = 3, lm = TRUE, jiggle = FALSE, rug = FALSE, cex.cor = .7, method = 'spearman', xlim = c(-2.5,2.5), ylim = c(-2.5, 2.5))
dev.off()

png('BaismanNexradPixPrecipScatterPlotMatrix_BufferPix_WithRainGauge.png', res = 300, height = 8, width = 8, units = 'in')
pairs.panels(x = log10(BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),-1][which((BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),6] > 0) & (BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),3] > 0) & !is.nan(BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),14])),]), 
             scale = FALSE, density = FALSE, ellipses = FALSE, smooth = FALSE, 
             digits = 3, lm = TRUE, jiggle = FALSE, rug = FALSE, cex.cor = .7, method = 'spearman', xlim = c(-2.5,2.5), ylim = c(-2.5, 2.5))
dev.off()

#With pixels from MD Sci and BWI gauges
png('BaismanNexradPixPrecipScatterPlotMatrix_GaugePixels_WithRainGauge.png', res = 300, height = 8, width = 8, units = 'in')
pairs.panels(x = log10(cbind(BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),c(6,7,8,11,14)][which((BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),6] > 0) & !is.nan(BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),14])),], 
                             BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat), -1][which((BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),6] > 0) & !is.nan(BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),14])),])+.0001), scale = FALSE, density = FALSE, ellipses = FALSE, smooth = FALSE, 
             digits = 3, lm = TRUE, jiggle = FALSE, rug = FALSE, cex.cor = .7, method = 'spearman', xlim = c(-2.5,2.5), ylim = c(-2.5, 2.5))
dev.off()

# Compare to the BWI rain gauge----
min(as.Date(MetStations$USW00093721$date))
min(as.Date(BaisNexGauges_Precip_d[[2]]$Date))
IndStart = which(as.Date(BaisNexGauges_Precip_d[[2]]$Date) == min(as.Date(BaisNexGauges_Precip_d[[2]]$Date)))
IndStart_Met = which(as.Date(MetStations$USW00093721$date) == min(as.Date(BaisNexGauges_Precip_d[[2]]$Date)))

BaisNexPrecipMat_Gauges$BWIGauge = NA 
BaisNexPrecipMat_Gauges$BWIGauge[IndStart:nrow(BaisNexPrecipMat_Gauges)] = MetStations$USW00093721$prcp[IndStart_Met:length(MetStations$USW00093721$prcp)][(as.Date(MetStations$USW00093721$date[IndStart_Met:length(MetStations$USW00093721$prcp)]) <= max(as.Date(BaisNexGauges_Precip_d[[1]]$Date)))]/10

#With pixels from MD Sci and BWI gauges
png('BaismanNexradPixPrecipScatterPlotMatrix_GaugePixels_BWIRainGauge.png', res = 300, height = 8, width = 8, units = 'in')
pairs.panels(x = log10(cbind(BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),c(6,7,8,11)][which((BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),6] > 0) & (!is.nan(BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges),3])) & (BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges),3] > 0) & (BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges), 4] > 0) & (BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges), 2] > 0)),], 
                             BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat), -1][which((BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),6] > 0) & (!is.nan(BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges),3])) & (BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges),3] > 0) & (BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges), 4] > 0) & (BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges), 2] > 0)),])), scale = FALSE, density = FALSE, ellipses = FALSE, smooth = FALSE, 
             digits = 3, lm = FALSE, jiggle = FALSE, rug = FALSE, cex.cor = .7, method = 'spearman', xlim = c(-2.5,2.5), ylim = c(-2.5, 2.5))
dev.off()

# Compare to the MD Science Center rain gauge----
min(as.Date(MetStations$USW00093784$date))
min(as.Date(BaisNexGauges_Precip_d[[2]]$Date))
IndStart = which(as.Date(BaisNexGauges_Precip_d[[1]]$Date) == min(as.Date(BaisNexGauges_Precip_d[[1]]$Date)))
IndStart_Met = which(as.Date(MetStations$USW00093784$date) == min(as.Date(BaisNexGauges_Precip_d[[1]]$Date)))

BaisNexPrecipMat_Gauges$MDSciGauge = NA 
BaisNexPrecipMat_Gauges$MDSciGauge[IndStart:nrow(BaisNexPrecipMat_Gauges)] = MetStations$USW00093784$prcp[IndStart_Met:length(MetStations$USW00093784$prcp)][(as.Date(MetStations$USW00093784$date[IndStart_Met:length(MetStations$USW00093784$prcp)]) <= max(as.Date(BaisNexGauges_Precip_d[[1]]$Date)))]/10

#With pixels from MD Sci and BWI gauges
png('BaismanNexradPixPrecipScatterPlotMatrix_GaugePixels_MDSciRainGauge.png', res = 300, height = 8, width = 8, units = 'in')
pairs.panels(x = log10(cbind(BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),c(6,7,8,11)][which((BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),6] > 0) & (!is.nan(BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges),3])) & (BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges),3] > 0) & (BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges), 5] > 0) & (BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges), 2] > 0)),], 
                             BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat), -c(1,4)][which((BaisNexPrecipMat[IndStart:nrow(BaisNexPrecipMat),6] > 0) & (!is.nan(BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges),3])) & (BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges),3] > 0) & (BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges), 5] > 0) & (BaisNexPrecipMat_Gauges[IndStart:nrow(BaisNexPrecipMat_Gauges), 2] > 0)),])), scale = FALSE, density = FALSE, ellipses = FALSE, smooth = FALSE, 
             digits = 3, lm = FALSE, jiggle = FALSE, rug = FALSE, cex.cor = .7, method = 'spearman', xlim = c(-2.5,2.5), ylim = c(-2.5, 2.5))
dev.off()

# Fixme: Fill in rasters with gauge values on days without radar information----
# Save the Nexrad data----
BaisNexPrecipMat$Date = as.Date(BaisNexPrecipMat$Date)
write.csv(BaisNexPrecipMat, file = f_BaismanNexradPixels_d, row.names = FALSE)
write.csv(NexradMat, file = f_BaismanNexradPixels_15min, row.names = FALSE)
writeOGR(BaisPix, dsn = getwd(), layer = f_BaismanNexradPixelsMap, driver = 'ESRI Shapefile')

BaisNexPrecipMat_Gauges$Date = as.Date(BaisNexPrecipMat_Gauges$Date)
write.csv(BaisNexPrecipMat_Gauges, file = f_BaismanNexradPixels_Gauges_d, row.names = FALSE)
write.csv(NexradMat_Gauges, file = f_BaismanNexradPixels_Gauges_15min, row.names = FALSE)
writeOGR(BaisPix_Gauges, dsn = getwd(), layer = f_BaismanNexradPixelsMap_Gauges, driver = 'ESRI Shapefile')

# Fixme: Convert to a raster for use in RHESSys patch precip assignment----
#gridded(Nexrad) = TRUE
#gridded(Nexrad2) = TRUE
#NexradRas = raster(Nexrad)
#NexradRas = projectRaster(NexradRas, crs = CRS(pCRS))

#NexradRas2 = raster(Nexrad2)
#NexradRas2 = projectRaster(NexradRas2, crs = CRS(pCRS))
