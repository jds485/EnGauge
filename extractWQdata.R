#Function to extract timeseries for each variable for a water quality site and write to a txt file.
extractWQdata = function(StationList, 
                         fName #Name to append to the start of the text file
                         ){
  for (i in 1:length(StationList)){
    #Find all of the unique CharacteristicName, ResultSampleFractionText combinations
    us = unique(StationList[[i]][,c("CharacteristicName", "ResultSampleFractionText")])
    #Write separate files for each combination
    for (u in 1:nrow(us)){
      #Check for NA values in us
      if(is.na(us[u,1])){
        u_Ind1 = which(is.na(StationList[[i]][, c("CharacteristicName", "ResultSampleFractionText")[1]]))
      }else{
        u_Ind1 = which(StationList[[i]][, c("CharacteristicName", "ResultSampleFractionText")[1]] == us[u,1])
      }
      if(is.na(us[u,2])){
        u_Ind2 = which(is.na(StationList[[i]][, c("CharacteristicName", "ResultSampleFractionText")[2]]))
      }else{
        u_Ind2 = which(StationList[[i]][, c("CharacteristicName", "ResultSampleFractionText")[2]] == us[u,2])
      }
      #Gather only the indices for this unique combination
      u_Ind = u_Ind1[u_Ind1 %in% u_Ind2]
      
      #Check that ActivityMediaName = Water for stream-only nitrogen data
      w_Ind = which(StationList[[i]][u_Ind,"ActivityMediaName"] != "Water")
      if (length(w_Ind) > 0){
        print(paste('Some of the samples are taken in media ', unique(StationList[[i]]$ActivityMediaName), ' for station ', StationList[[i]]$MonitoringLocationIdentifier[1]))
      }
      
      #place the timeseries in chronological order
      co = StationList[[i]][u_Ind,][order(StationList[[i]][u_Ind, 'ActivityStartDate']),]
      
      #Note: subbing ~ in for / because files will not save with / in the name.
      write.table(co, 
                  paste0(getwd(), '/', fName, '_', StationList[[i]]$MonitoringLocationIdentifier[1], "_cn", gsub(pattern = "/", x = us[u,1], replacement = "~", fixed = TRUE), '_rt', us[u,2], ".txt"), 
                  sep = "\t")
    }
  }
}

#Function to aggregate a specific water quality variable into a list containing elements (data.frame) for each sampling site 
makeWQStationList = function(pattern, #the grep pattern to search for in the files
                             wd,       #the directory to search for files
                             Sites #spatial dataset containing the sites
                           ){
  od = getwd()
  setwd(wd)
  #Gather the records for each gauge into a list of dataframes
  StationList = list()
  #Find all of the Nitrogen station file indices in directory
  Ind_f = list.files()[grep(pattern = pattern, x = list.files(), ignore.case = FALSE, fixed = TRUE)]
  #Check that the length is equal to the elements in the Sites file
  if(length(Ind_f) > nrow(Sites)){
    stop('Number of station data files is greater than number of sites')
  }
  for (i in 1:length(Ind_f)){
    #Read file
    f = read.table(Ind_f[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    #Add to list
    StationList = c(StationList, list(f))
  }
  setwd(od)
  return(StationList)
}

selectWQDataType = function(wd,       #the directory to search for files
                            charName, #the CharacteristicName keyword to grep for
                            resName  #the ResultSampleFractionText keyword to grep for
                            ){
  od = getwd()
  setwd(wd)
  #get all files in working directory
  lst = list.files()
  #gather all files with the CharacteristicName for charName measurements only
  cn = grep(x = lst, pattern = charName, ignore.case = TRUE)
  #gather all files with the ResultSampleFractionText = resName
  rn = grep(x = lst[cn], pattern = resName, ignore.case = TRUE)
  #grepped files only
  f_TN = lst[cn][rn]
  
  #Make a list of all of the gathered datasets
  TN = list()
  for (i in 1:length(f_TN)){
    f = read.table(f_TN[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    #Remove any NA dates. These entries are useless.
    f = f[!is.na(f$ActivityStartDate),]
    
    #Add the full date and time string as available in each record.
    # NOTE there is a field ActivityStartDateTime that already has this, but it has NAs when the time zone or the time are NA. That's undesireable.
    f$DateTimeNoNA = str_c(f$ActivityStartDate, ' ', ifelse(is.na(f$ActivityStartTime.Time), '', f$ActivityStartTime.Time), ' ', ifelse(is.na(f$ActivityStartTime.TimeZoneCode), '', f$ActivityStartTime.TimeZoneCode))
    
    #Average date and time for those days with an end date and time listed
    f$AvgDate = (as.Date(f$ActivityEndDate) - as.Date(f$ActivityStartDate))*.5 + as.Date(f$ActivityStartDate)
    f$AvgDateTime = (as.POSIXct(f$ActivityEndDateTime) - as.POSIXct(f$ActivityStartDateTime))*.5 + as.POSIXct(f$ActivityStartDateTime)
    
    #Make a new column for the sort date.
    f$SortDate = as.Date(f$ActivityStartDate)
    f$SortDateTime = as.POSIXct(f$ActivityStartDateTime)
    #Add any averaged dates to the sort date
    f$SortDate[!is.na(f$AvgDate)] = as.Date(f$AvgDate[!is.na(f$AvgDate)])
    f$SortDateTime[!is.na(f$AvgDateTime)] = as.POSIXct(f$AvgDateTime[!is.na(f$AvgDateTime)])
    
    #Place the measurements in chronological order by the sort date and time
    # First order by time (which may have some NA values), then order by date (no NA values)
    f = f[order(as.POSIXct(f$SortDateTime)),][order(f[order(as.POSIXct(f$SortDateTime)),]$SortDate),]
    
    #Check that the units are all the same for each measurement
    us = unique(f$ResultMeasure.MeasureUnitCode)
    #remove NAs if the results with NA units are all also NA
    if(all(is.na(f$ResultMeasureValue[is.na(f$ResultMeasure.MeasureUnitCode)]))){
      us = us[!is.na(us)]
    }
    if (length(us) > 1){
      print(paste('Warning: more than one measurement unit for station ', f$MonitoringLocationIdentifier[1], '. Fix manually.'))
    }
    
    #Check for detection limits
    dl = unique(f$DetectionQuantitationLimitMeasure.MeasureValue)
    #remove NAs
    dl = dl[!is.na(dl)]
    if(length(dl) != 0){
      print(paste0('Warning: detection limits for station ', f$MonitoringLocationIdentifier[1], '. Fix manually.'))
    }
    
    #add data to list
    TN = c(TN, list(f))
    names(TN)[length(TN)] = f$MonitoringLocationIdentifier[1]
  }
  
  setwd(od)
  return(TN)
}
