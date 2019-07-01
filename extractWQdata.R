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


