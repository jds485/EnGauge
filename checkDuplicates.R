#Function to check for and remove duplicate records, as recommended in USGS publications.
checkDuplicatesAndRemove = function(StationList, #List of station timeseries
                                    colNames = NA #an optional character vector of column names to check strictly for duplication
                                    ){
  for (i in 1:length(StationList)){
    #Check for number of unique records equal to the total length of records.
    if(nrow(unique(StationList[[i]])) < nrow(StationList[[i]])){
      #There are duplicate records. Get a dataset without duplicates.
      StationList[[i]] = StationList[[i]][-which(duplicated(StationList[[i]]) == TRUE),]
    }
    #Now there are no exact duplicate records. Check the specific columns that the user requested for duplicates.
    if (any(!is.na(colNames))){
      inds = which(duplicated(StationList[[i]][ , colNames]) == TRUE)
      if (length(inds) > 0){
        StationList[[i]] = StationList[[i]][-inds,]
      }
    }
  }
  
  return(StationList)
}