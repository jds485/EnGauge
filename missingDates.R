#Identify missing dates and fill them into the timeseries----
FillMissingDates = function(Dataset, #Where the missing date information will be added as counts 
                            StationList, #list of streamflow stations to check for missing dates
                            Var = 'X_00060_00003' #variable in StationList to check for missing dates
){
  #Check that Dataset and StationList are the same length
  if(length(Dataset) != length(StationList)){
    stop('Error: Dataset and StationList must be the same length and have the same stations.')
  }
  Dataset$MissingData = NA
  for (i in 1:length(StationList)){
    #Missing data that are reported as NA cells
    Dataset$MissingData[which(Dataset$site_no == StationList[[i]]$site_no[1])] = length(which(is.na(StationList[[i]][,Var])))
    #Add to those the missing data resulting from gaps > 1 day in the record
    # NOTE: this assumes data are daily measurements
    gaps = c(with(data = StationList[[i]], as.numeric(Date[-1]) - as.numeric(Date[-nrow(StationList[[i]])])))
    
    Dataset$MissingData[which(Dataset$site_no == StationList[[i]]$site_no[1])] = Dataset$MissingData[which(Dataset$site_no == StationList[[i]]$site_no[1])] + sum(gaps[which(gaps > 1)])
    
    #Fill in the missing data dates with NA values to have a complete time series for all records
    Inds = which(gaps > 1)
    if (length(Inds) > 0){
      for (j in 1:length(Inds)){
        NumNAs = as.numeric(StationList[[i]][(Inds[j]+1),]$Date - StationList[[i]][Inds[j],]$Date) - 1
        DateNAs = seq(StationList[[i]][Inds[j],]$Date+1, StationList[[i]][Inds[j]+1,]$Date-1, 1)
        #Assign NumNAs new rows to this dataframe with NA streamflow
        for (k in 1:NumNAs){
          r = cbind(StationList[[i]][1,1:2], DateNAs[k], NA, NA)
          colnames(r) = colnames(StationList[[i]]) 
          StationList[[i]] = rbind(StationList[[i]], r)
        }
      }
    }
    #Sort the streamflow series by date
    StationList[[i]] = StationList[[i]][order(StationList[[i]]$Date),]
  }
  
  #return the two datasets
  list(Dataset = Dataset, StationList = StationList)
}