#Identify missing dates and fill them into the timeseries----
#Fixme: make parallelization general for different operating systems

FillMissingDates = function(Dataset, #Where the missing date information will be added as counts 
                            StationList, #list of streamflow stations to check for missing dates
                            Var = 'X_00060_00003', #variable in StationList to check for missing dates
                            Date = 'Date',  #Name of the date variable in StationList
                            gapType = NA,  #length of the desired gap in the StationList dataset: d - day, m - month, a - annual
                            site_no_D = 'site_no', #column name containing the site number for Dataset
                            site_no_SL = 'site_no', #column name containing the site number for StationList
                            NoNAcols #character vector of column names that should not be made NA
){
  if (is.na(gapType)){
    stop('Error: gapType must be specified.')
  }
  #Check that Dataset and StationList are the same length
  if(length(Dataset) != length(StationList)){
    stop('Error: Dataset and StationList must be the same length and have the same stations.')
  }
  
  #Add a missing data column depending on the type of gap specified
  if (gapType == 'd'){
    Dataset$MissingData_d = NA
  }else if (gapType == 'm'){
    Dataset$MissingData_m = NA
  }else if (gapType == 'a'){
    Dataset$MissingData_a = NA
  }
  
  for (i in 1:length(StationList)){
    #Sort the streamflow series by date
    StationList[[i]] = StationList[[i]][order(StationList[[i]][ , Date]),]
    
    #Missing data that are reported as NA cells
    NAs = length(which(is.na(StationList[[i]][,Var])))
    
    #Add also those missing data resulting from gaps > gapType in the record
    # NOTE: this assumes data are provided in a format such that sequential records are already at the gap desired (e.g. 1 day for daily)
    if (gapType == 'd'){
      agaps = as.numeric(StationList[[i]][-1, Date]) - as.numeric(StationList[[i]][-nrow(StationList[[i]]), Date])
      Dataset$MissingData_d[which(Dataset[, site_no_D]@data == StationList[[i]][1, site_no_SL])] = NAs + sum(agaps[which(agaps > 1)]-1)
      
      #Fill in the missing data dates with NA values to have a complete time series for all records
      Inds = which(agaps > 1)
      if (length(Inds) > 0){
        for (j in 1:length(Inds)){
          NumNAs = as.numeric(StationList[[i]][(Inds[j]+1), Date] - StationList[[i]][Inds[j], Date]) - 1
          DateNAs = seq(StationList[[i]][Inds[j], Date]+1, StationList[[i]][Inds[j]+1, Date]-1, 1)
          #Assign NumNAs new rows to this dataframe with NA streamflow
          for (k in 1:NumNAs){
            r = StationList[[i]][1,]
            r[,Date] = DateNAs[k]
            #Set selected columns to NA
            r[,-which(colnames(r) %in% c(Date, site_no_SL, NoNAcols))] = NA
            StationList[[i]] = rbind(StationList[[i]], r)
          }
        }
      }
    }else if (gapType == 'm'){
      #Get theoretical gaps for each month, given the starting location
      tgaps = as.Date(StationList[[i]][, Date]) %m+% rep(months(1), length(as.numeric(as.Date(StationList[[i]][, Date])))) - as.Date(StationList[[i]][, Date])
      #Last entry doesn't have a gap. Remove.
      tgaps = tgaps[-length(tgaps)]
      #Actual gap
      agaps = as.numeric(as.Date(StationList[[i]][, Date])[-1]) - as.numeric(as.Date(StationList[[i]][, Date])[-nrow(StationList[[i]])])
      
      Dataset$MissingData_m[which(Dataset[, site_no_D]@data == StationList[[i]][1, site_no_SL])] = NAs + length(agaps[which(agaps > tgaps)])
      
      #Fill in the missing data dates with NA values to have a complete time series for all records
      Inds = which(agaps > tgaps)
      if (length(Inds) > 0){
        for (j in 1:length(Inds)){
          NumNAs = length(seq(from=StationList[[i]][Inds[j], Date], to=StationList[[i]][(Inds[j]+1), Date], by='month')) - 2
          DateNAs = seq(from=StationList[[i]][Inds[j], Date], to=StationList[[i]][(Inds[j]+1), Date], by='month')[-c(1,length(seq(from=StationList[[i]][Inds[j], Date], to=StationList[[i]][(Inds[j]+1), Date], by='month')))]
          #Assign NumNAs new rows to this dataframe with NA streamflow
          for (k in 1:NumNAs){
            r = StationList[[i]][1,]
            r[,Date] = DateNAs[k]
            #Set selected columns to NA
            r[,-which(colnames(r) %in% c(Date, site_no_SL, NoNAcols))] = NA
            StationList[[i]] = rbind(StationList[[i]], r)
          }
        }
      }
    }else if (gapType == 'a'){
      #Get theoretical gaps for each year, given the starting location
      tgaps = as.Date(StationList[[i]][, Date]) + years(1) - as.Date(StationList[[i]][, Date])
      #Last entry doesn't have a gap. Remove.
      tgaps = tgaps[-length(tgaps)]
      #Actual gap
      agaps = as.numeric(as.Date(StationList[[i]][, Date])[-1]) - as.numeric(as.Date(StationList[[i]][, Date])[-nrow(StationList[[i]])])
      
      Dataset$MissingData_a[which(Dataset[, site_no_D]@data == StationList[[i]][1, site_no_SL])] = NAs + length(agaps[which(agaps > tgaps)])
      
      #Fill in the missing data dates with NA values to have a complete time series for all records
      Inds = which(agaps > tgaps)
      if (length(Inds) > 0){
        for (j in 1:length(Inds)){
          NumNAs = length(seq(from=StationList[[i]][Inds[j], Date], to=StationList[[i]][(Inds[j]+1), Date], by='years')) - 2
          DateNAs = seq(from=StationList[[i]][Inds[j], Date], to=StationList[[i]][(Inds[j]+1), Date], by='years')[-c(1,length(seq(from=StationList[[i]][Inds[j], Date], to=StationList[[i]][(Inds[j]+1), Date], by='years')))]
          #Assign NumNAs new rows to this dataframe with NA streamflow
          for (k in 1:NumNAs){
            r = StationList[[i]][1,]
            r[,Date] = DateNAs[k]
            #Set selected columns to NA
            r[,-which(colnames(r) %in% c(Date, site_no_SL, NoNAcols))] = NA
            StationList[[i]] = rbind(StationList[[i]], r)
          }
        }
      }
    }
    #Sort the streamflow series by date again because NAs were added out of chronological order
    StationList[[i]] = StationList[[i]][order(StationList[[i]][ , Date]),]
  }
  
  #return the two datasets
  return(list(Dataset = Dataset, StationList = StationList))
}

#Same function as above, in parallel for Windows
FillMissingDates_par = function(Dataset, #Where the missing date information will be added as counts 
                            StationList, #list of streamflow stations to check for missing dates
                            Var = 'X_00060_00003', #variable in StationList to check for missing dates
                            Date = 'Date',  #Name of the date variable in StationList
                            gapType = NA,  #length of the desired gap in the StationList dataset: d - day, m - month, a - annual
                            site_no_D = 'site_no', #column name containing the site number for Dataset
                            site_no_SL = 'site_no', #column name containing the site number for StationList
                            NoNAcols #character vector of column names that should not be made NA
){
  if (is.na(gapType)){
    stop('Error: gapType must be specified.')
  }
  #Check that Dataset and StationList are the same length
  if(length(Dataset) != length(StationList)){
    stop('Error: Dataset and StationList must be the same length and have the same stations.')
  }
  
  #Add a missing data column depending on the type of gap specified
  if (gapType == 'd'){
    Dataset$MissingData_d = NA
  }else if (gapType == 'm'){
    Dataset$MissingData_m = NA
  }else if (gapType == 'a'){
    Dataset$MissingData_a = NA
  }
  
  cl = makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  lst = foreach (i = 1:length(StationList), .combine = 'c', .packages = 'lubridate') %dopar%{
    #Sort the streamflow series by date
    f = StationList[[i]][order(StationList[[i]][ , Date]),]
    
    #Missing data that are reported as NA cells
    NAs = length(which(is.na(f[,Var])))
    
    #Add also those missing data resulting from gaps > gapType in the record
    # NOTE: this assumes data are provided in a format such that sequential records are already at the gap desired (e.g. 1 day for daily)
    if (gapType == 'd'){
      agaps = as.numeric(f[-1, Date]) - as.numeric(f[-nrow(f), Date])
      Dataset$MissingData_d[which(Dataset[, site_no_D]@data == f[1, site_no_SL])] = NAs + sum(agaps[which(agaps > 1)]-1)
      
      #Fill in the missing data dates with NA values to have a complete time series for all records
      Inds = which(agaps > 1)
      if (length(Inds) > 0){
        for (j in 1:length(Inds)){
          NumNAs = as.numeric(f[(Inds[j]+1), Date] - f[Inds[j], Date]) - 1
          DateNAs = seq(f[Inds[j], Date]+1, f[Inds[j]+1, Date]-1, 1)
          #Assign NumNAs new rows to this dataframe with NA streamflow
          for (k in 1:NumNAs){
            r = f[1,]
            r[,Date] = DateNAs[k]
            #Set selected columns to NA
            r[,-which(colnames(r) %in% c(Date, site_no_SL, NoNAcols))] = NA
            f = rbind(f, r)
          }
        }
      }
      #Sort the streamflow series by date again because NAs were added out of chronological order
      list(f[order(f[ , Date]),], Dataset$MissingData_d[which(Dataset[, site_no_D]@data == f[1, site_no_SL])])
    }else if (gapType == 'm'){
      #Get theoretical gaps for each month, given the starting location
      tgaps = as.Date(f[, Date]) %m+% rep(months(1), length(as.numeric(as.Date(f[, Date])))) - as.Date(f[, Date])
      #Last entry doesn't have a gap. Remove.
      tgaps = tgaps[-length(tgaps)]
      #Actual gap
      agaps = as.numeric(as.Date(f[, Date])[-1]) - as.numeric(as.Date(f[, Date])[-nrow(f)])
      
      Dataset$MissingData_m[which(Dataset[, site_no_D]@data == f[1, site_no_SL])] = NAs + length(agaps[which(agaps > tgaps)])
      
      #Fill in the missing data dates with NA values to have a complete time series for all records
      Inds = which(agaps > tgaps)
      if (length(Inds) > 0){
        for (j in 1:length(Inds)){
          NumNAs = length(seq(from=f[Inds[j], Date], to=f[(Inds[j]+1), Date], by='month')) - 2
          DateNAs = seq(from=f[Inds[j], Date], to=f[(Inds[j]+1), Date], by='month')[-c(1,length(seq(from=f[Inds[j], Date], to=f[(Inds[j]+1), Date], by='month')))]
          #Assign NumNAs new rows to this dataframe with NA streamflow
          for (k in 1:NumNAs){
            r = f[1,]
            r[,Date] = DateNAs[k]
            #Set selected columns to NA
            r[,-which(colnames(r) %in% c(Date, site_no_SL, NoNAcols))] = NA
            f = rbind(f, r)
          }
        }
      }
      #Sort the streamflow series by date again because NAs were added out of chronological order
      list(f[order(f[ , Date]),], Dataset$MissingData_m[which(Dataset[, site_no_D]@data == f[1, site_no_SL])])
    }else if (gapType == 'a'){
      #Get theoretical gaps for each year, given the starting location
      tgaps = as.Date(f[, Date]) + years(1) - as.Date(f[, Date])
      #Last entry doesn't have a gap. Remove.
      tgaps = tgaps[-length(tgaps)]
      #Actual gap
      agaps = as.numeric(as.Date(f[, Date])[-1]) - as.numeric(as.Date(f[, Date])[-nrow(f)])
      
      Dataset$MissingData_a[which(Dataset[, site_no_D]@data == f[1, site_no_SL])] = NAs + length(agaps[which(agaps > tgaps)])
      
      #Fill in the missing data dates with NA values to have a complete time series for all records
      Inds = which(agaps > tgaps)
      if (length(Inds) > 0){
        for (j in 1:length(Inds)){
          NumNAs = length(seq(from=f[Inds[j], Date], to=f[(Inds[j]+1), Date], by='years')) - 2
          DateNAs = seq(from=f[Inds[j], Date], to=f[(Inds[j]+1), Date], by='years')[-c(1,length(seq(from=f[Inds[j], Date], to=f[(Inds[j]+1), Date], by='years')))]
          #Assign NumNAs new rows to this dataframe with NA streamflow
          for (k in 1:NumNAs){
            r = f[1,]
            r[,Date] = DateNAs[k]
            #Set selected columns to NA
            r[,-which(colnames(r) %in% c(Date, site_no_SL, NoNAcols))] = NA
            f = rbind(f, r)
          }
        }
      }
      #Sort the streamflow series by date again because NAs were added out of chronological order
      list(f[order(f[ , Date]),], Dataset$MissingData_a[which(Dataset[, site_no_D]@data == f[1, site_no_SL])])
    }
  }
  stopCluster(cl)
  
  #Extract the spatial data vector
  if (gapType == 'd'){
    Dataset$MissingData_d = unlist(lst[seq(2, length(lst), 2)])
  }else if (gapType == 'm'){
    Dataset$MissingData_m = unlist(lst[seq(2, length(lst), 2)])
  }else if (gapType == 'a'){
    Dataset$MissingData_a = unlist(lst[seq(2, length(lst), 2)])
  }
  lst = lst[seq(1, length(lst)-1, 2)]
  names(lst) = names(StationList)
  
  #return the two datasets
  return(list(Dataset = Dataset, StationList = lst))
}
