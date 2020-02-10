#Check for zeros and negative values in records

#Fixme: Add an option to interpolate for negative values

checkZerosNegs = function(StationList, #List of station records in dataframes
                          Var = 'X_00060_00003', #variable to check in the StationList
                          NegReplace = NA, #Value to replace negatives with. NA is default
                          ZeroReplace = 0  #Value to replace zeros with. 0 is the default (no replacement)
){
  for (i in 1:length(StationList)){
    #Add an indicator column for zeros and negative values
    StationList[[i]]$Neg = StationList[[i]]$Zero = NA
    
    #Check if there are records with values 0 or below
    if(any(StationList[[i]][,Var] <= 0, na.rm = TRUE)){
      #There are records with zero or negative values.
      #Find the negatives. These will be assumed NA values and replaced.
      IndNeg = which(StationList[[i]][,Var] < 0)
      #Add these values to the Neg column
      StationList[[i]]$Neg[IndNeg] = StationList[[i]][IndNeg,Var]
      
      #Replace with NegReplace in the StationList
      StationList[[i]][IndNeg,Var] = NegReplace
      
      #Find the zeros
      IndZero = which(StationList[[i]][,Var] == 0)
      #Add these values to the Zero column
      StationList[[i]]$Zero[IndZero] = StationList[[i]][IndZero,Var]
      
      #Replace with ZeroReplace in the StationList
      StationList[[i]][IndZero,Var] = ZeroReplace
    }
  }
  return(StationList)
}

#Add indicators for negatives and zeros to spatial dataset. Useful for plotting purposes.
addNegsToSpatialDataset = function(StationList, 
                                   SpatialDataset,
                                   site_D, #site ID name in Dataset
                                   site_SL #site ID name in StationList
                                   ){
  SpatialDataset$Neg = NA
  for (i in 1:length(StationList)){
    if(length(which(!is.na(StationList[[i]]$Neg))) > 0){
      print(paste('Number of negatives for station', StationList[[i]][1, site_SL], '=', length(which(!is.na(StationList[[i]]$Neg)))))
      SpatialDataset$Neg[SpatialDataset[, site_D]@data == StationList[[i]][1, site_SL]] = length(which(!is.na(StationList[[i]]$Neg)))
    }
  }
  return(SpatialDataset)
}

addZerosToSpatialDataset = function(StationList, 
                                    SpatialDataset,
                                    site_D, #site ID name in Dataset
                                    site_SL #site ID name in StationList
                                    ){
  SpatialDataset$Zero = NA
  for (i in 1:length(StationList)){
    if(length(which(!is.na(StationList[[i]]$Zero))) > 0){
      print(paste('Number of zeros for station', StationList[[i]][1, site_SL], '=', length(which(!is.na(StationList[[i]]$Zero)))))
      SpatialDataset$Zero[SpatialDataset[, site_D]@data == StationList[[i]][1, site_SL]] = length(which(!is.na(StationList[[i]]$Zero)))
    }
  }
  return(SpatialDataset)
}
