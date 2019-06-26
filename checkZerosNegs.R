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