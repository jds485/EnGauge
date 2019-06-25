#Function to add leading zeros to NWIS gauges
# NOTE: suppressing warnings for NAs introduced by coercion, which is intended.
#       Users should check that other warnings are not also being suppressed.
addZeros = function(Stations){
  StationStart = substr(Stations$GaugeNum, start = 1, stop = 1)
  for (i in 1:nrow(Stations)){
    if(Stations$Source[i] == 'NWIS'){
      # Some of the gauges are not numbers so only add 0 to number gauges
      if (suppressWarnings(is.na(as.numeric(StationStart[i])))){
        #This station starts with a character. Retain original name
        Stations$GaugeNum[i] = Stations$GaugeNum[i]
      }else{
        #Add a leading 0
        Stations$GaugeNum[i] = paste0("0", Stations$GaugeNum[i])
      }
    }
  }
  return(Stations)
}
