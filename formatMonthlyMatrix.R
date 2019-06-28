# Function for transforming a zoo monthly timeseries to a matrix with months (columns) by year (row) 
formatMonthlyMatrix = function(mts #zoo monthly timeseries data
                               #WaterYearStart #number of the month that the water year begins. Matrix will be returned with this month as the first column.
                               ){
  #Fill in the months before and after data are available. Full calendar years are needed.
  if(which(month.abb == format(time(mts)[1], '%b')) != 1){
    #Generate vector of NAs for the length needed to fill in
    NAvec = rep(NA, (which(month.abb == format(time(mts)[1], '%b')) - 1))
    #Create vector of the values that should be added to the start of the timeseries + rest of timeseries
    mts_plt = c(NAvec, coredata(mts))
    #Make a zoo object, and add months as indices for the zoo object.
    mts_plt = zoo(x = mts_plt, order.by = c(time(mts)[1] %m-% months(seq(length(NAvec),1,-1)), time(mts)))
  }else{
    mts_plt = mts
  }
  #Repeat for the end of the timeseries
  if(which(month.abb == format(time(mts_plt)[length(mts_plt)], '%b')) != 12){
    NAvec = rep(NA, (12 - which(month.abb == format(time(mts_plt)[length(mts_plt)], '%b'))))
    mts_plt2 = c(coredata(mts_plt), NAvec)
    mts_plt2 = zoo(x = mts_plt2, order.by = c(time(mts_plt), time(mts_plt)[length(mts_plt)] %m+% months(seq(1,length(NAvec),1))))
  }else{
    mts_plt2 = mts_plt
  }
  #Form into a matrix
  M <- matrix(mts_plt2, ncol=12, byrow=TRUE)
  #Add column names as months in the order they appear in the timeseries
  colnames(M) <- month.abb[c(which(month.abb == format(time(mts_plt2)[1], '%b')):12, 1:which(month.abb == format(time(mts_plt2)[12], '%b')))[1:12]]
  #Add row names as years in the order they appear in the timeseries
  rownames(M) <- unique(format(time(mts_plt2), "%Y"))
  
  #Fixme: Rearrange to have start of water year in the first column
  #if(WaterYearStart != 1){
    #Need to switch around NAs after the months are switched around, and then the years also have to be changed.
  #  M = M[,c(WaterYearStart:12, 1:(WaterYearStart-1))]
  #}
  
  return(M)
}
