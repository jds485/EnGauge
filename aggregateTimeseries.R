#Function to aggregate timeseries into average values.

#Fixme: make an option to supply any function for aggregating.
#Fixme: This is currently only for WQP data. Make a generic aggregation field name in the function call

aggregateTimesteps = function(StationList, #named list of station timeseries
                              aggVal, #text: d, m, a. Daily, Monthly, Annual
                              aggVar, #variable in StationList to aggregate by
                              date, #fieldname supplying the date
                              site  #fieldname supplying the site ID
){
  ld = list()
  lm = list()
  la = list()
  
  for (i in 1:length(StationList)){
    if ('d' %in% aggVal){
      #Average days with multiple measurements into a daily average concentration
      #make a new file with only the values by date
      daily = stats::aggregate(StationList[[i]][,aggVar] ~ StationList[[i]][,date], data = StationList[[i]], FUN = mean)
      colnames(daily)[colnames(daily) == 'StationList[[i]][, aggVar]'] = aggVar
      colnames(daily)[colnames(daily) == 'StationList[[i]][, date]'] = date
      daily[,site] = rep(StationList[[i]][1,site], nrow(daily))
      ld = c(ld, list(daily))
    }
    if ('m' %in% aggVal){
      #Month-year averaging:
      mthyr = stats::aggregate(StationList[[i]][,aggVar] ~ month(StationList[[i]][,date]) + year(StationList[[i]][,date]), data = StationList[[i]], FUN = mean)
      colnames(mthyr)[1] = 'month'
      colnames(mthyr)[2] = 'year'
      colnames(mthyr)[colnames(mthyr) == 'StationList[[i]][, aggVar]'] = aggVar
      mthyr$YrMthDy = as.Date(paste0(mthyr$year, '-', mthyr$month, '-01'))
      mthyr[,site] = rep(StationList[[i]][1,site], nrow(mthyr))
      lm = c(lm, list(mthyr))
    }
    if ('a' %in% aggVal){
      #Year averaging:
      ann = stats::aggregate(StationList[[i]][,aggVar] ~ year(StationList[[i]][,date]), data = StationList[[i]], FUN = mean)
      colnames(ann)[1] = 'year'
      colnames(ann)[colnames(ann) == 'StationList[[i]][, aggVar]'] = aggVar
      ann$YrMthDy = as.Date(paste0(ann$year, '-01-01'))
      ann[,site] = rep(StationList[[i]][1,site], nrow(ann))
      la = c(la, list(ann))
    }
  }
  
  if ('d' %in% aggVal){
    names(ld) = names(StationList)
  }
  if ('m' %in% aggVal){
    names(lm) = names(StationList)
  }
  if ('a' %in% aggVal){
    names(la) = names(StationList)
  }
  
  return(list(daily = ld, mthyr = lm, ann = la))
}
