#Function to aggregate timeseries into average values.

#Fixme: make an option to supply any function for aggregating.
#Fixme: This is currently only for WQP data. Make a generic aggregation field name in the function call

aggregateTimesteps = function(StationList, #named list of station timeseries
                              aggVal #text: d, m, a. Daily, Monthly, Annual
){
  ld = list()
  lm = list()
  la = list()
  
  for (i in 1:length(StationList)){
    if ('d' %in% aggVal){
      #Average days with multiple measurements into a daily average concentration
      #make a new file with only the values by date
      daily = stats::aggregate(ResultMeasureValue ~ SortDate, data = StationList[[i]], FUN = mean)
      daily$MonitoringLocationIdentifier = rep(StationList[[i]]$MonitoringLocationIdentifier[1], nrow(daily))
      ld = c(ld, list(daily))
    }
    if ('m' %in% aggVal){
      #Month-year averaging:
      mthyr = stats::aggregate(ResultMeasureValue ~ month(StationList[[i]]$SortDate) + year(StationList[[i]]$SortDate), data = StationList[[i]], FUN = mean)
      colnames(mthyr)[1] = 'month'
      colnames(mthyr)[2] = 'year'
      mthyr$YrMthDy = as.Date(paste0(mthyr$year, '-', mthyr$month, '-01'))
      mthyr$MonitoringLocationIdentifier = rep(StationList[[i]]$MonitoringLocationIdentifier[1], nrow(mthyr))
      lm = c(lm, list(mthyr))
    }
    if ('a' %in% aggVal){
      #Year averaging:
      ann = stats::aggregate(ResultMeasureValue ~ year(StationList[[i]]$SortDate), data = StationList[[i]], FUN = mean)
      colnames(ann)[1] = 'year'
      ann$YrMthDy = as.Date(paste0(ann$year, '-01-01'))
      ann$MonitoringLocationIdentifier = rep(StationList[[i]]$MonitoringLocationIdentifier[1], nrow(ann))
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
