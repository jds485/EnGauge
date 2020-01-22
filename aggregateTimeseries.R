#Function to aggregate timeseries into daily, monthly, or annual values.

aggregateTimesteps = function(StationList, #named list of station timeseries
                              aggVal, #character vector containing: d, m, and/or a. These stand for Daily, Monthly, and Annual
                              aggVar, #column name in StationList to aggregate by
                              date, #column name in StationList supplying the date
                              site,  #column name in StationList supplying the site ID
                              fun #function to use for aggregation (character)
){
  ld = list()
  lm = list()
  la = list()
  
  for (i in 1:length(StationList)){
    if ('d' %in% aggVal){
      #Average days with multiple measurements into a daily average concentration
      #make a new file with only the values by date
      daily = stats::aggregate(StationList[[i]][,aggVar] ~ StationList[[i]][,date], data = StationList[[i]], FUN = fun)
      colnames(daily)[colnames(daily) == 'StationList[[i]][, aggVar]'] = aggVar
      colnames(daily)[colnames(daily) == 'StationList[[i]][, date]'] = date
      daily[,site] = rep(StationList[[i]][1,site], nrow(daily))
      ld = c(ld, list(daily))
    }
    if ('m' %in% aggVal){
      #Month-year averaging:
      mthyr = stats::aggregate(StationList[[i]][,aggVar] ~ month(StationList[[i]][,date]) + year(StationList[[i]][,date]), data = StationList[[i]], FUN = fun)
      colnames(mthyr)[1] = 'month'
      colnames(mthyr)[2] = 'year'
      colnames(mthyr)[colnames(mthyr) == 'StationList[[i]][, aggVar]'] = aggVar
      mthyr$YrMthDy = as.Date(paste0(mthyr$year, '-', mthyr$month, '-01'))
      mthyr[,site] = rep(StationList[[i]][1,site], nrow(mthyr))
      lm = c(lm, list(mthyr))
    }
    if ('a' %in% aggVal){
      #Year averaging:
      ann = stats::aggregate(StationList[[i]][,aggVar] ~ year(StationList[[i]][,date]), data = StationList[[i]], FUN = fun)
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
