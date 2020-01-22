#Function for DEM mosaikinng from separate tiles
processDEM = function(dir_DEMs, #Character vector of directories for DEMs 
                      f_DEMs,    #Character vector of filenames of DEMs in those directories
                      pCRS      #Project coordinate system
){
  #Mosiac the tiles together
  for (d in 1:length(dir_DEMs)){
    if (d > 1){
      DEM2 = raster(x = paste0(dir_DEMs[d], "/", f_DEMs[d]))
      DEM = mosaic(DEM, DEM2, fun = mean)
    }else{
      DEM = raster(x = paste0(dir_DEMs[d], "/", f_DEMs[d]))
    }
  }
  #Project to project coordinate system
  DEM = projectRaster(DEM, crs = CRS(pCRS))
  return(DEM)
}