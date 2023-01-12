# construct masks for locations of plants in the early 20th century

library(terra)
getcropAreaYield <- function(cropName, dataType) { #dataType is area or yield
  # source of data
  tifZipUrl <-  "https://s3.us-east-2.amazonaws.com/earthstatdata/HarvestedAreaYield175Crops_Geotiff.zip"
  tifzipFile <- paste0("data-raw/crops/HarvestedAreaYield175Crops_Geotiff.zip")
  tifFileLoc <- "data-raw/geotiff/"
  tiffilecrop <- cropName
  if (dataType %in% "area") {
    tifcropFile <- paste0(tiffilecrop, "_HarvestedAreaHectares.tif")
  } else {tifcropFile <- paste0(tiffilecrop, "_YieldPerHectare.tif")
  }
  tifOut <- rast(paste0(tifFileLoc, tiffilecrop, "/", tifcropFile))
  return(tifOut)
}

crops <- c("almond", "apple", "cherry", "grape", "olive")
for (i in crops) {
  print(i)
  #  i <- "olive"
  rInArea_hires <- getcropAreaYield(i, "area") # returns a spatRaster
  # Earthstat data (using its harvest area) has a pixel size of 5 minute resolution. Need to convert to 1/2 degree to get to cmip6 cell size
  # 5 min = 0.0833333 degree
  # 30 min = 0.5 degree
  
  rInAreaAgg <- aggregate(rInArea_hires, fact = 6, fun = "sum")
 fileName_out <- paste0("data/crops/rasterMask_", i, ".tif")
   print(fileName_out)
  writeRaster(rInAreaAgg, fileName_out,  overwrite = TRUE)
}


