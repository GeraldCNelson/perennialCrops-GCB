# code to do various perennial crop calculations
# constants -----
library(terra)
sspChoices <- c("ssp126", "ssp585")
modelChoices <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
modelChoices_lower <- tolower(modelChoices)
startYearChoices <- c(2041, 2081)
hemispheres <- c("NH", "SH")
extent_NH <- ext(-180, 180, 0, 90)
extent_SH <- ext(-180, 180, -60, 0) #-60 gets rid of Antarctica for SH
extent_noAntarctica <- ext(-180, 180, -60, 90)
yearRange <- 19
speciesChoices <- c("almond_main", "apple_main", "cherry_main", "olive_main", "grape_main")
speciesNames <- gsub("_main", "", speciesChoices)
chillLevels <- c("_lo", "_main")
source("R/perennialFunctions.R")
source("R/perennialsPrep.R")
suitabilityLevel <- "good"
cropVals <- get(paste0("majorCropValues", "_main")) # used in several of the functions below
path_data <- "data/perennials/"
path_runs <- "data/runs/"
path_harvest_data <- "data-raw/geotiff/"
path_graphics <- "graphics/"

woptList <- list(gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL = 6", "NUM_THREADS=ALL_CPUS")) # reduces tiff file sizes
coastline <- vect("data-raw/regionInformation/ne_50m_coastline/ne_50m_coastline.shp")
coastline_cropped <- crop(coastline, extent_noAntarctica )


# the code in {} below runs all the raw data crunching - extreme cold, heat, frost, GDDs - except CPs which is done in chillPortions.R (the output from that script is available at https://doi.org/10.48364/ISIMIP.842396.1); the code below that does aggregation - means by model and ensemble means. 

# test values
k <- "ssp585"
l <- 2081
speciesName <- "apple"
hem <- "NH"
suitabilityLevel = "good"
modelChoice <- "GFDL-ESM4"
chillLevel <- "_lo"
# Note: that the difference between _lo and _main is just in the chill portions for the crops

{
  # extreme cold calcs -----
  # scenarios
  for (hem in hemispheres) {
    for (speciesName in speciesNames) {
      for (k in sspChoices) {
        for (l in startYearChoices) {
          f_extremeCold(k, l, speciesName, hem, modelChoices_lower, cropVals)
        }
      }
      # historical
      k <- "historical"
      l <- 1991
      f_extremeCold(k, l, speciesName, hem, modelChoices_lower, cropVals)
    }
  }
  
  # heat damage  -----
  # scenarios
  for (hem in hemispheres) {
    for (speciesName in speciesNames) {
      for (k in sspChoices) {
        for (l in startYearChoices) {
          f_heatDamage(k, l, modelChoices_lower, speciesName, hem, suitabilityLevel, cropVals)
        }
      }
      # historical
      k <- "historical"
      l <- 1991
      f_heatDamage(k, l, modelChoices_lower, speciesName, hem, suitabilityLevel, cropVals)
    }
  }
  
  # frost damage -----
  # scenarios
  for (hem in hemispheres) {
    for (speciesName in speciesNames) {
      for (k in sspChoices) {
        for (l in startYearChoices) {
          f_frostDamage(k, l, speciesName, hem, modelChoices_lower, suitabilityLevel, cropVals)
        }
      }
    }
    # historical
    k <- "historical"
    l <- 1991
    f_frostDamage(k, l, speciesName, hem, modelChoices_lower, suitabilityLevel, cropVals)
  }
  # gdd calcs ------
  # note that gdds are the same for all variety chill portions, so generate just with species name
  # scenarios
  for (speciesName in speciesNames) {
    speciesChoice <- paste0(speciesName, "_main") # needed to get a cropVals variable with the gdds, _lo would also work as well
    cropVals <- get(paste0("majorCropValues", "_main"))
    topt_min <- cropVals[cropName == speciesChoice, gddtb]
    topt_max <- cropVals[cropName == speciesChoice, GDD_opt]
    for (modelChoice in modelChoices) {
      for (k in sspChoices) {
        for (l in startYearChoices) {
          print(system.time(f_computeGDDs(k, l, speciesName, modelChoice, topt_min, topt_max)))
        }
      }
      # historical
      k <- "historical"
      l <- 1991
      print(system.time(f_computeGDDs(k, l, speciesName, modelChoice, topt_min, topt_max)))
    }
  }
}
# runs files -----
# uses runsParms created in perennialFunctions.R which is sourced above
# scenarios
for (modelChoice in modelChoices_lower) {
  for (k in sspChoices) {
    for (l in startYearChoices) {
      f_runsSetup(k, l, modelChoice, runsParms)
    }
  }
  # historical
  k <- "historical"
  l <- 1991
  f_runsSetup(k, l, modelChoice, runsParms)
}
# gdd sums  ------ this takes forever because it works year by year.
# scenarios
for (speciesName in speciesNames) {
  for (modelChoice in modelChoices) {
    for (hem in hemispheres) {
      for (k in sspChoices) {
        for (l in startYearChoices) {
          f_gddSums(k, l, speciesName, hem, runsParms, modelChoice)
        }
      }
      # historical
      k <- "historical"
      l <- 1991
      f_gddSums(k, l, speciesName, hem, runsParms, modelChoice)
    }
  }
}

# aggregation code -----
# runs only used in GDD calcs so no need to do means of runs files
{
  # GDD sum, means by model -----
  #  scenarios
  cropVals <- get(paste0("majorCropValues", "_main"))
  for (k in sspChoices) {
    for (l in startYearChoices) {
      for (modelChoice in modelChoices) {
        f_gddSum_mean(k, l, speciesNames, modelChoice, hemispheres)
      }
    }
    # historical -----
    k <- "historical"
    l <- 1991
    f_gddSum_mean(k, l, speciesNames, modelChoice, hemispheres)
  }
  
  # ensemble GDD sum -----
  # scenarios
  for (k in sspChoices) {
    for (l in startYearChoices) {
      f_ensemble_GDD_sum_mean(k, l, yearRange, speciesNames, cropVals)
    }
  }
  # historical
  k <- "historical"
  l <- 1991
  f_ensemble_GDD_sum_mean(k, l, yearRange, speciesNames, cropVals)
}

# gdds suitable-----
for (hem in hemispheres) {
  for (speciesName in speciesNames) {
    # scenarios
    for (k in sspChoices) {
      for (l in startYearChoices) {
        print(paste0("speciesName: ", speciesName, ", ssp choice: ", k, ", start year: ", l))
        f_gddsSuitability(k, l, speciesName, hem, cropVals)
      }
    }
    # historical
    k <- "historical"
    l <- 1991
    print(paste0("speciesName: ", speciesName, ", ssp choice: ", k, ", start year: ", l))
    f_gddsSuitability(k, l, speciesName, hem, cropVals)
  }
}

# combined damage, scenarios -----
# code to read in 1/0 metrics files and produce 1/0 tifs where the crop is potentially growable. The chill portions files are created in the chill_projections.R script
# chill_projections.R # uncomment this line to generate new chill portions files
# Important note: The chillPortions.R script must be run before combined damage whenever any chill portion value is changed.

 chillPortions.R 
{
  # combined damage -----
  for (chillLevel in chillLevels) {
    for (speciesName in speciesNames) {
      for (k in sspChoices) {
        for (l in startYearChoices) {
          print(paste0("speciesName: ", speciesName, ", chilllevel: ", chillLevel, ", ssp choice: ", k, ", start year: ", l))
          f_combinedDamage(k, l, speciesName, suitabilityLevel, chillLevel)
        }
      }
      # historical
      k <- "historical"
      l <- 1991
      print(paste0("speciesName: ", speciesName, ", chilllevel: ", chillLevel, ", ssp choice: ", k, ", start year: ", l))
      f_combinedDamage(k, l, speciesName, suitabilityLevel, chillLevel)
    }
  }
}

# suitable locations graphics, historical -----
for (hem in hemispheres) {
  for (speciesName in speciesNames) {
    k <- "historical"
    l <- 1991
    print(paste0("Suitability level:  ", suitabilityLevel, ", speciesName: ", speciesName, ", ssp choice: ", k, ", start year: ", l))
    f_suitableLocsGraphics(k, l, speciesName, suitabilityLevel)
    #scenarios
    for (k in sspChoices) {
      for (l in startYearChoices) {
        print(paste0("Suitability level:  ", suitabilityLevel, ", speciesName: ", speciesName, ", ssp choice: ", k, ", start year: ", l))
        f_suitableLocsGraphics(k, l, speciesName, suitabilityLevel)
      }
    }
  }
}
