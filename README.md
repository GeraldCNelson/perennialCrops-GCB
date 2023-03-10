[![DOI](https://zenodo.org/badge/588235215.svg)](https://zenodo.org/badge/latestdoi/588235215)

# "Assessing temperature-based adaptation limits to climate change of temperate perennial fruit crops" - R and c++ code and data

# Description of the data and file structure
The code and data needed to reproduce the results in Global Change Biology, "Assessing temperature-based adaptation limits to climate change of temperate perennial fruit crops" are available in the Dryad site (data, 79 GB) ([https://doi.org/10.5061/dryad.5dv41ns98](https://doi.org/10.5061/dryad.5dv41ns98)) and the Zenodo site (code and small data sets) ([https://doi.org/10.5281/zenodo.7539023](https://doi.org/10.5281/zenodo.7539023)).

To generate the results, the _Directory structure_ section below describes the directories needed (and should automatically be created when the zenodo files are downloaded) and what is contained in them.The _Order of operations_ section describes the order in which the R code files need to be run to generate the results.

## Directory structure

The directory structure for the code is described below. The data from the dryad download arrives as a zip file that contains a directory with 54 items. Of these, 52 are individual climate data geotiff files (e.g. gfdl-esm4_tas_cropped_ssp585_2041_2060.tif). All of the climate data _xxx.tif_ files need to be moved to the _climdata_ directory (see below). The other data are in a zip folder called _chillPortions.zip_. The contents of that folder need to be moved to the _data/chillPortions_ directory.

- climdata - climate data files. The ISIMIP project [https://www.isimip.org] prepares daily bias-corrected 1/2 degree resolution from five earth system models (ESMs - GFDL-ESM4, UKESM1-0-LL, MPI-ESM1-2-HR, MRI-ESM2-0, and IPSL-CM6A-LR). The paper uses the ISIMIP3b data from 
[https://doi.org/10.48364/ISIMIP.842396.1](https://doi.org/10.48364/ISIMIP.842396.1). The data are downloaded as netcdf files with 10 years of data. These files were combined to create 3 20-year periods for each ESM - 1991-2010, 2041-2060, and 208-2100 - and cropped to show data only for land areas. Each of the needed _.nc_ can be accessed using the urls in the csv file _data-raw/ISIMIPNCfiles.csv"_. A land-only mask is available - _data-raw/landseamask_no_antarctica.nc_.
- data-raw/geotiff - Contains geotiff files of early century crop areas for the five crops included in the paper - almond, apple, cherry, grape, and olive - created from data in [www.earthstat.org](http://www.earthstat.org). 
- data - processed data that are used to generate the final results. Generated by _perennialCalcs.R_.
  * chillPortions - The contents of the chillPortions.zip file should be put in this directory. They can also be generated with the R script _chillSpatial_projections.R_ but the process takes a long time (possibly 24 hours!) so we have also included a complete set of the resulting output. It is part of the dryad download and should be transferred to this directory.
  * crops - raster mask files for each crop created from the earthstat data.
  * perennials - an excel version of the table in the supplementary materials of the paper
  * growingDegreeDays - daily growing degree days by ESM and ensemble means across the models, locations where growing degree days are not   * regionInformation - files used to construct the boarders used in the graphics
limiting for one of the crops.
   * runs - files with the results of calculations on the number of days with a value below or above a limit, a run.
- graphics - destination for graphics results
- results - destination for _.csv_ and _.docx_ output files.

## Order of operations

All R scripts are in the R folder. The R script _perennialCalcs.R_ has code to run all of the calculations needed to produce the results in the paper. Each operation is done with a function that returns a 1 or 0 value for each location depending on whether it has met one of the plant requirements. For example the extreme cold function _f_extremeCold(k, l, speciesName, hem, modelChoices_lower, cropVals)_ take inputs for scenario (k), starting year (l), plant name (speciesName), hemisphere (hem), lower case version of the list of model names, and the cropvals table. 

The scripts need to be run in the order they appear in _perennialCalcs.R_. Some of them can take quite a while. These are identified in the code. 
When all the data crunching is done the individual outputs for the paper - tables and figures - can be run from the following scripts.
- GCBPerennials_Figure_1.R
- GCBPerennials_Figure_2.R
- GCBPerennials_Table_4.R
- GCBPerennials_Table_5.R
- GCBPerennials_Table_7.R

## Install R libraries

Here's some code to install any of the needed packages that are not already installed

`packages <- c("terra", "viridis", "data.table", "flextable", "officer", "crayon", "magrittr", "doParallel", "foreach", "Rcpp", "readxl")`
`installed_packages <- packages %in% rownames(installed.packages())`
`if (any(installed_packages == FALSE)) {  install.packages(packages[!installed_packages]) }`

The two figures for the paper use the pdf graphics driver built in to base R. The code used to produce Figure 1 and Figure 2 also uses the pdfcrop program to trim the white spaces around the edges of the pdf. 

For Mac users, run this code from a terminal to install pdfcrop. This uses the homebrew system (https://brew.sh). A plug - homebrew is a really nice way to manage all sorts of code that the MacOS should have but doesn't. It also works for linux users.
`brew install --cask mactex`

For windows users, these directions for installing pdfcrop might be  useful but haven't been verified.'

"You need first to install Perl support for Windows (for example https://www.activestate.com/activeperl) and after that use Miktex package manager to install pdfcrop."

More directions here [https://pdfcrop.sourceforge.net]

## Run scripts

- chillingCalcs.R - output the chilling portions file stored in data-raw. This will take a long time and the output files are already available. Run only if new versions are needed.
- perennialCalcs.R  - directions for running are included in the file. Read through the code before sourcing the file

The graphics for the paper are generated by the following files

- GCBPerennials_Figure_1.R
- GCBPerennials_Figure_2.R
- GCBPerennials_Table_4.R
- GCBPerennials_Table_5.R
- GCBPerennials_Table_7.R

Tables 1, 2, 6, and 6 are just word tables that are not generated in R.

## Sharing/Access information
The raw climate data files are available for download from ISIMIP using the urls in the csv file _data-raw/ISIMIPNCfiles.csv"_.
