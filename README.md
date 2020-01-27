# EnGauge README
This is a repository for gathering Environmental Gauge data, processing that data, and doing basic exploratory (spatial & temporal) data analysis.

The repository depends on several packages and functions developed by the USGS and NOAA to retreive data from government host servers. Significant processing is required to use some of these datasets for environmental modeling. This repository provides several functions to process datasets for this purpose.

It is likely that others have made similar functions to process these and other environmental gauge data. If you have scripts to process these datasets that provide additional functionality, please write to the authors of this repository and ask to contribute or link code repositories!

---
### Gauge datasets include:
1. USGS streamflow from the NWIS portal
2. EPA STORET, USGS, USDA and other water quality data via the water quality portal
3. NOAA ACIS, GHCN weather station data

---
### USGS and NOAA R libraries include:
1. dataRetrieval - main package for streamflow and water quality gauge data retrieval functions from government servers  
  Recommended resource - Package Readme: dataRetrieval.html (located somewhere on your computer after installing the R package)  
  [Function help](https://github.com/USGS-R/dataRetrieval)  
  [Blog](https://owi.usgs.gov/R/training-curriculum/usgs-packages/)  
  [Slideshow](https://owi.usgs.gov/R/dataRetrieval.html#1)

2. rnoaa - main package for NOAA weather station data retrieval.

---
### Other R Sources and Code Influences
[hydroTSM vignette](https://cran.r-project.org/web/packages/hydroTSM/vignettes/hydroTSM_Vignette-knitr.pdf)

[hydroTSM boxplots and flow duration curves](https://rstudio-pubs-static.s3.amazonaws.com/30850_e9f98501a0c3420897eb59de6abae6ab.html)

---
### Code Dependencies
1. [Geothermal_ESDA](https://github.com/jds485/Geothermal_ESDA) - the ColorFunctions.R script is used for plotting.

---
### Main script: USGSdataRetrieval.R
**Code Notes**
  
  This code was developed using R version 3.5.3, and all R packages were also downloaded for that version of R. You may try using other versions, but compatibility is not guaranteed.
  
  This code was developed on Windows 10 OS. The compatibility on other OS is untested. File an issue on GitHub if you experience issues with your OS.
  
  There are comment section and sub-section headers throughout this script that indicate blocks of code used for tasks and sub-tasks, respectively. Sub-headers are indented by one space more than the main header. In RStudio, these blocks can be opened and closed with arrows on the left side of the scripting window. This can be useful for navigating the script.

  Because of the wide variety of needs for different projects, this main script serves as examples of how downloading and processing data could be completed; you may want to modify for your project (e.g., if you do not want to use DEM processing, or have no water quality data needs).
  
  Input data for running this script are provided in the DataForExamples folder. Note that DEM data are too large to host on GitHub, but you can download the 2 DEM tiles used for the examples using [The National Map](https://viewer.nationalmap.gov/advanced-viewer/) and select the NED tiles corresponding to "1 arc-second n40w078 1 x 1 degree ArcGrid" and "1 arc-second n40w077 1 x 1 degree ArcGrid". You can also run the script without running the DEM sections, which all have "DEM" in their section header.

All examples in this script run independently of one another after specifying directory and file information, and code options.
   - Run the lines of code through the section entitled "Make Region of Interest (ROI) buffer"
   - Streamflow Example: Downloading and processing streamflow data using two methods.
   - Water Quality Example: Downloading and processing Total Nitrogen and Total Phosphorus water quality data using two methods. Outliers are not checked in this example because for this area the temporal data length is not great, and it would be difficult to detect and/or define an outlier.
   - NOAA Weather Station Example: Downloading weather station gauges and selecting sites with precipitation, maximum temperature, and minimum temperature. Duplicates, negatives, and zeros are not checked in this example, but they could be using the functions employed for streamflow and water quality. Aggreation of timeseries to monthly and annual is also not shown, but could be implemented.


**0. Preparing to Use the Code**

**_Install R Libraries_**

  Users may need to install several R libraries if using this code for the first time. These are listed in the "Load libraries and functions" heading within the script.


**_Prepare Input Data_**

  For each acquisition of streamflow, water quality, and weather station data, users should have either: 
  1. Method 1 (recommended, applies to streamflow, water quality, and weather stations): Scenario is that you know your region of interest but do not yet have data. Supply a polygon shapefile of the region of interest (ROI). The code will use this shapefile and an optional radial buffer around this shapefile to acquire data. There is an option to specify a different buffer for weather stations because it's common to need a larger region for these information. 
  OR
  2. Method 2 (streamflow and water quality only): Scenario is that you already have data (e.g. a gauge number). Supply a csv file of streamflow gauges (water quality gauges). The streamflow file must contain columns called Source, and GaugeNum. It may contain other columns. The water quality file must be the same as downloaded from the section below entitled "Water Quality Retrieval Method 2".


**_Set directory and file names_**

  The directories of input and output files, and the names of the files are specified at the top of the script.


**_Set the Project Coordinate System_**

  Users will need the [EPSG Code](https://spatialreference.org/ref/?page=2) for their coordinate system of choice for their project. All data will be projected into and saved in this coordinate system.


**_Set plot x- and y-axis limits_**

  If you do not want to use a particular limit, set equal to NULL. Other plot modifications will require editing of the script.


**_Set the quantile that defines a high outlier_**

  The default is 0.99, for the 99th percentile.


**_Obtain DEM Tiles (optional)_**

  Users who want to use the DEM tile merging/mosaicking functionality, processDEM, will have to download DEM raster tiles manually. These must have the same raster pixel resolution, be in the same coordinate system, and be spatially adjacent to each other. 



**1a. Streamflow Gauge Retrieval**

  There are two methods implemented to download streamflow gauge data, depending on what information you're starting with. Throughout the code there are labels for Method 1 only and Method 2 only processing steps. 


**_Streamflow Gauge Retrieval Method 1:_**  

  NOTE: This method seems more reliable than Method 2, and returned more gauges in several test regions.

The function `whatNWISsites()` will return gauge numbers with coordinates for a specified state (county, or various other criteria)  
  Search Criteria for `whatNWISsites()`: [Table 1 on waterqualitydata website](https://www.waterqualitydata.us/webservices_documentation/)  
   - References to "domain service" in Table 1 indicate to look at Table 2 on the same website.  
   - Likely the most common filtering method: [characteristicName options](https://www.waterqualitydata.us/public_srsnames/)  
     Example: whatNWISsites(statecode="MD", characteristicName="Phosphorus")  
   - Note: this function must use something like grep when characteristicName is specified. Searching for Phosphorus returned all instances of Phosphorus in the table. If you want a specific type of metric only, specify the parameterCd (parameter code) instead of the characteristicName.  
   - This function does not return full site information. For full site information, use `readNWISsite()` with the vector of gauge numbers that is returned using `whatNWISsites()`.


**_Streamflow Gauge Retrieval Method 2:_**  

You can find USGS gauges of any type in your region of interest [here](https://cida.usgs.gov/enddat/dataDiscovery.jsp)  
  1. Record the coordinates of the bounding box for your region of interest! They could be important for the next step. Select the data you want and copy the resulting table of values into a txt file. txt is important to retain leading zeros on gauge numbers. If you inhereted a csv or txt file without leading zeros, a script below can add them to NWIS gauges. Other gauge datasets seem to not require leading zeros.  
  2. Coordinates of gauges are not reported in that download :(  
     The USGS function `readNWISsite()` can look up the coordinates and all site info, given the gauge numbers :)  
     You can also find coordinates for your gauges if you provide a bounding box of coordinates [here](https://waterdata.usgs.gov/nwis/inventory?search_criteria=lat_long_bounding_box&submitted_form=introduction)


**1b. Streamflow Data Coordinate Reprojection**

  The downloaded gauges may be in several different coordinate systems. An automatic reprojection and merging is in development. For now, a manual method is implemented that requires users to look up EPSG codes for the projections of their downloaded datasets.
The reprojections are found under the headers "Make a spatial dataframe: Method X" where X is 1 or 2.



**2a. Water Quality Portal Data Retrieval**

  There are two methods implemented to download water quality data. Throughout the code there are labels for Method 1 only and Method 2 only processing steps.

  There are other water quality data that are not contained within the water quality portal. These datasets may also be processed by some functions in this repository, but may require the user to modify their dataset to be compatible with the functions.


**_Water Quality Portal Retrieval Method 1:_**  

  You can use the `whatWQPsites()`, using the same query criteria as the `whatNWISsites()` function, described in 1a. above. The `whatWQPsites()` function returns all site information, same as the waterqualitydata website described in Method 2 below. You will have to modify the function choices to meet your specific needs. These functions are below the header "Method 1: Read water quality station data using the USGS function"


**_Water Quality Portal Retrieval Method 2:_**  

You can find water quality sites [here](https://www.waterqualitydata.us/portal/)  
  Download "site data only" to receive a csv file with gauge/site information with coordinates. The MonitoringLocationIdentifier field is used to download data for each gauge/site.
  Note that the example in this script shows loading and processing data from a file, but makes plots only for the Method 1 approach. The sites selected may be different between the two methods, as shown on the maps that are plotted for each method's water quality sites.


**2b. Streamflow Data Coordinate Reprojection**

  Similarly to 1b above, the downloaded information may have several coordinate systems. The places where manual reprojections can be found are under the headers "Make a spatial dataframe: Method X" where X is 1 or 2.



**3. Weather Data Retrieval from NOAA**

  The `ghcnd_stations()` funuction grabs all available weather stations. The region of interest shapefile with radial buffer is used to select gauges.



**4. Comparing altitude of gauge vs. DEM**  

Sections for DEM processing are located at the end of the Streamflow and Weather Gauge sections under the headers "DEM for Streamflow" and "DEM for Weather Gauges" respectively.

Downloaded gauge/site coordinates may include the altitude of the gauge, which can be important vs. DEM elevation (e.g. if there's a cliff at the gauge vs. DEM mean elevation of the pixel).  
  - Ensure that the units and the vertical datum are the same for your DEM and all gauge altitudes. A common error in gauge altitudes is USGS data reported in m instead of in ft.  
  - DEMs in the US tend to be NAVD88 in m, whereas gauges tend to be referenced to NGVD29 in ft. Differences after unit conversion tend to be minor in the US, [except in the West](https://www.ngs.noaa.gov/TOOLS/Vertcon/vertcon.html)  
  - [Altitude datum codes](https://help.waterdata.usgs.gov/code/alt_datum_cd_query?fmt=html)  
  - [Altitude collection method codes](https://help.waterdata.usgs.gov/code/alt_meth_cd_query?fmt=html)



**5. Check Data Quality Codes**  

You can find data quality codes for USGS datasets [here](https://help.waterdata.usgs.gov/codes-and-parameters/codes#discharge_cd)  
  - You should always look at these quality codes, and process your data accordingly. This script colors streamflow timeseries by error code, but doesn't process further than that.  
  - [codes for streamflow 1](https://help.waterdata.usgs.gov/codes-and-parameters/daily-value-qualification-code-dv_rmk_cd)  
  - [codes for streamflow 2](https://help.waterdata.usgs.gov/codes-and-parameters/instantaneous-and-daily-value-status-codes)  
  Yes, there are 2 separate reference schemes for streamflow data.



**6. Process Data**

Functions in the code do the following checks and actions for downloaded products:
i.	Check for duplicate observations
ii.	Check for zeros and negative values
iii.	Fill date gaps (add NA values to dates that are missing from timeseries)
iv.	Aggregate to daily, monthly, and annual timesteps
v.	Spatio-temporal outlier detection (simple at the moment)



**7. Make Plots**

Plots for each downloaded dataset include:
Maps of data locations
Timeseries
Boxplots
Histograms
Hypsometric curves for DEMs

---
### Function files:

1. addZerosToGaugeNames.R - for adding leading zeros to NWIS streamflow gauge numbers.
2. extractWQdata.R - for extracting separate timeseries for each variable collected at sites downloaded from the water quality portal database. Also places timeseries in chronological order.
3. missingDates.R - for filling in missing dates in downloaded timeseries. Missing dates are assigned NA values. Also places timeseries in chronological order.
4. processDEM.R - for mosaicking DEM files into one DEM file.
5. checkDuplicates.R - checks timeseries for duplicate records, and keeps only one of them.
6. checkZerosNegs.R - several functions to check timeseries for zero and negative values, and replace with user-specified values.
7. fdcDefaultModification.R - a modified flow duration curve script from the [hydroTSM R package](https://www.rdocumentation.org/packages/hydroTSM/versions/0.5-1)
8. formatMonthlyMatrix.R - transforms a monthly timeseries into a matrix with 12 columns for the months and n rows for the years.
9. matplotDates.R - modified matplot() function that allows dates to be plotted on the x-axis in date format instead of numeric format.
10. aggregateTimeseries.R - function to summarize timeseries data by day, month, or year from its original format.
11. scatterHistCols.R - modified scatterHist() function from the psych package to accept colors, and provide density smoothing by color.

---
### Other Example Files - In development, data not all available here:

Example 1 is less in development than example 2. Both of these are based on earlier versions of the main script.
1. Method1Example_GFBR.R - example employing only the Method 1 downloads to the Baltimore Ecosystem Study Long Term Ecological Research Site in Gwynns Falls and Baisman Run (GFBR)
2. BES_GF_BR_SL.R - example of loading the processed data for further analysis, using site water quality data (as opposed to Water Quality Portal Data), downloading weather station data, among other tasks.

---
## Citation and Contact Information

Citation: Smith, J.D. et al. (2019). EnGauge. Online Github Repository. https://github.com/jds485/EnGauge

Contact Information:  
Jared D. Smith (js4yd@virginia.edu)  
Jonathan R. Lamontagne (Jonathan.Lamontagne@tufts.edu)  
Caitline A. Barber (Caitline.Barber@tufts.edu)  
Julianne D. Quinn (jdq6nn@virginia.edu)

---
## License Information

Please see the LICENSE file. You must cite this repository if you use it for your work.
