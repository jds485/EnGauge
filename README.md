# EnGauge README
Repository for gathering Environmental Gauge data, processing that data, and doing basic exploratory (spatial & temporal) data analysis.

The repository depends on several packages and functions developed by the USGS to retreive data from government host servers. Significant processing is required to use some of these datasets for environmental modeling. This repository provides several functions to process those datasets.

It is likely that others have made similar functions to process USGS and other environmental gauge data. If you have scripts to process these datasets, please write to the authors of this repository and ask to contribute or link code repositories!

**Full README In Development**
---
### USGS R libraries include:
1. dataRetrieval - main package for gauge data retrieval functions from government servers  
  Recommended resource - Package Readme: dataRetrieval.html (located somewhere on your computer after installing the package)  
  [Function help](https://github.com/USGS-R/dataRetrieval)  
  [Blog](https://owi.usgs.gov/R/training-curriculum/usgs-packages/)  
  [Slideshow](https://owi.usgs.gov/R/dataRetrieval.html#1)

2. EGRET - package for plotting water quality data  
  Resources ...

---
### Other R Sources and Code Influences
[hydroTSM vignette](https://cran.r-project.org/web/packages/hydroTSM/vignettes/hydroTSM_Vignette-knitr.pdf)
[hydroTSM boxplots and flow duration curves](https://rstudio-pubs-static.s3.amazonaws.com/30850_e9f98501a0c3420897eb59de6abae6ab.html)

---

### Gauge datasets include:
1. USGS streamflow from the NWIS portal
2. EPA STORET, USGS, USDA and other water quality data via the water quality portal
3. ...

---

### Main script: USGSdataRetrieval.R
**1. Streamflow Gauge Retrieval**
There are two methods implemented to download streamflow gauge data, depending on what information you're starting with. There are likely many other possible ways to setup data for download, but these two are common.

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
     The USGS function `readNWISsite()` can look up the coordinates and all site info, given the gauge numbers.  
     You can also find coordinates for your gauges if you provide a bounding box of coordinates [here](https://waterdata.usgs.gov/nwis/inventory?search_criteria=lat_long_bounding_box&submitted_form=introduction)

**2. Water Quality Data Retrieval**  
**_Water Quality Retrieval Method 1:_**  
You can use the `whatWQPsites()`, using the same query criteria as the `whatNWISsites()` function, described in 1. above. The `whatWQPsites()` function returns all site information, same as the waterqualitydata website described in Method 2.

**_Water Quality Retrieval Method 2:_**  
You can find water quality sites [here](https://www.waterqualitydata.us/portal/)  
  Download "site data only" to receive a csv file with gauge/site information with coordinates. The MonitoringLocationIdentifier field is used to download data for each gauge/site.

**3. Comparing altitude of gauge vs. DEM:**  
Downloaded gauge/site coordinates may include the altitude of the gauge, which can be important vs. DEM elevation (e.g. if there's a cliff at the gauge vs. DEM mean elevation of the pixel).  
  - Ensure that the units and the vertical datum are the same for your DEM and all gauge altitudes. A common error in gauge altitudes is USGS data reported in m instead of in ft.  
  - DEMs in the US tend to be NAVD88 in m, whereas gauges tend to be referenced to NGVD29 in ft. Differences after unit conversion tend to be minor in the US, [except in the West](https://www.ngs.noaa.gov/TOOLS/Vertcon/vertcon.html)  
  - [Altitude datum codes](https://help.waterdata.usgs.gov/code/alt_datum_cd_query?fmt=html)  
  - [Altitude collection method codes](https://help.waterdata.usgs.gov/code/alt_meth_cd_query?fmt=html)

**4. Check Data Quality Codes**  
You can find data quality codes for USGS datasets [here](https://help.waterdata.usgs.gov/codes-and-parameters/codes#discharge_cd)  
  - You should always look at these quality codes, and process your data accordingly. This script colors streamflow time series by error code, but doesn't process further than that.  
  - [codes for streamflow 1](https://help.waterdata.usgs.gov/codes-and-parameters/daily-value-qualification-code-dv_rmk_cd)  
  - [codes for streamflow 2](https://help.waterdata.usgs.gov/codes-and-parameters/instantaneous-and-daily-value-status-codes)  
  Yes, there are 2 separate reference schemes for streamflow data.

**5. Plots, EDA, Outlier Detection...

---
### Function files:
1. addZerosToGaugeNames.R - for adding leading zeros to NWIS streamflow gauge numbers.
2. extractWQdata.R - for extracting separate timeseries for each variable collected at sites downloaded from the water quality portal database. Also places timeseries in chronological order.
3. missingDates.R - for filling in missing dates in downloaded timeseries. Missing dates are assigned NA values. Also places timeseries in chronological order.
4. processDEM.R - for mosaicking DEM files into one DEM file.
5. checkDuplicates - checks timeseries for duplicate records, and keeps only one of them.
6. checkZerosNegs - several functions to check timeseries for zero and negative values, and replace with user-specified values.
7. fdcDefaultModification - a modified flow duration curve script from the [hydroTSM R package](https://www.rdocumentation.org/packages/hydroTSM/versions/0.5-1)
8. formatMonthlyMatrix - transforms a monthly timeseries into a matrix with 12 columns for the months and n rows for the years.
9. matplotDates - modified matplot() function that allows dates to be plotted on the x-axis in date format instead of numeric format.
10. aggregateTimeseries - function to summarize timeseries data by day, month, or year from its original format.
---
## Citation and Contact Information
Citation: Smith, J.D., ... (2019). EnGauge. Online Github Repository. https://github.com/jds485/EnGauge

Contact Information:  
Jared D. Smith (js4yd@virginia.edu)  
Jonathan R. Lamontagne (Jonathan.Lamontagne@tufts.edu)  
Caitline A. Barber (Caitline.Barber@tufts.edu)  
Julianne D. Quinn (jdq6nn@virginia.edu)

---
## License Information
You must cite this repository if you use it for your work.
...
