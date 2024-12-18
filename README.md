# Query of Hourly ERA5-Land Data and Aggregation of Daily Summary Heat Metrics
## Project Overview
Using the Python cdsapi package, ERA5-Land hourly temperature measures are downloaded from the Copernicus Climate Data Store (CDS). Heat index and humidex metrics are derived from ERA5-Land hourly 2-meter temperature and 2-meter dew point temperature. Daily mean, maximum, and minimum measures are assessed across administrative boundaries for 2-meter temperature, 2-meter dew point temperature, skin temperature, heat index and humidex. Kenya is used as a demonstration area for the Copernicus CDS API query and spatial aggregation, using administrative boundaries from the Database of Global Administrative Areas.

## Usage
This repository provides the building blocks for the query of any data from the Copernicus CDS. Code for spatial aggregation of ERA5-Land hourly netCDF files from their native raster format to daily measures averaged across administrative boundaries are also included, using the geopandas, shapely, xarray, and rioxarray packages. These spatial aggregation methods can be used to derive measures from ERA5-Land in analysis of sociodemographic disparities or epidemiologic studies of temperature and health outcomes. 

Comments are included within the code repository with descriptions of how to manipulate the ERA5 API language to query additional variables, time frames, or spatial extents. Please note that users will need to set up an account with the Copernicus CDS for the API query to function effectively. For more information on the cdsapi package and details about how to set up an account and access the necessary user ID and API key inputs for the API query please see: https://cds.climate.copernicus.eu/how-to-api#install-the-cds-api-client  

## Notes on Computation
- This code was run to generate the 24 years of heat measures across Kenya administrative boundaries using a shared computing cluster at Boston University. The use of a computing cluster allows for intensive computation to run more quickly and to run without the limits of conventional storage options on a single computer.
- The query of the data using cdsapi is done in monthly chunks due to restrictions in the CDS API.
- Each 1-month period of ERA5-Land data across Kenya with the three variables (2-m temp, dew point temp, skin temp) took about several minutes to download and is approximately 24MB for the raster file storage.
- For the 24 years, and 12 files each year, this results in just under 9 GB of storage for the raw raster data.
- The aggregation loop that computes the ward-level results takes about several hours per year to output the resulting .csv file on a local computer. In a .csv format, each year's data is ~162MB after including additional columns with geographic parameters for the larger geographies in which wards are nested.
- The computation efficiency and size will vary by use case. Aggregating to a smaller administrative boundary, or using a less spatially resolved raw product (such as ERA5, Non-Land with is a 25km grid) will reduce the needed computation resources.

## Data Sources
- [ERA5-Land hourly data from 1950 to present](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land)
- [Database of Global Administrative Areas](https://gadm.org/)

## Workflow
Three scripts are included in the code repository.    

1) 01_Query_ERA5_CDSAPI_Kenya_v1.py applies the cdsapi package to query hourly ERA5-Land data from the Copernicus CDS using the cdsapi package.  
   - The necessary code to set your API key in your Python environment is provided, as well as a description of the necessary inputs for an effective API query, including the dataset, spatial extent, time frame, and variables requested.
   - The query is broken down into queries for data every month. A for loop is used to query these one-month data cuts across the 24-year period of the full dataset. These cuts were set as the Copernicus CDS API restricts the amount of data accessible in a single pull.
   - *Please Note*: The API query can be time-consuming and the resulting output can be extremely large if the API query is not correctly specified. When first using the API for a new use case, consider requesting only a few hours or days of data for a small spatial area in order to get an estimate for the time and data size that will be queried with larger sets.
2) 02_Aggregate_ERA5_Kenya_SetFishnet_v1.py is the first in a two-step process for aggregating the queried rasters to the administrative boundaries.
   - A fishnet approach is used, where a grid of points is established and the merged with points representing the spatial extent of the administrative boundaries for the area of interest.
   - This is set as a separate file because it can be very time consuming with a highly spatially resolved raster. The ERA5-Land dataset is shared at a 9 km grid resolution, so this runs quickly in this application.
   - More details about the use of a fishnet for spatial aggregation of raster data can be found here: [https://github.com/Climate-CAFE/population_weighting_raster_data](https://github.com/Climate-CAFE/population_weighting_raster_data)
3) 03_Aggregate_ERA5_Kenya_UseFishnet_v1.py uses the fishnet grid set up in script 2, and merges the grid with centroids for the administrative boundaries in Kenya.
   - The raster data time series is reset from Coordinated Universal Time (UTC) to East Africa Time (EAT).
   - Hourly heat index and humidex are estimated from the 2-m temperature and 2-m dew point temperature.
   - Daily mean, minimum, and maximum are assessed across the raster, and then joined to the extraction points (the union of the fishnet grid and adminstrative centroids).
   - Spatially weighted averages across the administrative boundaries are computed from the extraction points.
   - The resulting output is a dataframe with Kenya administrative boundaries and a time series of daily mean, minimum, and maximum for 2-m temperature, 2-m dew point temperature, skin temperature, heat index, and humidex.
  
## Dependencies
Packages used in this repository include:

- import cdsapi for era5 data query
- import geopandas for spatial operations on geometric types
- import xarray for raster data
- import rioxarray for raster data
- import shapely for manipulation and analysis of geometric objects 
- from osgeo import gdal, ogr for creating a fishnet grid
- import os for creating/removing a file path
- import glob for findng all the pathnames matching a specified pattern 
- import numpy for scientific computing
- import pandas for manipulation of structured data 
- from functools import reduce for merging datasets
- import math for mathematical calculations
- import pytz for cross platform timezone calculations
- from datetime import datetime for manipulating dates and times
- import xvec for extending Xarray to work with vector data 
- import re for using regular expression
  
## References
Humidex estimation references:

- K.R. Spangler, S. Liang, and G.A. Wellenius. "Wet-Bulb Globe Temperature, Universal Thermal Climate Index, and Other Heat Metrics for US Counties, 2000-2020." Scientific Data (2022). doi: 10.1038/s41597-022-01405-3
- Lawrence, M. G. The relationship between relative humidity and the dewpoint temperature in moist air - A simple conversion and applications. B. Am. Meteorol. Soc. 86, 225–233, https://doi.org/10.1175/Bams-86-2-225 (2005).

Heat index (weathermetrics package) reference:

- Anderson GB, Bell ML, Peng RD. 2013. Methods to calculate the heat index as an exposure metric in environmental health research. Environmental Health Perspectives 121(10):1111-1119.

