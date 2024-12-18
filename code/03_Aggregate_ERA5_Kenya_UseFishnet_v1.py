# Date Created: 12/5/2024
# Version Number: v1
# Date Modified: 
# Modifications: 
# ************************************************************** #
# ~~~~~  ERA5 Re-Analysis Raster Processing Step 2        ~~~~~~ #
# ************************************************************** #
## Purpose: Process ERA5 rasters to Kenya administrative boundaries (wards). This
##    script is the first in a two-step raster processing process. In this 
##    a grid-based polygon will be derived from the raster grid of ERA5 data.
##    
## Overall Processing Steps:
##    Script: 02_Aggregate_ERA5_Kenya_SetFishnet_v1.py
## 1) Create Fishnet that can be used to extract ERA5 data from raster stack
##    including ERA5 hourly data (this file). 
##    This will allow for extraction from raster stack  
##    without the large computational burden of a terra::zonal loop (as below)
##
##    Script: 03_Aggregate_ERA5_Kenya_UseFishnet_v1.py
## 2) In this file, we will join the fishnet, which is a polygon grid with lines 
##    surrounding the grid of the ERA5 raster with the Ward geometries that have
##    been queried from the Database of Global Administrative Boundaries
##    (gadm.org). In merging the polygon grid with the ward polygon data, 
##    we ensure that every ward will be aligned with the relevant
##    ERA5 temperature metrics for the area.
## 3) Create extraction points from the union of the wards and fishnet. These 
##    are what we can use to extract values from the raster that overlaps with
##    with the points aligning to each wards (this file).
## 4) Estimate the ward-level exposure to ERA5, accounting for the availability
##    of data within the wards (this file).

import geopandas as gpd
import shapely
from shapely.validation import make_valid, explain_validity
import pandas as pd
import numpy as np
from functools import reduce
import glob
import math
import os 
import xarray
import rioxarray
import pytz
from pytz import timezone
from datetime import datetime
import xvec 
from functools import reduce
import re


# S2 is for computing distances, areas, etc. on a SPHERE (using
# geographic coordinates, i.e., lat/lon in decimal-degrees); no need for this
# extra computational processing time if using PROJECTED coordinates,
# since these are already mapped to a flat surface. Here, ERA5 data
# is indeed in geographic coordinates, but the scale of areas we are 
# interested in is very small, and hence the error introduced by 
# ignoring the Earth's curvature over these tiny areas is negligible and
# a reasonable trade off given the dramatic reduction in processing time. Moreover,
# the areas we calculate are not an integral part of the process
# and any error in that step would not materially impact the final output.
# Python GeoPandas package uses Shapely under the hood, which is a library
# for manipulating and analyzing planar not spherical geometric objects.

# Check package version numbers

if gpd.__version__ < "1.0.1" or xarray.__version__ < "2024.9.0" \
    or rioxarray.__version__ < "0.17.0" or shapely.__version__ < "2.0.6" \
    or np.__version__ >= "2.0.0" or xvec.__version__ < "0.3.0":

	print("WARNING: packages are outdated and may result in errors.")

# Set up directories to read in and output data
#

era_dir = "D:\\CAFE_DATA_MANAGEMENT\\ERA5_Python\\era5_daily_heat_aggregation\\ERA5_Out"
#era_dir = "YOUR LOCAL PATH TO CREATED FISHNET SHAPE FILE"
geo_dir = "D:\\CAFE_DATA_MANAGEMENT\\ERA5_Python\\era5_daily_heat_aggregation\\GEO"
#geo_dir = "YOUR LOCAL PATH"
outdir = "D:\\CAFE_DATA_MANAGEMENT\\ERA5_Python\\era5_daily_heat_aggregation\\OUT"
#outdir = "YOUR LOCAL PATH"

# %%%%%%%%%%%%%%%%%%%%%% LOAD WARDS SHAPEFILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#
# Read in wards shapefile
#

##### Go to https://gadm.org/download_country.html, download a country shapefile.
##### Downloading as a gpkg will embed multiple layers of shapefiles. A single
##### layer can be read in based on the "layer =" input. Here, we want to use 
##### the administrative boundaries at which the final metrics will be available.
##### We use the smallest boundaries available here (level 3) for Kenya

country_wards = gpd.read_file(geo_dir+"\\gadm41_KEN.gpkg", layer = "ADM_ADM_3") 

# Rename geom column to geometry if necessary, to align with fishnet
#
if "geometry" not in country_wards.columns: 
  country_wards = country_wards.rename(columns = {"geom": "geometry"})

# Read in the fishnet created previously
#
era_fishnet = gpd.read_file(era_dir+"\\era_fishnet.shp")

# Match the CRS of the wards shapefile to the fishnet and era data and confirm match
#
country_wards = country_wards.to_crs(era_fishnet.crs)

# Run check to ensure CRS are equal
#

assert country_wards.crs == era_fishnet.crs, "ERROR: CRS's don't match."

# %%%%%%%%%%%%%%%%%% CREATE UNION BETWEEN FISHNET AND WARDS %%%%%%%%% #
# Reference/credit: https://stackoverflow.com/a/68713743
#

def my_union_a(a, b):
	b_combined = b.dissolve(by=None, aggfunc="first")
	op1_geo = a['geometry'].difference(b_combined['geometry'][0]) 
	op1 = a
	op1['geometry']=op1_geo
	a_combined = a.dissolve(by=None, aggfunc="first")
	op2_geo = b['geometry'].difference(a_combined['geometry'][0])
	op2 = b
	op2['geometry'] = op2_geo
	op3 = gpd.overlay(b, a, how='intersection') 
	union = pd.concat([op1, op2, op3], ignore_index=True)
	return gpd.GeoDataFrame(union)

# Ensure geometries are valid
#
country_wards.geometry = country_wards.apply(lambda row: make_valid(row.geometry) if not explain_validity(row.geometry) == 'Valid Geometry' else row.geometry, axis=1)
era_fishnet.geometry = era_fishnet.apply(lambda row: make_valid(row.geometry) if not explain_validity(row.geometry) == 'Valid Geometry' else row.geometry, axis=1)

# Create the union between the fishnet and blocks layers
#
fishnetward = my_union_a(era_fishnet, country_wards)
fishnetward['UniqueID'] = list(range(1, fishnetward.shape[0]+1))
#fishnetward.to_csv('D:\\CAFE_DATA_MANAGEMENT\\ERA5_Python\\fishnetward.csv')

# Automated QC -- Check to see if the union has introduced any geometry errors
#                 and fix as appropriate
#
try:
	era_fishnet.apply(lambda row: make_valid(row.geometry) if not explain_validity(row.geometry) == 'Valid Geometry' else row.geometry, axis=1)
except Exception as e:
	print("An error occurred:", e)
	print("There is an issue with the geopandas object \n")
	print("..... Attempting fix \n")

	geo_types = fishnetward['geometry'].geom_type.unique()
	print("..... Geometry types in geopandas object 'fishnetward':", geo_types, "\n")

	for j in range(len(geo_types)):
		fishnetward_subset = fishnetward[fishnetward['geometry'].geom_type == geo_types[j]]
		if j == 0:
			updated_fishnetward = fishnetward_subset
			continue
		pd.concat([updated_fishnetward, fishnetward_subset], ignore_index=True)

	try:
		updated_fishnetward.apply(lambda row: make_valid(row.geometry) if not explain_validity(row.geometry) == 'Valid Geometry' else row.geometry, axis=1)
	except Exception as e:
		print("An error occurred:", e)
		print("..... ERROR NOT RESOLVED")
	else:
		print("..... :) issue has been fixed!")

		updated_fishnetward = updated_fishnetward.sort_values(by="UniqueID")
		result = np.isclose(updated_fishnetward['UniqueID'], fishnetward['UniqueID'])
		if all(result):
			print(":) unique ID's match. Reassigning 'updated_fishnetward' to 'fishnetward'")
			fishnetward <- updated_fishnetward    
		else:
			print("ERROR: Unique IDs do not match") 
else:
	print('Nothing went wrong')


# Automated QC check -- ensure that there is a variable identifying the geographies
#                       that will be used for aggregation (here Kenya wards). 
#                       Note that these variable names may change
#                       depending on the version and country used. 
# Specify name of the geographic identifier in your administrative boundary data
#
geo_name = "GID_3"

# Extract name from dataset
#
geo_id_var = fishnetward.filter(like=geo_name).columns.to_list()[0]

# Identify the polygons of the fishnet that do not intersect with the ward
# data; drop them.
#
before_dim = fishnetward.shape[0]
fishnetward = fishnetward[fishnetward[geo_id_var].notna()]
after_dim = fishnetward.shape[0]

print("Dropped", before_dim - after_dim, "polygons that do not intersect with census data")

# Some polygons formed in the union are incredibly small -- this adds unnecessary
# computation time without materially reducing error. Drop the small polygons.
# NOTE: Typically, when calculating areas of polygons, you would want to convert to
#       a projected CRS appropriate for your study domain. For the purpose of identifying
#       negligibly small areas to drop here, the error introduced by using geographic
#       coordinates for calculating area at this scale is negligible.
#
fishnetward = fishnetward.to_crs(3857)
fishnetward['Area_m2'] = fishnetward['geometry'].area.astype(float)

fishnetward = fishnetward[fishnetward['Area_m2'] > 10]
fishnetward = fishnetward.to_crs(4326)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%% CONVERT POLYGON TO POINTS %%%%%%%%%%%%%%%%%%%%% #
#
# The final step is to create the extraction points. This is a point shapefile
# that will enable us to extract ERA5 data from an entire stack of rasters rather
# than individually processing zonal statistics on each raster layer. 
#
# NOTE: This step throws a warning message related to using geographic coordinates 
#       rather than a projected CRS. This step is only placing a point inside the 
#       polygon to identify which ERA5 grid cell we need to extract from; as all
#       of the input data are on the same CRS and the spatial scale of the polygons
#       is extremely small, this does not introduce substantive error.
#
fishnetward['geometry'] = fishnetward['geometry'].representative_point()
extraction_pts = fishnetward

# %%%%%%%%%%%%%%%%%%%%% CALCULATE LAND-AREA WEIGHTED AVG BY WARD %%%%%%%%%%%% #

# Get the total area by ward to calculate the spatial weight value (typically 1.0)
#

# Define a function to calculate sums such that if all values are NA then it returns
# NA rather than 0.
#
def sumfun(x):
	if x.isnull().all():
		return None
	else:
		return x.sum(skipna=True)

# Calculate the sum area by ward
#
ptstotal = extraction_pts.groupby(geo_id_var, as_index=False).agg({'Area_m2': lambda x: sumfun(x)})
ptstotal = ptstotal.rename(columns={'Area_m2': 'Area_m2.sumfun'})

# Merge area and calculate spatial weight of points
#

data_frames = [extraction_pts, ptstotal]
extraction_pts = reduce(lambda  left,right: pd.merge(left,right,on=[geo_id_var], how='left'), data_frames)

extraction_pts['SpatWt'] = extraction_pts['Area_m2'] / extraction_pts['Area_m2.sumfun']

# %%%%%%%%%%%%%%%%%%%%%% EXTRACT ERA5 VALUES FOR EACH WARD POINT %%%%%%%%%%%% #
# 
# In this step, we will use the extraction points to extract the ERA5 grid cell
# underlying each portion of a geography (ward) across the entire raster stack of values.
#
# This process is set to run as a loop by year to cut down the size of the 
# file that will be produced by the processing. For each year, and in using 
# a cluster computing environment, the processing takes about 20 minutes and uses
# about 15GB of R memory on a computing cluster. 

# This can be further reduced by running as a loop for
# each month. The hourly data are transformed to a long format for processing of daily
# summary measures, and included the full 20+ year period would lead to a 
# very large dataframe (1400 polygons * 365 days * 24 hours * X years)
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%% READ IN THE ERA5 RASTER DATA %%%%%%%%%%%%%%%%%% #
#
era_files = glob.glob(os.path.join(era_dir, '*.nc'))

# Humidex function from heatmetrics package:
# Reference: K.R. Spangler, S. Liang, and G.A. Wellenius. "Wet-Bulb Globe 
# Temperature, Universal Thermal Climate Index, and Other Heat Metrics for US 
# Counties, 2000-2020." Scientific Data (2022). doi: 10.1038/s41597-022-01405-3
#
# Lawrence, M. G. The relationship between relative humidity and the dewpoint 
# temperature in moist air - A simple conversion and applications. B. Am. 
# Meteorol. Soc. 86, 225–233, https://doi.org/10.1175/Bams-86-2-225 (2005).
#

def humidex(t, td):
	hx = (t + (5/9) * ((6.1094 * math.exp((17.625 * td) / (243.04 + td))) - 10))
	return hx

def heatindex(t, td):
	f_t = (t * 9/5) + 32
	humidity = (6.112 * math.exp((17.67*td)/(td+243.5))) / (6.112 * math.exp((17.67*t)/(t+243.5))) * 100
	f_heat_index = -42.379 + 2.04901523*f_t + 10.14333127*humidity + (-0.22475541)*f_t*humidity + \
		(-6.83783*0.001)*f_t*f_t + (-5.481717*0.01)*humidity*humidity + \
		(1.22874*0.001)*f_t*f_t*humidity + (8.5282*0.0001)*f_t*humidity*humidity + \
		(-1.99*0.000001)*f_t*f_t*humidity*humidity
	heat_index = (f_heat_index - 32) * 5/9
	return heat_index

# Time zone specification. 
# The SpatRaster as downloaded from Copernicus will include hourly data based on
# the UTC time zone. When we calculate our daily summary statistics in the loop
# below, we want to make sure we are calculating statistics from midnight to
# midnight, using the local time zone. Specify the time zone relevant for your 
# data below. Note: If you are attempting to aggregate across multiple distinct
# time zones, additional processing is necessary. This capability is in development
# for an additional version of this process. To see a list of all time zones,
# run: OlsonNames()
#
tz_country = "Africa/Nairobi"

# Set year to process
#
years_to_agg = list(range(2000, 2001))

for year in years_to_agg: 
	print("Now processing {}\n".format(year))

	# Subset to year from all files 
	#
	era_files_yr = [file for file in era_files if str(year) or str(year-1) in file]
	
	# Stack all of the daily files by year
	#
	era_stack = xarray.open_mfdataset(era_files_yr, decode_coords="all")

	# Reset times to align with specified time zone
	#

	country_tz = timezone(tz_country)

	valid_time_list = list()

	for d in era_stack.valid_time.data:
		mytime = datetime.utcfromtimestamp(d.tolist()/1e9)
		utctz = pytz.utc
		myutctime = utctz.localize(mytime)
		mylocaltime = myutctime.astimezone(country_tz)
		valid_time_list.append(mylocaltime)
	valid_time_list = np.array(valid_time_list)

	era_stack = era_stack.assign_coords(valid_time=valid_time_list)

	# Subset stack to exclude the times that run past specified year due to 
	# time zone adjustment, and to exclude the previous year that was read in 
	# for time zone adjustment
	#
	end_basetime = datetime.strptime(str(year)+"-12-31"+" 21:00:00", "%Y-%m-%d %H:%M:%S")
	utctz = pytz.utc
	end_basetime_utc = utctz.localize(end_basetime)
	end_basetime_local = end_basetime_utc.astimezone(country_tz)

	start_basetime = datetime.strptime(str(year-1)+"-12-31"+" 21:00:00", "%Y-%m-%d %H:%M:%S")
	utctz = pytz.utc
	start_basetime_utc=utctz.localize(start_basetime)
	start_basetime_local = start_basetime_utc.astimezone(country_tz)

	mask = (era_stack.valid_time < pd.Timestamp(end_basetime_local)) & (era_stack.valid_time >= pd.Timestamp(start_basetime_local))
	era_stack = era_stack.where(mask, drop=True)

	# Our ERA stack includes three variables (2m dew point temperature,
	# skin temperature, and 2m temperature). Create subsets to perform
	# daily aggregation on
	#
	era_stack_d2m = era_stack['d2m']
	era_stack_t2m = era_stack['t2m']
	era_stack_skt = era_stack['skt']

	# Subset raster to each measure and convert Kelvin to Celsius
	#
	era_stack_d2m = era_stack_d2m - 273.15
	era_stack_t2m = era_stack_t2m - 273.15
	era_stack_skt = era_stack_skt - 273.15

	# Apply the heat index function to the SpatRaster objects
	#
	f1 = np.vectorize(heatindex)
	era_stack_hti = xarray.apply_ufunc(f1, era_stack_t2m, era_stack_d2m, dask="parallelized")
	
	print('I am here.')

	# Apply the humidex function to the SpatRaster objects
	#
	
	f2 = np.vectorize(humidex)
	era_stack_hum = xarray.apply_ufunc(f2, era_stack_t2m, era_stack_d2m, dask="parallelized")

	print('I am here.')

	# Assign names to new layers
	era_stack_hti = era_stack_hti.rename('hti')
	era_stack_hum = era_stack_hum.rename('hum')

	# Confirm all layers are same length
	#
	res = all(np.equal(era_stack_d2m.shape, era_stack_t2m.shape)) & \
		all(np.equal(era_stack_d2m.shape, era_stack_skt.shape)) & \
		all(np.equal(era_stack_d2m.shape, era_stack_hti.shape)) & \
		all(np.equal(era_stack_d2m.shape, era_stack_hum.shape))
	if res: 
		layer_shape = era_stack_d2m.shape
		print("Same number of layers in all stacks\n")
	else:
		print("Different number of layers, assess whether timing is consistent\n")
		break
 
	# Create a time sequence starting from January first of the year date
	#
	start_date = pd.to_datetime(str(year-1)+"-12-31 21:00", utc=True).tz_convert(country_tz)

	time_seq = pd.date_range(start = start_date, periods = 8784, freq = "h")

	# Convert our time sequence to a factor format. This will allow for use as 
	# a grouping variable in assessing daily level summary measures
	#
	time_seq_ind = xarray.DataArray(time_seq.date, dims=['valid_time'])

	# Set list of rasters. We will do the same processing of daily minimum, mean,
	# and maximum for the three variables we queried through ERA5, so we set
	# the lasters in a list and conduct the processing as below
	#
	list_rasters = [era_stack_d2m, era_stack_t2m, era_stack_skt, era_stack_hti, era_stack_hum]

	for i in range(len(list_rasters)):    
		# Tracker for viewing progress
		#
		print("Now processing raster {}\n ".format(i))
    
		# Aggregate to daily mean temperature
		#
		daily_mean = list_rasters[i].groupby(time_seq_ind, restore_coord_dims=True).apply(np.mean)
		daily_mean = daily_mean.expand_dims({'latitude':list_rasters[i]['latitude']})
		daily_mean = daily_mean.expand_dims({'longitude':list_rasters[i]['longitude']})

		print('I am here.')

		# Aggregate to daily maximum temperature
		#
		daily_max = list_rasters[i].groupby(time_seq_ind, restore_coord_dims=True).apply(np.max)
		daily_max = daily_max.expand_dims({'latitude':list_rasters[i]['latitude']})
		daily_max = daily_max.expand_dims({'longitude':list_rasters[i]['longitude']})

		print('I am here.')
		
		# Aggregate to daily minimum temperature
		#
		daily_min = list_rasters[i].groupby(time_seq_ind, restore_coord_dims=True).apply(np.min)
		daily_min = daily_min.expand_dims({'latitude':list_rasters[i]['latitude']})
		daily_min = daily_min.expand_dims({'longitude':list_rasters[i]['longitude']})

		print('I am here.')
		
		# Project points to WGS84 (coordinate system for ERA stack)
		#
		extraction_pts = extraction_pts.set_crs(epsg=4326)

		# Extract daily summaries to point-based grid
		#
		mean_pts = daily_mean.xvec.extract_points(extraction_pts['geometry'], x_coords="longitude", y_coords="latitude")
		max_pts = daily_max.xvec.extract_points(extraction_pts['geometry'], x_coords="longitude", y_coords="latitude")
		min_pts = daily_min.xvec.extract_points(extraction_pts['geometry'], x_coords="longitude", y_coords="latitude")

		print('I am here.') 

		# Join results with extraction points (includes geographic identifiers)
		#
		mean_pts = gpd.GeoDataFrame(mean_pts, geometry=mean_pts.coords['geometry'].values, columns=[date.strftime("%Y-%m-%d") for date in mean_pts.coords['group'].values])
		mean_pts = extraction_pts.merge(mean_pts, on='geometry')
		max_pts = gpd.GeoDataFrame(max_pts, geometry=max_pts.coords['geometry'].values, columns=[date.strftime("%Y-%m-%d") for date in max_pts.coords['group'].values])
		max_pts = extraction_pts.merge(max_pts, on='geometry')
		min_pts = gpd.GeoDataFrame(min_pts, geometry=min_pts.coords['geometry'].values, columns=[date.strftime("%Y-%m-%d") for date in min_pts.coords['group'].values])
		min_pts = extraction_pts.merge(min_pts, on='geometry')

		print('I am here.') 

		# Use package sf to join points on coast and not covered by ERA5 land to
		# nearby temperature measures
		#

		# Identify a date in year to assess missingness for points outside
		# land extent
		#
		varname = str(year)+'-01-01'

		# Subset missing data
		#
		mean_pts_missing = mean_pts[mean_pts[varname].isnull()]
		max_pts_missing = max_pts[max_pts[varname].isnull()]
		min_pts_missing = min_pts[min_pts[varname].isnull()]
    
		# Subset available data
		#
		mean_pts_avail = mean_pts[~mean_pts[varname].isnull()]
		max_pts_avail = max_pts[~max_pts[varname].isnull()]
		min_pts_avail = min_pts[~min_pts[varname].isnull()]

		# Join these with nearest features with available temp data
		#
		names_temp = [column for column in mean_pts_avail.columns if column not in extraction_pts.columns]    
		
		mean_pts_missing_join = mean_pts_missing[extraction_pts.columns]
		mean_pts_avail_join = mean_pts_avail[names_temp+['geometry']]

		try:
			mean_pts_missing = gpd.sjoin_nearest(mean_pts_missing_join, mean_pts_avail_join)
			mean_pts_missing = mean_pts_missing[mean_pts_avail.columns.to_list()]
		except:
			print('Try this 1')
			mean_pts_missing = gpd.GeoDataFrame(columns=mean_pts_avail.columns.to_list())
		
		print(mean_pts_missing.shape)

		max_pts_missing_join = max_pts_missing[extraction_pts.columns]
		max_pts_avail_join = max_pts_avail[names_temp+['geometry']]

		try:
			max_pts_missing = gpd.sjoin_nearest(max_pts_missing_join, max_pts_avail_join)
			max_pts_missing = max_pts_missing[max_pts_avail.columns.to_list()]
		except:
			print('Try this 2')
			max_pts_missing = gpd.GeoDataFrame(columns=max_pts_avail.columns.to_list())

		print(max_pts_missing.shape)

		min_pts_missing_join = min_pts_missing[extraction_pts.columns]
		min_pts_avail_join = min_pts_avail[names_temp+['geometry']]

		try:
			min_pts_missing = gpd.sjoin_nearest(min_pts_missing_join, min_pts_avail_join)
			min_pts_missing = min_pts_missing[min_pts_avail.columns.to_list()]
		except:
			print('Try this 3')
			min_pts_missing = gpd.GeoDataFrame(columns=min_pts_avail.columns.to_list())
		
		print(min_pts_missing.shape)

		# Rejoin updated data with those already including temp
		#
		mean_pts = pd.concat([mean_pts_avail, mean_pts_missing], axis=0)
		max_pts = pd.concat([max_pts_avail, max_pts_missing], axis=0)
		min_pts = pd.concat([min_pts_avail, min_pts_missing], axis=0)

		print(mean_pts.shape)
		print(max_pts.shape)
		print(min_pts.shape)

		# Convert the extracted data to a pandas data frame
		# Remove geometry column
		mean_pts_df = pd.DataFrame(mean_pts.drop(columns='geometry'))
		mean_pts_df['ID'] = range(1, len(mean_pts_df) + 1)
		max_pts_df = pd.DataFrame(max_pts.drop(columns='geometry'))
		max_pts_df['ID'] = range(1, len(max_pts_df) + 1)
		min_pts_df = pd.DataFrame(min_pts.drop(columns='geometry'))
		min_pts_df['ID'] = range(1, len(min_pts_df) + 1)

		# Extract columns relevant to ERA5 data
		era5_cols = [column for column in mean_pts_df.columns if column not in extraction_pts.columns]

		# Set names for ERA5 variables based on input raster naming
		#
		mean_name = list_rasters[i].name + "_mean"
		max_name = list_rasters[i].name + "_max"
		min_name = list_rasters[i].name + "_min"

		id_vars_list = [column for column in extraction_pts.columns.to_list() if column not in ['geometry']]+['ID']
		
		# Transpose the data frame to get time series format (long-form)
		#
		mean_long = pd.melt(mean_pts_df, id_vars = id_vars_list, var_name="date", value_name=mean_name).sort_values(by=[geo_name, 'date']).reset_index(drop=True)
		max_long = pd.melt(max_pts_df, id_vars = id_vars_list, var_name="date", value_name=max_name).sort_values(by=[geo_name, 'date']).reset_index(drop=True)
		max_long = max_long[['UniqueID', 'date', max_name]]
		min_long = pd.melt(min_pts_df, id_vars = id_vars_list, var_name="date", value_name=min_name).sort_values(by=[geo_name, 'date']).reset_index(drop=True)
		min_long = min_long[['UniqueID', 'date', min_name]]

		# Combine data into single dataframe
		#
		data_frames = [mean_long, max_long, min_long]
		era5_long = reduce(lambda  left,right: pd.merge(left,right,on=["UniqueID", "date"], how='left'), data_frames)

		# Join together all measures
		#
		if i == 0:
			era5_full = era5_long
		elif i != 0:
			era5_long = era5_long[["UniqueID", "date", mean_name, max_name, min_name]]
			era5_full = reduce(lambda  left,right: pd.merge(left,right,on=["UniqueID", "date"], how='left'), [era5_full, era5_long])
    
	# In this step, we will use the extraction points to extract the ERA5 grid cell
	# underlying each portion of a ward across the entire raster stack of values.
	# We again will follow the same process for each individual variable, and use 
	# a loop to conduct the processing
	#
	varnames = [column for column in era5_full.columns if column not in extraction_pts.columns.to_list() + ["date", "ID"]]
	
	for i in range(len(varnames)):
		print("Now processing {} \n".format(varnames[i]))
    
		# Before we calculate the final weighted average of the ERA5 measure, we need to check for missing data.
		# If a value is NA on one of the polygons, then it will give an underestimate of the
		# temperature since the weights will no longer add to 1. Example: there are two
		# polygons, each with 50% area. If Tmax is 30 C in one and NA in the other, then
		# the area weighted average (which removes NA values) would give: (30 * 0.5) + (NA * 0.5) = 15 C.
		# Therefore, we need to re-weight the weights based on the availability of data.
		#
		varname = varnames[i]
		era5_full = era5_full[~era5_full[varname].isnull()]
		avail = era5_full.groupby(['GID_3', 'date'], as_index=False).agg({'SpatWt': lambda x: sumfun(x)})
		avail = avail.rename(columns={'SpatWt': 'SpatWt.sumfun'})

		# Merge this value back into the longform ERA5 data
		#
		era5_full = reduce(lambda  left,right: pd.merge(left,right,on=["GID_3", "date"], how='left'), [era5_full, avail])
		
		# Re-weight the area weight by dividing by total available weight
		#
		era5_full['SpatWt'] = era5_full['SpatWt'] / era5_full['SpatWt.sumfun']

		# QC: check that the weights of *available data* all add to 1
		#
		era5_full = era5_full[~era5_full[varname].isnull()]
		check = era5_full.groupby(['GID_3', 'date'], as_index=False).agg({'SpatWt': lambda x: sumfun(x)})
		check = check.rename(columns={'SpatWt': 'SpatWt.sumfun'})


		if len(check[round(check['SpatWt.sumfun'], 4) != 1]) > 0: 
			print("ERROR: weights do not sum to 1\n") 
			break 
		else:
			print(":) weights sum to 1\n")
			era5_full=era5_full.drop('SpatWt.sumfun', axis=1)

		# Multiply the variable of interest (here "newvarname") by the weighting value and then
		# sum up the resultant values within admin boundaries. This is an area-weighted average.
		#
		tempvar = varname+"_Wt"
		era5_full[tempvar] = era5_full[varname] * era5_full["SpatWt"]

		final = era5_full.groupby(['GID_3', 'date'], as_index=False).agg({tempvar: lambda x: sumfun(x)})
		
		# Automated QC to confirm that the dimensions are correct
		#
		if len(extraction_pts['GID_3'].unique())  * len(era5_full['date'].unique()) != final.shape[0]:
			print("ERROR: incorrect dimensions of final df\n")
			break
		else:
			print(":) dimensions of final df are as expected\n")
    

		# Set name for output
		#
		matched_names = [x for x in final.columns if bool(re.search(varname, x))]
		for name in matched_names:
			final = final.rename(columns={name: varname})

		if i == 0:
			finaloutput = final
		else:
			finaloutput = reduce(lambda  left,right: pd.merge(left,right,on=["GID_3", "date"], how='left'), [finaloutput, final])

	print("The final output has {} rows.\n".format(finaloutput.shape[0]))
	print("The first few lines of the output are: \n")
	print(finaloutput.head(10))

	# Automated QC: missing data
	#
	missing_t2m_max = finaloutput[finaloutput['t2m_max'].isnull()]
	missing_t2m_min = finaloutput[finaloutput['t2m_min'].isnull()]
	missing_t2m_mean = finaloutput[finaloutput['t2m_mean'].isnull()]
	missing_d2m_max = finaloutput[finaloutput['d2m_max'].isnull()]
	missing_d2m_min = finaloutput[finaloutput['d2m_min'].isnull()]
	missing_d2m_mean = finaloutput[finaloutput['d2m_mean'].isnull()]
	missing_skt_max = finaloutput[finaloutput['skt_max'].isnull()]
	missing_skt_min = finaloutput[finaloutput['skt_min'].isnull()]
	missing_skt_mean = finaloutput[finaloutput['skt_mean'].isnull()]
	missing_hti_max = finaloutput[finaloutput['hti_max'].isnull()]
	missing_hti_min = finaloutput[finaloutput['hti_min'].isnull()]
	missing_hti_mean = finaloutput[finaloutput['hti_mean'].isnull()]
	missing_hum_max = finaloutput[finaloutput['hum_max'].isnull()]
	missing_hum_min = finaloutput[finaloutput['hum_min'].isnull()]
	missing_hum_mean = finaloutput[finaloutput['hum_mean'].isnull()]

	if len(missing_t2m_max) > 0 | len(missing_d2m_max) > 0 | len(missing_skt_max) > 0 \
      | len(missing_hti_max) > 0 | len(missing_hum_max) > 0:
		print("WARNING: Note the number of missing ward-days by variable: \n")
		print("Dew Max: {}\n".format(len(missing_d2m_max)))
		print("Temp Max: {}\n".format(len(missing_t2m_max)))
		print("Skin Temp Max: {}\n".format(len(missing_skt_max)))
		print("Heat Index: {}\n".format(len(missing_hti_max)))
		print("Humidex: {}\n".format(len(missing_hum_max)))
    
		print("The first few lines of missing T2M_max (if any) are printed below: \n")
		print(missing_t2m_max.head())

		print("The first few lines of missing D2M_max (if any) are printed below: \n")
		print(missing_d2m_max.head())

		print("The first few lines of missing SKT_max (if any) are printed below: \n")
		print(missing_skt_max.head())

		print("The first few lines of missing HTI_max (if any) are printed below: \n")
		print(missing_hti_max.head())

		print("The first few lines of missing HUM_max (if any) are printed below: \n")
		print(missing_hum_max.head())
	else: 
		print(":) No missing temperature values! \n") 

	# Automated QC: impossible temperature values
	#
	num_temp_errors_t2m = len(finaloutput[(finaloutput['t2m_max'] < finaloutput['t2m_mean']) | 
									   (finaloutput['t2m_max'] < finaloutput['t2m_min']) | 
									   (finaloutput['t2m_min'] > finaloutput['t2m_mean']) | 
									   (finaloutput['t2m_min'] > finaloutput['t2m_max'])])

	num_temp_errors_d2m = len(finaloutput[(finaloutput['d2m_max'] < finaloutput['d2m_mean']) | 
									   (finaloutput['d2m_max'] < finaloutput['d2m_min']) | 
									   (finaloutput['d2m_min'] > finaloutput['d2m_mean']) | 
									   (finaloutput['d2m_min'] > finaloutput['d2m_max'])])
  
	num_temp_errors_skt = len(finaloutput[(finaloutput['skt_max'] < finaloutput['skt_mean']) | 
									   (finaloutput['skt_max'] < finaloutput['skt_min']) | 
									   (finaloutput['skt_min'] > finaloutput['skt_mean']) |
									   (finaloutput['skt_min'] > finaloutput['skt_max'])])
  
	num_temp_errors_hti = len(finaloutput[(finaloutput['hti_max'] < finaloutput['hti_mean']) | 
									   (finaloutput['hti_max'] < finaloutput['hti_min']) |
									   (finaloutput['hti_min'] > finaloutput['hti_mean']) |
									   (finaloutput['hti_min'] > finaloutput['hti_max'])])
  
	num_temp_errors_hum = len(finaloutput[(finaloutput['hum_max'] < finaloutput['hum_mean']) |
									   (finaloutput['hum_max'] < finaloutput['hum_min']) |
									   (finaloutput['hum_min'] > finaloutput['hum_mean']) |
									   (finaloutput['hum_min'] > finaloutput['hum_max'])])

	if num_temp_errors_t2m > 0 | num_temp_errors_d2m > 0 | num_temp_errors_skt > 0 | num_temp_errors_hti > 0  | num_temp_errors_hum > 0:
		print("ERROR: impossible temperature values. Applicable rows printed below:")
		print(finaloutput[(finaloutput['t2m_max'] < finaloutput['t2m_mean']) |
					(finaloutput['t2m_max'] < finaloutput['t2m_min']) |
					(finaloutput['t2m_min'] > finaloutput['t2m_mean']) |
					(finaloutput['t2m_min'] > finaloutput['t2m_max'])])

		print(finaloutput[(finaloutput['d2m_max'] < finaloutput['d2m_mean']) |
					(finaloutput['d2m_max'] < finaloutput['d2m_min']) |
					(finaloutput['d2m_min'] > finaloutput['d2m_mean']) |
					(finaloutput['d2m_min'] > finaloutput['d2m_max'])])

		print(finaloutput[(finaloutput['skt_max'] < finaloutput['skt_mean']) |
					(finaloutput['skt_max'] < finaloutput['skt_min']) |
					(finaloutput['skt_min'] > finaloutput['skt_mean']) |
					(finaloutput['skt_min'] > finaloutput['skt_max'])])

		print(finaloutput[(finaloutput['hti_max'] < finaloutput['hti_mean']) |
					(finaloutput['hti_max'] < finaloutput['hti_min']) |
					(finaloutput['hti_min'] > finaloutput['hti_mean']) |
					(finaloutput['hti_min'] > finaloutput['hti_max'])])

		print(finaloutput[(finaloutput['hum_max'] < finaloutput['hum_mean']) |
					(finaloutput['hum_max'] < finaloutput['hum_min']) |
					(finaloutput['hum_min'] > finaloutput['hum_mean']) |
					(finaloutput['hum_min'] > finaloutput['hum_max'])])
	else: 
		print(":) all temperature values are of correct *relative* magnitude") 

	# Output results by year to output directory
	#
	finaloutput.to_csv(os.path.join(outdir, "country_agg_era5_"+str(year)+"_d2m_t2m_skt_hti_hum.csv"))
  
  
  
  
  
  
 