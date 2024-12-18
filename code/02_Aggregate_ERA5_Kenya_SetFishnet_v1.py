# Date Created: 10/15/2024
# Version Number: v1
# Date Modified: 
# Modifications:
# ************************************************************** #
# ~~~~~~~  ERA5 Re-Analysis Raster Processing Step 1     ~~~~~~~ #
# ************************************************************** #
## Purpose: Process ERA5 rasters to administrative boundaries. This
##    script is the first in a two-step raster processing process. In this 
##    a grid-based polygon will be derived from the raster grid of ERA5 data.
##    An example is provided for Kenya data.
##    
## Overall Processing Steps:
##    Script: 02_Aggregate_ERA5_Kenya_SetFishnet_v1.py
## 1) Create Fishnet that can be used to extract ERA5 data from raster stack
##    including ERA5 hourly data (this file). 
##    This will allow for extraction from raster stack  
##    without the large computational burden of a loop (as below)
##
##    Script: 03_Aggregate_ERA5_Kenya_UseFishnet_v1.py
## 2) Load administrative boundaries
## 3) Create extraction points from the union of the block and fishnet. These 
##    are what we can use to extract values from the raster that overlaps with
##    the points aligning to each block (next file).
## 4) Estimate the ward-level exposure to ERA5, accounting for the availability
##    of data within the block (next file).

# Note need to install the latest version of packages to meet version requirement.
# To install GDAL for Python on Windows, first download binaries of gdal and python
# bindings from https://github.com/cgohlke/geospatial-wheels/releases. Then type:
# pip install C:\Users\liuji\Downloads\GDAL-3.9.2-cp310-cp310-win_amd64.whl, for
# example for GDAL version 3.9.2, python version 3.10, and 64 bits on Windows.

import geopandas as gpd
import os 
import glob
import xarray
import rioxarray
import shapely
from shapely.geometry import Polygon
import numpy
from osgeo import gdal, ogr

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
    or numpy.__version__ >= "2.0.0" or gdal.__version__ < "3.9.2":

    print("WARNING: packages are outdated and may result in errors.")

# %%%%%%%%%%%%%%%%%%%%%%% USER-DEFINED FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# Create a fishnet grid using GDAL
def make_fishnet(outputGridfn,xmin,xmax,ymin,ymax,rows,cols):
    # Calculate grid parameters
    xmin = float(xmin)
    xmax = float(xmax)
    ymin = float(ymin)
    ymax = float(ymax)
    gridWidth = float((xmax-xmin) / cols)
    gridHeight = float((ymax-ymin) / rows)

    # Start grid cell envelope
    ringXleftOrigin = xmin
    ringXrightOrigin = xmin + gridWidth
    ringYtopOrigin = ymax
    ringYbottomOrigin = ymax-gridHeight

    # Create the output shapefile
    outDriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(outputGridfn):
        os.remove(outputGridfn)
    outDataSource = outDriver.CreateDataSource(outputGridfn)
    outLayer = outDataSource.CreateLayer(outputGridfn, geom_type=ogr.wkbPolygon)
    # Add fields to the layer
    featureDefn = outLayer.GetLayerDefn()

    # Create grid cells
    countcols = 0
    while countcols < cols:
        countcols += 1

        # Reset envelope for rows
        ringYtop = ringYtopOrigin
        ringYbottom =ringYbottomOrigin
        countrows = 0

        while countrows < rows:
            countrows += 1
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)

            # Add new geom to layer
            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(poly)
            outLayer.CreateFeature(outFeature)
            outFeature = None

            # New envelope for next poly
            ringYtop = ringYtop - gridHeight
            ringYbottom = ringYbottom - gridHeight

        # New envelope for next poly
        ringXleftOrigin = ringXleftOrigin + gridWidth
        ringXrightOrigin = ringXrightOrigin + gridWidth

    # Save and close DataSources
    outDataSource = None


# %%%%%%%%%%%%%%%%%%%%%%% USER-DEFINED PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# Set up directories to read in and output data

#era_dir = "D:\\CAFE_DATA_MANAGEMENT\\ERA5_Python\\era5_daily_heat_aggregation\\ERA5_Out"
era_dir = "YOUR LOCAL PATH TO DOWNLOADED .NC FILES"

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%% CREATE ERA5 EXTRACTION POINTS  %%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#
# %%%%%%%%%%%%%%%%%%%%%%%%% READ IN THE ERA5 RASTER DATA %%%%%%%%%%%%%%%%%%%%% #

era_files = glob.glob(os.path.join(era_dir, '*.nc'))

# Stack all of the hourly files by year

era_stack = xarray.open_mfdataset(era_files, decode_coords="all")

# When the file has latitude and longitude coordinates, the grid_mapping parameter is None.
# To assign a CRS (Coordinate Reference System) to the "grid_mapping" attribute of longitude
# and latitude data, you typically need to specify the CRS code, usually using the EPSG code 
# for "WGS84" (EPSG:4326), within the relevant metadata associated with your data, often using
# a "crs" variable that references the grid_mapping attribute in your data format (like NetCDF) 
# according to CF conventions; essentially, you're telling the system that your longitude and 
# latitude values are based on the WGS84 standard. 

era_stack.rio.write_crs("WGS 84", inplace=True)

# %%%%%%%%%%%%%%%%%%%% CREATE A FISHNET GRID OF THE RASTER EXTENT %%%%%%%%%%%% #
#
# Here, we are making a shapefile that is a fishnet grid of the raster extent.
# It will essentially be a polygon of lines surrounding each ERA5 cell.
#
# Reference/credit: https://gis.stackexchange.com/a/243585

era_extent = era_stack.rio.bounds()

xmin = era_extent[0]
xmax = era_extent[2]
ymin = era_extent[1]
ymax = era_extent[3]

height = era_stack.rio.height
width = era_stack.rio.width

era_coords = [(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)]
era_polygon = Polygon(era_coords)

# Create fishnet of the ERA5 polygon. This takes some time.
# And write fishnet to output for use in later scripts at the same time.

make_fishnet(os.path.join(era_dir, 'era_fishnet.shp'), xmin,xmax,ymin,ymax,height,width)

# Read in the created fishnet shape file and check for its crs. In this case, crs is None.

ogr_shp = gpd.read_file(os.path.join(era_dir, 'era_fishnet.shp'))
print(ogr_shp.crs)

# Set the CRS for the created fishnet shape file as the same one from the stacked raster files. 
new_crs = era_stack.rio.crs.data  # Replace with the desired CRS
ogr_shp = ogr_shp.set_crs(new_crs)

# Save the shapefile with the new CRS
ogr_shp.to_file(os.path.join(era_dir, 'era_fishnet.shp'))



