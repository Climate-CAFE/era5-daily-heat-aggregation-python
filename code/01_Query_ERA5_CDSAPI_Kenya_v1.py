# Date Created: 10/8/2024
# Version Number: v1
# Overview:
#     This code uses the CDS API to download ERA5 temperature measures
#     from The Copernicus Climate Data Store. More details on the CDS API
#     setup, installation and usage are accessible at:
#     https://cds.climate.copernicus.eu/how-to-api#install-the-cds-api-client
#     Please explore the website provided and note the provided instructions
#     for how to sign up for a Copernicus account and to access the relevant
#     user and key inputs below required to query data.
#
# Load required packages
#
# Note: need to install the latest version of packages to meet version requirement. 
# Version 0.7.2 or higher for the cdsapi package is required in order to be able to 
# use the new data stores.

import cdsapi 
import geopandas as gpd
import os

# Set directory. Establishing this at the start of the script is useful in case
# future adjustments are made to the path. There are subdirectories within this
# directory for reading in Global Administrative boundaries and for outputting
# the ERA5 rasters queried in this script.

#ecmw_dir =  "D:\\CAFE_DATA_MANAGEMENT\\ERA5_Python\\era5_daily_heat_aggregation"
ecmw_dir = "YOUR LOCAL PATH"

# You don't need to set key for web API in the program. This is becuase the web API
# url and your key have been saved to the file C:\Users\YOURUSERNAME\.cdsapirc. When
# you make a request, the program can automatically find those information from that
# file. Note: your user ID and API key should never be shared externally.

# Do these steps to create account, get API key, and agree to the Term of Use:
# 1. Create an ECMWF account by self-registering. Go to: https://www.ecmwf.int/. Click Log in.
# 2. Once you log in, go to Climate Data Store webpage: https://cds.climate.copernicus.eu/.
# Click your log-in icon, then click "Your profile" to get user ID and key. 
# 3. Visit user profile to accept the terms and conditions in the profile page. 
# 4. One must agree to the Terms of Use of a dataset before downloading any data out of it. 
# This step must be done manually from the dataset page (at the bottom of the download form).
# For example, Go to ERA-5 land data Download page 
# https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=download.
# Go to "Terms of use" block to accept the data licence to use products.
# 5. Visit user profile page again to double check that Dataset licences to use Copernicus 
# products shows up there and has been accepted. 

# Read in Kenya boundaries geopackage. These data were downloaded from GADM
# https://gadm.org/download_country.html. The layer specification "ADM_ADM_0"
# indicates that the level 0 boundaries should be read in, representing the 
# national boundary. We read in the boundaries here to establish the geographic
# extent that should be queried from ERA5, so ward-level boundaries are not 
# needed.

# Note: need to first create a subfolder "Kenya_GADM" on your path, and put 
# gadm41_KEN.gpkg file in that subfolder.

kenya_shape =  gpd.read_file(os.path.join(ecmw_dir, "Kenya_GADM\\gadm41_KEN.gpkg"), layer = "ADM_ADM_0")

# Assess bounding box. The bounding box represents the coordinates of the 
# extent of the shapefile, and will be used to specify the area we would like
# to query from Copernicus Climate Data Store. The API will allow any bounding 
# parameters; however, values that deviate from the original model grid scale
# will be interpolated onto a new grid. Therefore, it’s recommended that for 
# ERA5-Land (which is 0.1˚ resolution) the bounding coordinates be divisible by 
# 0.1 (e.g., 49.5˚N, -66.8˚E, etc.), and that coordinates for ERA5 be divisible
# by 0.25 (e.g., 49.25˚N, -66.75˚E, etc.)

kenya_bbox = kenya_shape.total_bounds

# Add a small buffer around the bounding box to ensure the whole region 
# is queried, and round the parameters to a 0.1 resolution. A 0.1 resolution
# is applied because the resolution of netCDF ERA5 data is .25x.25
# https://confluence.ecmwf.int/display/CKB/ERA5%3A+What+is+the+spatial+reference

kenya_bbox[0] = round(kenya_bbox[0], 1) - 0.1
kenya_bbox[1] = round(kenya_bbox[1], 1) - 0.1
kenya_bbox[2] = round(kenya_bbox[2], 1) + 0.1
kenya_bbox[3] = round(kenya_bbox[3], 1) + 0.1

# The ERA5 area query requires the area to be formatted specifically as below,
# as a list of xmin, ymin, xmax, ymax. 

query_area = [kenya_bbox[0], kenya_bbox[1], kenya_bbox[2], kenya_bbox[3]]

# The set of inputs below specify the range of years to request, and set each
# month throughout the year with which API requests will be cut. This is because
# CDS API only allows to request a single month of data at most, since there is 
# a limit on the data size that can be downloaded in a given request. So we loop
# each year, and within each year we loop each month to ensure the request
# is completed.

query_years = list(range(2000, 2024))
query_years_str = [str(x) for x in query_years]

query_months = list(range(1, 13))
query_months_str = [str(x).zfill(2) for x in query_months]

# Use for loop to query in each month blocks, by year

for year_str in query_years_str:
    # Track progress
    print("Now processing year ", year_str, "\n")

    # For each year, the query is divided into each month sections. 
    # If a request is too large, it will not be accepted by the CDS servers, 
    # so this division of requests is required.

    for month_str in query_months_str:
        # Track progress
        print("Now processing month ", month_str, "\n")

        # The below is the formatted API request language. All of the inputs
        # specified below in proper formatting can be identified by forming a 
        # request using the Copernicus CDS point-and-click interface for data
        # requests. https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form
        # Select the variables, timing, and netcdf as the output format, and then 
        # select "Show API Request" at the bottom of the screen. 
    
        # Note that the argument in the download() function is the file path and 
        # file name that data will be exported to and stored at. If using a loop, 
        # ensure that the unique features of each request are noted in the output.

        # Note: need to create "ERA5_Out" subfolder on your path
 
        dataset = "reanalysis-era5-land"
        request = {
                    "product_type": "reanalysis",
                    "variable": ["2m_dewpoint_temperature",
                                 "2m_temperature",
                                 "skin_temperature"], 
                    "year": year_str,
                    "month": month_str,
                    "day": [  
                            "01", "02", "03",
                            "04", "05", "06",
                            "07", "08", "09",
                            "10", "11", "12",
                            "13", "14", "15",
                            "16", "17", "18",
                            "19", "20", "21",
                            "22", "23", "24",
                            "25", "26", "27",
                            "28", "29", "30",
                            "31"],
                    "time": [
                             "00:00", "01:00", "02:00",
                             "03:00", "04:00", "05:00",
                             "06:00", "07:00", "08:00",
                             "09:00", "10:00", "11:00",
                             "12:00", "13:00", "14:00",
                             "15:00", "16:00", "17:00",
                             "18:00", "19:00", "20:00",
                             "21:00", "22:00", "23:00"],
                    "data_format": "netcdf",
                    "download_format": "unarchived",
                    "area": query_area
        }

        client = cdsapi.Client()
        client.retrieve(dataset, request).download(os.path.join(ecmw_dir,
                                                                "ERA5_Out", 
                                                                "{}_{}.nc".format(year_str, month_str)))

# The ERA5 data is distributed in UTC. We want to calculate our daily measures
# based on Kenya local time. To accommodate this, we will query the last three 
# hours of 1999 as Kenya is 3 hours ahead of UTC.

dataset = "reanalysis-era5-land"
request_1999 = {
                "product_type": "reanalysis",
                "variable": ["2m_dewpoint_temperature",
                             "2m_temperature",
                             "skin_temperature"], 
                "year": "1999",
                "month": "12",
                "day": ["31"],
                "time": ["21:00", "22:00", "23:00"],
                "data_format": "netcdf",
                "download_format": "unarchived",
                "area": query_area
}


client = cdsapi.Client()
client.retrieve(dataset, request_1999).download(os.path.join(ecmw_dir,
                                                             "ERA5_Out", 
                                                             "1999_12-31.nc"))

