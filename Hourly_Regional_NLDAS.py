#!/usr/bin/env python
# coding: utf-8

# ![HPC Masthead](https://kyrill.ias.sdsmt.edu/wjc/eduresources/AES_519_Masthead.png)
# 
# # Extracting Hourly NLDAS 
# 
# This is a brief overview of how to retrieve gridded data from NASA's EarthAccess Services.
# 
# ---
# ## Get your EarthAccess Account and your Token.
# 
# [https://urs.earthdata.nasa.gov](https://urs.earthdata.nasa.gov)
# 
# ## Create Files for Automated File Logins.
# 
# (do the "manual option" it's easier)
# 
# [https://github.com/nasa/gesdisc-tutorials/blob/main/notebooks/How_to_Generate_Earthdata_Prerequisite_Files.ipynb](https://github.com/nasa/gesdisc-tutorials/blob/main/notebooks/How_to_Generate_Earthdata_Prerequisite_Files.ipynb)
# 
# ---
# ## Step 1: Crack the Libraries
# 
# You need the EarthAccess package, which should now be on both Kernels on the Jupyter Hub.
# 
# [https://earthaccess.readthedocs.io/en/latest/](https://earthaccess.readthedocs.io/en/latest/)

# In[ ]:


#################################################
#
# Libraries
#

import numpy             as np
import xarray            as xr
import matplotlib.pyplot as plt
import pandas            as pd
import h5netcdf          as h5netcdf

import earthaccess       as earthaccess

import platform          as platform
import os                as os

from functools import partial


#
#################################################

#################################################
#
# For clip xarray.open_mfdataset read.
#

def _preprocess(x, lat_clip, lon_clip):
    return x.sel(lon = slice(*lon_clip), 
                 lat = slice(*lat_clip))

#
#################################################


# ---
# ## Step 2: Go shopping!  
# 
# Here's the URL.  
# 
# [https://disc.gsfc.nasa.gov/datasets/](https://disc.gsfc.nasa.gov/datasets/)
# 
# ## Get the data "doi" for the Earthdata catalog.
# 
# This example is the monthly averaged NLDAS2 Forcing Data
# 
# [https://disc.gsfc.nasa.gov/datasets/NLDAS_FORA0125_H_2.0/summary](https://disc.gsfc.nasa.gov/datasets/NLDAS_FORA0125_H_2.0/summary)
# 
# You need the DOI value: It's in the bottom material.

# In[ ]:


#################################################
#
# The Earthdata Catalog DOI Code
#

doi_nldas_2_forA = "10.5067/THUF4J1RLSYG"

# clip_polygon 41.66°, -106.393° : 48.099°, -95.14°
#

min_lat = 41.66 # degrees north
max_lat = 48.99 # degrees north

min_lon = -106.393 # degrees east
max_lon =   -95.14 # degrees east

#
#

lat_clip = (min_lat, max_lat)
lon_clip = (min_lon, max_lon)

#
# Apply xarray multifile dataset
#


partial_func = partial(_preprocess, 
                       lon_clip = lon_clip, 
                       lat_clip = lat_clip)

#
# Name for Subset for Archiving
#

subset_name = "SODAK"

#
# Name for Subset for Archiving
#

if (platform.system() == "Linux"):
    root_dir = "/data/DATASETS/NLDAS/hourly/"+subset_name + "/"
else:
    root_dir = "./data/"

output_directory = root_dir
#
#################################################


# ---
# ## Get Metadata Lookup Table
# 
# NASA doesn't use CF Metadata Standard Data

# In[ ]:


#################################################
#
# Variable Lookup Table to Convert NASA Names
#   to Standard Names
#

df_varlut = pd.read_csv(filepath_or_buffer = "./nldas_lut_lookup.csv",
                        index_col          = "INDEX")





encoding = {"time"      :{"units":"hours since 1970-01-01 00:00:00",
                          "dtype":np.float64},
            "time_bnds" :{"units":"hours since 1970-01-01 00:00:00",
                          "dtype":np.float64}}

variable_list = ["air_pressure",
                 "air_temperature",
                 "specific_humidity",
                 "eastward_wind",
                 "northward_wind",
                 "surface_downward_longwave_flux",
                 "surface_downward_shortwave_flux",
                 "atmosphere_convective_available_potential_energy",
                 "precipitation_amount",
                 "convective_precipitation_fraction",
                 "water_potential_evaporation_flux"]


for variable in variable_list:
    encoding[variable] = dict(zlib      =       True,
                              complevel =          7, 
                              dtype     = np.float32)


latitude_longitude = xr.DataArray(name  = "latitude_longitude",
                                  data  = np.int32(0),
                                  attrs = {"grid_mapping_name":"latitude_longitude",
                                           "longitude_of_prime_meridian":0.0 ,
                                           "semi_major_axis":6378137.0 ,
                                           "inverse_flattening":298.257223563})

print(latitude_longitude)

#
#################################################


# ---
# ## Step 5: Use your EarthData credentials to log in.
# 
# To open your session, you will need to "log on" to NASA using [earthaccess.login()](https://earthaccess.readthedocs.io/en/latest/user-reference/api/api/#earthaccess.api.login).
# 
# It will prompt you for a login, but make sure .netrc and .edl_token are in your "~" [home] directory so you don't have to!
# 

# In[ ]:


#################################################
#
# Authentication for EarthAccess
#

auth = earthaccess.login()

print("       Account: ",auth.username)
print("        System: ",auth.system)
print(" Authenticated: ",auth.authenticated)

print("")
print("---------------------------")

#
#################################################


# ---
# ## Step 6: Get Date Range
# 
# Don't get greedy; the actual pulling of the files can take a while.  You can always use loops.  The longer it tries, the harder it is to complete the request.
# 
# A monthly period extraction is used in this script.  For a daily example, the code block below will work. 
# 
# ```
# #################################################
# #
# # Date Range (Daily)
# #
# 
# start_date = np.datetime64('1980-01-01')
# end_date   = np.datetime64('1980-01-01')
# 
# 
# 
# time_series_to_end = np.arange(start = start_date, 
#                                stop  = end_date+1, 
#                                step  = np.timedelta64(1, 'D'))
# 
# print(time_series_to_end)
# 
# #
# #################################################
# ```

# In[ ]:


#################################################
#
# Date Range (Monthly)
#

start_date = np.datetime64('2009-01-01')
end_date   = np.datetime64('2025-01-01')



date_range = np.arange(start = start_date, 
                       stop  = end_date+1, 
                       step  = np.timedelta64(1, 'D'))

print(date_range)

#
#################################################


# ---
# ## The master retrieval Loop
# 
# 

# In[ ]:


for working_date in date_range:

    print("============================")


    #################################################
    #
    # Search for Available Data
    #

    daily_date       = pd.to_datetime(working_date)
    output_directory = root_dir + daily_date.strftime('%Y/%m/')
    filedate         = daily_date.strftime('%Y-%m-%d')


    output_filename = "NLDAS2_FORA_DAILY_"+subset_name+"_"+filedate+".nc"

    print(output_directory + output_filename)

    results_nldas  = earthaccess.search_data(doi      = doi_nldas_2_forA,
                                             temporal = (str(filedate) + ' 00:00:00',
                                                         str(filedate) + ' 23:59:00')) 
    # print(results_nldas)

    print("")
    print("---------------------------")

    #
    #################################################

    #################################################
    #
    # Retrieve the Requested EarthAccess Objects
    #

    print("Touching Requested NLDAS Objects")

    fs_nldas = earthaccess.open(results_nldas) # Extracts URLs from the results variable    

    #print(fs_nldas)

    print("")
    print("---------------------------")

    #
    #################################################

    #################################################
    #
    # Import the Dataframes from NASA
    #


    print("")
    print("---------------------------")
    print("Opening Pulling NLDAS FOR-A")
    ds_nldas  = xr.open_mfdataset(paths      = fs_nldas,
                                  preprocess = partial_func)





    ds_nldas["latitude_longitude"] = latitude_longitude

    print("")
    print("---------------------------")

    #
    #################################################

    #################################################
    #
    # Patch Attributes and Names
    #

    for index, row in df_varlut.iterrows():
        ds_nldas = ds_nldas.rename_vars(name_dict = {index:row["new_name"]})
        ds_nldas[row["new_name"]].attrs["long_name"]     = row["long_name"]
        ds_nldas[row["new_name"]].attrs["description"]   = row["long_name"]
        ds_nldas[row["new_name"]].attrs["standard_name"] = row["standard_name"]
        ds_nldas[row["new_name"]].attrs["units"]         = row["units_in"]
        if ("vmax" in ds_nldas[row["new_name"]].attrs):
            del ds_nldas[row["new_name"]].attrs["vmax"]
        if ("vmin" in ds_nldas[row["new_name"]].attrs):
            del ds_nldas[row["new_name"]].attrs["vmin"]
        ds_nldas[row["new_name"]].attrs["grid_mapping"] = "latitude_longitude"


    del ds_nldas["lat"].attrs["vmax"]
    del ds_nldas["lat"].attrs["vmin"]
    del ds_nldas["lon"].attrs["vmax"]
    del ds_nldas["lon"].attrs["vmin"]

    del ds_nldas["time"].attrs["time_increment"]
    del ds_nldas["time"].attrs["begin_date"]
    del ds_nldas["time"].attrs["end_date"]
    del ds_nldas["time"].attrs["begin_time"]
    del ds_nldas["time"].attrs["end_time"]
    ds_nldas["time"].attrs["standard_name"] = "time"
    del ds_nldas["time_bnds"].attrs["grid_mapping"]
    del ds_nldas["time_bnds"].attrs["units"]

    #
    #################################################

    #################################################
    #
    # Import the Dataframes from NASA
    #

    try:
        os.makedirs(output_directory)
    except FileExistsError:
        print(output_directory + " Exists")
    print("writing "+output_directory + output_filename)    
    output_fileloc = output_directory + output_filename

    ds_nldas.to_netcdf(path           = output_fileloc,
                       encoding       = encoding,
                       unlimited_dims = ['time'])

    #
    #################################################


# In[ ]:




