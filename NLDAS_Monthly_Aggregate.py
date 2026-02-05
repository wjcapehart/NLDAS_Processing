#!/usr/bin/env python
# coding: utf-8

# # Creating Monthly  NLDAS Data

# In[ ]:





# In[ ]:


#################################################
#
# Library 
#

import numpy             as np
import xarray            as xr
import matplotlib        as mpl
import matplotlib.pyplot as plt
import datetime          as datetime
import pandas            as pd
import os                as os
import h5netcdf          as h5netcdf
import platform          as platform

import cartopy.crs       as ccrs
import cartopy.feature   as cfeature

#
#################################################


# In[ ]:





# In[ ]:


#################################################
#
# Getting Spatial Metadata 
#

xf_metadata = xr.open_dataset(filename_or_obj = "/data/DATASETS/NLDAS/THREDDS/NLDAS_METADATA.nc",
                              engine          = "h5netcdf")

nx   = xf_metadata["lon"].values.size
ny   = xf_metadata["lat"].values.size
nt_h = 24

xf_metadata
#
#################################################


# In[ ]:





# In[ ]:


#################################################
#
# Getting Variable Lookup Tables 
#

df_varlut = pd.read_csv(filepath_or_buffer = "./metadata_lookup.csv", index_col="INDEX")
df_varlut

#
#################################################


# In[ ]:





# In[ ]:


#################################################
#
# Data Sets and EarthAccess 
#

input_daily_root_dir = "/data/DATASETS/NLDAS/THREDDS/DAILY/CONUS/"
out_monthly_dir      = "/data/DATASETS/NLDAS/THREDDS/MONTHLY"

#
#################################################


# In[ ]:


#################################################
#
# Date Range
#

start_date = np.datetime64('2026-01')
end_date   = np.datetime64('2026-01')
date_range = np.arange(start_date, end_date + np.timedelta64(1, 'M'))

for working_date in date_range:
    print(working_date)


#
#################################################


# In[ ]:


#################################################
#
# Compression Encoding
#

variable_list = ["mean_air_temperature", 
                 "specific_humidity", 
                 "air_pressure", 
                 "eastward_wind", 
                 "northward_wind", 
                 "atmosphere_convective_available_potential_energy", 
                 "surface_downward_shortwave_flux", 
                 "surface_downward_longwave_flux", 
                 "surface_net_downward_shortwave_flux", 
                 "surface_net_downward_longwave_flux", 
                 "surface_upward_latent_heat_flux", 
                 "surface_sensible_heat_flux", 
                 "downward_heat_flux_at_ground_level_in_soil", 
                 "surface_snow_and_ice_refreezing_flux", 
                 "snowfall_amount", 
                 "precipitation_amount", 
                 "water_evapotranspiration_amount", 
                 "surface_runoff_amount", 
                 "subsurface_runoff_amount", 
                 "surface_snow_melt_amount", 
                 "surface_temperature", 
                 "surface_albedo", 
                 "liquid_water_content_of_surface_snow", 
                 "surface_snow_thickness", 
                 "surface_snow_area_fraction", 
                 "mass_content_of_water_in_soil_layer_defined_by_root_depth", 
                 "surface_upward_potential_latent_heat_flux", 
                 "upward_latent_heat_flux_into_air_due_to_evaporation_of_intercepted_precipitation", 
                 "upward_latent_heat_flux_into_air_due_to_transpiration", 
                 "upward_latent_heat_flux_into_air_due_to_evaporation_from_soil", 
                 "surface_snow_sublimation_heat_flux", 
                 "canopy_water_amount", 
                 "leaf_area_index", 
                 "vegetation_area_fraction", 
                 "water_volume_transport_in_river_channel", 
                 "liquid_water_content_of_soil_layer",
                 "mass_content_of_water_in_soil_layer",
                 "soil_temperature",
                 "maximum_air_temperature", 
                 "minimum_air_temperature",
                 "maximum_maximum_air_temperature", 
                 "minimum_minimum_air_temperature",
                 "maximum_precipitation_amount",
                 "maximum_water_volume_transport_in_river_channel"]


encoding = {"time"       :{"units":"seconds since 1970-01-01 00:00:00",
                           "dtype":np.float64},
           "time_bnds"   :{"units":"seconds since 1970-01-01 00:00:00",
                           "dtype":np.float64}}

for variable in variable_list:
    encoding[variable] = dict(zlib      =       True,
                              complevel =          7, 
                              dtype     = np.float32)

#
#################################################


# In[ ]:


for date in date_range:

    #
    # File Names and Locations
    #

    date_slash          = datetime.datetime.strftime(pd.to_datetime(date), "%Y/%m")
    date_dash           = datetime.datetime.strftime(pd.to_datetime(date), "%Y-%m")

    daily_files         = input_daily_root_dir + date_slash + "/NLDAS_NOAH_DAILY_"   + date_dash + "-??.nc"
    monthly_output_file = out_monthly_dir                   + "/NLDAS_NOAH_MONTHLY_" + date_dash +    ".nc"

    start_month = np.datetime64(str(date)+                          "-01 00:00:00")
    end_month   = np.datetime64(str(date + np.timedelta64(1,'M')) + "-01 00:00:00")

    time_bnds   = xr.DataArray(data  = np.array([[start_month, end_month]]),
                               dims  = {"time":1,
                                        "bnds":2},
                               attrs = {"long_name": "time bounds",
                                      "description":"time bounds",
                                    "standard_name":"time"})

    print("Processing " + monthly_output_file)

    #################################################
    #
    # Pull Full Month
    #

    xf_noah_daily = xr.open_mfdataset(paths     = daily_files,
                                      data_vars = 'all')

    #
    #################################################

    #################################################
    #
    # Aggregate Daily to Monthly Data
    #

    #print("  Aggregating Daily Data")

    xf_noah_monthly = xf_noah_daily.coarsen(time   = len(xf_noah_daily.time)). \
                                        sum(skipna = False)

    del xf_noah_monthly["soil_layer_thickness"]
    del xf_noah_monthly["soil_layer_depth_bnds"]

    df_varlut_mean = df_varlut[ df_varlut["Aggregation"] == "MEAN" ]["new_name"].to_list()

    df_varlut_mean.append("maximum_air_temperature")
    df_varlut_mean.append("minimum_air_temperature")

    for variable in df_varlut_mean:
        # print("      - meaning for : " + variable)
        xf_noah_monthly[variable].values  = xf_noah_monthly[variable].values / 24.

    #
    #################################################


    #################################################
    #
    # Aggregating Max and Min Series
    #    

    #
    # Monthly Maximum Daily Air Temperature 
    #

    xf_noah_monthly["maximum_maximum_air_temperature"] = xf_noah_daily["maximum_air_temperature"]. \
                                                            coarsen(time   = len(xf_noah_daily.time)). \
                                                            max(skipna = False)
    xf_noah_monthly["maximum_maximum_air_temperature"].attrs["long_name"]     = "2-meter monthly maximum air temperature"
    xf_noah_monthly["maximum_maximum_air_temperature"].attrs["description"]   = "2-meter monthly maximum air temperature"
    xf_noah_monthly["maximum_maximum_air_temperature"].attrs["standard_name"] = "air_temperature"
    xf_noah_monthly["maximum_maximum_air_temperature"].attrs["units"]         = "K"

    #
    # Monthly Minimum Daily Air Temperature 
    #

    xf_noah_monthly["minimum_minimum_air_temperature"] = xf_noah_daily["minimum_air_temperature"]. \
                                                            coarsen(time   = len(xf_noah_daily.time)). \
                                                            min(skipna = False)
    xf_noah_monthly["minimum_minimum_air_temperature"].attrs["long_name"]     = "2-meter monthly minimum air temperature"
    xf_noah_monthly["minimum_minimum_air_temperature"].attrs["description"]   = "2-meter monthly minimum air temperature"
    xf_noah_monthly["minimum_minimum_air_temperature"].attrs["standard_name"] = "air_temperature"
    xf_noah_monthly["minimum_minimum_air_temperature"].attrs["units"]         = "K"

    #
    # Monthly Maximum Daily Precipitation
    #

    xf_noah_monthly["maximum_precipitation_amount"] = xf_noah_daily["precipitation_amount"].           \
                                                            coarsen(time   = len(xf_noah_daily.time)). \
                                                            max(skipna = False)
    xf_noah_monthly["maximum_precipitation_amount"].attrs["long_name"]     = "Monthly Maximum Precipitation Amount"
    xf_noah_monthly["maximum_precipitation_amount"].attrs["description"]   = "Monthly Maximum Precipitation Amount"
    xf_noah_monthly["maximum_precipitation_amount"].attrs["standard_name"] = "precipitation_amount"
    xf_noah_monthly["maximum_precipitation_amount"].attrs["units"]         = "kg m-2"

    #
    # Monthly Maximum Streamflow
    #

    xf_noah_monthly["maximum_water_volume_transport_in_river_channel"] = xf_noah_daily["water_volume_transport_in_river_channel"].           \
                                                            coarsen(time   = len(xf_noah_daily.time)). \
                                                            max(skipna = False)
    xf_noah_monthly["maximum_water_volume_transport_in_river_channel"].attrs["long_name"]     = "Monthly Maximum Streamflow"
    xf_noah_monthly["maximum_water_volume_transport_in_river_channel"].attrs["description"]   = "Monthly Maximum Streamflow"
    xf_noah_monthly["maximum_water_volume_transport_in_river_channel"].attrs["standard_name"] = "water_volume_transport_in_river_channel"
    xf_noah_monthly["maximum_water_volume_transport_in_river_channel"].attrs["units"]         = "m3"



    xf_noah_monthly["time_bnds"]              = time_bnds
    xf_noah_monthly["lon_bnds"]               = xf_metadata["lon_bnds"]
    xf_noah_monthly["lat_bnds"]               = xf_metadata["lat_bnds"]
    xf_noah_monthly["soil_layer_thickness"]   = xf_metadata["soil_layer_thickness"]
    xf_noah_monthly["soil_depth_bnds"]        = xf_metadata["soil_depth_bnds"]
    #
    #################################################

    del xf_noah_daily["time"].attrs["time_increment"]
    del xf_noah_daily["time"].attrs["begin_date"]
    del xf_noah_daily["time"].attrs["begin_time"]
    del xf_noah_daily["time"].attrs["end_date"]
    del xf_noah_daily["time"].attrs["end_time"]



    #################################################
    #
    # Drop NetCDF Files
    #

    xf_noah_monthly.to_netcdf(path         = monthly_output_file,
                            unlimited_dims = "time",
                            engine         = "h5netcdf",
                            encoding       = encoding)

    #
    #################################################


# In[ ]:





# In[ ]:





# In[ ]:




