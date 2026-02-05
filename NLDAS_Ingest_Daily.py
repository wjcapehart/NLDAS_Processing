#!/usr/bin/env python
# coding: utf-8

# # Processing DODS Stored NLDAS Data

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
import earthaccess       as earthaccess
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

xf_metadata = xr.open_dataset(filename_or_obj = "./NLDAS_METADATA.nc")
xf_metadata

nx   = xf_metadata["lon"].values.size
ny   = xf_metadata["lat"].values.size
nt_h = 24

#
#################################################


# In[ ]:


#################################################
#
# Getting Variable Lookup Tables 
#

df_varlut = pd.read_csv(filepath_or_buffer = "./metadata_lookup.csv", index_col="INDEX")

#
#################################################


# In[ ]:


#################################################
#
# Data Sets and EarthAccess 
#

if (platform.system() == "Linux"):
    root_dir = "/data/DATASETS/NLDAS/THREDDS/DAILY/CONUS/"
else:
    root_dir = "./data/"


doi_nldas_noah  = "10.5067/T4OW83T8EXDO"
doi_nldas_force = "10.5067/THUF4J1RLSYG"

auth            = earthaccess.login()

dt              = 3600.

#
#################################################


# In[ ]:


#################################################
#
# Date Range
#

start_date = np.datetime64('2026-02-01')
end_date   = np.datetime64('2026-02-04')
date_range = np.arange(start_date, end_date + np.timedelta64(1, 'D'))

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
                 "minimum_air_temperature"]


encoding = {"time" :{"units":"seconds since 1970-01-01 00:00:00",
                     "dtype":np.float64}}

for variable in variable_list:
    encoding[variable] = dict(zlib      =       True,
                              complevel =          7, 
                              dtype     = np.float32)

#
#################################################


# In[ ]:


for working_date in date_range:

    print("============================")
    
    
    #################################################
    #
    # Extract One Day of Data 
    #
    
    daily_date       = pd.to_datetime(working_date)
    output_directory = root_dir + daily_date.strftime('%Y/%m/')
    filedate         = daily_date.strftime('%Y-%m-%d')
    
    fileout          = "NLDAS_NOAH_DAILY_"+filedate+".nc"
    print("Processing "+filedate)
    # This will work if Earthdata prerequisite files have already been generated
    
    results_noah  = earthaccess.search_data(doi      = doi_nldas_noah,
                                            temporal = (filedate + ' 00:00:00', 
                                                        filedate + ' 23:59:00')) 
    
    results_force = earthaccess.search_data(doi      = doi_nldas_force,
                                            temporal = (filedate + ' 00:00:00', 
                                                        filedate + ' 23:59:00')) 
    
    print("Opening EA-Noah")
    fs_noah  = earthaccess.open(results_noah)  # Extracts URLs from the results variable
    
    print("")
    print("---------------------------")
    print("Opening EA-Forcings")
    fs_force = earthaccess.open(results_force) # Extracts URLs from the results variable
    
    print("")
    print("---------------------------")
    print("Opening XR-Noah")
    ds_noah  = xr.open_mfdataset(paths = fs_noah)
    
    print("")
    print("---------------------------")
    print("Opening XR-Forcings")
    ds_force = xr.open_mfdataset(paths = fs_force)
    
    #
    #################################################
    
    #################################################
    #
    # Clean and Merge the NOAH and Forcing datasets
    #
    
    ds_noah  =  ds_noah.drop_vars(["SMAvail_0_100cm",
                                   "SMAvail_0_200cm",
                                   "SoilM_0_100cm",
                                   "SoilM_0_200cm",
                                   "ACond",
                                   "CCond",
                                   "RCS",
                                   "RCT",
                                   "RCQ",
                                   "RCSOL",
                                   "RSmin",
                                   "RSMacr"])
    
    ds_force = ds_force.drop_vars(["SWdown",
                                   "LWdown",
                                   "PotEvap",
                                   "Rainf",
                                   "CRainf_frac"])
    
    xf_noah_hourly = xr.merge([ds_force, 
                               ds_noah])
    
    del ds_noah
    del ds_force
    
    
    #
    #################################################
    
    #################################################
    #
    # Reshape the Soils Parameters
    #
    
    
    xf_noah_hourly["soil_depth"]      = xf_metadata["soil_depth"]
    xf_noah_hourly["soil_thickness"]  = xf_metadata["soil_thickness"]
    xf_noah_hourly["soil_depth_bnds"] = xf_metadata["soil_depth_bnds"]
    
    soil_temp = xr.DataArray(data   = np.zeros(shape = [nt_h, 4, ny, nx], 
                                               dtype = np.float32),
                             dims   = ["time","soil_depth","lat","lon"],
                             coords = {"time":xf_noah_hourly["time"],
                                       "soil_depth":xf_metadata["soil_depth"],
                                       "lat":xf_metadata["lat"],
                                       "lon":xf_metadata["lon"]},
                             attrs  = {"description":"Soil Temperature",
                                       "long_name":"Soil Temperature",
                                       "standard_name":"soil_temperature",
                                       "units":"K"})
    
    soil_moist = xr.DataArray(data   = np.zeros(shape = [nt_h, 4, ny, nx], 
                                               dtype = np.float32),
                             dims   = ["time","soil_depth","lat","lon"],
                             coords = {"time":xf_noah_hourly["time"],
                                       "soil_depth":xf_metadata["soil_depth"],
                                       "lat":xf_metadata["lat"],
                                       "lon":xf_metadata["lon"]},
                             attrs  = {"description":"Soil Liquid Water Content",
                                       "long_name":"Soil Liquid Water Content",
                                       "standard_name":"mass_content_of_water_in_soil_layer",
                                       "units":"kg m-2"})
    
    soil_liq  = xr.DataArray(data   = np.zeros(shape = [nt_h, 4, ny, nx], 
                                               dtype = np.float32),
                             dims   = ["time","soil_depth","lat","lon"],
                             coords = {"time":xf_noah_hourly["time"],
                                       "soil_depth":xf_metadata["soil_depth"],
                                       "lat":xf_metadata["lat"],
                                       "lon":xf_metadata["lon"]},
                             attrs  = {"description":"Soil Liquid Water Content",
                                       "long_name":"Soil Liquid Water Content",
                                       "standard_name":"liquid_water_content_of_soil_layer",
                                       "units":"kg m-2"})
    
    xf_noah_hourly["liquid_water_content_of_soil_layer"]  = soil_liq
    xf_noah_hourly["mass_content_of_water_in_soil_layer"] = soil_moist
    xf_noah_hourly["soil_temperature"]                    = soil_temp
    
    soil_liq.values[:,0,:,:] = xf_noah_hourly["SMLiq_0_10cm"].values
    soil_liq.values[:,1,:,:] = xf_noah_hourly["SMLiq_10_40cm"].values
    soil_liq.values[:,2,:,:] = xf_noah_hourly["SMLiq_40_100cm"].values
    soil_liq.values[:,3,:,:] = xf_noah_hourly["SMLiq_100_200cm"].values 
    
    soil_moist.values[:,0,:,:] = xf_noah_hourly["SoilM_0_10cm"].values
    soil_moist.values[:,1,:,:] = xf_noah_hourly["SoilM_10_40cm"].values
    soil_moist.values[:,2,:,:] = xf_noah_hourly["SoilM_40_100cm"].values
    soil_moist.values[:,3,:,:] = xf_noah_hourly["SoilM_100_200cm"].values 
    
    soil_temp.values[:,0,:,:] = xf_noah_hourly["SoilT_0_10cm"].values
    soil_temp.values[:,1,:,:] = xf_noah_hourly["SoilT_10_40cm"].values
    soil_temp.values[:,2,:,:] = xf_noah_hourly["SoilT_40_100cm"].values
    soil_temp.values[:,3,:,:] = xf_noah_hourly["SoilT_100_200cm"].values 
    
    del xf_noah_hourly["SoilT_0_10cm"]
    del xf_noah_hourly["SoilT_10_40cm"]
    del xf_noah_hourly["SoilT_40_100cm"]
    del xf_noah_hourly["SoilT_100_200cm"]
    
    del xf_noah_hourly["SoilM_0_10cm"]
    del xf_noah_hourly["SoilM_10_40cm"]
    del xf_noah_hourly["SoilM_40_100cm"]
    del xf_noah_hourly["SoilM_100_200cm"]
    
    del xf_noah_hourly["SMLiq_0_10cm"]
    del xf_noah_hourly["SMLiq_10_40cm"]
    del xf_noah_hourly["SMLiq_40_100cm"]
    del xf_noah_hourly["SMLiq_100_200cm"]
    
    #
    #################################################
    
    #################################################
    #
    # Change Variable Metadata for Hourly Data
    #
    
    for index, row in df_varlut.iterrows():
        xf_noah_hourly  = xf_noah_hourly.rename_vars(name_dict = {index:row["new_name"]})
        xf_noah_hourly[row["new_name"]].attrs["long_name"]     = row["long_name"]
        xf_noah_hourly[row["new_name"]].attrs["description"]   = row["long_name"]
        xf_noah_hourly[row["new_name"]].attrs["standard_name"] = row["standard_name"]
        xf_noah_hourly[row["new_name"]].attrs["units"]         = row["units_in"]
    
        if ("vmax" in xf_noah_hourly[row["new_name"]].attrs):
            del xf_noah_hourly[row["new_name"]].attrs["vmax"]
        if ("vmin" in xf_noah_hourly[row["new_name"]].attrs):
            del xf_noah_hourly[row["new_name"]].attrs["vmin"]
    print("--------")
    
    #
    #################################################
    
    #################################################
    #
    # Aggregate Hourly to Daily Data
    #
    
    print("Aggregating Daily Data")
    
    xf_noah_daily = xf_noah_hourly.drop_vars("time_bounds").coarsen(time=24).sum(skipna=False)
    for index, row in df_varlut.drop(axis="columns", index=["time","time_bnds"]).iterrows():
        xf_noah_daily[row["new_name"]].attrs["long_name"]     = row["long_name"]
        xf_noah_daily[row["new_name"]].attrs["description"]   = row["long_name"]
        xf_noah_daily[row["new_name"]].attrs["standard_name"] = row["standard_name"]
        xf_noah_daily[row["new_name"]].attrs["units"]         = row["units_in"]
    
    df_varlut_mean = df_varlut[ df_varlut["Aggregation"] ==     "MEAN" ]["new_name"].to_list()
    df_varlut_sum  = df_varlut[ df_varlut["Aggregation"] ==      "SUM" ]["new_name"].to_list()
    df_varlut_int  = df_varlut[ df_varlut["Aggregation"] == "INTEGRAL" ]["new_name"].to_list()
    
    for variable in df_varlut_int:
        xf_noah_daily[variable].values  = xf_noah_daily[variable].values * dt
    
    for variable in df_varlut_mean:
        xf_noah_daily[variable].values  = xf_noah_daily[variable].values / 24.
    
    xf_noah_daily["maximum_air_temperature"] = xf_noah_hourly["mean_air_temperature"].coarsen(time=24).max(skipna=False)
    xf_noah_daily["maximum_air_temperature"].attrs["long_name"]     = "2-meter maximum air temperature"
    xf_noah_daily["maximum_air_temperature"].attrs["description"]   = "2-meter maximum air temperature"
    xf_noah_daily["maximum_air_temperature"].attrs["standard_name"] = "air_temperature"
    xf_noah_daily["maximum_air_temperature"].attrs["units"]         = "K"
    
    xf_noah_daily["minimum_air_temperature"] = xf_noah_hourly["mean_air_temperature"].coarsen(time=24).min(skipna=False)
    xf_noah_daily["minimum_air_temperature"].attrs["long_name"]     = "2-meter minimum air temperature"
    xf_noah_daily["minimum_air_temperature"].attrs["description"]   = "2-meter minimum air temperature"
    xf_noah_daily["minimum_air_temperature"].attrs["standard_name"] = "air_temperature"
    xf_noah_daily["minimum_air_temperature"].attrs["units"]         = "K"
    
    
    #
    #########################################################
    
    
    #########################################################
    #
    # Write to File
    #
    
    
    del xf_noah_daily["time"].attrs["units"]
    
    try:
        os.makedirs(output_directory)
    except FileExistsError:
        print(output_directory + " Exists")
    print("writing "+output_directory + fileout)
    xf_noah_daily.to_netcdf(path           = output_directory + fileout,
                            unlimited_dims = "time",
                            engine         = "h5netcdf",
                            encoding       = encoding)
    
    #
    #########################################################

    


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




