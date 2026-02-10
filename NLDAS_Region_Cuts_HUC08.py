#!/usr/bin/env python
# coding: utf-8

# In[ ]:


##############################################
#
# Library
#

import numpy                as np
import matplotlib.pyplot    as plt
import mpl_toolkits         as mpl

import cartopy.crs          as ccrs  
import cartopy.feature      as cfeature
import cartopy.io.img_tiles as cimgt
import cartopy.mpl          as cmpl

import geopandas            as geopandas
import pandas               as pd

import xarray               as xr

import rasterio             as rio

import pathlib              as pathlib

import metpy                as metpy

try:
    __IPYTHON__
    print("Running in IPython/Jupyter environment")
    iPython_On = True
except NameError:
    print("Not running in IPython environment (running in standard Python interpreter)")
    iPython_On = False

#
##############################################


# In[ ]:


##############################################
#
# Import Data
#


web_nldas_url  = "https://thredds.ias.sdsmt.edu:8443/thredds/dodsC/NLDAS/"
fs_nldas_url   = "/data/DATASETS/NLDAS/THREDDS/"


if pathlib.Path("/ias_raid/").is_dir():
    root_nldas_url = fs_nldas_url
else:
    root_nldas_url = web_nldas_url


NLDAS_CONUS_File = root_nldas_url + "MONTHLY/NLDAS_NOAH_MONTHLY_1979-02_to_2026-01.nc"
MASK_CONUS_File  = root_nldas_url + "NLDAS_METADATA.nc"

ds_CONUS = xr.open_dataset(filename_or_obj = NLDAS_CONUS_File)
ds_MASK  = xr.open_dataset(filename_or_obj =  MASK_CONUS_File)



if (root_nldas_url == web_nldas_url):

    #
    # String Conversions.
    #

    ds_MASK["climdiv_name"]       = ds_MASK["climdiv_name"].str.decode('utf-8')
    ds_MASK["climdiv_state_abrv"] = ds_MASK["climdiv_state_abrv"].str.decode('utf-8')
    ds_MASK["climdiv_state_name"] = ds_MASK["climdiv_state_name"].str.decode('utf-8')
    ds_MASK["huc08_name"]         = ds_MASK["huc08_name"].str.decode('utf-8')
    ds_MASK["huc06_name"]         = ds_MASK["huc06_name"].str.decode('utf-8')
    ds_MASK["huc04_name"]         = ds_MASK["huc04_name"].str.decode('utf-8')
    ds_MASK["huc02_name"]         = ds_MASK["huc02_name"].str.decode('utf-8')
    ds_MASK["huc08_state"]        = ds_MASK["huc08_state"].str.decode('utf-8')

#
##############################################


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



encoding = {"time"       :{"units":"seconds since 1970-01-01 00:00:00",
                           "dtype":np.float64} }

for variable in variable_list:
    encoding[variable] = dict(zlib      =       True,
                              complevel =          7, 
                              dtype     = np.float32)

#
#################################################


# In[ ]:


##############################################
#
# Mapping Coordinates
#


ccrs_laee = ccrs.LambertAzimuthalEqualArea(central_latitude  = ds_MASK["lat"].mean().values,
                                           central_longitude = ds_MASK["lon"].mean().values)


#
##############################################


# In[ ]:


##############################################
#
# Pull Climdiv tabular data 
# 



variables_with_dim = [var_name for var_name, var in ds_MASK.data_vars.items() if "climdiv" in var.dims]
ds_climdiv_table   = ds_MASK[variables_with_dim]
df_climdiv_table   = ds_MASK[variables_with_dim].to_pandas()

variables_with_dim = [var_name for var_name, var in ds_MASK.data_vars.items() if "huc08" in var.dims]
ds_huc08_table   = ds_MASK[variables_with_dim]
df_huc08_table   = ds_MASK[variables_with_dim].to_pandas()
df_huc08_table
#
##############################################


# In[ ]:


##############################################
#
# Pull Climdiv tabular data 
# 



localle = "nCLIMDIV-"
localle = "HUC08-"

climdiv_state_abrv        = "SD"

df_climdiv_table_regional = df_climdiv_table[df_climdiv_table["climdiv_state_abrv"]  == climdiv_state_abrv]

df_climdiv_table

#
##############################################


# In[ ]:


##############################################
#
# Create Master Mask
# 

nldas_mask    = ds_MASK["CONUS_mask"].where(ds_MASK["CONUS_mask"] == 1)

regional_mask = ds_MASK["CLIMDIV"] * nldas_mask
regional_mask = ds_MASK["HUC_08"] * nldas_mask



#
##############################################


# In[ ]:


#regional_mask.plot()


# In[ ]:


##############################################
#
# Loop Through Regions
# 




for my_region in df_huc08_table.index:


    if (localle == "nCLIMDIV-"):
        localle_places = 4
        file_prefix    = root_nldas_url +  "MONTHLY/CLIMDIV/" + localle + str(my_region).zfill(localle_places)

    else:
        localle_places = 8
        file_prefix    = root_nldas_url +  "MONTHLY/HUC08/" + localle + str(my_region).zfill(localle_places)

    print("-> Processing Region : ", str(my_region).zfill(localle_places))

    #
    # Create Mask
    #

    print("   -> Selecting Mask or Region", str(my_region).zfill(localle_places))

    local_mask = regional_mask.where(regional_mask == my_region)
    local_mask = local_mask/local_mask

    local_clip = local_mask.dropna(dim="lat", how="all").\
                            dropna(dim="lon", how="all")

    #
    # Resolvable Boundaries?
    #

    if (0 in (local_clip.shape)):
        print("   -> Empty Box ", str(my_region).zfill(localle_places))
    else:
        lon_max    = local_clip.coords["lon"].max().values
        lat_max    = local_clip.coords["lat"].max().values

        lon_min    = local_clip.coords["lon"].min().values
        lat_min    = local_clip.coords["lat"].min().values

        print("   -> Apply Box Clip: ",lon_min, lon_max,lat_min, lat_max)    

        ds_CONUS_clipped =   ds_CONUS.sel(lon=slice(lon_min, lon_max), 
                                          lat=slice(lat_min, lat_max))

        local_clip       = local_mask.sel(lon=slice(lon_min, lon_max), 
                                          lat=slice(lat_min, lat_max))

        print("   -> Apply Mask")

        #
        # Mask
        #

        df_nldas_regional       = ds_CONUS_clipped.where(local_clip>0)

        print("   -> Aggregate Mask")

        #
        # Aggregate
        #

        df_nldas_regional_mean  = df_nldas_regional.mean(dim        = ['lat', 'lon'],
                                                         keep_attrs = True)





        #################################################
        #
        # Drop NetCDF Files
        #

        aggregate_output_file = file_prefix+"_NLDAS_NOAH_MONTHLY_1979-02_to_2026-01.nc"

        print("   -> Export NETCDF to ", aggregate_output_file)


        df_nldas_regional_mean.to_netcdf(path           = aggregate_output_file,
                                         unlimited_dims = "time",
                                         engine         = "h5netcdf",
                                         encoding       = encoding)

        #
        #################################################

        print("-----------------")

        #
        ##############################################


# In[ ]:





# In[ ]:




