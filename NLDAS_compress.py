#!/usr/bin/env python
# coding: utf-8

# # Processing DODS Stored NLDAS Data

# In[ ]:





# In[4]:


#################################################
#
# Library 
#

import numpy             as np
import xarray            as xr
import datetime          as datetime
import pandas            as pd
import h5netcdf          as h5netcdf
import platform          as platform


#
#################################################


# In[ ]:


#################################################
#
# Date Range
#

root_dir = "/data/DATASETS/NLDAS/netcdf/"

start_date = np.datetime64('2024-01-02')
end_date   = np.datetime64('2024-01-04')
date_range = np.arange(start_date, end_date + np.timedelta64(1, 'D'))

for working_date in date_range:
    print(working_date)


#
#################################################


# In[6]:


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

    #
    #################################################
        
    #########################################################
    #
    # Open File
    #

     

    print("  reading "+output_directory + fileout)
    xf_noah_daily = xr.open_dataset(output_directory + fileout)
    
    #
    #########################################################


    
    #########################################################
    #
    # Write to File
    #
    
    
    del xf_noah_daily["time"].attrs["units"]
    

    print("  writing "+output_directory + fileout)
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




