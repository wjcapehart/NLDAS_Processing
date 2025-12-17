#!/usr/bin/env python
# coding: utf-8

# # Processing DODS Stored MERRA2 Data

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

xf_metadata = xr.open_dataset(filename_or_obj = "./MERRA2_101.const_2d_asm_Nx.00000000.nc4.ncml.nc")

nx   = xf_metadata["lon"].values.size
ny   = xf_metadata["lat"].values.size
nt_h = 24
nz   = 42

g =  9.807 # m/s²

sfc_z = xf_metadata["PHIS"] 

sfc_z = sfc_z/g

print(sfc_z)

sfc_z.plot()



#
#################################################


# In[ ]:


#################################################
#
# Getting Variable Lookup Tables 
#

#
#################################################


# In[ ]:


#################################################
#
# Data Sets and EarthAccess 
#

if (platform.system() == "Linux"):
    root_dir = "/data/DATASETS/MERRA2/netcdf/"
else:
    root_dir = "./data/"

output_directory = root_dir


doi_merra2_3danl   = "10.5067/V92O8XZ30XBI"
doi_merra2_3dassim = "10.5067/2E096JV59PK7"
doi_merra2_static  = "10.5067/ME5QX6Q5IGGU"

auth            = earthaccess.login()

auth.token

#
#################################################


# In[ ]:


#################################################
#
# Date Range
#

start_year = 1980
end_year   = 2025


start_date = np.datetime64('1980-01')
end_date   = np.datetime64('2024-12')
date_range = np.arange(start_date, end_date + np.timedelta64(1, 'M'))
year_range = np.arange(start_year, end_year + .1, dtype=int)

print(date_range)

print(year_range)

date_range.size

#
#################################################


# 

# In[ ]:


#################################################
#
# Compression Encoding
#

variable_list = ["H", 
                 "OMEGA", 
                 "PS", 
                 "QI", 
                 "QL", 
                 "QV", 
                 "SLP", 
                 "T", 
                 "U", 
                 "V", 
                 "Var_H", 
                 "Var_OMEGA", 
                 "Var_PS", 
                 "Var_QI", 
                 "Var_QL", 
                 "Var_QV", 
                 "Var_RH", 
                 "Var_SLP", 
                 "Var_T", 
                 "Var_U", 
                 "Var_V"]


encoding = {"time" :{"units":"seconds since 1970-01-01 00:00:00",
                     "dtype":np.float64}}

#for variable in variable_list:
#    encoding[variable] = dict(zlib      =       True,
#                              complevel =          7, 
#                              dtype     = np.float32)

variable_list = ["time","lev","lat","lon",
                 "H", 
                 "OMEGA", 
                 "PS", 
                 "QI", 
                 "QL", 
                 "QV", 
                 "SLP", 
                 "T", 
                 "U", 
                 "V", 
                 "Var_H", 
                 "Var_OMEGA", 
                 "Var_PS", 
                 "Var_QI", 
                 "Var_QL", 
                 "Var_QV", 
                 "Var_RH", 
                 "Var_SLP", 
                 "Var_T", 
                 "Var_U", 
                 "Var_V"]

print(encoding)
#
#################################################


# In[ ]:


for year in year_range:

    str_year = str(year)
    fileout          = "MERRA_3Da_MONTHLY_"+str_year+".nc"


    results_merra = earthaccess.search_data(doi      = doi_merra2_3dassim,
                                            temporal = (str_year + '-01-01 00:00:00', 
                                                        str_year + '-12-31 23:59:59')) 
    print("")
    print("---------------------------")
    print("Opening EA-MERRA")
    fs_merra = earthaccess.open(results_merra) # Extracts URLs from the results variable

    print("")
    print("---------------------------")
    print("Opening XR-MERRA")
    ds_merra  = xr.open_mfdataset(paths = fs_merra)[variable_list]




    #########################################################
    #
    # Write to File
    #

    print("writing file")
    #del ds_merra["time"].attrs["units"]
    ds_merra.to_netcdf(path           = output_directory + fileout,
                       unlimited_dims = "time",
                       engine         = "h5netcdf",
                       encoding       = encoding)

    print(fileout, "written")
    print(fileout, "============================")

    #
    #########################################################


# In[ ]:





# In[ ]:




