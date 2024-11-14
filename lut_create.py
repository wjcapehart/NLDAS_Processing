#################################################
#
# Create Soil Depths 
#

soil_depth = xr.DataArray(np.array([0.05, 0.25, 0.70, 1.50], 
                                   dtype = np.float32),
                          dims  = ["soil_depth"],
                          attrs = {"long_names":"Soil Depth",
                                   "standard_name":"depth",
                                    "bounds":"soil_depth_bnds",
                                   "units":"m"})

soil_thickness = xr.DataArray(data  = np.array([0.10, 0.40, 0.50, 1.00], 
                                               dtype = np.float32),
                              dims  = ["soil_depth"],
                              attrs = {"long_names":"Soil Depth Bounds",
                                        "standard_name":"depth",
                                        "bounds":"soil_depth_bnds",
                                       "units":"m"})

soil_depth_bnds      = np.zeros(shape = [4, 2],  
                                dtype = np.float32)
soil_depth_bnds[:,0] = soil_depth.values - soil_thickness.values/2
soil_depth_bnds[:,1] = soil_depth.values + soil_thickness.values/2

soil_depth_bnds = xr.DataArray(data  = soil_depth_bnds,
                               dims  = ["soil_depth",
                                        "bnds"],
                               attrs = {"long_name":"soil depth bounds",
                                            "units":"m"})
df_metadata["soil_depth"] = soil_depth
df_metadata["soil_depth_bnds"] = soil_depth_bnds
df_metadata["soil_thickness"] = soil_thickness

df_metadata.to_netcdf("./NLDAS_METADATA.nc4")

#################################################
#

variables =  list(xf_noah_hourly.coords.keys())  + list(xf_noah_hourly.data_vars.keys()) 
variables_df = pd.DataFrame(index=variables,
                            columns=["long_name","standard_name","units"])

for var in variables:
    print("Processing :", var)
    try:
        long_name = xf_noah_hourly[var].attrs["long_name"]
        print("   -> long_name : ", long_name)
        variables_df.loc[var, "long_name"] = long_name
    except KeyError:
        print("   -> long_name : NA")
    try:
        standard_name = xf_noah_hourly[var].attrs["standard_name"]
        print("   -> standard_name : ", standard_name)
        variables_df.loc[var, "standard_name"] = standard_name
    except KeyError:
        print("   -> standard_name : NA")
    try:
        units = xf_noah_hourly[var].attrs["units"]
        print("   -> units : ", units)
        variables_df.loc[var, "units"] = units
    except KeyError:
        print("   -> units : NA")

#
#################################################

