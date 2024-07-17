#!/usr/bin/env python3

import xarray as xr
import numpy as np

def convert_coords_dtype(ds, dtype):
    new_coords = {coord: ds.coords[coord].astype(dtype) for coord in ["Z", "XC", "YC"]}
    return ds.assign_coords(**new_coords)

root_dir = '/Users/zhihua/Documents/Work/Research/Projects/SBtracer/data_original/'

cnames = ['ctrl', 'GM2X', 'Redi2X', 'Kappa2X']
for cname in cnames:
    ds = xr.open_dataset(root_dir + 'single_basin_'+cname+'_compressed.nc')
    ds.close()
    
    dsm = ds[['Ttave', 'Stave', 'maskC', 'PTRtave01']].mean('time', keep_attrs=True).isel(Z=slice(1,None))
    dsm = dsm.rename_vars({'Ttave':'THETA', 'Stave':'SALT', 'PTRtave01':'AGE'})
    dsm = convert_coords_dtype(dsm, np.float64)
    # dsm['SALT']  = dsm.SALT.where(dsm.Z!=dsm.Z[0], dsm.SALT.isel(Z=1))
    # dsm['THETA'] = dsm.THETA.where(dsm.Z!=dsm.Z[0], dsm.THETA.isel(Z=1))
    # dsm['AGE']   = dsm.AGE.where(dsm.Z!=dsm.Z[0], dsm.AGE.isel(Z=1))
    dsm.to_netcdf('/Users/zhihua/Documents/Work/Research/Projects/SBtracer/data/SB_'+cname+'.nc')