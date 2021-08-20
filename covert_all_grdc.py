# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 18:52:20 2021

#!/usr/bin/python3.8



@author:  Matthias Büchner, buechner@pik-potsdam.de, Emmanuel Nyenah 
"""


import netCDF4 as nc
import numpy as np
import xarray as xr
ref_year = '1700'
ds=xr.open_dataset('GRDC-Monthly.nc' ,decode_times=False)
file = 'GRDC-Monthly.nc'

dataset = nc.Dataset(file)

n_stations = dataset.dimensions['id'].size
n_timesteps = dataset.dimensions['time'].size


ds_lat = dataset.variables['geo_y'][:]
ds_lon = dataset.variables['geo_x'][:]


ds_discharge = dataset.variables['runoff_mean'][:]
ds_time = dataset.variables['time'][:]


ncout = nc.Dataset( 'affresh.nc', 'w',format='NETCDF4_CLASSIC')
ncout.createDimension('time', n_timesteps)
ncout.createDimension('lon', 720)
ncout.createDimension('lat', 360)

time = ncout.createVariable('time',np.dtype('float32').char,('time'))
lat = ncout.createVariable('lat',np.dtype('float32').char,('lat'))
lon = ncout.createVariable('lon',np.dtype('float32').char,('lon'))
var = ncout.createVariable('dis',np.dtype('float32').char,('time','lat','lon'), zlib=True, complevel=5, fill_value=1e+20)

time.standard_name = 'time'
time.long_name = 'Time'
time.axis = 'T'
time.calendar = 'standard'
time.units = 'days since ' + ref_year + '-01-01 00:00:00'
time[:] = ds_time

#=========================================================
res    = 0.5
latbnd = np.asarray([np.arange(- 90    , 90     ,res),
                     np.arange(- 90+res, 90+0.01,res)]).T
lonbnd = np.asarray([np.arange(-180    ,180     ,res),
                     np.arange(-180+res,180+0.01,res)]).T

lat_r   = latbnd.mean(axis=1)
lon_r  = lonbnd.mean(axis=1)
#=========================================================

lat.long_name = 'Latitude'
lat.standard_name = 'latitude'
lat.units = 'degrees_north'
lat.axis = 'Y'
lat[:] = lat_r

lon.long_name = 'Longitude'
lon.standard_name = 'longitude'
lon.units = 'degrees_east'
lon.axis = 'X'
lon[:] = lon_r

#=============================================

lat_fix=ds['geo_y'][:]
lat_fix=lat_fix.values

lon_fix=ds['geo_x'][:]
lon_fix=lon_fix.values
runoff=ds['runoff_mean']
runoff=runoff.values
c=np.ones((len(runoff),360,720))
c[:]=np.nan
for j in range(0,len(runoff)):
    for i in range( 0,len(lat_fix)):
        lat_ts=lat_fix[i]
        lon_ts=lon_fix[i]
      
        sq_dis_lat=(lat-lat_ts)**2
        sq_dis_lon=(lon-lon_ts)**2
      
        #getting index
        min_lat_index=sq_dis_lat.argmin()
        min_lon_index=sq_dis_lon.argmin()
        print(j,min_lat_index,min_lon_index)
        print(runoff[j,i])
        c[j,min_lat_index,min_lon_index]=runoff[j,i]

#==============================================
var.comment = 'GRDC calculated from daily data'
var.units = 'm3 s-1'
var[:] = c


ncout.title = 'Mean daily discharge (Q)'
ncout.references = 'grdc.bafg.de'
ncout.institution = 'GRDC'
ncout.history = 'Download from GRDC Database, 29/07/2021'
