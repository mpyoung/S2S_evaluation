# Import libraries
import netCDF4 as nc4
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap, shiftgrid

from netCDF4 import Dataset


'''
1. Grab NCEP total rainfall over a certain region
'''

def grab_ncep_region(date,inputfile,lonlim,latlim):
  #print 'loading tamsat data'
  nc_fid = Dataset(inputfile, 'r')
  lat = np.array(nc_fid.variables['latitude'][:])  # extract/copy the data
  lon = np.array(nc_fid.variables['longitude'][:])
  rfe = np.array(nc_fid.variables['tp'][:])
  t_time = np.array(nc_fid.variables['t'][:])
  #t_units = str(nc_fid.variables['t'].units)

  rfe,lon = shiftgrid(180.,rfe,lon,start=False)
  # t_time = nc_fid.variables['time'][:]
  # # convert time to date using netCDF4 function
  units='hours since '+date+' 00:00:00'
  all_dates = nc4.num2date(t_time,units) # need to check this works
  # all_dates = t_time

  # split dates into year, month and day columns
  final_times = np.zeros((len(all_dates),4))
  for t in np.arange(0,len(all_dates)):
    curr_time = all_dates[t]
    final_times[t,0] = curr_time.year
    final_times[t,1] = curr_time.month
    final_times[t,2] = curr_time.day
    final_times[t,3] = curr_time.hour


  # grab rfe at specific latitude and longitude (inlat,inlon)
  rg_lon = np.where((lon >= lonlim[0]) & (lon <= lonlim[1]))
  rg_lat = np.where((lat >= latlim[0]) & (lat <= latlim[1]))

  match_lon = lon[rg_lon[0]]
  match_lat = lat[rg_lat[0]]

  rfe_rg = rfe[:,:,rg_lat[0],:][:,:,:,rg_lon[0]].squeeze()  # shape is time, lat, lon as shown above
  nc_fid.close()
  #rfe_rg[rfe_rg < 0] = 0

  return (rfe_rg,match_lon,match_lat,final_times)
