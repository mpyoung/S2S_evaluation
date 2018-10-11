# Import libraries
import netCDF4 as nc4
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap, shiftgrid

from netCDF4 import Dataset


'''
1. Grab GPCP rainfall over a certain region
'''

def grab_gpcp_region(inputfile,lonlim,latlim):
  #print 'loading tamsat data'
  nc_fid = Dataset(inputfile, 'r')
  lat = np.array(nc_fid.variables['latitude'][:])  # extract/copy the data
  lon = np.array(nc_fid.variables['longitude'][:])
  rfe = np.array(nc_fid.variables['precip'][:])

  rfe,lon = shiftgrid(180.,rfe,lon,start=False)
  # t_time = nc_fid.variables['time'][:]
  # # convert time to date using netCDF4 function
  # units='pentads since 2016-01-01'
  # #all_dates = nc4.num2date(time,units)
  # all_dates = t_time

  # grab rfe at specific latitude and longitude (inlat,inlon)
  rg_lon = np.where((lon >= lonlim[0]) & (lon <= lonlim[1]))
  rg_lat = np.where((lat >= latlim[0]) & (lat <= latlim[1]))

  match_lon = lon[rg_lon[0]]
  match_lat = lat[rg_lat[0]]

  rfe_rg = rfe[:,rg_lat[0],:][:,:,rg_lon[0]].squeeze()  # shape is time, lat, lon as shown above
  nc_fid.close()
  #rfe_rg[rfe_rg < 0] = 0

  return (rfe_rg,match_lon,match_lat)

'''
2. Grab GPCP at a specific point
'''
def grab_gpcp_point(inputfile,inlon,inlat):
  # print 'loading tamsat data'
  nc_fid = Dataset(inputfile, 'r')
  lat = np.array(nc_fid.variables['lat'][:])  # extract/copy the data
  lon = np.array(nc_fid.variables['lon'][:])
  t_time = nc_fid.variables['time'][:]
  rfe = np.array(nc_fid.variables['rfe'][:],dtype='f')  # shape is time, lat, lon as shown above
  nc_fid.close()

  # convert time to date using netCDF4 function
  units='days since 1901-01-01 00:00:00'
  all_dates = nc4.num2date(t_time,units)
  all_dates = t_time

  # grab rfe at specific latitude and longitude (inlat,inlon)
  diff_lon = abs(lon - inlon)
  diff_lat = abs(lat - inlat)

  match_lon = lon[diff_lon.argmin()]
  match_lat = lat[diff_lat.argmin()]

  rfe_point = rfe[:,diff_lat.argmin(),diff_lon.argmin()]

  return (rfe_point,match_lon,match_lat,all_dates)
