'''
Functions to open and spatially subset
data files from:

GPCP 1 degree precipitation estimates
S2S forecasts from UKMO, NCEP & ECMWF

All this data can be found in:
/group_workspaces/jasmin2/klingaman/datasets/

M.Young 29/06/2018
'''

# Import libraries
import netCDF4 as nc4
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset


'''
Grab GPCP rainfall over a certain region
'''
def grab_gpcp_region(inputfile,lonlim,latlim):
  #print 'loading tamsat data'
  nc_fid	= Dataset(inputfile, 'r')
  lat		= np.array(nc_fid.variables['latitude'][:])		# extract/copy the data
  lon		= np.array(nc_fid.variables['longitude'][:])
  rfe		= np.array(nc_fid.variables['precip'][:])

  # shiftgrid: shift global lat/lon grid east or west
  # start, if True, 180. represents the starting of the new grid; if False, 180. represent the ending longitude.
  rfe,lon	= shiftgrid(180.,rfe,lon,start=False)
  # t_time = nc_fid.variables['time'][:]
  # # convert time to date using netCDF4 function
  # units='pentads since 2016-01-01'
  # #all_dates = nc4.num2date(time,units)
  # all_dates = t_time

  # grab rfe at specific latitude and longitude (inlat,inlon)
  rg_lon	= np.where( (lon >= lonlim[0]) & (lon <= lonlim[1]) )
  rg_lat	= np.where( (lat >= latlim[0]) & (lat <= latlim[1]) )

  match_lon	= lon[rg_lon[0]]
  match_lat	= lat[rg_lat[0]]

  rfe_rg	= rfe[:,rg_lat[0],:][:,:,rg_lon[0]].squeeze()		# shape is time, lat, lon as shown above
  nc_fid.close()
  #rfe_rg[rfe_rg < 0] = 0

  return(rfe_rg,match_lon,match_lat)


'''
2. Grab UKMO GloSea5-GC2 rainfall over a certain region
'''
def grab_ukmo_region(inputfile,lonlim,latlim):
  nc_fid	= Dataset(inputfile,'r')
  lat		= np.array(nc_fid.variables['latitude'][:] )
  lon		= np.array(nc_fid.variables['longitude'][:])
  rfe		= np.array(nc_fid.variables['tp'][:]       )

  t_time	= np.array(nc_fid.variables['t'][:])
  t_units	= str(nc_fid.variables['t'].units)

  rfe,lon	= shiftgrid(180.,rfe,lon,start=False)
  # return datetime objects given numeric time values.
  all_dates	= nc4.num2date(t_time,t_units)

  # split dates into year, month and day columns
  final_times	= np.zeros((len(all_dates),3))
  for t in np.arange(0,len(all_dates)):
    curr_time		= all_dates[t]
    final_times[t,0]	= curr_time.year
    final_times[t,1]	= curr_time.month
    final_times[t,2]	= curr_time.day

  # grab rfe at specific latitude and longitude (inlat,inlon)
  # this is interesting: if only condition is given for np.where, return condition.nonzero(), which is a tuple!
  rg_lon	= np.where( (lon>=lonlim[0]) & (lon<=lonlim[1]) )
  rg_lat	= np.where( (lat>=latlim[0]) & (lat<=latlim[1]) )

  # [0] here does nothing
  match_lon	= lon[rg_lon[0]]
  match_lat	= lat[rg_lat[0]]

  # indexing arrays cannot be broadcast togather
  # squeeze() here does nothing
  rfe_rg	= rfe[:,:,rg_lat[0],:][:,:,:,rg_lon[0]].squeeze()
  nc_fid.close()
  #rfe_rg[rfe_rg < 0] = 0

  return (rfe_rg,match_lon,match_lat,final_times)


'''
3. Grab NCEP total rainfall over a certain region
'''

def grab_ncep_region(date,inputfile,lonlim,latlim):
  nc_fid	= Dataset(inputfile,'r')
  lat		= np.array(nc_fid.variables['latitude'][:] )
  lon		= np.array(nc_fid.variables['longitude'][:])
  rfe		= np.array(nc_fid.variables['tp'][:]       )

  t_time	= np.array(nc_fid.variables['t'][:])
# t_units	= str(nc_fid.variables['t'].units)			# this units is wrong!
  units		='hours since '+date+' 00:00:00'

  rfe,lon	= shiftgrid(180.,rfe,lon,start=False)
  # t_time	= nc_fid.variables['time'][:]
# # convert time to date using netCDF4 function
  all_dates	= nc4.num2date(t_time,units)				# need to check this works
  # all_dates = t_time

  # split dates into year, month and day columns
  final_times	= np.zeros((len(all_dates),4))
  for t in np.arange(0,len(all_dates)):
    curr_time		= all_dates[t]
    final_times[t,0]	= curr_time.year
    final_times[t,1]	= curr_time.month
    final_times[t,2]	= curr_time.day
    final_times[t,3]	= curr_time.hour				# why hour is needed?

  # grab rfe at specific latitude and longitude (inlat,inlon)
  rg_lon	= np.where( (lon>=lonlim[0]) & (lon<=lonlim[1]) )
  rg_lat	= np.where( (lat>=latlim[0]) & (lat<=latlim[1]) )

  match_lon	= lon[rg_lon[0]]
  match_lat	= lat[rg_lat[0]]

  rfe_rg	= rfe[:,:,rg_lat[0],:][:,:,:,rg_lon[0]].squeeze()	# shape is time, lat, lon as shown above
  nc_fid.close()
  #rfe_rg[rfe_rg < 0] = 0

  return (rfe_rg,match_lon,match_lat,final_times)


'''
4. Grab ECMWF rainfall over a certain region

# e.g. filepath:
/group_workspaces/jasmin2/klingaman/datasets/S2S/ECMWF/
ENSO/hindcasts_1120_2017/ecmwf_precip_ens_19981120.nc
'''

def grab_ecmf_region(date,inputfile,lonlim,latlim):
  #print 'loading tamsat data'
  nc_fid= Dataset(inputfile, 'r')
  lat	= np.array(nc_fid.variables['latitude'][:])
  lon	= np.array(nc_fid.variables['longitude'][:])
  #t_units = str(nc_fid.variables['t'].units) # time units are wrong
  t_time= np.array(nc_fid.variables['t'][:])
   # convert time to date using netCDF4 function
  units='hours since '+date+' 00:00:00'
  all_dates = nc4.num2date(t_time,units)
  # all_dates = t_time

  # split dates into year, month and day columns
  final_times = np.zeros((len(all_dates),4))
  for t in np.arange(0,len(all_dates)):
    curr_time = all_dates[t]
    final_times[t,0] = curr_time.year
    final_times[t,1] = curr_time.month
    final_times[t,2] = curr_time.day
    final_times[t,3] = curr_time.hour

  rfe = np.array(nc_fid.variables['tp'][:])
  nc_fid.close()

  rfe,lon = shiftgrid(180.,rfe,lon,start=False)

  # grab rfe at specific latitude and longitude (inlat,inlon)
  rg_lon = np.where((lon >= lonlim[0]) & (lon <= lonlim[1]))
  rg_lat = np.where((lat >= latlim[0]) & (lat <= latlim[1]))

  match_lon = lon[rg_lon[0]]
  match_lat = lat[rg_lat[0]]

  rfe_rg = rfe[:,:,rg_lat[0],:][:,:,:,rg_lon[0]].squeeze()

  return (rfe_rg,match_lon,match_lat,final_times)


'''
##############################################################################
############## Functions for opening temperature data ########################
##############################################################################
'''
'''
Grab ERA-I t2m over a certain region for a list of dates
dates format []
'''
def grab_t2m_region(inputfile,lonlim,latlim,dates):
  #print 'loading tamsat data'
  nc_fid	= Dataset(inputfile, 'r')
  lat		= np.array(nc_fid.variables['latitude'][:])		# extract/copy the data
  lon		= np.array(nc_fid.variables['longitude'][:])
  t_time = np.array(nc_fid.variables['t'][:])
  t_units	= str(nc_fid.variables['t'].units)			# this units is wrong!
  all_dates = nc4.num2date(t_time,t_units)

  # split dates into year, month and day columns
  final_times	= np.zeros((len(all_dates),4))
  for t in np.arange(0,len(all_dates)):
    curr_time		= all_dates[t]
    final_times[t,0]	= curr_time.year
    final_times[t,1]	= curr_time.month
    final_times[t,2]	= curr_time.day
    final_times[t,3]	= curr_time.hour				# why hour is needed?

  y1_id	= np.where(final_times[:,0]==dates[ 0,0])[0]
  y2_id	= np.where(final_times[:,0]==dates[-1,0])[0]
  m1_id	= np.where(final_times[:,1]==dates[ 0,1])[0]
  m2_id	= np.where(final_times[:,1]==dates[-1,1])[0]
  d1_id	= np.where(final_times[:,2]==dates[ 0,2])[0]
  d2_id	= np.where(final_times[:,2]==dates[-1,2])[0]

  # always start accumlation from hour 0 (assuming this is the accumulation between 0z-6z)
  h1_id	= np.where(final_times[:,3]==0)[0]
  ym1_id	= np.intersect1d(y1_id,m1_id)
  ym2_id	= np.intersect1d(y2_id,m2_id)
  ymd1_id	= np.intersect1d(ym1_id,d1_id)
  ymdh1_id	= np.intersect1d(ymd1_id,h1_id)
  ymd2_id	= np.intersect1d(ym2_id,d2_id)
  # have to add an extra index onto the end of time id to
  # rainfall accumulated upto zero Z
  time_id	= np.arange(ymdh1_id[0],ymd2_id[len(ymd2_id)-1]+1)	# value is 28, 6-hourly data for a week, and shift one day backward

  t2	= np.mean(np.array(nc_fid.variables['T2'][time_id,:,:,:]),axis=0)

  # shiftgrid: shift global lat/lon grid east or west
  # start, if True, 180. represents the starting of the new grid; if False, 180. represent the ending longitude.
  t2,lon	= shiftgrid(180.,t2,lon,start=False)

  # grab rfe at specific latitude and longitude (inlat,inlon)
  rg_lon	= np.where( (lon >= lonlim[0]) & (lon <= lonlim[1]) )
  rg_lat	= np.where( (lat >= latlim[0]) & (lat <= latlim[1]) )

  match_lon	= lon[rg_lon[0]]
  match_lat	= lat[rg_lat[0]]

  t2_rg	= t2[:,rg_lat[0],:][:,:,rg_lon[0]].squeeze()		# shape is time, lat, lon as shown above
  nc_fid.close()
  #rfe_rg[rfe_rg < 0] = 0

  return(t2_rg,match_lon,match_lat)

'''
2. Grab UKMO GloSea5-GC2 2m temp over a certain region
'''
def grab_ukmo_region_t2m(date,inputfile,lonlim,latlim):
  nc_fid	= Dataset(inputfile,'r')
  lat		= np.array(nc_fid.variables['latitude'][:] )
  lon		= np.array(nc_fid.variables['longitude'][:])
  rfe		= np.array(nc_fid.variables['t2m'][:])

  t_time	= np.array(nc_fid.variables['t'][:])
  # t_units	= str(nc_fid.variables['t'].units)# these units are wrong!
  t_units		='hours since '+date+' 00:00:00'

  rfe,lon	= shiftgrid(180.,rfe,lon,start=False)
  # return datetime objects given numeric time values.
  all_dates	= nc4.num2date(t_time,t_units)

  # split dates into year, month and day columns
  final_times	= np.zeros((len(all_dates),3))
  for t in np.arange(0,len(all_dates)):
    curr_time		= all_dates[t]
    final_times[t,0]	= curr_time.year
    final_times[t,1]	= curr_time.month
    final_times[t,2]	= curr_time.day

  # grab rfe at specific latitude and longitude (inlat,inlon)
  # this is interesting: if only condition is given for np.where, return condition.nonzero(), which is a tuple!
  rg_lon	= np.where( (lon>=lonlim[0]) & (lon<=lonlim[1]) )
  rg_lat	= np.where( (lat>=latlim[0]) & (lat<=latlim[1]) )

  # [0] here does nothing
  match_lon	= lon[rg_lon[0]]
  match_lat	= lat[rg_lat[0]]

  # indexing arrays cannot be broadcast togather
  # squeeze() here does nothing
  rfe_rg	= rfe[:,:,rg_lat[0],:][:,:,:,rg_lon[0]].squeeze()
  nc_fid.close()
  #rfe_rg[rfe_rg < 0] = 0

  return (rfe_rg,match_lon,match_lat,final_times)



'''
3. Grab NCEP t2m over a certain region
'''

def grab_ncep_region_t2m(date,inputfile,lonlim,latlim):
  nc_fid	= Dataset(inputfile,'r')
  lat		= np.array(nc_fid.variables['latitude'][:])
  lon		= np.array(nc_fid.variables['longitude'][:])
  rfe		= np.array(nc_fid.variables['t2m'][:])

  t_time	= np.array(nc_fid.variables['t'][:])
# t_units	= str(nc_fid.variables['t'].units)			# this units is wrong!
  units		='hours since '+date+' 00:00:00'

  rfe,lon	= shiftgrid(180.,rfe,lon,start=False)
  # t_time	= nc_fid.variables['time'][:]
# # convert time to date using netCDF4 function
  all_dates	= nc4.num2date(t_time,units)				# need to check this works
  # all_dates = t_time

  # split dates into year, month and day columns
  final_times	= np.zeros((len(all_dates),4))
  for t in np.arange(0,len(all_dates)):
    curr_time		= all_dates[t]
    final_times[t,0]	= curr_time.year
    final_times[t,1]	= curr_time.month
    final_times[t,2]	= curr_time.day
    final_times[t,3]	= curr_time.hour				# why hour is needed?

  # grab rfe at specific latitude and longitude (inlat,inlon)
  rg_lon	= np.where( (lon>=lonlim[0]) & (lon<=lonlim[1]) )
  rg_lat	= np.where( (lat>=latlim[0]) & (lat<=latlim[1]) )

  match_lon	= lon[rg_lon[0]]
  match_lat	= lat[rg_lat[0]]

  rfe_rg	= rfe[:,:,rg_lat[0],:][:,:,:,rg_lon[0]].squeeze()	# shape is time, lat, lon as shown above
  nc_fid.close()
  #rfe_rg[rfe_rg < 0] = 0

  return (rfe_rg,match_lon,match_lat,final_times)


'''
4. Grab ECMWF t2m over a certain region

# e.g. filepath:
/group_workspaces/jasmin2/klingaman/datasets/S2S/ECMWF/
ENSO/hindcasts_1120_2017/ecmwf_precip_ens_19981120.nc
'''

def grab_ecmf_region_t2m(date,inputfile,lonlim,latlim):
  #print 'loading tamsat data'
  nc_fid= Dataset(inputfile, 'r')
  lat	= np.array(nc_fid.variables['latitude'][:])
  lon	= np.array(nc_fid.variables['longitude'][:])
  #t_units = str(nc_fid.variables['t'].units) # time units are wrong
  t_time= np.array(nc_fid.variables['t'][:])
   # convert time to date using netCDF4 function
  units='hours since '+date+' 00:00:00'
  all_dates = nc4.num2date(t_time,units)
  # all_dates = t_time

  # split dates into year, month and day columns
  final_times = np.zeros((len(all_dates),4))
  for t in np.arange(0,len(all_dates)):
    curr_time = all_dates[t]
    final_times[t,0] = curr_time.year
    final_times[t,1] = curr_time.month
    final_times[t,2] = curr_time.day
    final_times[t,3] = curr_time.hour

  rfe = np.array(nc_fid.variables['t2m'][:])
  nc_fid.close()

  rfe,lon = shiftgrid(180.,rfe,lon,start=False)

  # grab rfe at specific latitude and longitude (inlat,inlon)
  rg_lon = np.where((lon >= lonlim[0]) & (lon <= lonlim[1]))
  rg_lat = np.where((lat >= latlim[0]) & (lat <= latlim[1]))

  match_lon = lon[rg_lon[0]]
  match_lat = lat[rg_lat[0]]

  rfe_rg = rfe[:,:,rg_lat[0],:][:,:,:,rg_lon[0]].squeeze()

  return (rfe_rg,match_lon,match_lat,final_times)
