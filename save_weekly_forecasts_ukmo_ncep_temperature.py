'''
Save S2S precipitation forecasts over a given region
from UKMO (initialised weekly) and NCEP
(initialised daily), w.r.t weekly forecasts
of UKMO (intialised once per week)

Saves lagged ensembles of NCEP forecasts
i.e. saves the 7 NCEP forecasts intialised
on and before each UKMO initialisation date

UKMO data are 1 degree resolution
NCEP data are 1.5 degree resolution

M. Young 29/06/2018
---
UKMO data has been changed since 16/08/2018
Inpot time record switchs from daily to 6-hourly.

L. Guo
'''
from __future__ import division
import glob
import numpy as np
from netCDF4 import Dataset
import time as tt
from datetime import datetime, timedelta
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import basemap
import os.path

# execfile() is a python built-in function
# it is different from the import statement in that it does not use the module administation - it reads the file unconditionally and does not create a new module.
execfile('date_str.py')
execfile('grab_data.py')

# input directories
dir_ukmo = '/group_workspaces/jasmin2/klingaman/datasets/S2S/UKMO/ENSO/'
dir_ncep = '/group_workspaces/jasmin2/klingaman/datasets/S2S/NCEP/ENSO/'
dir_ecmf = '/group_workspaces/jasmin2/klingaman/datasets/S2S/ECMWF/ENSO/'
dir_era= '/group_workspaces/jasmin2/klingaman/datasets/ERA-INTERIM/T2/'
era_t2m_file = dir_era+'t2m_1998-2017.nc'

# output directories to save forecast data
dir_out = '/group_workspaces/jasmin2/ncas_climate/users/myoung02/datasets/S2S_forecasts/weekly_dubstep_style/'

nleads = 5 # number of lead times (weeks)
nmembers_ukmo	= 7 # number of UKMO ensemble members
nmembers_ncep	= 4 # number of NCEP ensemble members
nlags_ncep	= 7 # number of lags for NCEP - for lagged ensemble means
# Years over which to
# Years available, UKMO 1996-2017, NCEP 1999-2010, ECMWF 1998-2016
years	= np.arange(2000,2010+1,1) #years to read and save


# Define region for analysis over Brazil
# region = 'Brazil'
# latlim = [-40,20]
# lonlim = [-90,-20]
region = 'Africa'
latlim = [-40,40]
lonlim = [-30,60]

#ukmo_lon = np.arange(lonlim[0],lonlim[1]+1,1)
#ukmo_lat = np.arange(latlim[0],latlim[1]+1,1)
#lons_ukmo_rg,lats_ukmo_rg = np.meshgrid(ukmo_lon,ukmo_lat)

# Note, following variables are (string) list
# list start dates for ukmo (just for DJF)
# change months for another seasons
month	= ['01','02','03','04','05','06','07','08','09','10','11','12']#['01','02','03','04','05','06','07','08','09','10','11','12']
starts_ukmo	= ['01','09','17','25']
# Use version 2017 for 25 March to 25 December and version 18 for 1Jan-17Mar
version_list1 = ['0801','0809','0817','0825',
                 '0901','0909','0917','0925',
                 '1001','1009','1017','1025',
                 '1101','1109','1117','1125',
                 '1201','1209','1217','1225']
version_list2 = ['0101','0109','0117','0125',
                 '0201','0209','0217','0225',
                 '0301','0309','0317','0325',
                 '0401','0409','0417','0425',
                 '0501','0509','0517','0525',
                 '0601','0609','0617','0625',
                 '0701','0709','0717','0725']

yc = 0
for y in years:
  yc = yc + 1
  print y

  # loop through UKMO start dates
  for m in np.arange(0,len(month)):
    for s in np.arange(0,len(starts_ukmo)):

      print '===================================='
      print str(y)+'-'+month[m]+'-'+starts_ukmo[s]
      curr_start= month[m]+starts_ukmo[s]	# get current start date
      curr_m	= int(month[m])			# month from start date
      curr_d	= int(starts_ukmo[s])		# day from start date

      # Version of UKMO to use
      if curr_start in version_list1:
        v_year = '2017'
      elif curr_start in version_list2:
        v_year = '2018'

      # create date time array for current date
      # datetime.strptime() creates a datetime object from a string representing a date and time and a corresponding format string.
      curr_datetime	= datetime.strptime(str(y)+curr_start,'%Y%m%d')
      curr_dir_ukmo	= dir_ukmo+'hindcasts_'+curr_start+'_'+v_year+'_6hr/'
      ukmo_fstamp	= 'ukmo_t2m_ens_'+str(y)+str(curr_start)+'.nc'
      f_ls = glob.glob(curr_dir_ukmo+ukmo_fstamp)
      print '------------'
      print 'opening ukmo'
      print '------------'
      print f_ls[0]+"	"+str(os.path.exists(f_ls[0]))
      k_date = curr_datetime.strftime('%Y-%m-%d')
      (ukmo_rfe,ukmo_lon,ukmo_lat,ukmo_dates)	= grab_ukmo_region_t2m(k_date,f_ls[0],lonlim,latlim)

      # create an empty array to store all the data
      if (s == 0):
        all_ukmo = np.zeros((nleads,len(starts_ukmo),nmembers_ukmo,len(ukmo_lat),len(ukmo_lon)))		# (5, 4, 7, 40, 47)

      # loop through the lead times, 6-hourly data
      for lead in np.arange(0,5):

        if lead == 0: # week 1
          week_id = np.arange(1,7+1)
        elif lead == 1: # week 2
          week_id = np.arange(8,14+1)
        elif lead == 2: # week 3
          week_id = np.arange(15,21+1)
        elif lead == 3: # week 4
          week_id = np.arange(22,28+1)
        elif lead == 4: # week 5
          week_id = np.arange(29,35+1)

        curr_dates	= ukmo_dates[week_id-1,:]
        sub_ukmo_rfe	= ukmo_rfe[:,week_id-1,:,:]

        # compute average of ukmo data for the current week, kg/m2 per day?
        all_ukmo[lead,s,:,:,:] = np.mean(sub_ukmo_rfe,axis=1)
        #
        # if lead == 0:		# only open the ncep, ecmwf and gpcp if lead == 0
        #   '''
        #   Open NCEP & compute -6d lagged ensemble mean
        #   '''
        #   print '------------'
        #   print 'opening ncep'
        #   print '------------'
        #   # get lagged dates from day 0 to 6 days before
        #   for l in np.arange(0,nlags_ncep):
        #     print 'reading lead '+str(l+1)+' of '+str(nlags_ncep)
        #     # subtract 'l' days from current date to get start
        #     # date to open lagged NCEP initialised forecasts
        #     lag_date		= curr_datetime - timedelta(days=l)
        #     # there is no 29th Feb in the leap year
        #     # the follow treatment will aount the 29th Feb twice. I pressuem this does not cause large bias.
        #     if lag_date.strftime('%m%d') == '0229':
        #       lag_date		= lag_date - timedelta(days=1)
        #     str_lag_start	= lag_date.strftime('%m%d')
        #     str_lag_year	= lag_date.strftime('%Y')
        #     curr_dir_ncep	= dir_ncep+'hindcasts_'+str_lag_start+'_6hr/'
        #     ncep_fstamp	= curr_dir_ncep+'ncep_t2m_ens_'+str_lag_year+str_lag_start+'.nc'
        #     date = lag_date.strftime('%Y-%m-%d')
        #     print ncep_fstamp+'	'+str(os.path.exists(ncep_fstamp))
        #     (ncep_rfe,ncep_lon,ncep_lat,ncep_dates) = grab_ncep_region_t2m(date,ncep_fstamp,lonlim,latlim)
        #
        #     # create an empty array to store all the data
        #     if (s == 0) & (l == 0):	# note, lead == 0
        #       all_ncep = np.zeros((nleads,len(starts_ukmo),nlags_ncep,nmembers_ncep,len(ncep_lat),len(ncep_lon)))	# (5, 4, 7, 4, 40, 47)
        #       # return coordinate matrices from coordinate vectors.
        #       lons_rg,lats_rg = np.meshgrid(ncep_lon,ncep_lat)		# both are 2d arrays (40, 47)
        #
        #     if l == 0:
        #       tmp_lens = np.zeros((nlags_ncep,nmembers_ncep,len(ncep_dates),len(ncep_lat),len(ncep_lon)))	# ( 7, 4, 176, 40, 47)
        #       tmp_dates	= np.zeros((nlags_ncep,len(ncep_dates),4))						# ( 7, 176, 4)
        #     # store lags and lagged dates
        #     tmp_lens[l,:,:,:,:]	= np.copy(ncep_rfe)
        #     tmp_dates[l,:,:] = np.copy(ncep_dates)
        #
        # '''
        # Find lagged ensembles for the current week
        # '''
        # # find lagged ensemble for current lead time date NCEP
        # for l in np.arange(0,nlags_ncep):
        #   ncep_y1_id	= np.where(tmp_dates[l,:,0]==curr_dates[ 0,0])[0]
        #   ncep_y2_id	= np.where(tmp_dates[l,:,0]==curr_dates[-1,0])[0]
        #   ncep_m1_id	= np.where(tmp_dates[l,:,1]==curr_dates[ 0,1])[0]
        #   ncep_m2_id	= np.where(tmp_dates[l,:,1]==curr_dates[-1,1])[0]
        #   ncep_d1_id	= np.where(tmp_dates[l,:,2]==curr_dates[ 0,2])[0]
        #   ncep_d2_id	= np.where(tmp_dates[l,:,2]==curr_dates[-1,2])[0]
        #   ncep_ym1_id	= np.intersect1d(ncep_y1_id,ncep_m1_id)
        #   ncep_ym2_id	= np.intersect1d(ncep_y2_id,ncep_m2_id)
        #   ncep_ymd1_id	= np.intersect1d(ncep_ym1_id ,ncep_d1_id)
        #   ncep_ymd2_id	= np.intersect1d(ncep_ym2_id ,ncep_d2_id)
        #   # have to add an extra index onto the end of time id to
        #   # rainfall accumulated upto zero Z
        #   ncep_time_id	= np.arange(ncep_ymd1_id[0],ncep_ymd2_id[0]+1,1)	# value is 28, 6-hourly data for a week, and shift one day backward
        #
        #   # store lagged ensemble
        #   #print tmp_lens.shape				# (7, 4, 176, 40, 47)
        #   #print tmp_lens[l,:,:,ncep_time_id,:].shape	# (28, 4, 40, 47)
        #   # why the time dimention is permutated to the left?
        #   all_ncep[lead,s,l,:,:,:] = np.sum(tmp_lens[l,:,ncep_time_id,:,:],axis=0)/7		# units: kg/m2 per day?

        '''
        Open ERA-I 2m temperature data
        '''
        # print '----------------------------'
        # print 'opening ERA-I 2m temperature'
        # print '----------------------------'
        #
        #
        # # open mean ERA-I t2m data for the week
        # t2 = []
        # t2,era_lon,era_lat = grab_t2m_region(era_t2m_file,lonlim,latlim,curr_dates)
        #
        # # regrid to ncep grid
        # t2_rg	= []	# it is a list, but what does this line do?
        # # interpolate data on a rectilinear grid (with x = gpcp_lon, y = gpcp_lat) to a grid with x =lons_rg, y = lats_rg.
        # # gpcp_p is a rank-2 array with 1st dimension corresponding to y, 2nd dimension x.
        # # gpcp_lon, gpcp_lat are rank-1 arrays containing x and y of gpcp_p grid in INCREASING order.
        # # lons_rg, lats_rg are rand-2 arrays containing x and y of desired output grid.
        # t2_rg	= basemap.interp(t2.squeeze(),era_lon,era_lat[::-1],lons_rg,lats_rg[::-1,:],order=1)	# (40, 47)
        #
        # if (lead == 0) & (s == 0):
        #   t2_era_rg = np.zeros((nleads,len(starts_ukmo),len(ncep_lat),len(ncep_lon)))			# (5, 4, 40, 47)
        #
        # t2_era_rg[lead,s,:,:] = np.copy(t2_rg)				# units: mm/week, turn latitude to [20,-40]



    # at end of each month
    # save UKMO as netcdf
    nlats	= len(ukmo_lat)
    nlons	= len(ukmo_lon)
    fname_ukmo= dir_out+region+'_UKMO_ensemble_weekly_t2m_hindcasts_lead1-5_'+str(y)+'_'+month[m]+'.nc'
    dataset	= Dataset(fname_ukmo,'w',format='NETCDF3_CLASSIC')
    lat		= dataset.createDimension('lat',nlats )
    lon		= dataset.createDimension('lon',nlons)
    lead_nc	= dataset.createDimension('lead',nleads)
    week_nc	= dataset.createDimension('week',None)
    ens_nc	= dataset.createDimension('ensemble',nmembers_ukmo)

    # create variables
    latitudes	= dataset.createVariable('latitude','f',('lat',))
    longitudes	= dataset.createVariable('longitude','f',('lon',))
    weeks	= dataset.createVariable('week',np.float64,('week',))
    leads	= dataset.createVariable('lead',np.float64,('lead',))
    ensembles	= dataset.createVariable('ensemble',np.float64,('ensemble',))
    p_mean	= dataset.createVariable('week_t2m','d',('week','lead','ensemble','lat','lon'))

    # Global Attributes (will need modified accordingly)
    dataset.description	= 'UKMO GloSea5-GC2 ensembles weekly mean 2m temperature'
    dataset.history	= 'Created '+tt.ctime(tt.time())
    dataset.source	= 'Subset by M. Young'
    # Variable Attributes
    latitudes.units	= 'degrees_north'
    longitudes.units	= 'degrees_east'
    weeks.units		= 'weeks since '+str(y)+'-'+month[m]+'-01'
    weeks.calendar	= 'none'
    p_mean.units	= 'degrees Kelvin (K)'
    p_mean.long_name	= 'Weekly mean 2m temperature'

    # Fill variables with data
    latitudes[:]	= ukmo_lat
    longitudes[:]	= ukmo_lon
    p_mean[:]	= np.swapaxes(all_ukmo,0,1)
    weeks[:] = np.arange(1,len(starts_ukmo)+1)
    leads[:] = np.arange(1,nleads+1)
    ensembles[:] = np.arange(1,nmembers_ukmo+1)
    dataset.close()
    print ''
    print 'saved ukmo file for '+str(y)+'-'+month[m]

    '''
    # save NCEP as netcdf
    nlats	= len(ncep_lat)
    nlons	= len(ncep_lon)
    fname_ncep= dir_out+region+'_NCEP_ensemble_weekly_t2m_hindcasts_lead1-5_'+str(y)+'_'+month[m]+'.nc'
    dataset	= Dataset(fname_ncep,'w',format='NETCDF3_CLASSIC')
    lat	= dataset.createDimension('lat',nlats)
    lon	= dataset.createDimension('lon',nlons)
    lead_nc	= dataset.createDimension('lead',nleads)
    week_nc	= dataset.createDimension('week',None)
    lag_nc	= dataset.createDimension('lag',nlags_ncep)
    ens_nc	= dataset.createDimension('ensemble',nmembers_ncep)

    # create variables
    latitudes	= dataset.createVariable('latitude','f',('lat',))
    longitudes= dataset.createVariable('longitude','f',('lon',))
    leads	= dataset.createVariable('lead',np.float64,('lead',))
    weeks	= dataset.createVariable('week',np.float64,('week',))
    lags	= dataset.createVariable('lag',np.float64,('lag',))
    ensembles	= dataset.createVariable('ensemble',np.float64,('ensemble',))
    p_mean	= dataset.createVariable('week_t2m','d',('week','lead','lag','ensemble','lat','lon'))

    # Global Attributes (will need modified accordingly)
    dataset.description	= 'NCEP lagged ensembles weekly mean 2m temperature'
    dataset.history	= 'Created ' + tt.ctime(tt.time())
    dataset.source	= 'Subset by M. Young'
    # Variable Attributes
    latitudes.units	= 'degrees_north'
    longitudes.units	= 'degrees_east'
    weeks.units		= 'weeks since '+str(y)+'-'+month[m]+'-01'
    weeks.calendar	= 'none'
    p_mean.units	= 'degrees Kelvin (K)'
    p_mean.long_name	= 'Weekly mean 2m temperature'

    # Fill variables with data
    latitudes[:]	= ncep_lat
    longitudes[:]	= ncep_lon
    p_mean[:]		= np.swapaxes(all_ncep,0,1)
    leads[:]		= np.arange(1,nleads+1)
    weeks[:]		= np.arange(1,len(starts_ukmo)+1)
    lags[:]		= np.arange(0,nlags_ncep)
    ensembles[:]	= np.arange(1,nmembers_ncep+1)
    dataset.close()
    print ''
    print 'saved ncep file for '+str(y)+'-'+month[m]


    # save GPCP  as netcdf
    nlats	= len(ncep_lat)
    nlons	= len(ncep_lon)
    fname_gpcp	= dir_out+region+'_ERA-I_weekly_t2m_'+str(y)+'_'+month[m]+'.nc'
    dataset	= Dataset(fname_gpcp,'w',format='NETCDF3_CLASSIC')
    lat		= dataset.createDimension('lat' ,nlats)
    lon		= dataset.createDimension('lon' ,nlons)
    lead_nc	= dataset.createDimension('lead',nleads)
    week_nc	= dataset.createDimension('week',None )

    # create variables
    latitudes	= dataset.createVariable('latitude','f',('lat',))
    longitudes	= dataset.createVariable('longitude','f',('lon',))
    leads	= dataset.createVariable('lead',np.float64,('lead',))
    weeks	= dataset.createVariable('week',np.float64,('week',))
    p_mean	= dataset.createVariable('week_t2m','d',('week','lead','lat','lon'))

    # Global Attributes (will need modified accordingly)
    dataset.description	= 'ERA-I weekly mean 2m temperature at forecast res, based on UKMO intialisation dates'
    dataset.history	= 'Created ' + tt.ctime(tt.time())
    dataset.source	= 'Subset by M. Young'
    # Variable Attributes
    latitudes.units	= 'degrees_north'
    longitudes.units	= 'degrees_east'
    weeks.units		= 'weeks since '+str(y)+'-'+month[m]+'-01'
    weeks.calendar	= 'none'
    p_mean.units	= 'degrees Kelvin (K)'
    p_mean.long_name	= 'Weekly mean 2m temperature'

    # Fill variables with data
    latitudes[:]	= ncep_lat
    longitudes[:]	= ncep_lon
    p_mean[:]		= np.swapaxes(t2_era_rg,0,1)			# convert weekly mean of mm d rather than weekly total
    leads[:]		= np.arange(1,nleads+1)
    weeks[:]		= np.arange(1,len(starts_ukmo)+1)
    dataset.close()
    print ''
    print 'saved ERA-I file for '+str(y)+'-'+month[m]
    '''
# # save gpcp data as netcdf
# nlats = len(gpcp_lat)
# nlons = len(gpcp_lon)
# fname_gpcp = dir_out+region+'_GPCP_weekly_DJF_UKMO_times_'+str(y)+'.nc'
# dataset = Dataset(fname_gpcp,'w',format='NETCDF3_CLASSIC')
# lat = dataset.createDimension('lat', nlats) # create lat (dims depend on region)
# lon = dataset.createDimension('lon', nlons) # create lon
# lead_nc = dataset.createDimension('lead',nleads) # create time
# week_nc = dataset.createDimension('week',len(starts_ukmo)) # create time
# # create variables
# latitudes = dataset.createVariable('latitude','f',('lat',))
# longitudes = dataset.createVariable('longitude','f',('lon',))
# weeks = dataset.createVariable('week',np.float64,('week',))
# leads = dataset.createVariable('lead',np.float64,('lead',))
# p_mean = dataset.createVariable('week_precip','d',('lead','week','lat','lon'))

# # Global Attributes (will need modified accordingly)
# dataset.description = 'GPCP weekly mean precipitation at 1 degree res, based on UKMO intialisation dates'
# dataset.history = 'Created ' + tt.ctime(tt.time())
# dataset.source = 'Subset by M. Young'
# # Variable Attributes
# latitudes.units = 'degrees_north'
# longitudes.units = 'degrees_east'
# f_week = starts_ukmo[0]
# weeks.units = 'weeks since'+str(y-1)+'-'+f_week[0:2]+'-'+f_week[2:4]
# weeks.calendar = 'gregorian'
# p_mean.units = 'mm d-1'
# p_mean.long_name = 'Weekly mean precipitation'

# # Fill variables with data
# latitudes[:] = gpcp_lat
# longitudes[:] = gpcp_lon
# p_mean[:] = gpcp4ukmo/7 # convert weekly mean of mm d rather than weekly total
# weeks[:] = np.arange(1,len(starts_ukmo)+1)
# leads[:] = np.arange(1,nleads+1)
# dataset.close()
# print 'saved gpcp file for '+str(y)
