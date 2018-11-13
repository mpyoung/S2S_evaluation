'''
Analyse skill of weekly S2S 2m temperature forecasts
from UKMO, NCEP and ECMWF over a given region

1) Reads in data saved over the given region from
'save_weekly_forecasts_ukmo_ncep_temperature.py'
'save_weekly_forecasts_ecmf_temperature.py'

Also reads in probabilistic percentiles (e.g. tercile) forecasts saved from:
'save_weekly_forecast_percentile_probabilities_temperature.py'

2) Computes mean temperature, bias, anomaly correlation coefficient
and brier & ranked probability skill scores

3) Plots results as spatial maps

M. Young 29/06/2018
'''
from __future__ import division
import glob
import numpy as np
from netCDF4 import Dataset
from netCDF4 import MFDataset
import time as tt
from datetime import datetime, timedelta
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import calendar
import os.path
execfile('date_str.py')
execfile('grab_data.py')

dir_in	= '/group_workspaces/jasmin2/ncas_climate/users/myoung02/datasets/S2S_forecasts/weekly_dubstep_style/'
dir_out	= '/home/users/myoung02/S2S_evaluation_results/seasonal/'

# number of ensemble members for each forecast
ukmo_nmembers = 7
ncep_nmembers = 4
ecmf_nmembers = 11
# number of lagged ensembles
ncep_lags = 7
ecmf_lags = 3

# Define region for analysis over Brazil
# region	= 'Brazil'
# latlim	= [-40,20]
# lonlim	= [-90,-20]
region = 'Africa'
latlim = [-40,40]
lonlim = [-30,60]

plot_all = 0 # option to plot bias, rmse etc (probably leave)
# define percentile categories to analyse
percentiles = [100./3.,200./3.]
p_cat = 0.333 # probability of falling into tercile
p_ptile = np.array([1.0/3.0,2.0/3.0,1]) # cumulative probability for RPSS
p_name = 'tercile' # give categories an appropriate name
p_name_ls = ['Below normal','Normal','Above normal']

percentiles = [100./3.,200./3.]
p_cat = 0.333 # probability of falling into tercile
p_ptile = np.array([1.0/3.0,2.0/3.0,1]) # cumulative probability for RPSS
p_name = 'BC_tercile' # give categories an appropriate name
p_name_ls = ['Below normal','Normal','Above normal']

# percentiles = [20.0,40.0,60.0,80.0]
# p_cat = 0.2 # probability of falling into category
# p_ptile = np.array([0.2,0.4,0.6,0.8,1.0]) # cumulative probability for RPSS
# p_name = 'quintile' # give categories an appropriate name
# p_name_ls = ['Below 20','20 to 40','40 to 60','60 to 80','Above 80']


years	= np.arange(2000,2010+1,1)	# december will always correspond to year-1
nleads	= 5				# number of lead times (in weeks) in the data
nweeks = 20 # number of weeks in season
#season	= 'JJA'
seasons = ['JJA']#['DJF','MAM','JJA','SON']

'''
Function: 'mask_percentiles':
Function takes in a 2d spatial forecast array and 2d spatial array of corresponding percentiles
and masks the forecast array according to input arguments, 'arg':
'between': mask the forecast between 2 percentiles
'below': mask the forecast below a percentile
'above': mask the forecast above a percentile
'''
def mask_percentiles(data,percentile,arg):
  data_mask= np.copy(data)

  if arg == 'between': # find between 2 percentiles
    # the method of masking below doesn't work
    #data_mask[(data_mask >= percentile[0,:,:]) & (data_mask <= percentile[1,:,:])] = np.nan
    idx = np.where((data_mask >= percentile[0,:,:]) & (data_mask <= percentile[1,:,:]))
    data_mask[idx] = np.nan
    data_mask[np.isnan(data_mask)==False] = 0
    data_mask[np.isnan(data_mask)==True] = 1

  elif arg == 'below': # find data below percentile
    data_mask[data_mask < percentile] = np.nan
    data_mask[data_mask >= percentile] = 0
    data_mask[np.isnan(data_mask)==True] = 1

  elif arg == 'above': # find data above percentile
    data_mask[data_mask > percentile] = np.nan
    data_mask[data_mask <= percentile] = 0
    data_mask[np.isnan(data_mask)==True] = 1

  return data_mask


for season in seasons:

  if season == 'DJF':
    # for DJF, prepend five weeks starting 1025
    starts_ukmo = ['1001','1009','1017','1025',
                   '1101','1109','1117','1125',
                   '1201','1209','1217','1225',
                   '0101','0109','0117','0125',
                   '0201','0209','0217','0225']
    starts_mon = ['10','11','12','01','02']
    starts_day = ['01','09','17','25']
    keep_mon = [12,1,2]

  if season == 'MAM':
    starts_ukmo = ['0101','0109','0117','0125',
                   '0201','0209','0217','0225',
                   '0301','0309','0317','0325',
                   '0401','0409','0417','0425',
                   '0501','0509','0517','0525']
    starts_mon = ['01','02','03','04','05']
    starts_day = ['01','09','17','25']
    keep_mon = [3,4,5]

  if season == 'JJA':
    starts_ukmo = ['0401','0409','0417','0425',
                   '0501','0509','0517','0525',
                   '0601','0609','0617','0625',
                   '0701','0709','0717','0725',
                   '0801','0809','0817','0825']
    starts_mon = ['04','05','06','07','08']
    starts_day = ['01','09','17','25']
    keep_mon = [6,7,8]

  if season == 'SON':
    starts_ukmo = ['0701','0709','0717','0725',
                   '0801','0809','0817','0825',
                   '0901','0909','0917','0925',
                   '1001','1009','1017','1025',
                   '1101','1109','1117','1125']
    starts_mon = ['07','08','09','10','11']
    starts_day = ['01','09','17','25']
    keep_mon = [9,10,11]

  '''
  Read data for one season and for all years
  Useful arrays:
  	all_ukmo_ens
  	all_necp_ens
  	all_ecmf_ens
  	all_gpcp
  '''
  week_mid_dates = np.zeros((len(years),nleads,20,3))			# (11, 5, 20, 3)

  yc = 0
  for y in years:

    print y

    yc = yc + 1

    # create a mask to mask out all dates not within DJF
    ukmo_file	= []
    ncep_file	= []
    ecmf_file	= []
    gpcp_file	= []

    i = 0
    for mo in np.arange(0,len(starts_mon)):

      if season == 'DJF' and starts_mon[mo] in ['10','11','12']:
        y_ec	= y-1
      else:
        y_ec	= np.copy(y)

      ukmo_file.append( dir_in+region+'_UKMO_ensemble_weekly_t2m_hindcasts_lead1-5_'+str(y_ec)+'_'+starts_mon[mo]+'.nc' )
      print ukmo_file[mo]+'	'+str(os.path.exists(ukmo_file[mo]))
      ncep_file.append( dir_in+region+'_NCEP_ensemble_weekly_t2m_hindcasts_lead1-5_'+str(y_ec)+'_'+starts_mon[mo]+'.nc' )
      print ncep_file[mo]+'	'+str(os.path.exists(ncep_file[mo]))
      ecmf_file.append( dir_in+region+'_ECMWF_ensemble_weekly_t2m_hindcasts_lead1-5_'+str(y_ec)+'_'+starts_mon[mo]+'.nc' )
      print ecmf_file[mo]+'	'+str(os.path.exists(ecmf_file[mo]))
      gpcp_file.append( dir_in+region+'_ERA-I_weekly_t2m_'+str(y_ec)+'_'+starts_mon[mo]+'.nc' )
      print gpcp_file[mo]+'	'+str(os.path.exists(gpcp_file[mo]))

      for dd in starts_day:

        # Get start date as date object
        curr_ec_datetime = datetime.strptime(str(y_ec)+starts_mon[mo]+dd,'%Y%m%d')
        print curr_ec_datetime

        for l in np.arange(0,nleads):
          # get dates of the 4th day of each week (this keeps samples consistent at each lead time)
          # Note:
          #  1) if you use the date of the first day in the week this excludes forecasts
          #     with start-dates e.g. 30th November, causing different samples at each lead time.
          #  2) if you use the date of the last day in the week this excludes forecasts where
          #     the end of the week stretches over into the next month.
          week_mid_dates[yc-1,l,i,0] = (curr_ec_datetime + timedelta(days=((l)*7)+3)).day
          week_mid_dates[yc-1,l,i,1] = (curr_ec_datetime + timedelta(days=((l)*7)+3)).month
          week_mid_dates[yc-1,l,i,2] = (curr_ec_datetime + timedelta(days=((l)*7)+3)).year

  #       curr_mid = str(int(week_mid_dates[yc-1,l,i,0])).zfill(2)+'/'+str(int(week_mid_dates[yc-1,l,i,1])).zfill(2)+'/'+str(int(week_mid_dates[yc-1,l,i,2]))
  #       print 'lead time week '+str(l+1)+' '+curr_mid

        i = i + 1

    # UKMO Data, units is kg/m2/day, which needs converted to mm/d. However, nothing need to be done.
    # UKMO has 7 ensemble members
    nc_fid = MFDataset(ukmo_file)
    if y == years[0]:
      ukmo_lat = np.array(nc_fid.variables['latitude'][:])
      ukmo_lon = np.array(nc_fid.variables['longitude'][:])
      all_ukmo_ens = np.zeros((len(years),nleads,20,ukmo_nmembers,len(ukmo_lat),len(ukmo_lon)))			# (11, 5, 20, 7,   , 40, 47)
    all_ukmo_ens[yc-1,:,:,:,:,:] = np.swapaxes(np.array(nc_fid.variables['week_t2m'][:]),0,1)-273.15
    nc_fid.close()

    # NCEP units don't need converted. NCEP has 7 lags and 4 ensemble members
    nc_fid = MFDataset(ncep_file)
    if y == years[0]:
      ncep_lat = np.array(nc_fid.variables['latitude'][:])
      ncep_lon = np.array(nc_fid.variables['longitude'][:])
      all_ncep_ens = np.zeros((len(years),nleads,20,ncep_lags,ncep_nmembers,len(ncep_lat),len(ncep_lon)))		# (11, 5, 20, 7,  4, 40, 47)
    all_ncep_ens[yc-1,:,:,:,:,:] = np.swapaxes(np.array(nc_fid.variables['week_t2m'][:]),0,1)-273.15
    nc_fid.close()

    # ECMWF units don't need converted. There are 3 lags and 11 ensemble members
    nc_fid = MFDataset(ecmf_file)
    if y == years[0]:
      ecmf_lat = np.array(nc_fid.variables['latitude'][:])
      ecmf_lon = np.array(nc_fid.variables['longitude'][:])
      all_ecmf_ens = np.zeros((len(years),nleads,20,ecmf_lags,ecmf_nmembers,len(ncep_lat),len(ncep_lon)))
    all_ecmf_ens[yc-1,:,:,:,:,:] = np.swapaxes(np.array(nc_fid.variables['week_t2m'][:]),0,1)-273.15			# (11, 5, 20, 3, 11, 40, 47)
    nc_fid.close()

    # GPCP
    nc_fid = MFDataset(gpcp_file)
    if y == years[0]:
      gpcp_lat = np.array(nc_fid.variables['latitude'][:])
      gpcp_lon = np.array(nc_fid.variables['longitude'][:])
      all_gpcp = np.zeros((len(years),nleads,20,len(gpcp_lat),len(gpcp_lon)))					# (11, 5, 20,  ,   , 40, 47)
    all_gpcp[yc-1,:,:,:,:] = np.swapaxes(np.array(nc_fid.variables['week_t2m'][:]),0,1)-273.15
    nc_fid.close()


  '''
  Mask data for a season, compute lag/ensemble means, and create a dry mask
  '''
  # for DJF, any other month is masked
  # the size is (440, 3), becasuse date_mask is tuple, therefore, it needs 3 dimensions to identify one position
  # the total is 1,100, 8 out of 20 months have been masked out.
  if season == 'DJF':
    date_mask = np.where((week_mid_dates[:,:,:,1] >  2) & (week_mid_dates[:,:,:,1] < 12))
  if season == 'MAM':
    date_mask = np.where((week_mid_dates[:,:,:,1] >  5) & (week_mid_dates[:,:,:,1] <  3))
  if season == 'JJA':
    date_mask = np.where((week_mid_dates[:,:,:,1] >  8) & (week_mid_dates[:,:,:,1] <  6))
  if season == 'SON':
    date_mask = np.where((week_mid_dates[:,:,:,1] > 11) & (week_mid_dates[:,:,:,1] <  9))
  all_ukmo_ens[date_mask[0],date_mask[1],date_mask[2],:,:,:] = np.nan
  all_ncep_ens[date_mask[0],date_mask[1],date_mask[2],:,:,:,:] = np.nan
  all_ecmf_ens[date_mask[0],date_mask[1],date_mask[2],:,:,:,:] = np.nan
  all_gpcp[date_mask[0],date_mask[1],date_mask[2],:,:] = np.nan

  # test the mask
  for y in np.arange(0,11):
    for l in np.arange(0,5):
      for i in np.arange(0,20):
        curr_mid = str(int(week_mid_dates[y,l,i,0])).zfill(2)+'/'+str(int(week_mid_dates[y,l,i,1])).zfill(2)+'/'+str(int(week_mid_dates[y,l,i,2]))
        print str(years[y])+'	lead week '+str(l)+'	'+starts_ukmo[i]+'	'+curr_mid+'	'+str(all_gpcp[y,l,i,20,23])

  # Forecast ensemble and lagged means
  # compute the arithmetic mean along the specific axis, ignoring NaNs
  all_ukmo = np.nanmean(all_ukmo_ens,axis=3)								# (11, 5, 20, 40, 47)
  all_ncep = np.nanmean(all_ncep_ens,axis=(3,4))
  all_ecmf = np.nanmean(all_ecmf_ens,axis=(3,4))

  # average forecasts over lags only (leaving ensembles)
  lag_ncep = np.nanmean(all_ncep_ens,axis=3)								# (11, 5, 20, 4, 40, 47)
  lag_ecmf = np.nanmean(all_ecmf_ens,axis=3)								# (11, 5, 20, 11, 40, 47)

  # Forecast means over years and ensembles (& lags for ncep/ecmf)
  week_ukmo_mean = np.nanmean(all_ukmo_ens,axis=(0,3))							# (5, 20, 40, 47)
  week_ncep_mean = np.nanmean(all_ncep_ens,axis=(0,3,4))
  week_ecmf_mean = np.nanmean(all_ecmf_ens,axis=(0,3,4))

  # GPCP mean for each week
  week_gpcp_mean = np.nanmean(all_gpcp,axis=0)								# (5, 20, 40, 47)

  # Create a dry mask (GPCP seasonal clim < 1 mmd)
  # this is used to mask Brier Skill Score later
  # as BSS results are dodgy over dry regions
  dry_mask = np.nanmean(np.copy(week_gpcp_mean),axis=(0,1))
  dry_mask[dry_mask <  1] = np.nan
  dry_mask[dry_mask >= 1] = 1

  '''
  Compute the Anomaly correlation coefficient (ACC), bias and RMSE
  '''
  # set up empty arrays
  ukmo_acc_pre = np.zeros((3,nleads,len(ukmo_lat),len(ukmo_lon)))						# (3, 5, 40, 47)
  ncep_acc_pre = np.zeros((3,nleads,len(ncep_lat),len(ncep_lon)))
  ecmf_acc_pre = np.zeros((3,nleads,len(ecmf_lat),len(ecmf_lon)))

  ukmo_bias = np.zeros((nleads,len(ukmo_lat),len(ukmo_lon)))						# (5, 40, 47)
  ncep_bias = np.zeros((nleads,len(ncep_lat),len(ncep_lon)))
  ecmf_bias = np.zeros((nleads,len(ecmf_lat),len(ecmf_lon)))
  ukmo_week_bias = np.zeros((nleads,20,len(ukmo_lat),len(ukmo_lon)))					# (5, 20, 40, 47)
  ncep_week_bias = np.zeros((nleads,20,len(ncep_lat),len(ncep_lon)))
  ecmf_week_bias = np.zeros((nleads,20,len(ecmf_lat),len(ecmf_lon)))

  ukmo_rmse = np.zeros((nleads,len(ukmo_lat),len(ukmo_lon)))						# (5, 40, 47)
  ncep_rmse = np.zeros((nleads,len(ncep_lat),len(ncep_lon)))
  ecmf_rmse = np.zeros((nleads,len(ecmf_lat),len(ecmf_lon)))

  nsamp = np.zeros((3,nleads))										# count number of samples at each lead time
  nsamp_week = np.zeros((3,nleads,20))
  nsamp_week_only = np.zeros((nleads))

  for y in np.arange(0,len(years)):

    curr_n = np.zeros((nleads))		# count the number of samples for each lead time in each year
    # compute climatologies independant of current year
    tmp_gpcp = np.delete(all_gpcp,y,axis=0) # (11,5, 20, 40, 47)
    week_gpcp_mean = np.nanmean(tmp_gpcp,axis=0)
    tmp_ukmo_ens = np.delete(all_ukmo_ens,y,axis=0)# (11,5, 20, 40, 47)
    week_ukmo_mean = np.nanmean(tmp_ukmo_ens,axis=(0,3))
    tmp_ncep_ens = np.delete(all_ncep_ens,y,axis=0)
    week_ncep_mean = np.nanmean(tmp_ncep_ens,axis=(0,3,4))
    tmp_ecmf_ens = np.delete(all_ecmf_ens,y,axis=0)
    week_ecmf_mean = np.nanmean(tmp_ecmf_ens,axis=(0,3,4))
    for w in np.arange(0,20):		# loop through starts
      for l in np.arange(0,nleads):	# loop through lead times

        # Only use forecasts in DJF since data is masked with NaN's for different months:
        if week_mid_dates[y,l,w,1] in keep_mon:
          print 'lead time week '+str(l+1)+' mid-week date: '+str(int(week_mid_dates[y,l,w,0])).zfill(2)+'/'+str(int(week_mid_dates[y,l,w,1])).zfill(2)+'/'+str(int(week_mid_dates[y,l,w,2]))

          nsamp_week_only[l] = nsamp_week_only[l] + 1	# count total number of weeks in season at each lead
          curr_n[l] = curr_n[l] + 1

          # Compute current anomalies for GPCP
          curr_gpcp_anom = all_gpcp[y,l,w,:,:] - week_gpcp_mean[l,w,:,:]

          # overall bias throughout whole season, weekly bias and overall rmse
          for e in np.arange(0,ecmf_nmembers):
            if e < ukmo_nmembers:
              nsamp[0,l] = nsamp[0,l] + 1		# counts the number of samples at each week lead time
              nsamp_week[0,l,w] = nsamp_week[0,l,w] + 1
              ukmo_bias[l,:,:] = ukmo_bias[l,:,:] + (all_ukmo_ens[y,l,w,e,:,:] - all_gpcp[y,l,w,:,:])
              ukmo_week_bias[l,w,:,:] = ukmo_week_bias[l,w,:,:] + (all_ukmo_ens[y,l,w,e,:,:] - all_gpcp[y,l,w,:,:])
              ukmo_rmse[l,:,:] = ukmo_rmse[l,:,:] + (all_ukmo_ens[y,l,w,e,:,:] - all_gpcp[y,l,w,:,:])**2

            if e < ncep_nmembers:
              for g in np.arange(0,ncep_lags):	# loop through NCEP lags (7)
                nsamp[1,l] = nsamp[1,l] + 1
                nsamp_week[1,l,w] = nsamp_week[1,l,w] + 1
                ncep_bias[l,:,:] = ncep_bias[l,:,:] + (all_ncep_ens[y,l,w,g,e,:,:] - all_gpcp[y,l,w,:,:])
                ncep_week_bias[l,w,:,:] = ncep_week_bias[l,w,:,:] + (all_ncep_ens[y,l,w,g,e,:,:] - all_gpcp[y,l,w,:,:])
                ncep_rmse[l,:,:] = ncep_rmse[l,:,:] + (all_ncep_ens[y,l,w,g,e,:,:] - all_gpcp[y,l,w,:,:])**2

            for g in np.arange(0,ecmf_lags):	# loop through ECMWF lags (3)
              nsamp[2,l] = nsamp[2,l] + 1
              nsamp_week[2,l,w] = nsamp_week[2,l,w] + 1
              ecmf_bias[l,:,:] = ecmf_bias[l,:,:] + (all_ecmf_ens[y,l,w,g,e,:,:] - all_gpcp[y,l,w,:,:])
              ecmf_week_bias[l,w,:,:] = ecmf_week_bias[l,w,:,:] + (all_ecmf_ens[y,l,w,g,e,:,:] - all_gpcp[y,l,w,:,:])
              ecmf_rmse[l,:,:] = ecmf_rmse[l,:,:] + (all_ecmf_ens[y,l,w,g,e,:,:] - all_gpcp[y,l,w,:,:])**2

          # for anomaly corelation coefficent just compare lagged, ensemble means
          curr_ukmo_anom = []
          curr_ukmo_anom = all_ukmo[y,l,w,:,:] - week_ukmo_mean[l,w,:,:]
          ukmo_acc_pre[0,l,:,:] = ukmo_acc_pre[0,l,:,:] + (curr_ukmo_anom*curr_gpcp_anom)
          ukmo_acc_pre[1,l,:,:] = ukmo_acc_pre[1,l,:,:] + curr_ukmo_anom**2
          ukmo_acc_pre[2,l,:,:] = ukmo_acc_pre[2,l,:,:] + curr_gpcp_anom**2

          curr_ncep_anom = []
          curr_ncep_anom = all_ncep[y,l,w,:,:] - week_ncep_mean[l,w,:,:]
          ncep_acc_pre[0,l,:,:] = ncep_acc_pre[0,l,:,:] + (curr_ncep_anom*curr_gpcp_anom)
          ncep_acc_pre[1,l,:,:] = ncep_acc_pre[1,l,:,:] + curr_ncep_anom**2
          ncep_acc_pre[2,l,:,:] = ncep_acc_pre[2,l,:,:] + curr_gpcp_anom**2

          curr_ecmf_anom = []
          curr_ecmf_anom = all_ecmf[y,l,w,:,:] - week_ecmf_mean[l,w,:,:]
          ecmf_acc_pre[0,l,:,:] = ecmf_acc_pre[0,l,:,:] + (curr_ecmf_anom*curr_gpcp_anom)
          ecmf_acc_pre[1,l,:,:] = ecmf_acc_pre[1,l,:,:] + curr_ecmf_anom**2
          ecmf_acc_pre[2,l,:,:] = ecmf_acc_pre[2,l,:,:] + curr_gpcp_anom**2

        else: # if not DJF
          ukmo_week_bias[l,w,:,:] = np.nan
          ncep_week_bias[l,w,:,:] = np.nan
          ncep_week_bias[l,w,:,:] = np.nan
          nsamp_week[:,l,w] = np.nan

  # Compute ACC
  ukmo_acc = ukmo_acc_pre[0,:,:,:]/np.sqrt(ukmo_acc_pre[1,:,:,:]*ukmo_acc_pre[2,:,:,:])
  ncep_acc = ncep_acc_pre[0,:,:,:]/np.sqrt(ncep_acc_pre[1,:,:,:]*ncep_acc_pre[2,:,:,:])
  ecmf_acc = ecmf_acc_pre[0,:,:,:]/np.sqrt(ecmf_acc_pre[1,:,:,:]*ecmf_acc_pre[2,:,:,:])

  for l in np.arange(0,nleads):
    ukmo_bias[l,:,:]	= ukmo_bias[l,:,:]/nsamp[0,l]
    ncep_bias[l,:,:]	= ncep_bias[l,:,:]/nsamp[1,l]
    ecmf_bias[l,:,:]	= ecmf_bias[l,:,:]/nsamp[2,l]
    ukmo_rmse[l,:,:]	= np.sqrt(ukmo_rmse[l,:,:]/nsamp[0,l])
    ncep_rmse[l,:,:]	= np.sqrt(ncep_rmse[l,:,:]/nsamp[1,l])
    ecmf_rmse[l,:,:]	= np.sqrt(ecmf_rmse[l,:,:]/nsamp[2,l])
    for w in np.arange(0,20):
      ukmo_week_bias[l,w,:,:]	= ukmo_week_bias[l,w,:,:]/nsamp_week[0,l,w]
      ncep_week_bias[l,w,:,:]	= ncep_week_bias[l,w,:,:]/nsamp_week[0,l,w]
      ecmf_week_bias[l,w,:,:]	= ecmf_week_bias[l,w,:,:]/nsamp_week[0,l,w]

  # alternatively compute bias using lagged ensemble means
  # (there is a difference to the above method but the difference is v. small)
  ukmo_bias_s	= np.nanmean(all_ukmo - all_gpcp,axis=(0,2))
  ncep_bias_s	= np.nanmean(all_ncep - all_gpcp,axis=(0,2))
  ecmf_bias_s	= np.nanmean(all_ecmf - all_gpcp,axis=(0,2))


  '''
  Forecast probabilities
  '''
  fname_p_gpcp = dir_in+'stats/probabilities/'+region+'_ERA-I_t2m_'+p_name.replace('BC_','')+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'
  fname_p_ukmo = dir_in+'stats/probabilities/'+region+'_UKMO_t2m_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'
  fname_p_ncep = dir_in+'stats/probabilities/'+region+'_NCEP_t2m_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'
  fname_p_ecmf = dir_in+'stats/probabilities/'+region+'_ECMWF_t2m_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'

  if os.path.exists(fname_p_ukmo) == True:

    tmp_p_gpcp = np.load(fname_p_gpcp)
    tmp_p_ukmo = np.load(fname_p_ukmo)
    tmp_p_ncep = np.load(fname_p_ncep)
    tmp_p_ecmf = np.load(fname_p_ecmf)

    p_gpcp = np.copy(tmp_p_gpcp)
    p_ukmo = np.copy(tmp_p_ukmo)
    p_ncep = np.copy(tmp_p_ncep)
    p_ecmf = np.copy(tmp_p_ecmf)

    n_week_p = np.zeros((nleads))
    nan_week_p = np.zeros((nleads))

    for y in np.arange(0,len(years)):
      for w in np.arange(0,20):# loop through starts
        for l in np.arange(0,nleads):
          if week_mid_dates[y,l,w,1] in keep_mon:
            n_week_p[l] = n_week_p[l] +1
            print 'lead time week '+str(l+1)+' mid-week date:'
          else:
            nan_week_p[l] = nan_week_p[l] + 1
            p_gpcp[:,l,w,y,:,:] = np.nan
            p_ukmo[:,:,l,w,y,:,:] = np.nan
            p_ncep[:,:,l,w,y,:,:] = np.nan
            p_ecmf[:,:,l,w,y,:,:] = np.nan

    #make reliability diagram
    # take each tercile and lead time and find obs frequency for
    # each forecast probability bin
    p_bin = np.arange(0,1+0.2,0.2)
    reliability_fcst = np.zeros((3,len(percentiles)+1,nleads,2,len(p_bin)))
    n_fcst = np.zeros((3,len(percentiles)+1,nleads,2,len(p_bin)))
    for ter in np.arange(0,len(percentiles)+1):
      for l in np.arange(0,nleads):
        curr_gpcp = np.copy(p_gpcp[ter,l,:,:,:,:].flatten())
        for r in [0,1]:
          curr_ukmo = np.copy(p_ukmo[ter,r,l,:,:,:,:].flatten())
          curr_ncep = np.copy(p_ncep[ter,r,l,:,:,:,:].flatten())
          curr_ecmf = np.copy(p_ecmf[ter,r,l,:,:,:,:].flatten())

          for p in np.arange(0,len(p_bin)):
            if p == 0:
              p_id_ukmo = np.where((curr_ukmo == p_bin[p]))[0]
              p_id_ncep = np.where((curr_ncep == p_bin[p]))[0]
              p_id_ecmf = np.where((curr_ecmf == p_bin[p]))[0]
            elif 0 < p < len(p_bin)-1:
              p_id_ukmo = np.where((curr_ukmo >= p_bin[p]) & (curr_ukmo < p_bin[p+1]))[0]
              p_id_ncep = np.where((curr_ncep >= p_bin[p]) & (curr_ncep < p_bin[p+1]))[0]
              p_id_ecmf = np.where((curr_ecmf >= p_bin[p]) & (curr_ecmf < p_bin[p+1]))[0]
            elif p == len(p_bin)-1:
              p_id_ukmo = np.where((curr_ukmo == p_bin[p]))[0]
              p_id_ncep = np.where((curr_ncep == p_bin[p]))[0]
              p_id_ecmf = np.where((curr_ecmf == p_bin[p]))[0]

            reliability_fcst[0,ter,l,r,p] = np.sum(curr_gpcp[p_id_ukmo])/float(len(p_id_ukmo))
            reliability_fcst[1,ter,l,r,p] = np.sum(curr_gpcp[p_id_ncep])/float(len(p_id_ncep))
            reliability_fcst[2,ter,l,r,p] = np.sum(curr_gpcp[p_id_ecmf])/float(len(p_id_ecmf))
            n_fcst[0,ter,l,r,p] = float(len(p_id_ukmo))
            n_fcst[1,ter,l,r,p] = float(len(p_id_ncep))
            n_fcst[2,ter,l,r,p] = float(len(p_id_ecmf))

    # compute the climatological frequency
    f_clim = np.zeros((len(percentiles)+1,nleads))
    for ter in np.arange(0,len(percentiles)+1):
      for l in np.arange(0,nleads):
        f_clim[ter,l] = np.nanmean(p_gpcp[ter,l,:,:,:])


    '''
    fname_plot = dir_out+region+'_S2S_weekly_Rel_Diagram_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
    size1 = [8,7]	# width, height in inches
    fig = plt.figure(figsize=(size1[0],size1[1]))
    cnt = 0
    for tr in np.arange(0,len(percentiles)+1):
      for f in [0,1,2]:
        cnt = cnt + 1
        ax1 = plt.subplot(3,3,cnt)
        plt.plot(p_bin,p_bin,'-k',linewidth=1,label='_nolegend_') # one-to-one line
        plt.plot(np.repeat(f_clim[tr,0],len(p_bin)),p_bin,'--k',color='grey',label='_nolegend_') # Pf = Pc
        plt.plot(p_bin,np.repeat(f_clim[tr,0],len(p_bin)),'--k',color='grey',label='_nolegend_') # Pobs = Pc
        plt.plot(p_bin,(p_bin+f_clim[tr,0])/2,'--k',color='grey',label='_nolegend_') # no skill line
        # y1 = (p_bin+f_clim[tr,0])/2
        # y2 = p_bin
        # plt.fill_between(p_bin,y1,y2 ,color='grey', alpha='0.5')
        # plt.fill_between(p_bin,(p_bin+f_clim[tr,0])/2,p_bin,color='grey', alpha='0.5')
        for l in np.arange(0,nleads):
          plt.plot(p_bin,reliability_fcst[f,tr,l,1,:],'-ok',color=lead_col[l],ms=5,linewidth=1.5,alpha=0.8,label=lead_ls[l])
        plt.xlim([0,1])
        plt.ylim([0,1])

        if tr == len(percentiles):
          plt.xlabel('$P_{forecast}$')
        if f == 0:
          plt.ylabel('$P_{observations}$')
        if tr == 0:
          plt.title(model_ls[f])
        if (tr == len(percentiles)) & (f == 1):
          plt.legend(lead_ls,loc=0,prop={'size':8},ncol=2)

    plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=1.3)
    # plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
    plt.show()
    plt.close()
    '''

    model_ls= ['UKMO','NCEP','ECMWF']
    fmt_ls	= ['-o','-x','-+']
    ms_ls	= [8,10,10]
    mw_ls	= [1,2,2]
    col_ls	= ['firebrick','dodgerblue','forestgreen']
    leads	= np.arange(1,nleads+1,1)
    lead_col = ['blue','green','yellow','orange','red']
    lead_ls = ['W-1','W-2','W-3','W-4','W-5']
    fname_plot = dir_out+region+'_S2S_weekly_t2m_Rel_Diagram_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])

    size1 = [9,2.7*(len(percentiles)+1)]
    ff,axarr = plt.subplots(len(percentiles)+1,3,figsize=(size1[0],size1[1]))
    cnt = 0
    for tr in np.arange(0,len(percentiles)+1):
      for f in [0,1,2]:
        cnt = cnt + 1
        axarr[tr,f].plot(p_bin,p_bin,'-k',linewidth=1,label='_nolegend_') # one-to-one line
        axarr[tr,f].plot(np.repeat(f_clim[tr,0],len(p_bin)),p_bin,'--k',color='grey',label='_nolegend_') # Pf = Pc
        axarr[tr,f].plot(p_bin,np.repeat(f_clim[tr,0],len(p_bin)),'--k',color='grey',label='_nolegend_') # Pobs = Pc
        axarr[tr,f].plot(p_bin,(p_bin+f_clim[tr,0])/2,'--k',color='grey',label='_nolegend_') # no skill line
        # y1 = (p_bin+f_clim[tr,0])/2
        # y2 = p_bin
        # plt.fill_between(p_bin,y1,y2 ,color='grey', alpha='0.5')
        # plt.fill_between(p_bin,(p_bin+f_clim[tr,0])/2,p_bin,color='grey', alpha='0.5')
        for l in np.arange(0,nleads):
          axarr[tr,f].plot(p_bin,reliability_fcst[f,tr,l,1,:],'-ok',color=lead_col[l],ms=5,linewidth=1.5,alpha=0.8,label=lead_ls[l])
        axarr[tr,f].set_xlim([0,1])
        axarr[tr,f].set_ylim([0,1])

        if tr == len(percentiles):
          axarr[tr,f].set_xlabel('$P_{forecast}$')
        if f == 0:
          axarr[tr,f].set_ylabel(p_name_ls[tr]+'\n'+'$P_{observations}$',fontweight='bold')
        if tr == 0:
          axarr[tr,f].set_title(model_ls[f])
        if (tr == len(percentiles)) & (f == 2):
          axarr[tr,f].legend(lead_ls,loc=4,prop={'size':8},ncol=2)

        # this is an inset axes over the main axes
        ax_ins = inset_axes(axarr[tr,f],width="50%",height=0.7,loc=2,bbox_to_anchor=(0.06,1-0.48,1,.5), bbox_transform=axarr[tr,f].transAxes) # bbox (x, y, width, height) t
        t_width = 0
        for l in np.arange(0,nleads):
          ax_ins.bar(p_bin+t_width,n_fcst[f,tr,l,1,:],width=0.2/6.0,align='edge',color=lead_col[l],label=lead_ls[l])
          t_width = t_width+(0.2/6.0)
        ax_ins.set_xticks(p_bin[::2])
        ax_ins.set_xticklabels(p_bin[::2].round(2).astype('S'))
        ax_ins.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax_ins.tick_params(axis='both', labelsize=8)
    plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
    #plt.show()
    plt.close()


    # BSS
    nw = n_week_p[0] # number of data samples
    gpcp_bs_clim = np.zeros((len(percentiles)+1,nleads,len(ukmo_lat),len(ukmo_lon)))
    ukmo_bs = np.zeros((len(percentiles)+1,2,nleads,len(ukmo_lat),len(ukmo_lon)))
    ncep_bs = np.zeros((len(percentiles)+1,2,nleads,len(ukmo_lat),len(ukmo_lon)))
    ecmf_bs = np.zeros((len(percentiles)+1,2,nleads,len(ukmo_lat),len(ukmo_lon)))
    # brier score
    for l in np.arange(0,nleads):
      for ter in np.arange(0,len(percentiles)+1):
        gpcp_bs_clim[ter,l,:,:] = np.nansum((p_cat - p_gpcp[ter,l,:,:,:,:])**2,axis=(0,1))/nw
        for r in [0,1]:
          ukmo_bs[ter,r,l,:,:] = np.nansum((p_ukmo[ter,r,l,:,:,:,:] - p_gpcp[ter,l,:,:,:,:])**2,axis=(0,1))/nw # brier score sum
          ncep_bs[ter,r,l,:,:] = np.nansum((p_ncep[ter,r,l,:,:,:,:] - p_gpcp[ter,l,:,:,:,:])**2,axis=(0,1))/nw # brier score sum
          ecmf_bs[ter,r,l,:,:] = np.nansum((p_ecmf[ter,r,l,:,:,:,:] - p_gpcp[ter,l,:,:,:,:])**2,axis=(0,1))/nw # brier score sum

    # gpcp_bs_clim = gpcp_bs_clim/nw	# nw = 220
    ukmo_bss = np.zeros(ukmo_bs.shape)
    ncep_bss = np.zeros(ncep_bs.shape)
    ecmf_bss = np.zeros(ecmf_bs.shape)
    for q in [0,1]:
      ukmo_bss[:,q,:,:,:] = (gpcp_bs_clim - (ukmo_bs[:,q,:,:,:])) / gpcp_bs_clim
      ncep_bss[:,q,:,:,:] = (gpcp_bs_clim - (ncep_bs[:,q,:,:,:])) / gpcp_bs_clim
      ecmf_bss[:,q,:,:,:] = (gpcp_bs_clim - (ecmf_bs[:,q,:,:,:])) / gpcp_bs_clim

    rps_clim = np.zeros((2,nleads,nweeks,len(years),len(ukmo_lat),len(ukmo_lon)))
    rps_ukmo = np.zeros((2,nleads,nweeks,len(years),len(ukmo_lat),len(ukmo_lon)))
    rps_ncep = np.zeros((2,nleads,nweeks,len(years),len(ukmo_lat),len(ukmo_lon)))
    rps_ecmf = np.zeros((2,nleads,nweeks,len(years),len(ukmo_lat),len(ukmo_lon)))
    rpss_ukmo = []
    rpss_ncep = []
    rpss_ecmf = []
    # Compute ranked probability skill score
    for l in np.arange(nleads):
      for r in [0,1]:
        clim_tot = np.zeros((nweeks,len(years),len(ukmo_lat),len(ukmo_lon)))

        for pp in np.arange(0,len(p_ptile)):
          rps_clim[r,l,:,:,:,:] = rps_clim[r,l,:,:,:,:] + (p_ptile[pp] - np.sum(p_gpcp[0:pp+1,l,:,:,:,:],axis=0))**2

          rps_ukmo[r,l,:,:,:,:] = rps_ukmo[r,l,:,:,:,:] + (np.sum(p_ukmo[0:pp+1,r,l,:,:,:,:],axis=0) - np.sum(p_gpcp[0:pp+1,l,:,:,:,:],axis=0))**2

          rps_ncep[r,l,:,:,:,:] = rps_ncep[r,l,:,:,:,:] + (np.sum(p_ncep[0:pp+1,r,l,:,:,:,:],axis=0) - np.sum(p_gpcp[0:pp+1,l,:,:,:,:],axis=0))**2

          rps_ecmf[r,l,:,:,:,:] = rps_ecmf[r,l,:,:,:,:] + (np.sum(p_ecmf[0:pp+1,r,l,:,:,:,:],axis=0) - np.sum(p_gpcp[0:pp+1,l,:,:,:,:],axis=0))**2

        # rps_ukmo[r,l,:,:,:,:] = (p_ukmo[0,r,l,:,:,:,:] - p_gpcp[0,l,:,:,:,:])**2 +
        # ((p_ukmo[0,r,l,:,:,:,:] + p_ukmo[1,r,l,:,:,:,:]) - (p_gpcp[0,l,:,:,:,:] + p_gpcp[1,l,:,:,:,:]))**2 + ((p_ukmo[0,r,l,:,:,:,:] + p_ukmo[1,r,l,:,:,:,:] + p_ukmo[2,r,l,:,:,:,:]) - (p_gpcp[0,l,:,:,:,:] + p_gpcp[1,l,:,:,:,:] + p_gpcp[2,l,:,:,:,:]))**2
        #
        # rps_ncep[r,l,:,:,:,:] = (p_ncep[0,r,l,:,:,:,:] - p_gpcp[0,l,:,:,:,:])**2 + ((p_ncep[0,r,l,:,:,:,:] + p_ncep[1,r,l,:,:,:,:]) - (p_gpcp[0,l,:,:,:,:] + p_gpcp[1,l,:,:,:,:]))**2 + ((p_ncep[0,r,l,:,:,:,:] + p_ncep[1,r,l,:,:,:,:] + p_ncep[2,r,l,:,:,:,:]) - (p_gpcp[0,l,:,:,:,:] + p_gpcp[1,l,:,:,:,:] + p_gpcp[2,l,:,:,:,:]))**2
        #
        # rps_ecmf[r,l,:,:,:,:] = (p_ecmf[0,r,l,:,:,:,:] - p_gpcp[0,l,:,:,:,:])**2 + ((p_ecmf[0,r,l,:,:,:,:] + p_ecmf[1,r,l,:,:,:,:]) - (p_gpcp[0,l,:,:,:,:] + p_gpcp[1,l,:,:,:,:]))**2 + ((p_ecmf[0,r,l,:,:,:,:] + p_ecmf[1,r,l,:,:,:,:] + p_ecmf[2,r,l,:,:,:,:]) - (p_gpcp[0,l,:,:,:,:] + p_gpcp[1,l,:,:,:,:] + p_gpcp[2,l,:,:,:,:]))**2

    rpss_ukmo = 1 - (np.nanmean(rps_ukmo,axis=(2,3))/np.nanmean(rps_clim,axis=(2,3)))
    rpss_ncep = 1 - (np.nanmean(rps_ncep,axis=(2,3))/np.nanmean(rps_clim,axis=(2,3)))
    rpss_ecmf = 1 - (np.nanmean(rps_ecmf,axis=(2,3))/np.nanmean(rps_clim,axis=(2,3)))

  '''
  Plot metrics
  - regional averages at each lead time
  - spatial maps of mean precipitation, model biases, ACC, BSS
  '''
  print 'Plotting figs...'

  if plot_all == 1:
    # Average ACC, RMSE and BIAS over a sub-region at each lead time
    sublon = [-60,-35]
    sublat = [-25,-5]
    subregion_ukmo_lat_id	= np.where((ukmo_lat >= sublat[0]) & (ukmo_lat <= sublat[1]))[0]
    subregion_ukmo_lon_id	= np.where((ukmo_lon >= sublon[0]) & (ukmo_lon <= sublon[1]))[0]
    subregion_ncep_lat_id	= np.where((ncep_lat >= sublat[0]) & (ncep_lat <= sublat[1]))[0]
    subregion_ncep_lon_id	= np.where((ncep_lon >= sublon[0]) & (ncep_lon <= sublon[1]))[0]

    acc_region	= np.zeros((3,nleads))
    rmse_region	= np.zeros((3,nleads))
    bias_region	= np.zeros((3,nleads))
    acc_region[0,:]	= np.nanmean(ukmo_acc[:,subregion_ukmo_lat_id,:][:,:,subregion_ukmo_lon_id],axis=(1,2))
    acc_region[1,:]	= np.nanmean(ncep_acc[:,subregion_ncep_lat_id,:][:,:,subregion_ncep_lon_id],axis=(1,2))
    acc_region[2,:]	= np.nanmean(ecmf_acc[:,subregion_ncep_lat_id,:][:,:,subregion_ncep_lon_id],axis=(1,2))

    bias_region[0,:]= np.nanmean(ukmo_bias[:,subregion_ukmo_lat_id,:][:,:,subregion_ukmo_lon_id],axis=(1,2))
    bias_region[1,:]= np.nanmean(ncep_bias[:,subregion_ncep_lat_id,:][:,:,subregion_ncep_lon_id],axis=(1,2))
    bias_region[2,:]= np.nanmean(ecmf_bias[:,subregion_ncep_lat_id,:][:,:,subregion_ncep_lon_id],axis=(1,2))

    rmse_region[0,:]= np.nanmean(ukmo_rmse[:,subregion_ukmo_lat_id,:][:,:,subregion_ukmo_lon_id],axis=(1,2))
    rmse_region[1,:]= np.nanmean(ncep_rmse[:,subregion_ncep_lat_id,:][:,:,subregion_ncep_lon_id],axis=(1,2))
    rmse_region[2,:]= np.nanmean(ecmf_rmse[:,subregion_ncep_lat_id,:][:,:,subregion_ncep_lon_id],axis=(1,2))

    model_ls= ['UKMO','NCEP','ECMWF']
    fmt_ls	= ['-o','-x','-+']
    ms_ls	= [8,10,10]
    mw_ls	= [1,2,2]
    col_ls	= ['firebrick','dodgerblue','forestgreen']
    leads	= np.arange(1,nleads+1,1)
    #
    #
    #
    # fname_plot = dir_out+region+'_S2S_weekly_ACC_'+season+'_region_average_'+str(years[0])+'_'+str(years[len(years)-1])
    # size1 = [5,3.5]	# width, height in inches
    # fig = plt.figure(figsize=(size1[0],size1[1]))
    # for p in [0,1,2]:
    #   # ms - markersize
    #   # mew - markeredgewidth
    #   # alpha - 0.0 transparent; 1.0 opaque
    #   plt.plot(leads,acc_region[p,:],fmt_ls[p],color=col_ls[p],linewidth=2,ms=ms_ls[p],alpha=0.9,mew=mw_ls[p])
    # plt.ylim([0,0.62])
    # plt.xlim([0.8,5.2])
    # plt.xticks(np.arange(1,nleads+1,1))
    # plt.ylabel('ACC')
    # plt.xlabel('Lead time (weeks)')
    # plt.legend(model_ls,loc=0)	# loc=0/'best'
    # plt.savefig(fname_plot+'.pdf',bbox_inches='tight')	# not sure that I understand 'bbox'
    # #plt.show()
    # plt.close()
    #
    #
    #
    # fname_plot = dir_out+region+'_S2S_weekly_RMSE_'+season+'_region_average_'+str(years[0])+'_'+str(years[len(years)-1])
    # size1 =[5,3.5]
    # fig = plt.figure(figsize=(size1[0],size1[1]))
    # for p in [0,1,2]:
    #   plt.plot(leads,rmse_region[p,:],fmt_ls[p],color=col_ls[p],linewidth=2,ms=ms_ls[p],alpha=0.9,mew=mw_ls[p])
    # plt.ylim([0.9*np.min(rmse_region),1.1*np.max(rmse_region)])
    # plt.xlim([0.8,5.2])
    # plt.xticks(np.arange(1,nleads+1,1))
    # plt.xticks(leads)
    # plt.ylabel('RMSE (mm d$^{-1}$)')
    # plt.xlabel('Lead time (weeks)')
    # plt.legend(model_ls,loc=0)
    # plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
    # #plt.show()
    # plt.close()
    #
    #
    #
    # fname_plot = dir_out+region+'_S2S_weekly_Bias_'+season+'_region_average_'+str(years[0])+'_'+str(years[len(years)-1])
    # size1 =[5,3.5]
    # fig = plt.figure(figsize=(size1[0],size1[1]))
    # for p in [0,1,2]:
    #   plt.plot(leads,bias_region[p,:],fmt_ls[p],color=col_ls[p],linewidth=2,ms=ms_ls[p],alpha=0.9,mew=mw_ls[p])
    # plt.plot(leads,np.zeros(len(leads)),'-k')
    # plt.ylim([1.1*np.min(bias_region),1.1*np.max(bias_region)])
    # plt.xlim([0.8,5.2])
    # plt.xticks(np.arange(1,nleads+1,1))
    # plt.ylabel('Bias (mm d$^{-1}$)')
    # plt.xlabel('Lead time (weeks)')
    # plt.legend(model_ls,loc=0,prop={'size': 10})
    # plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
    # #plt.show()
    # plt.close()



    # Plot seasonal mean precipitation
    # web colours
    # these are hexadecimal colour codes
    # a hex triplet is a six-digit, three-byte hexadecimal number used in many computing applications to represent colours.
    # bytes represent the red, green and blue components of the colour.
    # one byte represent a number in the range 00 to FF (in hexadecimal notation) or 0 to 255 in decimal notation.
    #precip_colors	= ["#fdfdfd","#f2f2f2","#bfbfbf","#04e9e7","#019ff4","#0300f4","#02fd02","#01c501","#008e00","#fdf802","#e5bc00","#fd9500","#fd0000", "#d40000","#bc0000","#f800fd","#9854c6"]
    #precip_colormap	= matplotlib.colors.ListedColormap(precip_colors)
    #cols		= precip_colormap
    # a sequential colormap
    # more details can be found https://matplotlib.org/tutorials/colors/colormaps.html
    cols = 'OrRd'
    cmin = 15
    cmax = 35
    cspc = 2.5
    clevs = np.arange(cmin,cmax+cspc,cspc)
    clabel1 = 'Mean 2m temperature ($^{\circ}C$)'
    # matplotlib.colors.BoundaryNorm generates a colormap index based on discrete intervals.
    # BoundaryNorm maps values to integers.
    # BoundaryNorm defines the edges of bins, and data falling within a bin is mapped to the color with the same index.
    # if the number of bins doesn't equal ncolors, the color is chosen by linear interpolation of the bin number onto color numbers.
    norm = BoundaryNorm(boundaries=clevs,ncolors=256)
    lw = 1
    gl = 20

    '''
    size1 = [5,3.5]
    fname_plot = dir_out+region+'_GPCC_weekly_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
    fig = plt.figure(figsize=(size1[0],size1[1]))
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
    mymap.drawcoastlines(linewidth=lw)
    mymap.drawcountries(linewidth=lw)
    # labels list of 4 values that control whether parallels are labelled where they intersect the left, right, top or bottom of the plot.
    # +/-, north and south latitudes are labelled with "+" and "-", otherwise they are labelled with "N" and "S".
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x,y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
    # create a psesudocolor plot with a non-regular rectuanglar grid
    # norm scales the data values to the canonocal colormap range [0,1] for mapping to colors
    uncal = mymap.pcolormesh(x,y,np.nanmean(week_gpcp_mean,axis=(0,1)),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.title('GPCP '+season)
    plt.colorbar(uncal,label=clabel1,extend='max')
    plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
    plt.show()
    plt.close()
    '''

    nrow = 4
    ncol = nleads
    size1 =[15,12]
    fname_plot = dir_out+region+'_S2S_weekly_t2m_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
    fig = plt.figure(figsize=(size1[0],size1[1]))
    for n in np.arange(0,nleads):
      if n == 2:
        plt.subplot(nrow,ncol,3)
        mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
        mymap.drawcoastlines(linewidth=lw)
        #mymap.drawcountries(linewidth=lw)
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
        mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
        x,y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
    #   uncal = mymap.pcolormesh(x,y,np.nanmean(week_gpcp_mean,axis=(0,1)),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
        uncal = mymap.pcolormesh(x,y,np.nanmean(week_gpcp_mean[0,:,:,:],axis=0),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
        plt.title('ERA-I')

      plt.subplot(nrow,ncol,n+6)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x,y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
    # uncal = mymap.pcolormesh(x,y,np.nanmean(week_ukmo_mean,axis=(0,1)),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      uncal = mymap.pcolormesh(x,y,np.nanmean(week_ukmo_mean[n,:,:,:],axis=0),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('UKMO Week '+str(n+1))

      plt.subplot(nrow,ncol,n+11)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x,y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
    # uncal = mymap.pcolormesh(x,y,np.nanmean(week_ncep_mean,axis=(0,1)),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      uncal = mymap.pcolormesh(x,y,np.nanmean(week_ncep_mean[n,:,:,:],axis=0),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('NCEP Week '+str(n+1))

      plt.subplot(nrow,ncol,n+16)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x,y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
    # uncal = mymap.pcolormesh(x,y,np.nanmean(week_ecmf_mean,axis=(0,1)),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      uncal = mymap.pcolormesh(x,y,np.nanmean(week_ecmf_mean[n,:,:,:],axis=0),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('ECMWF Week '+str(n+1))

    plt.suptitle(season,fontweight='bold',fontsize=14,horizontalalignment='left')
    plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=1.3)
    fig.subplots_adjust(right=0.90)
    cbar_pos = [0.92, 0.10, 0.015, 0.35] #[left, bottom, width, height]
    cbar_ax = fig.add_axes(cbar_pos)
    cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='both')
    plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
    #plt.show()
    plt.close()



    '''
    # plot ACC just at lead times Week 1, Week 3 and Week 5 (save on panels)
    cols = 'RdYlBu'
    cmin = -1
    cmax = 1
    cspc = 0.2
    clevs = np.arange(cmin,cmax+cspc,cspc)
    clabel1 = 'ACC'
    norm = BoundaryNorm(boundaries=clevs, ncolors=256)

    lw = 1
    nrow = 3
    ncol = 3
    size1=[10,9.5]
    gl = 20
    fname_plot = dir_out+region+'_S2S_weekly_sub_ACC_DJF_'+str(years[0])+'_'+str(years[len(years)-1])
    fig = plt.figure(figsize=(size1[0],size1[1]))
    nc = 0
    for n in [0,2,4]:
      nc = nc + 1
      plt.subplot(nrow,ncol,nc)
      mymap = Basemap(projection='cyl',resolution='l',\
              llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
              llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if nc in [1,4,7]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
      uncal = mymap.pcolormesh(x,y,ukmo_acc[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)

      plt.title('UKMO Week '+str(n+1))

      nc = nc + 1
      plt.subplot(nrow,ncol,nc)
      mymap = Basemap(projection='cyl',resolution='l',\
              llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
              llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if nc in [1,4,7]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
      uncal = mymap.pcolormesh(x,y,ncep_acc[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('NCEP Week '+str(n+1))

      nc = nc + 1
      plt.subplot(nrow,ncol,nc)
      mymap = Basemap(projection='cyl',resolution='l',\
              llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
              llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if nc in [1,4,7]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
      uncal = mymap.pcolormesh(x,y,ecmf_acc[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('ECMWF Week '+str(n+1))

    plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.6)
    fig.subplots_adjust(right=0.90)
    cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]
    cbar_ax = fig.add_axes(cbar_pos)
    cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='both')
    plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
    plt.show()
    plt.close()
    '''



    # plot ACC
    cols = 'RdYlBu'
    cmin = -1
    cmax = 1
    cspc = 0.2
    clevs = np.arange(cmin,cmax+cspc,cspc)
    clabel1 = 'ACC'
    norm = BoundaryNorm(boundaries=clevs, ncolors=256)

    lw = 1
    nrow = 3
    ncol = nleads
    size1=[15,9]
    gl = 20
    fname_plot = dir_out+region+'_S2S_weekly_t2m_ACC_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
    fig = plt.figure(figsize=(size1[0],size1[1]))
    for n in np.arange(0,nleads):
      plt.subplot(nrow,ncol,n+1)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      #mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
      uncal = mymap.pcolormesh(x,y,ukmo_acc[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)

      plt.title('UKMO Week '+str(n+1))

      plt.subplot(nrow,ncol,n+6)
      mymap = Basemap(projection='cyl',resolution='l',\
              llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
              llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      #mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
      uncal = mymap.pcolormesh(x,y,ncep_acc[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('NCEP Week '+str(n+1))

      plt.subplot(nrow,ncol,n+11)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x,y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
      uncal = mymap.pcolormesh(x,y,ecmf_acc[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('ECMWF Week '+str(n+1))

    plt.suptitle(season,fontweight='bold',fontsize=14,horizontalalignment='left')
    plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.5)
    fig.subplots_adjust(right=0.90)
    cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]
    cbar_ax = fig.add_axes(cbar_pos)
    cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='both')
    plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
    #plt.show()
    plt.close()



    # plot Bias
    cols = 'RdBu_r'
    cmin = -5
    cmax = 5
    cspc = 1
    clevs = np.arange(cmin,cmax+cspc,cspc)
    clabel1 = 'Bias ($^{\circ}C$)'
    norm = BoundaryNorm(boundaries=clevs, ncolors=256)

    lw = 1
    nrow = 3
    ncol = nleads
    gl = 20
    fname_plot = dir_out+region+'_S2S_weekly_t2m_Bias_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
    fig = plt.figure(figsize=(15,9))
    for n in np.arange(0,nleads):
      plt.subplot(nrow,ncol,n+1)#indexes goes from 1 to nrows * ncols
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      #mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x,y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
      tmp_bias = []
      tmp_bias = np.copy(ukmo_bias[n,:,:])
      tmp_bias[tmp_bias == 0] = np.nan
      uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('UKMO Week '+str(n+1))

      plt.subplot(nrow,ncol,n+6)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      #mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x,y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
      tmp_bias = []
      tmp_bias = np.copy(ncep_bias[n,:,:])
      tmp_bias[tmp_bias == 0] = np.nan
      uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('NCEP Week '+str(n+1))

      plt.subplot(nrow,ncol,n+11)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
      tmp_bias = []
      tmp_bias = np.copy(ecmf_bias[n,:,:])
      tmp_bias[tmp_bias == 0] = np.nan
      uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('ECMWF Week '+str(n+1))

    plt.suptitle(season,fontweight='bold',fontsize=14,horizontalalignment='left')
    plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.5)
    fig.subplots_adjust(right=0.90)
    cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]
    cbar_ax = fig.add_axes(cbar_pos)
    cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='both')
    plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
    #plt.show()
    plt.close()



    '''
    # plot bias for lead times week 1 and 5 only
    nrow = 2
    ncol = 3
    size1=[10,7]
    week_want = [0,4] # index for weeks want
    fname_plot = dir_out+region+'_S2S_weekly_1&5_Bias_DJF_'+str(years[0])+'_'+str(years[len(years)-1])
    fig = plt.figure(figsize=(size1[0],size1[1]))
    nc = 0
    for n in week_want:
      nc = nc + 1
      plt.subplot(nrow,ncol,nc)
      mymap = Basemap(projection='cyl',resolution='l',\
              llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
              llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if nc in [1,4]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
      uncal = mymap.pcolormesh(x,y,ukmo_bias[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('UKMO Week '+str(n+1))
      nc = nc + 1

      plt.subplot(nrow,ncol,nc)
      mymap = Basemap(projection='cyl',resolution='l',\
              llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
              llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      mymap.drawcountries(linewidth=lw)

      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if nc in [1,4]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
      #uncal = mymap.contourf(x,y,ukmo_acc[n,:,:],clevs,cmap=cols,extend='both')
      uncal = mymap.pcolormesh(x,y,ncep_bias[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('NCEP Week '+str(n+1))

      nc = nc + 1
      plt.subplot(nrow,ncol,nc)
      mymap = Basemap(projection='cyl',resolution='l',\
              llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
              llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if nc in [1,4]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
      uncal = mymap.pcolormesh(x,y,ecmf_bias[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('ECMWF Week '+str(n+1))

    plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.4)
    fig.subplots_adjust(right=0.90)
    cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]
    cbar_ax = fig.add_axes(cbar_pos)
    cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='both')
    plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
    #plt.show()
    plt.close()
    '''



    # plot RMSE
    cols = 'YlOrRd'
    cmin = 0
    cmax = 8
    cspc = 1
    clevs = np.arange(cmin,cmax+cspc,cspc)
    clabel1 = 'RMSE ($^{\circ}C$)'
    norm = BoundaryNorm(boundaries=clevs, ncolors=256)

    lw = 1
    nrow = 3
    ncol = nleads
    size1=[15,9]
    gl = 20
    fname_plot = dir_out+region+'_S2S_weekly_t2m_RMSE_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
    fig = plt.figure(figsize=(size1[0],size1[1]))
    for n in np.arange(0,nleads):
      plt.subplot(nrow,ncol,n+1)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      #mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
      tmp_bias = []
      tmp_bias = np.copy(ukmo_rmse[n,:,:])
      tmp_bias[tmp_bias == 0] = np.nan
      uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('UKMO Week '+str(n+1))

      plt.subplot(nrow,ncol,n+6)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      #mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
      tmp_bias = []
      tmp_bias = np.copy(ncep_rmse[n,:,:])
      tmp_bias[tmp_bias == 0] = np.nan
      uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('NCEP Week '+str(n+1))

      plt.subplot(nrow,ncol,n+11)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
      tmp_bias = []
      tmp_bias = np.copy(ecmf_rmse[n,:,:])
      tmp_bias[tmp_bias == 0] = np.nan
      uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('ECMWF Week '+str(n+1))

    plt.suptitle(season,fontweight='bold',fontsize=14,horizontalalignment='left')
    plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.4)
    fig.subplots_adjust(right=0.90)
    cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]
    cbar_ax = fig.add_axes(cbar_pos)
    cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='max')
    plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
    #plt.show()
    plt.close()



  # plot BSS
  cols = 'RdYlBu'
  clabel1 = 'BSS'
  lw = 1
  nrow = 3
  ncol = 5
  size1 = [15,9]
  gl = 20

  cmin = -0.5
  cmax = 0.5
  cspc = 0.1
  clevs = np.arange(cmin,cmax+cspc,cspc)
  norm = BoundaryNorm(boundaries=clevs,ncolors=256)

  for vv in [0,1]: # loop for own brier score
    if vv == 0:
      flag = 'comp_gpcc_clim'
    elif vv == 1:
      flag = 'comp_own_clim'
    for tc in np.arange(0,len(percentiles)+1):
      fname_plot = dir_out+region+'_S2S_weekly_t2m_BSS_'+p_name+'_'+flag+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'_'+p_name_ls[tc].replace(" ","")+'_v1'
      fig = plt.figure(figsize=(size1[0],size1[1]))

      for n in np.arange(0,nleads):
        plt.subplot(nrow,ncol,n+1)
        mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
        mymap.drawcoastlines(linewidth=lw)
        #mymap.drawcountries(linewidth=lw)
        mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
        mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
        if n in [0]:
          mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
        #mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
        x,y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
        mask_array = ukmo_bss[tc,vv,n,:,:]#np.ma.array(ukmo_bss[tc,vv,n,:,:]*dry_mask, mask=np.isnan(dry_mask))
        uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
        plt.title('UKMO Week '+str(n+1))

        plt.subplot(nrow,ncol,n+6)
        mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
        mymap.drawcoastlines(linewidth=lw)
        #mymap.drawcountries(linewidth=lw)
        mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
        mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
        if n in [0]:
          mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
        #mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
        x,y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
        mask_array = ncep_bss[tc,vv,n,:,:]#np.ma.array(ncep_bss[tc,vv,n,:,:]*dry_mask, mask=np.isnan(dry_mask))
        uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
        plt.title('NCEP Week '+str(n+1))

        plt.subplot(nrow,ncol,n+11)
        mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
        mymap.drawcoastlines(linewidth=lw)
        #mymap.drawcountries(linewidth=lw)
        mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
        mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
        if n in [0]:
          mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
        mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
        x,y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
        mask_array = ecmf_bss[tc,vv,n,:,:]#np.ma.array(ecmf_bss[tc,vv,n,:,:]*dry_mask, mask=np.isnan(dry_mask))
        uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
        plt.title('ECMWF Week '+str(n+1))

      plt.suptitle(season,fontweight='bold',fontsize=14,horizontalalignment='left')
      plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.4)
      fig.subplots_adjust(right=0.90)
      cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]

      cbar_ax = fig.add_axes(cbar_pos)
      cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='both')
      plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
    # plt.show()
      plt.close()



  # plot RPSS
  cols = 'YlGnBu'
  clabel1 = 'RPSS'
  lw = 1
  nrow = 3
  ncol = 5
  size1 = [15,9]
  gl = 20

  cmin = 0
  cmax = 0.5
  cspc = 0.05
  clevs = np.arange(cmin,cmax+cspc,cspc)
  norm = BoundaryNorm(boundaries=clevs,ncolors=256)

  for vv in [0,1]: # loop for own brier score
    if vv == 0:
      flag = 'comp_gpcc_clim'
    elif vv == 1:
      flag = 'comp_own_clim'
    fname_plot = dir_out+region+'_S2S_weekly_t2m_RPSS_'+p_name+'_'+flag+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
    fig = plt.figure(figsize=(size1[0],size1[1]))
    for n in np.arange(0,nleads):
      plt.subplot(nrow,ncol,n+1)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      #mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x,y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
      mask_array = rpss_ukmo[vv,n,:,:]#np.ma.array(rpss_ukmo[vv,n,:,:]*dry_mask, mask=np.isnan(dry_mask))
      uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('UKMO Week '+str(n+1))

      plt.subplot(nrow,ncol,n+6)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      #mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x,y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
      mask_array = rpss_ncep[vv,n,:,:]#np.ma.array(rpss_ncep[vv,n,:,:]*dry_mask, mask=np.isnan(dry_mask))
      uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('NCEP Week '+str(n+1))

      plt.subplot(nrow,ncol,n+11)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x,y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
      mask_array = rpss_ecmf[vv,n,:,:]#np.ma.array(rpss_ecmf[vv,n,:,:]*dry_mask, mask=np.isnan(dry_mask))
      uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('ECMWF Week '+str(n+1))

    plt.suptitle(season,fontweight='bold',fontsize=14,horizontalalignment='left')
    plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.4)
    fig.subplots_adjust(right=0.90)
    cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]

    cbar_ax = fig.add_axes(cbar_pos)
    cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='max')
    plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
  # plt.show()
    plt.close()
