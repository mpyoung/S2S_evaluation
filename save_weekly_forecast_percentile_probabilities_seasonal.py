'''
Analyse skill of weekly S2S precipitation forecasts
from UKMO, NCEP and ECMWF over a given region

1) Reads in data saved over the given region from
'save_weekly_forecasts_ukmo_ncep.py'
'save_weekly_forecasts_ecmf.py'

2) Computes mean precipitation, bias, anomaly correlation coefficient
and brier skill scores

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

# define percentile categories to analyse
percentiles = [100./3.,200./3.]
p_name = 'tercile' # give categories an appropriate name

years	= np.arange(2000,2010+1,1)	# december will always correspond to year-1
nleads	= 5				# number of lead times (in weeks) in the data
nweeks = 20 # number of weeks in season
#season	= 'JJA'
seasons = ['SON']#['JJA']#['DJF','MAM','JJA','SON']


def mask_percentiles(data,percentile,arg):
  '''
  Function: 'mask_percentiles':
  Function takes in a 2d spatial forecast array and 2d spatial array of corresponding percentiles
  and masks the forecast array according to input arguments, 'arg':
  'between': mask the forecast between 2 percentiles
  'below': mask the forecast below a percentile
  'above': mask the forecast above a percentile
  '''
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

      ukmo_file.append(dir_in+region+'_UKMO_ensemble_weekly_hindcasts_lead1-5_'+str(y_ec)+'_'+starts_mon[mo]+'.nc')
      print ukmo_file[mo]+'	'+str(os.path.exists(ukmo_file[mo]))
      ncep_file.append(dir_in+region+'_NCEP_ensemble_weekly_hindcasts_lead1-5_'+str(y_ec)+'_'+starts_mon[mo]+'.nc')
      print ncep_file[mo]+'	'+str(os.path.exists(ncep_file[mo]))
      ecmf_file.append(dir_in+region+'_ECMWF_ensemble_weekly_hindcasts_lead1-5_'+str(y_ec)+'_'+starts_mon[mo]+'.nc')
      print ecmf_file[mo]+'	'+str(os.path.exists(ecmf_file[mo]))
      gpcp_file.append(dir_in+region+'_GPCP_weekly_'+str(y_ec)+'_'+starts_mon[mo]+'.nc')
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
      all_ukmo_ens = np.zeros((len(years),nleads,nweeks,ukmo_nmembers,len(ukmo_lat),len(ukmo_lon)))			# (11, 5, 20, 7,   , 40, 47)
    all_ukmo_ens[yc-1,:,:,:,:,:] = np.swapaxes(np.array(nc_fid.variables['week_precip'][:]),0,1)
    nc_fid.close()

    # NCEP units don't need converted. NCEP has 7 lags and 4 ensemble members
    nc_fid = MFDataset(ncep_file)
    if y == years[0]:
      ncep_lat = np.array(nc_fid.variables['latitude'][:])
      ncep_lon = np.array(nc_fid.variables['longitude'][:])
      all_ncep_ens = np.zeros((len(years),nleads,nweeks,ncep_lags,ncep_nmembers,len(ncep_lat),len(ncep_lon)))		# (11, 5, 20, 7,  4, 40, 47)
    all_ncep_ens[yc-1,:,:,:,:,:] = np.swapaxes(np.array(nc_fid.variables['week_precip'][:]),0,1)
    nc_fid.close()

    # ECMWF units don't need converted. There are 3 lags and 11 ensemble members
    nc_fid = MFDataset(ecmf_file)
    if y == years[0]:
      ecmf_lat = np.array(nc_fid.variables['latitude'][:])
      ecmf_lon = np.array(nc_fid.variables['longitude'][:])
      all_ecmf_ens = np.zeros((len(years),nleads,nweeks,ecmf_lags,ecmf_nmembers,len(ncep_lat),len(ncep_lon)))
    all_ecmf_ens[yc-1,:,:,:,:,:] = np.swapaxes(np.array(nc_fid.variables['week_precip'][:]),0,1)			# (11, 5, 20, 3, 11, 40, 47)
    nc_fid.close()

    # GPCP
    nc_fid = MFDataset(gpcp_file)
    if y == years[0]:
      gpcp_lat = np.array(nc_fid.variables['latitude'][:])
      gpcp_lon = np.array(nc_fid.variables['longitude'][:])
      all_gpcp = np.zeros((len(years),nleads,nweeks,len(gpcp_lat),len(gpcp_lon)))					# (11, 5, 20,  ,   , 40, 47)
    all_gpcp[yc-1,:,:,:,:] = np.swapaxes(np.array(nc_fid.variables['week_precip'][:]),0,1)
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

  # average forecasts over lags only (leaving ensembles)
  lag_ncep = np.nanmean(all_ncep_ens,axis=3) # (11, 5, 20, 4, 40, 47)
  lag_ecmf = np.nanmean(all_ecmf_ens,axis=3) # (11, 5, 20, 11, 40, 47)

  '''
  Compute forecast probabilities
  '''
  fname_p_gpcp = dir_in+'stats/probabilities/'+region+'_GPCP_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'
  fname_p_ukmo = dir_in+'stats/probabilities/'+region+'_UKMO_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'
  fname_p_ncep = dir_in+'stats/probabilities/'+region+'_NCEP_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'
  fname_p_ecmf = dir_in+'stats/probabilities/'+region+'_ECMWF_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'

  if os.path.exists(fname_p_ukmo) == False:

    p_gpcp = np.zeros((len(percentiles)+1,nleads,nweeks,len(years),len(ukmo_lat),len(ukmo_lon)))
    p_ukmo = np.zeros((len(percentiles)+1,2,nleads,nweeks,len(years),len(ukmo_lat),len(ukmo_lon)))
    p_ncep = np.zeros((len(percentiles)+1,2,nleads,nweeks,len(years),len(ukmo_lat),len(ukmo_lon)))
    p_ecmf = np.zeros((len(percentiles)+1,2,nleads,nweeks,len(years),len(ukmo_lat),len(ukmo_lon)))

    for l in np.arange(nleads):
      for y in np.arange(0,len(years)):
        nw = 0
        for w in np.arange(0,20):
          print 'week '+str(w+1)+' of 20'
          if l == 0:
            nw = nw + 1
          if week_mid_dates[y,l,w,1] in keep_mon:
            # compute the qth percentile of the data along the specified axis, while ignoring nan values.
            # here, it means that find 33% and 66% of gpcp weekly mean precipitation value at each grid cell.

            # for BSS analysis, we should be removing the current year being assessed
            # from the calculation of percentiles
            tmp_gpcp = []
            tmp_gpcp = np.delete(all_gpcp[:,l,w,:,:],y,axis=0)
            gpcp_percentiles = np.nanpercentile(tmp_gpcp,percentiles,axis=0)
            curr_gpcp = np.copy(all_gpcp[y,l,w,:,:])

            # 1. Compute observed probability (will be 1 or 0)
            for nn in np.arange(0,len(percentiles)+1):
              if nn == 0:
                p_gpcp[nn,l,w,y,:,:] = mask_percentiles(curr_gpcp,gpcp_percentiles[nn,:,:],'below')
              elif nn == len(percentiles):
                p_gpcp[nn,l,w,y,:,:] = mask_percentiles(curr_gpcp,gpcp_percentiles[nn-1,:,:],'above')
              else:
                p_gpcp[nn,l,w,y,:,:] = mask_percentiles(curr_gpcp,gpcp_percentiles[nn-1:nn+1,:,:],'between')

            # compute the probability that the ensemble was
            # above/below each tercile
            for e in np.arange(0,ecmf_nmembers):
              if e < 7: # 1) do UKMO
                #get percentiles for UKMO
                tmp_ukmo = []
                tmp_ukmo = np.delete(all_ukmo_ens[:,l,w,e,:,:],y,axis=0)
                ukmo_percentiles = np.nanpercentile(tmp_ukmo,percentiles,axis=0)# all_gpcp in (11, 5, 20, 40, 47); gpcp_percentiles in (2, 40, 47)
                curr_ukmo = np.copy(all_ukmo_ens[y,l,w,e,:,:].squeeze())

                for nn in np.arange(0,len(percentiles)+1):
                  if nn == 0:
                    p_ukmo[nn,0,l,w,y,:,:] = p_ukmo[nn,0,l,w,y,:,:] + mask_percentiles(curr_ukmo,gpcp_percentiles[nn,:,:],'below')
                    p_ukmo[nn,1,l,w,y,:,:] = p_ukmo[nn,1,l,w,y,:,:] + mask_percentiles(curr_ukmo,ukmo_percentiles[nn,:,:],'below')
                  elif nn == len(percentiles):
                    p_ukmo[nn,0,l,w,y,:,:] = p_ukmo[nn,0,l,w,y,:,:] + mask_percentiles(curr_ukmo,gpcp_percentiles[nn-1,:,:],'above')
                    p_ukmo[nn,1,l,w,y,:,:] = p_ukmo[nn,1,l,w,y,:,:] + mask_percentiles(curr_ukmo,ukmo_percentiles[nn-1,:,:],'above')
                  else:
                    p_ukmo[nn,0,l,w,y,:,:] = p_ukmo[nn,0,l,w,y,:,:] + mask_percentiles(curr_ukmo,gpcp_percentiles[nn-1:nn+1,:,:],'between')
                    p_ukmo[nn,1,l,w,y,:,:] = p_ukmo[nn,1,l,w,y,:,:] + mask_percentiles(curr_ukmo,ukmo_percentiles[nn-1:nn+1,:,:],'between')

              if e < 4:
                #get percentiles for NCEP
                tmp_ncep = []
                tmp_ncep = np.delete(lag_ncep[:,l,w,e,:,:],y,axis=0)
                ncep_percentiles = np.nanpercentile(tmp_ncep,percentiles,axis=0)
                curr_ncep = np.copy(lag_ncep[y,l,w,e,:,:])
                for nn in np.arange(0,len(percentiles)+1):
                  if nn == 0:
                    p_ncep[nn,0,l,w,y,:,:] = p_ncep[nn,0,l,w,y,:,:] + mask_percentiles(curr_ncep,gpcp_percentiles[nn,:,:],'below')
                    p_ncep[nn,1,l,w,y,:,:] = p_ncep[nn,1,l,w,y,:,:] + mask_percentiles(curr_ncep,ncep_percentiles[nn,:,:],'below')
                  elif nn == len(percentiles):
                    p_ncep[nn,0,l,w,y,:,:] = p_ncep[nn,0,l,w,y,:,:] + mask_percentiles(curr_ncep,gpcp_percentiles[nn-1,:,:],'above')
                    p_ncep[nn,1,l,w,y,:,:] = p_ncep[nn,1,l,w,y,:,:] + mask_percentiles(curr_ncep,ncep_percentiles[nn-1,:,:],'above')
                  else:
                    p_ncep[nn,0,l,w,y,:,:] = p_ncep[nn,0,l,w,y,:,:] + mask_percentiles(curr_ncep,gpcp_percentiles[nn-1:nn+1,:,:],'between')
                    p_ncep[nn,1,l,w,y,:,:] = p_ncep[nn,1,l,w,y,:,:] + mask_percentiles(curr_ncep,ncep_percentiles[nn-1:nn+1,:,:],'between')

              #get percentiles for ECMWF
              tmp_ecmf = []
              tmp_ecmf = np.delete(lag_ecmf[:,l,w,e,:,:],y,axis=0)
              ecmf_percentiles = np.nanpercentile(tmp_ecmf,percentiles,axis=0)
              curr_ecmf = np.copy(lag_ecmf[y,l,w,e,:,:])
              for nn in np.arange(0,len(percentiles)+1):
                if nn == 0:
                  p_ecmf[nn,0,l,w,y,:,:] = p_ecmf[nn,0,l,w,y,:,:] + mask_percentiles(curr_ecmf,gpcp_percentiles[nn,:,:],'below')
                  p_ecmf[nn,1,l,w,y,:,:] = p_ecmf[nn,1,l,w,y,:,:] + mask_percentiles(curr_ecmf,ecmf_percentiles[nn,:,:],'below')
                elif nn == len(percentiles):
                  p_ecmf[nn,0,l,w,y,:,:] = p_ecmf[nn,0,l,w,y,:,:] + mask_percentiles(curr_ecmf,gpcp_percentiles[nn-1,:,:],'above')
                  p_ecmf[nn,1,l,w,y,:,:] = p_ecmf[nn,1,l,w,y,:,:] + mask_percentiles(curr_ecmf,ecmf_percentiles[nn-1,:,:],'above')
                else:
                  p_ecmf[nn,0,l,w,y,:,:] = p_ecmf[nn,0,l,w,y,:,:] + mask_percentiles(curr_ecmf,gpcp_percentiles[nn-1:nn+1,:,:],'between')
                  p_ecmf[nn,1,l,w,y,:,:] = p_ecmf[nn,1,l,w,y,:,:] + mask_percentiles(curr_ecmf,ecmf_percentiles[nn-1:nn+1,:,:],'between')


    # At end compute actual Probability
    pf_ukmo = p_ukmo/ukmo_nmembers # probability
    pf_ncep = p_ncep/ncep_nmembers # probability, normalise by N ensemble members
    pf_ecmf = p_ecmf/ecmf_nmembers # probability

    # save probabilities
    np.save(fname_p_gpcp,p_gpcp)
    np.save(fname_p_ukmo,pf_ukmo)
    np.save(fname_p_ncep,pf_ncep)
    np.save(fname_p_ecmf,pf_ecmf)
