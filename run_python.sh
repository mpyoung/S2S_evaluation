#nohup ipython save_download_gsod_isd.py > outfile_download_gsod_africa.txt < /dev/null &
# nohup ipython get_ecmwf_hindcasts.py > outfile_download_ecmwf.txt < /dev/null &
# nohup ipython matts_organise_and_convert.py > outfile_convert_ecmwf.txt < /dev/null &
# nohup ipython save_ncep_ltime_week_ens.py > outfile_save_ncep.txt < /dev/null &

#nohup python2.7 save_weekly_forecasts_ukmo_ncep.py > outfile_save_Africa_monthly_UKMO_NCEP.txt < /dev/null &
# nohup python2.7 save_weekly_forecasts_ukmo_ncep_temperature.py > outfile_save_Africa_monthly_UKMO_NCEP_temp.txt < /dev/null &
#nohup python2.7 save_weekly_forecasts_ecmf_temperature.py > outfile_save_Africa_monthly_ecmf_temp.txt < /dev/null &
# nohup python2.7 analyse_weekly_forecasts_seasonal.py > outfile_analyse_Africa_seasonal_SON.txt < /dev/null &
#nohup python2.7 save_weekly_forecast_percentile_probabilities_seasonal.py > outfile_save_Africa_p_probs_SON.txt < /dev/null &
nohup python2.7 save_weekly_forecast_percentile_probabilities_temperature_seasonal.py > outfile_save_Africa_t_probs_JJA.txt < /dev/null &
