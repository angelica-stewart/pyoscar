import numpy as np
from .computation.interp_and_grad import *
from .computation.compute_currents import *
from .utils.data_utils import *
from .config.setup import *  # contains: compute_currents, plot_currents, validate_currents, start_date, end_date, etc.
from .analysis.plot_currents import *
from .analysis.drifter_func import *
import pandas as pd
from .utils.download import *
import warnings
from collections import defaultdict
from .analysis.prepare_pdf import *
import os, glob

warnings.filterwarnings("ignore")

def run_compute_currents(dates_to_process, oscar_mode, save_file):
    for date in dates_to_process:
        
        # SSH
       
        ssh_ds = load_ds(date, oscar_mode, var='ssh')
        ssh_ds = interpolate_dataset(ssh_ds)
        ssh_ds = calculate_gradient(ssh_ds, 'ssh')
        #plot_interp_and_grad(ssh_ds, 'test', 'inter.png', '/Users/stewarta/Desktop/oscarpy/','sshy', 'ssh')

        # SST
        sst_ds = load_ds(date, oscar_mode, var='sst')
        sst_ds = interpolate_dataset(sst_ds)
        sst_ds = calculate_gradient(sst_ds, 'sst')

        # Wind
        wind_ds = load_ds(date, oscar_mode, var='wind')
        wind_ds = interpolate_dataset(wind_ds)

        # Compute and write OSCAR
        ref_ds = ssh_ds.drop_vars(['sshx', 'sshy'])
        print(f"The currents computation for {date} is starting now...")
        Ug, Uw, Ub = compute_surface_currents(ssh_ds, wind_ds, sst_ds, do_eq)

        write_oscar_to_podaac(ref_ds, Ug, Uw, Ub, save_file, podaacdir,
                              ssh_long_desc, wind_long_desc, sst_long_desc,
                              oscar_long_desc, oscar_summary, oscar_id,
                              doi, SSH_MODE)
        print(f"Finished computing currents for {date}")
    print(f"ALL CURRENTS COMPUTATION FOR {dates_to_process[0]} - {dates_to_process[-1]} IS COMPLETED AND SAVED FOR PODAAC")

def run_plotting(dates, oscar_mode):

    for date in dates:
        ds_path= get_file_path(date, oscar_mode)
        title = get_title(date)
        ds = xr.open_dataset(ds_path)
        u, v, lon, lat, speed, lon2d, lat2d = extract_data(ds)
        plot_currents(ds_path, lon2d, lat2d, u, v, speed, title, SSH_MODE)

def run_validation(dates):
    #1) create list of dates 

    validation_mode_cmems = determine_validation_mode(dates, 'cmems')
    validation_mode_neurost = determine_validation_mode(dates, 'neurost')

   #if the paths dont exist at all and the DATADIR/CURRENTS/validation_mode_cmems.upper()/first[0:4]/first[4:6]/oscar_currents_{validation_mode_cmems}{first[0].replace('-','').nc
    master_path_cmems_pdf = os.path.join(VALIDATION_DIR, f'validation_report_cmems_{validation_mode_cmems}*.pdf')
    master_path_neurost_pdf = os.path.join(VALIDATION_DIR, f'validation_report_neurost_{validation_mode_neurost}*.pdf')
     
 

    #2) create the paths for the CMEMS and NEUROST currents depending on the validation mode

    if validation_mode_cmems != None:
        ok_cmems_date = is_master_consecutive(master_path_cmems_pdf, dates)
        
    else:
        ok_cmems_date = False
 
   

    if validation_mode_neurost != None:
        ok_neurost_date = is_master_consecutive(master_path_neurost_pdf, dates)
    else:
        ok_neurost_date = False


    if not ok_cmems_date and  not ok_neurost_date:
        raise ValueError (' None of the dates for the OSCAR inputs are consecutive, the validation is designed to only validate consecutive files')

    by_month = defaultdict(list)


    
    #2)Setting up the CMEMS and NEUROST directory depending on the Validation Mode: 'final' or 'interim'
        # 2a) loop through the list of dates and set up the folders for validation down the line. 
    for date in dates:
        year, month, day_str = date.split("-")

        if ok_cmems_date:
            if validation_mode_cmems =='final':
                cmems_dir = os.path.join(currents_cmems_final, year, month)
                cmems_pattern = f'{cmems_dir}/oscar_currents_{validation_mode_cmems}{year}{month}{day_str}.nc'
            if validation_mode_cmems == "interim":
                cmems_dir = os.path.join(currents_cmems_interim, year, month)
                cmems_pattern = f'{cmems_dir}/oscar_currents_{validation_mode_cmems}{year}{month}{day_str}.nc'
            cmems_ok   = bool(glob.glob(os.path.join(cmems_dir, cmems_pattern)))
        
        if ok_neurost_date:
            if validation_mode_neurost == 'final':
                neurost_dir = os.path.join(currents_neurost_final, year, month)
                neurost_pattern = f'{neurost_dir}/oscar_currents_{validation_mode_neurost}{year}{month}{day_str}.nc'
            if validation_mode_neurost == 'interim':
                neurost_dir = os.path.join(currents_neurost_interim, year, month)
                neurost_pattern = f'{neurost_dir}/oscar_currents_{validation_mode_neurost}{year}{month}{day_str}.nc'
            neurost_ok = bool(glob.glob(os.path.join(neurost_dir, neurost_pattern)))
 
            
    

        
      
        drifter_dir = os.path.join(DRIFTER_SRC_DIR, year, month)
        drifter_pattern = f'drifter_6hour_qc_{year}_{month}.nc'
        drifter_ok = bool(glob.glob(os.path.join(drifter_dir, drifter_pattern)))


        if ok_cmems_date:
            if drifter_ok and cmems_ok:

                by_month[(year, month)].append(date)
            else: 
                print(f"CMEMS files not found for {date}: pattern {cmems_pattern} in {cmems_dir}")
        
        if ok_neurost_date:
            if drifter_ok and neurost_ok:
                by_month[(year, month)].append(date)
            else:
                print(f"NEUROST files not found for {date}: pattern {neurost_pattern} in '{neurost_dir}")
    
    
        if not drifter_ok:
            raise FileNotFoundError(
                f"DRIFTER files not found for {date}: pattern '{drifter_pattern}' in '{drifter_dir}'"
                        )
            
          
  

    # 3)We are looping by year_month, so grab the dates in each month 
    for (year, month), dates in by_month.items():

        # 3a) creating the directories we wiill need 
        month_root = os.path.join(DRIFTER_SRC_DIR, year, month)
        unique_drifter_month_dir = os.path.join(month_root, "drifter_by_id")
        daily_avg_dir = os.path.join(month_root, "daily_avg")
        os.makedirs(unique_drifter_month_dir, exist_ok=True)
        os.makedirs(daily_avg_dir, exist_ok=True)

        # 3b) opens the raw drifter file for that month so for the year 2020 
        # the file that will be opened is drifter_6hour_qc_2020_01.nc
        monthly_nc = get_drifter_monthly_file(DRIFTER_SRC_DIR, year, month)
      
        

        # 3c) check if the drifter_by_id folder for that month is empty, if it has files in it, 
              #we are assuminng the monthly average is already computed an no need to recompute 
              #(we can out a flag here to override). if it is empty the daily average is not computed so we do 3d)
        if not os.listdir(unique_drifter_month_dir):
            sentinel = os.path.join(unique_drifter_month_dir, ".monthly_done")
            if not os.path.exists(sentinel):
                with xr.open_dataset(monthly_nc) as ds_month:

                    #3d) we create the drifter_by_id folder, the get_uniq_drifter_daily_avg_folder creates a dataset 
                    # each representing one drifter and ech drifter has the daily average for that day. 
                    get_unique_drifter_daily_avg(ds_month, unique_drifter_month_dir)
                    
                open(sentinel, "w").close()

        

        # now that we have the daily averages for the month for each drifter, 
        # we want daily files of all the drifter so we can match it to OSCAR. 
        # note to coder: we are still in the year_month loop 

        # 4) For each day in that month, extract the daily data from  each drifter in the drifter_by_id folder
            #we are going to save each validation file by month. 
        #month_of_validation_data = []
        drifter_vs_neurost_list = []
        drifter_vs_cmems_list = []
        for d in dates:
            day_str = d[8:10]
    
            expected = os.path.join(daily_avg_dir, f"drifters_{str(d).replace('-', '_')}.nc")
            if not os.path.isfile(expected):
                # 3a) this function creates the daily_avg folder which has for each day, all the drifters, their lat, lon, ve, vn etc. 
                get_daily_avg_all_drifters(unique_drifter_month_dir, daily_avg_dir, d)
            
            #Now we have the daily averages to compare against OSCAR, however, the coordinates arent aligned, so we need to interpolate OSCAR to drifter coordinates. 
        

            #this block is her twice -- can improve to avaoid redundancy but right now a necessary redundancy

            if ok_cmems_date:
                if validation_mode_cmems =='final':
                    cmems_dir = os.path.join(currents_cmems_final, year, month)
                elif validation_mode_cmems == 'interim':
                    cmems_dir = os.path.join(currents_cmems_interim, year, month)
               
            if ok_neurost_date:
                if validation_mode_neurost == 'final':
                    neurost_dir = os.path.join(currents_neurost_final, year, month)
                elif validation_mode_neurost == "interim":
                    neurost_dir = os.path.join(currents_neurost_interim, year, month)
    
            
        
            ds_drifter = xr.open_dataset(f'{DRIFTER_SRC_DIR}/{year}/{month}/daily_avg/drifters_{year}_{month}_{day_str}.nc')
            print(ds_drifter.time.max())

            if ok_cmems_date:
                ds_cmems = xr.open_dataset(f'{cmems_dir}/oscar_currents_{validation_mode_cmems}{year}{month}{day_str}.nc')
                ds_cmems = ds_cmems.rename({'longitude': 'lon', 'latitude': 'lat'})
                ds_cmems = ds_cmems.set_coords(['lat', 'lon'])
                u_cmems = ds_cmems['u'].transpose('time', 'lat', 'lon')
                u_cmems = u_cmems.assign_coords(lat=ds_cmems['lat'], lon=ds_cmems['lon'])
                drifter_vs_cmems_ds = interpolate_oscar_to_drifters(ds_cmems, ds_drifter, u_var='ve', v_var='vn', prefix='cmems')
                drifter_vs_cmems_list.append(drifter_vs_cmems_ds)

            if ok_neurost_date:
                ds_neurost = xr.open_dataset(f'{neurost_dir}/oscar_currents_{validation_mode_neurost}{year}{month}{day_str}.nc')
                ds_neurost = ds_neurost.rename({'longitude': 'lon', 'latitude': 'lat'})
                ds_neurost = ds_neurost.set_coords(['lat', 'lon'])
                u_neurost = ds_neurost['u'].transpose('time', 'lat', 'lon')
                u_neurost = u_neurost.assign_coords(lat=ds_neurost['lat'], lon=ds_neurost['lon'])
                drifter_vs_neurost_ds = interpolate_oscar_to_drifters(ds_neurost, ds_drifter, u_var='ve', v_var='vn', prefix='neurost')
                print(drifter_vs_neurost_ds.time.max())
                drifter_vs_neurost_list.append(drifter_vs_neurost_ds)
            
    


            if ok_cmems_date:
                validation_path_cmems = os.path.join(VALIDATION_DIR, 'CMEMS', validation_mode_cmems, year, f'validation_{year}_{month}.nc')
                if not os.path.exists(validation_path_cmems):
                    if drifter_vs_cmems_list:
                        #validation_ds = xr.open_dataset(validation_path)
                        drifter_vs_cmems_month = xr.concat(drifter_vs_cmems_list, dim='drifter')
                        print(drifter_vs_cmems_month.time.max())
                        os.makedirs(os.path.dirname(validation_path_cmems), exist_ok=True)
                        drifter_vs_cmems_month.to_netcdf(validation_path_cmems, mode = 'w')


                #append_to_master(year, month, 'cmems', validation_mode_cmems)
                drifter_vs_cmems_list.clear()

                master_path_cmems = os.path.join(VALIDATION_DIR, f'validation_cmems_{validation_mode_cmems}_master_*.nc')
                cmems_file = glob.glob(master_path_cmems)
                if len(cmems_file) == 0:
                    df_cmems_start_date = dates[0].replace('-', '')
                    df_cmems_end_date = dates[-1].replace('-', '')
                    append_to_master(year, month, df_cmems_start_date, df_cmems_end_date, 'cmems', validation_mode_cmems)
                    new_cmems_file = os.path.join(VALIDATION_DIR, f'validation_cmems_{validation_mode_cmems}_master_{df_cmems_start_date}_{df_cmems_end_date}.nc')
                    ds_master_cmems = xr.open_dataset(new_cmems_file)
                    df_cmems = clean_validation_df(ds_master_cmems, 'cmems')
                
                else:
                    ds_master_cmems = xr.open_dataset(cmems_file[0])
                    df_cmems = clean_validation_df(ds_master_cmems, 'cmems')
                    df_cmems_start_date = pd.to_datetime(df_cmems['time'].iloc[0]).strftime("%Y%m%d")
                    df_cmems_end_date   = pd.to_datetime(df_cmems['time'].iloc[-1]).strftime("%Y%m%d")
                    append_to_master(year, month, df_cmems_start_date, df_cmems_end_date, 'cmems', validation_mode_cmems)

                corr_map_cmems = compute_binned_correlations(
                df_cmems,
                models=['cmems'],
                components=('ve','vn'),
                lat_range=(-60, 60),
                lon_range=(0, 360),
                bin_size=2,
                min_samples=5)
                
            if ok_neurost_date:
        
                validation_path_neurost = os.path.join(VALIDATION_DIR, 'NEUROST', validation_mode_neurost, year, f'validation_{year}_{month}.nc')
            

                print(validation_path_neurost)
                if not os.path.exists(validation_path_neurost):
                    if drifter_vs_neurost_list:
                        print(drifter_vs_neurost_list)
                        drifter_vs_neurost_month = xr.concat(drifter_vs_neurost_list, dim='drifter')
                        print(drifter_vs_neurost_month.time.min())
                        print(drifter_vs_neurost_month.time.max())
                        os.makedirs(os.path.dirname(validation_path_neurost), exist_ok=True)
                        drifter_vs_neurost_month.to_netcdf(validation_path_neurost, mode = 'w')
                if os.path.exists(validation_path_neurost):
                    validation_ds = xr.open_dataset(validation_path_neurost)
                    drifter_vs_cmems_list.insert(0, validation_ds)
                    drifter_vs_neurost_month = xr.concat(drifter_vs_neurost_list, dim='drifter')
                    print(drifter_vs_neurost_month.time.min())
                    print(drifter_vs_neurost_month.time.max())
                    validation_ds.close()
                    os.makedirs(os.path.dirname(validation_path_neurost), exist_ok=True)
                    drifter_vs_neurost_month.to_netcdf(validation_path_neurost, mode = 'w')

                drifter_vs_neurost_list.clear()
                master_path_neurost = os.path.join(VALIDATION_DIR, f'validation_neurost_{validation_mode_neurost}_master_*.nc')
                neurost_file = glob.glob(master_path_neurost)
                if len(neurost_file) == 0:
                    df_neurost_start_date = dates[0].replace('-', '')
                    df_neurost_end_date   = dates[-1].replace('-', '')
                    append_to_master(year, month, df_neurost_start_date, df_neurost_end_date,
                                    'neurost', validation_mode_neurost)
                    new_neurost_file = os.path.join(VALIDATION_DIR, f'validation_neurost_{validation_mode_neurost}_master_{df_neurost_start_date}_{df_neurost_end_date}.nc')
                    ds_master_neurost = xr.open_dataset(new_neurost_file)
                    df_neurost = clean_validation_df(ds_master_neurost, 'neurost')

                else:
                    ds_master_neurost = xr.open_dataset(neurost_file[0])
                    print(ds_master_neurost.time.max())
                    df_neurost = clean_validation_df(ds_master_neurost, 'neurost')
                    df_neurost_start_date = pd.to_datetime(df_neurost['time'].iloc[0]).strftime("%Y%m%d")
                    df_neurost_end_date   = pd.to_datetime(df_neurost['time'].iloc[-1]).strftime("%Y%m%d")
                    append_to_master(year, month, df_neurost_start_date, df_neurost_end_date,
                                    'neurost', validation_mode_neurost)



                corr_map_neurost = compute_binned_correlations(
                    df_neurost,
                    models=['neurost'],
                    components=('ve','vn'),
                    lat_range=(-60, 60),
                    lon_range=(0, 360),
                    bin_size=2,
                    min_samples=5
                )
            
            if ok_cmems_date:
                temporal_range_cmems = summarize_time_jumps(df_cmems, "time")

                figs_cmems = []

                figs_cmems.append(plot_velocity_comparison_scatter(df_cmems, temporal_range_cmems, 'cmems'))
                figs_cmems.append(plot_binned_correlation_map(temporal_range_cmems, corr_map_cmems[('cmems','ve')], 'cmems', 've', zoom_extent=None, title_prefix="Jan 2020"))
                figs_cmems.append(plot_binned_correlation_map(temporal_range_cmems, corr_map_cmems[('cmems','vn')], 'cmems', 'vn', zoom_extent=None, title_prefix="Jan 2020"))
                figs_cmems.append(plot_density(df_cmems, temporal_range_cmems, 'cmems', 've', 'CMEMS', save_path=None))
                figs_cmems.append(plot_density(df_cmems, temporal_range_cmems, 'cmems', 'vn', 'CMEMS', save_path=None))
                figs_cmems.append(make_slope_map(df_cmems, 've', 'cmems_ve', f'CMEMS Zonal Slope {temporal_range_cmems}' ))
                figs_cmems.append(make_slope_map(df_cmems, 'vn', 'cmems_vn', f'CMEMS Meridional Slope {temporal_range_cmems}'))
                figs_cmems.append(make_rms_map(df_cmems, 'cmems', 've', f'CMEMS Zonal RMSD {temporal_range_cmems}'))
                figs_cmems.append(make_rms_map(df_cmems, 'cmems', 'vn', f'CMEMS Meridional RMSD {temporal_range_cmems}'))
                figs_cmems.append(make_residual_corr_map(df_cmems, 'cmems', 've', f'CMEMS Zonal Residual–Truth Corr {temporal_range_cmems}'))
                figs_cmems.append(make_residual_corr_map(df_cmems, 'cmems', 'vn', f'CMEMS Meridional Residual–Truth Corr {temporal_range_cmems}'))
                save_plots_to_pdf(figs_cmems, explanations,
                                pdf_path=VALIDATION_DIR + f"/validation_report_cmems_{validation_mode_cmems}{df_cmems_start_date}_{df_cmems_end_date}.pdf",
                                metadata={"Title":"Validation Report","Author":"OSCAR"})

            

            if ok_neurost_date:
                temporal_range_neurost = summarize_time_jumps(df_neurost, "time")

                figs_neurost = []
                figs_neurost.append(plot_velocity_comparison_scatter(df_neurost, temporal_range_neurost, 'neurost'))
                figs_neurost.append(plot_binned_correlation_map(temporal_range_neurost, corr_map_neurost[('neurost','ve')], 'neurost', 've', zoom_extent=None, title_prefix="Jan 2020"))
                figs_neurost.append(plot_binned_correlation_map(temporal_range_neurost, corr_map_neurost[('neurost','vn')], 'neurost', 'vn', zoom_extent=None, title_prefix="Jan 2020"))
                figs_neurost.append(plot_density(df_neurost, temporal_range_neurost, "neurost", "ve", "NEUROST", save_path=None))
                figs_neurost.append(plot_density(df_neurost, temporal_range_neurost, "neurost", "vn", "NEUROST", save_path=None))
                figs_neurost.append(make_slope_map(df_neurost, 've', 'neurost_ve', f'NEUROST Zonal Slope  {temporal_range_neurost}'))
                figs_neurost.append(make_slope_map(df_neurost, 'vn', 'neurost_vn', f'NEUROST Meridional Slope {temporal_range_neurost}'))
                figs_neurost.append(make_rms_map(df_neurost, 'neurost', 've', f'NEUROST Zonal RMSD {temporal_range_neurost}'))
                figs_neurost.append(make_rms_map(df_neurost, 'neurost', 'vn', f'NEUROST Meridional RMSD {temporal_range_neurost}'))
                figs_neurost.append(make_residual_corr_map(df_neurost, 'neurost', 've', f'NEUROST Zonal Residual–Truth Cor {temporal_range_neurost}'))
                figs_neurost.append(make_residual_corr_map(df_neurost, 'neurost', 'vn', f'NEUROST Meridional Residual–Truth Corr {temporal_range_neurost}'))
                save_plots_to_pdf(figs_neurost, explanations,
                                pdf_path=VALIDATION_DIR + f"/validation_report_neurost_{validation_mode_neurost}{df_neurost_start_date}_{df_neurost_end_date}.pdf",
                                metadata={"Title":"Validation Report NEUROST","Author":"OSCAR"})
        
                
        
        print('print out paths')
    return

def run_download(dates):

        
    download_ssh_cmems(dates)

    download_ssh_neurost(dates)

    download_wind_era5(dates)

    download_sst_cmc(dates)

    download_drifter_data(dates)

    return

def main():
    np.seterr(divide='ignore', invalid='ignore', over='ignore')
    start_dt = np.datetime64(START_DATE)
    end_dt = np.datetime64(END_DATE)
    dates_to_process = get_dates_to_process(start_dt, end_dt, output_dir, override)

    # if DO_DOWNLOAD:
    print("STARTING DOWNLOAD...")

    if not OVERRIDE_DOWNLOAD:
        run_download(dates_to_process)

    
    if len(dates_to_process) != 0:
        print("STARTING CURRENTS")  
        np.seterr(divide='ignore', invalid='ignore', over='ignore')
        start_dt = np.datetime64(START_DATE)
        end_dt = np.datetime64(END_DATE)
        dates_to_process = get_dates_to_process(start_dt, end_dt, output_dir, override)
        oscar_mode = determine_oscar_mode(dates_to_process, SSH_MODE)
        save_file = f"oscar_currents_{oscar_mode}"
        run_compute_currents(dates_to_process, oscar_mode, save_file)
   

    


        if PLOT_CURRENTS:
            print("STARTING PLOTTING CURRENTS...")
            run_plotting(dates_to_process, oscar_mode)

        if DO_VALIDATION:
            print("STARTING VALIDATION...")
            run_validation(dates_to_process)
      


        

if __name__ == "__main__":
    main()
