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

def run_compute_currents(dates_to_process, oscar_mode, save_file, ssh_pattern, ssh_mode_list_val = None):
    print(ssh_mode_list_val)
    for date in dates_to_process:
        
        # SSH
        print(ssh_pattern)
        ssh_ds = load_ds(date, oscar_mode, var='ssh', ssh_pattern = ssh_pattern, ssh_mode_list_val = ssh_mode_list_val)
        ssh_ds = interpolate_dataset(ssh_ds)
        ssh_ds = calculate_gradient(ssh_ds, 'ssh')

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

    validation_mode = determine_validation_mode(dates)
    

    by_month = defaultdict(list)


    
    #2)Setting up the CMEMS and NEUROST directory depending on the Validation Mode: 'final' or 'interim'
        # 2a) loop through the list of dates and set up the folders for validation down the line. 
    for date in dates:
        year, month, day_str = date.split("-")
        if validation_mode =='final':
            cmems_dir = os.path.join(currents_cmems_final, year, month)
            cmems_pattern = f'{cmems_dir}/oscar_currents_{validation_mode}{year}{month}{day_str}.nc'
            neurost_dir = os.path.join(currents_neurost_final, year, month)
            neurost_pattern = f'{neurost_dir}/oscar_currents_{validation_mode}{year}{month}{day_str}.nc'
        elif validation_mode == "interim":
            cmems_dir = os.path.join(currents_cmems_interim, year, month)
            cmems_pattern = f'{cmems_dir}/oscar_currents_{validation_mode}{year}{month}{day_str}.nc'
            neurost_dir = os.path.join(currents_neurost_interim, year, month)
            cmems_pattern = f'{neurost_dir}/oscar_currents_{validation_mode}{year}{month}{day_str}.nc'
        else:
            raise ValueError ("please specific if you want to validate interim or final")
        
      
        drifter_dir = os.path.join(drifter_src_dir, year, month)
        drifter_pattern = f'drifter_6hour_qc_{year}_{month}.nc'
        cmems_ok   = bool(glob.glob(os.path.join(cmems_dir, cmems_pattern)))
        neurost_ok = bool(glob.glob(os.path.join(neurost_dir, neurost_pattern)))
        drifter_ok = bool(glob.glob(os.path.join(drifter_dir, drifter_pattern)))
        
        if cmems_ok and neurost_ok and drifter_ok:
            #2b) if the files exist to do the validation add the dates that belongs to the year month
            by_month[(year, month)].append(date)
        else: 
            if not cmems_ok: 
                raise FileNotFoundError(
                        f"CMEMS files not found for {date}: pattern '{cmems_pattern}' in '{cmems_dir}'"
                                )
            if not neurost_ok: 
                raise FileNotFoundError(
                            f"NEUROST files not found for {date}: pattern '{neurost_pattern}' in '{neurost_dir}'"
                                    )
            if not drifter_ok:
                raise FileNotFoundError(
                    f"DRIFTER files not found for {date}: pattern '{drifter_pattern}' in '{drifter_dir}'"
                        )

  

    # 3)We are looping by year_month, so grab the dates in each month 
    for (year, month), dates in by_month.items():

        # 3a) creating the directories we wiill need 
        month_root = os.path.join(drifter_src_dir, year, month)
        unique_drifter_month_dir = os.path.join(month_root, "drifter_by_id")
        daily_avg_dir = os.path.join(month_root, "daily_avg")
        os.makedirs(unique_drifter_month_dir, exist_ok=True)
        os.makedirs(daily_avg_dir, exist_ok=True)

        # 3b) opens the raw drifter file for that month so for the year 2020 
        # the file that will be opened is drifter_6hour_qc_2020_01.nc
        monthly_nc = get_drifter_monthly_file(drifter_src_dir, year, month)
      
        

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
        month_of_validation_data = []
        for d in dates:
            day_str = d[8:10]
            
            expected = os.path.join(daily_avg_dir, f"drifters_{str(d).replace("-", "_")}.nc")
            if not os.path.isfile(expected):
                # 3a) this function creates the daily_avg folder which has for each day, all the drifters, their lat, lon, ve, vn etc. 
                get_daily_avg_all_drifters(unique_drifter_month_dir, daily_avg_dir, d)
            
            #Now we have the daily averages to compare against OSCAR, however, the coordinates arent aligned, so we need to interpolate OSCAR to drifter coordinates. 
        

            #this block is her twice -- can improve to avaoid redundancy but right now a necessary redundancy

            if validation_mode =='final':
                cmems_dir = os.path.join(currents_cmems_final, year, month)
                neurost_dir = os.path.join(currents_neurost_final, year, month)
            elif validation_mode == "interim":
                cmems_dir = os.path.join(currents_cmems_interim, year, month)
                neurost_dir = os.path.join(currents_neurost_interim, year, month)
            else:
                raise ValueError ("please specific if you want to validate interim or final")
            
            ds_cmems = xr.open_dataset(f'{cmems_dir}/oscar_currents_{validation_mode}{year}{month}{day_str}.nc')
            ds_neurost = xr.open_dataset(f'{neurost_dir}/oscar_currents_{validation_mode}{year}{month}{day_str}.nc')
            ds_drifter = xr.open_dataset(f'{drifter_src_dir}/{year}/{month}/daily_avg/drifters_{year}_{month}_{day_str}.nc')


            ds_cmems = ds_cmems.rename({'longitude': 'lon', 'latitude': 'lat'})
            ds_cmems = ds_cmems.set_coords(['lat', 'lon'])
            u_cmems = ds_cmems['u'].transpose('time', 'lat', 'lon')
            u_cmems = u_cmems.assign_coords(lat=ds_cmems['lat'], lon=ds_cmems['lon'])


            
            ds_neurost = ds_neurost.rename({'longitude': 'lon', 'latitude': 'lat'})
            ds_neurost = ds_neurost.set_coords(['lat', 'lon'])
            u_neurost = ds_neurost['u'].transpose('time', 'lat', 'lon')
            u_neurost = u_neurost.assign_coords(lat=ds_neurost['lat'], lon=ds_neurost['lon'])

            # Interpolate CMEMS
            ds_drifter = interpolate_oscar_to_drifters(ds_cmems, ds_drifter, u_var='ve', v_var='vn', prefix='cmems')

     

            # # Interpolate NEUROST
            ds_drifter = interpolate_oscar_to_drifters(ds_neurost, ds_drifter, u_var='ve', v_var='vn', prefix='neurost')

    
            # we are adding each interpolated day to a list, so we should have a month of
            #drifter_d, drifter u, drifer v, cmems u, cmems v, neurost u, neurost v, lat lon, time
          

            month_of_validation_data.append(ds_drifter)
            
    


        validation_path = os.path.join(validation_dir, year, f'validation_{year}_{month}.nc')

        # we can implement override functionality here 
        if not os.path.exists(validation_path):
            if month_of_validation_data:
                #validation_ds = xr.open_dataset(validation_path)
                ds_monthly_validation = xr.concat(month_of_validation_data, dim='drifter')
               
                #validation_ds.close()
                

                os.makedirs(os.path.dirname(validation_path), exist_ok=True)
                ds_monthly_validation.to_netcdf(validation_path, mode = 'w')
                #now for the above, we have one month of validation data, yayy, time to plot
        #after the above block of code runs we now have one .nc file with a month of validation data
        #validation data is the interpplated, neurost and cmems 
        #now we need to appedn that to the master folder 
  
        #ACCOUNT FOR DUPLICATION
        
        append_to_master(year, month)

        month_of_validation_data.clear()

    master_path = os.path.join(validation_dir, 'validation_master.nc')
    ds = xr.open_dataset(master_path)





    df = ds[['latitude', 'longitude', 've', 'vn', 'cmems_ve', 'cmems_vn', 'neurost_ve', 'neurost_vn', 'time']].to_dataframe().dropna()
    df = df[
    (df['ve'] >= -2) & (df['ve'] <= 2) &
    (df['vn'] >= -2) & (df['vn'] <= 2)
    ]
    

    corr_maps = compute_binned_correlations(
        df,
        models=('cmems','neurost'),
        components=('ve','vn'),
        lat_range=(-60, 60),
        lon_range=(0, 360),
        bin_size=2,
        min_samples=5
    )



    temporal_range = summarize_time_jumps(df, "time")

    figs = []
    figs.append(plot_velocity_comparison_scatter(df, temporal_range))
    figs.append(plot_validation_metrics(df, temporal_range, models=['cmems', 'neurost'], components=['ve', 'vn']))
    figs.append(plot_binned_correlation_map(temporal_range,corr_maps[('cmems','ve')], 'cmems', 've', zoom_extent=None, title_prefix="Jan 2020"))
    figs.append(plot_binned_correlation_map(temporal_range,corr_maps[('neurost','ve')], 'neurost', 've', zoom_extent=None, title_prefix="Jan 2020"))
    figs.append(plot_binned_correlation_map(temporal_range,corr_maps[('cmems','vn')], 'cmems', 'vn', zoom_extent=None, title_prefix="Jan 2020"))
    figs.append(plot_binned_correlation_map(temporal_range,corr_maps[('neurost','vn')], 'neurost', 'vn', zoom_extent=None, title_prefix="Jan 2020"))
    figs.append(plot_density(df, temporal_range, 'cmems', 've', 'CMEMS', save_path=None))
    figs.append(plot_density(df, temporal_range, "cmems", "vn", "CMEMS", save_path = None))
    figs.append(plot_density(df, temporal_range, "neurost", "ve", "NEUROST", save_path = None))
    figs.append(plot_density(df, temporal_range, "neurost", "vn", "NEUROST", save_path = None))
    figs.append(make_slope_map(df, 've', 'cmems_ve', f'CMEMS Zonal Slope {temporal_range}' ))
    figs.append(make_slope_map(df, 'vn', 'cmems_vn', f'CMEMS Meridional Slope {temporal_range}'))
    figs.append(make_slope_map(df, 've', 'neurost_ve', f'NEUROST Zonal Slope  {temporal_range}'))
    figs.append(make_slope_map(df, 'vn', 'neurost_vn', f'NEUROST Meridional Slope {temporal_range}'))
    figs.append(make_rms_map(df, 'cmems',   've', f'CMEMS Zonal RMSD {temporal_range}'))
    figs.append(make_rms_map(df, 'cmems',   'vn', f'CMEMS Meridional RMSD {temporal_range}'))
    figs.append(make_rms_map(df, 'neurost', 've', f'NEUROST Zonal RMSD {temporal_range}'))
    figs.append(make_rms_map(df, 'neurost', 'vn', f'NEUROST Meridional RMSD {temporal_range}'))
    figs.append(make_residual_corr_map(df, 'cmems',   've', f'CMEMS Zonal Residual–Truth Corr {temporal_range}'))
    figs.append(make_residual_corr_map(df, 'cmems',   'vn', f'CMEMS Meridional Residual–Truth Corr {temporal_range}'))
    figs.append(make_residual_corr_map(df, 'neurost', 've', f'NEUROST Zonal Residual–Truth Cor {temporal_range}'))
    figs.append(make_residual_corr_map(df, 'neurost', 'vn', f'NEUROST Meridional Residual–Truth Corr {temporal_range}'))


    
    save_plots_to_pdf(
            figs,
            explanations,
            pdf_path="validation_report.pdf",
            metadata={"Title":"Validation Report","Author":"OSCAR"})
            
    
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

    if DO_DOWNLOAD:
        print("STARTING DOWNLOAD...")
        run_download(dates_to_process)

    
    if len(dates_to_process) != 0:
        if DO_CURRENTS:
            print("STARTING CURRENTS")
            ssh_mode_list = ['cmems', 'neurost']
            if SSH_MODE == 'both':
                
                for ssh in ssh_mode_list:
                    if ssh == 'cmems':
                        SSH_PATTERN = "ssh_*.nc"
                    elif ssh == 'neurost':
                        SSH_PATTERN = "NeurOST_SSH-SST_*.nc"
                    print(f"STARTING CURRENTS for {ssh}")
            
                    np.seterr(divide='ignore', invalid='ignore', over='ignore')
                    start_dt = np.datetime64(START_DATE)
                    end_dt = np.datetime64(END_DATE)
                    dates_to_process = get_dates_to_process(start_dt, end_dt, output_dir, override)
                    oscar_mode = determine_oscar_mode(dates_to_process, ssh)
                    save_file = f"oscar_currents_{oscar_mode}"
                    run_compute_currents(dates_to_process, oscar_mode, save_file, SSH_PATTERN, ssh)
                    run_plotting(dates_to_process, oscar_mode)
                    
            else:
                if SSH_MODE == 'cmems':
                    SSH_PATTERN = "ssh_*.nc"
                elif SSH_MODE == 'neurost':
                    SSH_PATTERN = "NeurOST_SSH-SST_*.nc"
                np.seterr(divide='ignore', invalid='ignore', over='ignore')
                start_dt = np.datetime64(START_DATE)
                end_dt = np.datetime64(END_DATE)
                dates_to_process = get_dates_to_process(start_dt, end_dt, output_dir, override)
                oscar_mode = determine_oscar_mode(dates_to_process, SSH_MODE)
                save_file = f"oscar_currents_{oscar_mode}"
                run_compute_currents(dates_to_process, oscar_mode, save_file, SSH_PATTERN)
                run_plotting(dates_to_process, oscar_mode)
        
    


    # if PLOT_CURRENTS:
    #     print("STARTING PLOTTING CURRENTS...")
    #     run_plotting(dates_to_process)

    if DO_VALIDATION:
        print("STARTING VALIDATION...")
        run_validation(dates_to_process)
      


        

if __name__ == "__main__":
    main()
