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

warnings.filterwarnings("ignore")

def run_compute_currents(dates_to_process):
    for date in dates_to_process:
        # SSH
        ssh_ds = load_ds(date, var=ssh_var)
        ssh_ds = interpolate_dataset(ssh_ds)
        ssh_ds = calculate_gradient(ssh_ds, ssh_var)

        # SST
        sst_ds = load_ds(date, var=sst_var)
        sst_ds = interpolate_dataset(sst_ds)
        sst_ds = calculate_gradient(sst_ds, sst_var)

        # Wind
        wind_ds = load_ds(date, var='wind')
        wind_ds = interpolate_dataset(wind_ds)

        # Compute and write OSCAR
        ref_ds = ssh_ds.drop_vars(['sshx', 'sshy'])
        print(f"The currents computation for {date} is starting now")
        Ug, Uw, Ub = compute_surface_currents(ssh_ds, wind_ds, sst_ds, do_eq)

        write_oscar_to_podaac(ref_ds, Ug, Uw, Ub, podaacfile, podaacdir,
                              ssh_long_desc, wind_long_desc, sst_long_desc,
                              oscar_long_desc, oscar_summary, oscar_id,
                              doi, ssh_mode)
        print("Currents computation done\n")

def run_plotting(dates):
    if len(dates) == 1:
        date_list = [dates[0]]
    elif len(dates) == 2:
        date_list = get_date_range(dates)
    else:
        raise ValueError("date list must contain 1 or 2 dates (start and end).")

    for date in date_list:
        ds_path= get_file_path(date)
        print(ds_path)
        title = get_title(date)
        ds = xr.open_dataset(ds_path)
        u, v, lon, lat, speed, lon2d, lat2d = extract_data(ds)
        plot_currents(ds_path, lon2d, lat2d, u, v, speed, title, ssh_mode_plots)

def run_validation(dates):
    if len(dates) == 1:
        date_list = [dates[0]]
    elif len(dates) == 2:
        date_list = get_date_range(dates)
    else:
        raise ValueError("date list must contain 1 or 2 dates (start and end).")

    by_month = defaultdict(list)
    for date in date_list:
        year, month, day= date.split("-")
        by_month[(year, month)].append(date)
  

    
    for (year, month), dates in by_month.items():
        month_root = os.path.join(drifter_src_dir, year, month)
        unique_drifter_month_dir = os.path.join(month_root, "drifter_by_id")
        daily_avg_dir = os.path.join(month_root, "daily_avg")
        os.makedirs(unique_drifter_month_dir, exist_ok=True)
        os.makedirs(daily_avg_dir, exist_ok=True)

        # 1) opens the raw drifter file for that month so for the year 2020 the first file that will be opened is drifter_6hour_qc_2020_01.nc
        monthly_nc = get_drifter_monthly_file(drifter_src_dir, year, month)
        

        # 2) check if the drifter_by_id folder for that month is empty, if it has files in it, 
        #    we are assuminng the monthly average is already computed an no need to recompute (we can out a flag here to override).
        #    if it is empty then e do 2a)
        if not os.listdir(unique_drifter_month_dir):
            sentinel = os.path.join(unique_drifter_month_dir, ".monthly_done")
            if not os.path.exists(sentinel):
                with xr.open_dataset(monthly_nc) as ds_month:

                    #2a) we create the drifter_by_id folder, the get_uniq_drifter_daily_avg_folder returns a list of dataset, 
                    # each representing one drifter and ech drifter has the daily average for that day. 
                    get_unique_drifter_daily_avg(ds_month, unique_drifter_month_dir)
                open(sentinel, "w").close()

        #now that we have the daily averages for the month for each drifter, we want daily files of all the drifter so we can match it to OSCAR. 

        # 3) For each day in that month, extract the daily data from  each drifter in the drifter_by_id folder
        month_of_validation_data = []
        for d in dates:
            day_str = d[8:10]
            
            if not os.listdir(daily_avg_dir):
                # 3a) this function creates the daily_avg folder which has for each day, all the drifters, their lat, lon, ve, vn etc. 
                get_daily_avg_all_drifters(unique_drifter_month_dir, daily_avg_dir, d)
            
            #Now we have the daily averages to compare against OSCAR, however, the coordinates arent aligned, so we need to interpolate OSCAR to drifter coordinates. 
        

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

           
            # we are adding each interpolated day to a list, so we should have a monthe of
            #drifter_d, drifter u, drifer v, cmems u, cmems v, neurost u, neurost v, lat lon, time
          

            month_of_validation_data.append(ds_drifter)
    


        print(month_of_validation_data[0])
        print(month_of_validation_data[0].dims)
        validation_path = os.path.join(validation_dir, year, f'validation_{year}_{month}.nc')
        if not os.path.exists(os.path.join(validation_dir, year, f'validation_{year}_{month}.nc')):
            if month_of_validation_data:
                #validation_ds = xr.open_dataset(validation_path)
                ds_monthly_validation = xr.concat(month_of_validation_data, dim='drifter')
                print(ds_monthly_validation)
                #validation_ds.close()
                ds_monthly_validation.to_netcdf(validation_path)
        #after the above block of code runs we now have one .nc file with a month of validation data
        #validation data is the interpplated, neurost and cmems 
        #now we need to appedn that to the master folder 
  
        #ACCOUNT FOR DUPLICATION
        
        append_to_master(year, month)

        month_of_validation_data.clear()

        master_path = os.path.join(validation_dir, 'validation_master.nc')
        ds = xr.open_dataset(master_path)



        df = ds[['latitude', 'longitude', 've', 'vn', 'cmems_ve', 'cmems_vn', 'neurost_ve', 'neurost_vn']].to_dataframe().dropna()
        df = df[
        (df['ve'] >= -2) & (df['ve'] <= 2) &
        (df['vn'] >= -2) & (df['vn'] <= 2)
        ]
        print(df)
        

        corr_maps = compute_binned_correlations(
            df,
            models=('cmems','neurost'),
            components=('ve','vn'),
            lat_range=(-60, 60),
            lon_range=(0, 360),
            bin_size=2,
            min_samples=5
        )


        figs = []
        figs.append(plot_velocity_comparison_scatter(df))
        figs.append(plot_validation_metrics(df, models=['cmems', 'neurost'], components=['ve', 'vn']))
        figs.append(plot_binned_correlation_map(corr_maps[('cmems','ve')], 'cmems', 've', zoom_extent=None, title_prefix="Jan 2020"))
        figs.append(plot_binned_correlation_map(corr_maps[('neurost','ve')], 'neurost', 've', zoom_extent=None, title_prefix="Jan 2020"))
        figs.append(plot_binned_correlation_map(corr_maps[('cmems','vn')], 'cmems', 'vn', zoom_extent=None, title_prefix="Jan 2020"))
        figs.append(plot_binned_correlation_map(corr_maps[('neurost','vn')], 'neurost', 'vn', zoom_extent=None, title_prefix="Jan 2020"))
        figs.append(plot_density(df, 'cmems', 've', 'CMEMS', save_path=None))
        figs.append(plot_density(df, "cmems", "vn", "CMEMS", save_path = None))
        figs.append(plot_density(df, "neurost", "ve", "NEUROST", save_path = None))
        figs.append(plot_density(df, "neurost", "vn", "NEUROST", save_path = None))
        figs.append(make_slope_map(df, 've', 'cmems_ve', 'CMEMS Zonal Slope' ))
        figs.append(make_slope_map(df, 'vn', 'cmems_vn', 'CMEMS Meridional Slope'))
        figs.append(make_slope_map(df, 've', 'neurost_ve', 'NEUROST Zonal Slope'))
        figs.append(make_slope_map(df, 'vn', 'neurost_vn', 'NEUROST Meridional Slope'))
        figs.append(make_rms_map(df, 'cmems',   've', 'CMEMS Zonal RMSD'))
        figs.append(make_rms_map(df, 'cmems',   'vn', 'CMEMS Meridional RMSD'))
        figs.append(make_rms_map(df, 'neurost', 've', 'NEUROST Zonal RMSD'))
        figs.append(make_rms_map(df, 'neurost', 'vn', 'NEUROST Meridional RMSD'))
        figs.append(make_residual_corr_map(df, 'cmems',   've', 'CMEMS Zonal Residual–Truth Corr'))
        figs.append(make_residual_corr_map(df, 'cmems',   'vn', 'CMEMS Meridional Residual–Truth Corr'))
        figs.append(make_residual_corr_map(df, 'neurost', 've', 'NEUROST Zonal Residual–Truth Corr'))
        figs.append(make_residual_corr_map(df, 'neurost', 'vn', 'NEUROST Meridional Residual–Truth Corr'))


        
        save_plots_to_pdf(
            figs,
            explanations,
            pdf_path="validation_report.pdf",
            metadata={"Title":"Validation Report","Author":"OSCAR"})
                
        
        
        

    return

def run_download():
    
    if download_cmems:
        dates = parse_dates(download_dates_cmems)
        download_ssh_cmems(dates)
    if download_neurost:
        dates = parse_dates(download_dates_neurost)
        download_ssh_neurost(dates)


    if download_wind:
        dates = parse_dates(download_dates_wind)
        download_wind_era5(dates)
    
    if download_sst:
        dates = parse_dates(download_dates_sst)
        download_sst_cmc(dates)

    if download_drifter:
        download_drifter_data(drifter_download_date)

    return

def main():
    np.seterr(divide='ignore', invalid='ignore', over='ignore')
    start_dt = np.datetime64(start_date)
    end_dt = np.datetime64(end_date)
    dates_to_process = get_dates_to_process(start_dt, end_dt, output_dir, override)

    if do_download:
        run_download()

    if ready_for_oscar_mode(dates_to_process, oscar_mode, ssh_mode):
        if do_currents:
            np.seterr(divide='ignore', invalid='ignore', over='ignore')
            start_dt = np.datetime64(start_date)
            end_dt = np.datetime64(end_date)
            dates_to_process = get_dates_to_process(start_dt, end_dt, output_dir, override)
            run_compute_currents(dates_to_process)
    else:
        print("you have missing files")


    if do_plots:
        run_plotting(dates_to_plot)

    if do_validation:
        run_validation(dates_to_validate)
      


        

if __name__ == "__main__":
    main()
