import numpy as np
from .computation.interp_and_grad import *
from .computation.compute_currents import *
from .utils.data_utils import *
from .config.setup import *  # contains: compute_currents, plot_currents, validate_currents, start_date, end_date, etc.
from .analysis.plot_currents import *
from .analysis.drifter_func import *
import pandas as pd
from .utils.download import *

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
    
    for date in date_list:

        #will create a "get_drifter_path_function"
        drifter_path = "/Users/stewarta/Desktop/oscarpy/datasets/DRIFTERS/drifter_6hour_qc_e065_93e9_5ff9_U1752698601011.nc"
        ds_drifter = xr.open_dataset(drifter_path)
        
        plot_drifter_map(ds_drifter)
        valid_points = np.isfinite(ds_drifter['ve'].values)
        print(np.sum(valid_points))
        #plot_drifter_currents(ds_drifter)
        
        cmems_path = "/Users/stewarta/Desktop/NASA_OSCAR/datasets/CURRENTS/INTERIM/PODAAC/2020/oscar_currents_interim_20200101.nc"

        #cmems_path = get_file_path(date, 'cmems')
        #neurost_path = get_file_path(date, 'neurost')

        neurost_path = "/Users/stewarta/Desktop/oscarpy/datasets/CURRENTS/INTERIM/PODAAC/NEUROST/2020/01/oscar_currents_interim20200101.nc"

        ds_oscar_cmems = xr.open_dataset(cmems_path)
        ds_oscar_neurost = xr.open_dataset(neurost_path)
        print("ve attr range:", ds_drifter['ve'].min().item())
        print("vn attr range:", ds_drifter['vn'].max().item())


        ds_drifter_cmems = bin_drifter_data(ds_drifter, 0.25)
        plot_binned_locations(ds_drifter_cmems)
        valid_points = np.isfinite(ds_drifter_cmems['ve'].values)
        print(np.sum(valid_points))
        ds_drifter_neurost = bin_drifter_data(ds_drifter, 0.1)
        plot_binned_locations(ds_drifter_neurost)
        valid_points = np.isfinite(ds_drifter_neurost['ve'].values)
        print(np.sum(valid_points))
        # # plot_drifter_component(ds_drifter_cmems, component="ve", region = 'gulf_stream')
        # # plot_drifter_component(ds_drifter_neurost, component="ve", region = 'gulf_stream')



        neurost_metrics = compare_velocity_components(ds_oscar_neurost, ds_drifter_neurost)
        cmems_metrics = compare_velocity_components(ds_oscar_cmems, ds_drifter_cmems)
        print(neurost_metrics)
        print(cmems_metrics)
        plot_velocity_comparison(ds_drifter_neurost['ve'].values, ds_oscar_neurost['u'].values)
        plot_log_histogram(ds_drifter_neurost['ve'].values, ds_oscar_neurost['u'].values, component='u')

        

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
