from .computation.interp_and_grad import *
from .utils.data_utils import *
from pyoscar.analysis.plot_currents import *
from .computation.compute_currents import *
from .config.setup import *


def main():
    np.seterr(divide='ignore', invalid='ignore', over='ignore')
    start_dt = np.datetime64(start_date)
    end_dt = np.datetime64(end_date)

    dates_to_process = get_dates_to_process(start_dt, end_dt, output_dir, override)

    for date in dates_to_process:
   

      ssh_ds= load_ds(date, var=ssh_var)
      ssh_ds = interpolate_dataset(ssh_ds)
      plot_interp_and_grad(ssh_ds[ssh_var],title = "Interpolated {var} on {day}",
                            filename = filename_interp, save_dir = save_dir_fig, var = ssh_var, var_mode = ssh_mode)
      ssh_ds = calculate_gradient(ssh_ds, ssh_var)
      plot_interp_and_grad(ssh_ds[ssh_var + 'x'], title = "{var} on {day}", 
                            filename = filename_gradient, save_dir = save_dir_fig,
                              var = ssh_var +'x', var_mode = ssh_mode)
      plot_interp_and_grad(ssh_ds[ssh_var + 'y'], title = "{var} on {day}", 
                            filename = filename_gradient, save_dir = save_dir_fig,
                              var = ssh_var + 'y', var_mode = ssh_mode)

      sst_ds = load_ds(date, var=sst_var)
      sst_ds = interpolate_dataset(sst_ds)
      plot_interp_and_grad(sst_ds[sst_var], title = "Interpolated {var} on {day}",
                            filename = filename_interp, save_dir = save_dir_fig,
                              var = sst_var, var_mode = sst_mode)
      sst_ds = calculate_gradient(sst_ds, sst_var)
      plot_interp_and_grad(sst_ds[sst_var + 'x'], title = "{var} on {day}", 
                            filename = filename_gradient, save_dir = save_dir_fig, 
                            var = sst_var +'x', var_mode = sst_mode)
      plot_interp_and_grad(sst_ds[sst_var + 'y'], title = "{var} on {day}", 
                            filename = filename_gradient, save_dir = save_dir_fig, 
                            var = sst_var +'y', var_mode = sst_mode)


      wind_ds = load_ds(date, var = 'wind')
      wind_ds = interpolate_dataset(wind_ds)
      plot_interp_and_grad(wind_ds[wind_var_u], 
                            title = "Interpolated {var} on {day}", 
                            filename = filename_interp, save_dir = save_dir_fig, 
                            var = wind_var_u, var_mode=wind_mode)
      plot_interp_and_grad(wind_ds[wind_var_u],
                            title = "Interpolated {var} on {day}", 
                            filename = filename_interp, save_dir = save_dir_fig, 
                            var = wind_var_v, var_mode =wind_mode)
      

      #part 2
      ref_ds = ssh_ds.drop_vars(['sshx', 'sshy'])

      print(f"The currents computation for {date} is starting now")
      Ug, Uw, Ub = compute_surface_currents(ssh_ds, wind_ds, sst_ds, do_eq)



      #part 3
      write_oscar_to_podaac(ref_ds, Ug, Uw, Ub, podaacfile, podaacdir,
                              ssh_long_desc, wind_long_desc, sst_long_desc, 
                              oscar_long_desc, oscar_summary, oscar_id, 
                              doi, ssh_mode)
      

      print("All done")


main()