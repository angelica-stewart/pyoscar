import os
import matplotlib.pyplot as plt
import numpy as np 
from ..utils.data_utils import extract_date, write_netcdf_file
from ..config.setup import REGION, FIG_ROOT



def plot_interp_and_grad(data_array, title, filename, save_dir, var, var_mode):

    if 'time' in data_array.dims:
        for i in range(data_array.sizes['time']):
 
            day = str(data_array['time'].values[i])[:10]
            year = str(day)[:4]
            month = day[5:7]

            plot_title = title.format(var=var.upper(), day=day)
            plot_filename = filename.format(var=var.lower(), day=day)

            folder_path = os.path.join(save_dir, var.lower(), var_mode, year, month)
            os.makedirs(folder_path, exist_ok=True)

            plt.figure(figsize=(10, 6))
            data_array.isel(time=i).scatter()
            plt.title(plot_title)
                    
            save_path = os.path.join(folder_path, plot_filename)
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            plt.close()

            print(f"Saved figure {var} to path: {save_path}")
            write_netcdf_file(data_array, folder_path, plot_filename)
         
    else:

        plt.figure(figsize=(10, 6))
        data_array.plot()
        plt.title(plot_title)
                
        save_path = os.path.join(save_dir, plot_filename)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"Saved figure {var} to path: {save_path}")
        write_netcdf_file(data_array, save_path, title)

def extract_data(ds):
    u = ds['u'].isel(time=0).where(ds['u'] != -999.)
    v = ds['v'].isel(time=0).where(ds['v'] != -999.)
    lon = ds['lon']
    lat = ds['lat']
    
    lon_adj = lon.where(lon <= 180, lon - 360)  
    #lon_adj = ((lon + 180) % 360) - 180
    if REGION == 'global':
        lon_adj = lon_adj.sortby(lon_adj)
    
    region_bounds = {
        'global':      {'lon_min': -180, 'lon_max': 180,  'lat_min': -90,  'lat_max': 90},
        'gulf_stream': {'lon_min': -85,  'lon_max': -60,  'lat_min': 25,   'lat_max': 45},
        'caribbean':  {'lon_min': -90,  'lon_max': -55,  'lat_min': 9,    'lat_max': 23},
        'indian_ocean':     {'lon_min': 20,   'lon_max': 120,  'lat_min': -50,  'lat_max': 30},
        'australia_nz': {'lon_min': 110,'lon_max': 180,'lat_min': -50,'lat_max': -10}
    }

    if REGION not in region_bounds:
        raise ValueError(f"Unknown region '{REGION}'. Choose from: {list(region_bounds.keys())}")

    bounds = region_bounds[REGION]

    lon_mask = (lon_adj >= bounds['lon_min']) & (lon_adj <= bounds['lon_max'])
    lat_mask = (lat >= bounds['lat_min']) & (lat <= bounds['lat_max'])



    u_sub = u.where(lon_mask & lat_mask, drop=True).squeeze()
    v_sub = v.where(lon_mask & lat_mask, drop=True).squeeze()
    lon_sub = lon_adj.where(lon_mask, drop=True)
   
    lat_sub = lat.where(lat_mask, drop=True)
    speed_sub = np.sqrt(u_sub**2 + v_sub**2)
    lon2d, lat2d = np.meshgrid(lon_sub, lat_sub)

  
    return u_sub, v_sub, lon_sub, lat_sub, speed_sub, lon2d, lat2d

def plot_currents(ds_path, lon2d, lat2d, u, v, speed, title, ssh_mode, vmin=0, vmax=1.6):
    year, month = extract_date(ds_path)
    save_dir = os.path.join(FIG_ROOT, ssh_mode, REGION, year, month)
    os.makedirs(save_dir, exist_ok=True)

    # File path
    filename = os.path.basename(ds_path).replace(".nc", ".png")
    fig_path = os.path.join(save_dir, filename)


    
    plt.figure(figsize=(12, 8))

    if REGION == 'global':
        skip = 30
    else:
        skip =3

    # Plot speed as background
    p = plt.pcolormesh(lon2d, lat2d, speed.T, shading='auto', cmap='turbo', vmin=vmin, vmax=vmax)
    plt.colorbar(p, label='m/s')



    # Plot current vectors (quiver)
    plt.quiver(
        lon2d[::skip, ::skip],
        lat2d[::skip, ::skip],
        u.T[::skip, ::skip],
        v.T[::skip, ::skip],
        scale=25,
        color='black',
        width=0.002
    )

    plt.title(title)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved plot to: {fig_path}")