import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from ..config.setup import *
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
import os
from glob import glob
from scipy.interpolate import RegularGridInterpolator
import pandas as pd
from matplotlib.colors import LogNorm
from scipy.stats import linregress
   


def bin_drifter_data(ds_drifter, res):
    '''
    -This function should take in drifter data (lat, lon, ve and vn)
    -Create an oscar grid 
    '''
    # Wrap longitude
    ds_drifter['longitude'] = ds_drifter['longitude'] % 360
    

    lat = ds_drifter['latitude'].values
    lon = ds_drifter['longitude'].values
    ve = ds_drifter['ve'].values
    vn = ds_drifter['vn'].values

    lat_bins = np.arange(-89.75, 90 + res, res)
    lon_bins = np.arange(0, 360 + res, res)


    #set up numpy array
    ve_sum = np.zeros((len(lat_bins) - 1, len(lon_bins) - 1))
    vn_sum = np.zeros((len(lat_bins) - 1, len(lon_bins) - 1))
    counts = np.zeros((len(lat_bins) - 1, len(lon_bins) - 1))


    # --- Bin the Data ---
    # Digitize returns the index of the bin each data point belongs to.
    # We subtract 1 because digitize is 1-based.
    lat_indices = np.digitize(lat, lat_bins) - 1 
    lon_indices = np.digitize(lon, lon_bins) - 1
    
    # print(f"Point {0}: lat = {lat[0]}, lon = {lon[0]}")
    # print(lat_indices[25])
    # print(lon_indices[25])

    # --- Loop Through Drifter Points and Accumulate ---
    for i in range(len(lat)):
        # Check for valid indices and non-NaN velocity values
        lat_idx = lat_indices[i]
        lon_idx = lon_indices[i]

        
        if (0 <= lat_idx < ve_sum.shape[0] and
            0 <= lon_idx < ve_sum.shape[1] and
            not np.isnan(ve[i]) and not np.isnan(vn[i])):

            
            ve_sum[lat_idx, lon_idx] += ve[i]
            vn_sum[lat_idx, lon_idx] += vn[i]
            counts[lat_idx, lon_idx] += 1

    
    # print(ve_sum[748, 616]) 
    # print(vn_sum[748, 616])
    # print(counts[748, 616])

    #     # Get bin edges
    # lat_edge_min = lat_bins[748]
    # lat_edge_max = lat_bins[749]
    # lon_edge_min = lon_bins[616]
    # lon_edge_max = lon_bins[617]

    # # Get bin center (optional — often used for plotting)
    # lat_center = 0.5 * (lat_edge_min + lat_edge_max)
    # lon_center = 0.5 * (lon_edge_min + lon_edge_max)

    # print(f"Lat bin index 231: [{lat_edge_min}, {lat_edge_max}] → center: {lat_center}")
    # print(f"Lon bin index 1040: [{lon_edge_min}, {lon_edge_max}] → center: {lon_center}")

    with np.errstate(divide='ignore', invalid='ignore'):
        ve_avg = np.where(counts > 0, ve_sum / counts, np.nan)
        vn_avg = np.where(counts > 0, vn_sum / counts, np.nan)

    lat_coords = (lat_bins[:-1]) 
    lon_coords = (lon_bins[:-1]) 


    output_ds = xr.Dataset(
        {
            've': (('latitude', 'longitude'), ve_avg),
            'vn': (('latitude', 'longitude'), vn_avg),
            'point_counts': (('latitude', 'longitude'), counts)
        },
        coords={
            'latitude': ('latitude', lat_coords),
            'longitude': ('longitude', lon_coords)
        }
    )
    return output_ds



def compare_velocity_components(ds_oscar, ds_drifter):
    """
    Compare velocity components (u,v) from OSCAR with (ve,vn) from drifters.
    Returns error metrics for both components.
    """

    # Align datasets to overlapping lat/lon grid
    ds_o, ds_d = xr.align(ds_oscar, ds_drifter, join="inner")

    results = {}

    for model_var, drifter_var, label in zip(['u', 'v'], ['ve', 'vn'], ['u', 'v']):
        model_vals = ds_o[model_var].values.flatten()
        obs_vals   = ds_d[drifter_var].values.flatten()

        # Apply mask to ignore NaNs
        mask = ~np.isnan(model_vals) & ~np.isnan(obs_vals)
        model_clean = model_vals[mask]
        obs_clean = obs_vals[mask]

        if len(model_clean) == 0:
            results[label] = {
                'R2': np.nan,
                'RMSD': np.nan,
                'Difference': np.nan,
                'Pearson': np.nan,
                'N_points': 0
            }
            continue

        # Compute metrics
        rmsd = mean_squared_error(obs_clean, model_clean)
        r2 = r2_score(obs_clean, model_clean)
        difference = np.mean(model_clean - obs_clean)
        corr, _ = pearsonr(model_clean, obs_clean)

        results[label] = {
            'R2': r2,
            'RMSD': rmsd,
            'Difference': difference,
            'Pearson': corr,
            'N_points': len(obs_clean)
        }

    return results

def plot_velocity_comparison(obs, model, component='u', limits=(-1, 1)):
    """
    Scatter plot comparing drifter (obs) vs model (OSCAR) velocity component.
    Computes and displays R², RMSD, Difference, and Pearson correlation.

    Parameters:
        obs (np.ndarray): Drifter data (ve or vn)
        model (np.ndarray): Model data (u or v)
        component (str): 'u' or 'v' to label axes and title
        limits (tuple): Axis limits for plot
    """
    # Flatten and clean NaNs
    obs = obs.flatten()
    model = model.flatten()
    mask = ~np.isnan(obs) & ~np.isnan(model)
    obs_clean = obs[mask]
    model_clean = model[mask]

    # Metrics
    r2 = r2_score(obs_clean, model_clean)
    rmsd = np.sqrt(mean_squared_error(obs_clean, model_clean))
    difference = np.mean(model_clean - obs_clean)
    corr, _ = pearsonr(obs_clean, model_clean)

    # Plot
    plt.figure(figsize=(6, 6))
    plt.scatter(obs_clean, model_clean, alpha=0.4, label="Data", s=15)
    plt.plot(limits, limits, 'k--', label='1:1 Line')

    plt.xlabel(f"Drifter {component} (m/s)")
    plt.ylabel(f"Model {component} (m/s)")
    plt.title(f"{component.upper()} Component Velocity Comparison")
    plt.xlim(limits)
    plt.ylim(limits)
    plt.grid(True)
    plt.legend()

    # Metrics annotation box
    text = (
        f"R²: {r2:.3f}\n"
        f"RMSD: {rmsd:.3f} m/s\n"
        f"Difference: {difference:.3f} m/s\n"
        f"Pearson: {corr:.3f}\n"
        f"N: {len(obs_clean)}"
    )
    plt.text(
        limits[0] + 0.05, limits[1] - 0.15,
        text,
        bbox=dict(facecolor='white', alpha=0.8),
        fontsize=10
    )

    plt.tight_layout()
    plt.show()

def plot_log_histogram(obs, model, component='u', bins=100):
    """
    Create a log-scaled 2D histogram of model vs obs velocity.
    """
    # Flatten and mask NaNs
    obs = obs.flatten()
    model = model.flatten()
    mask = ~np.isnan(obs) & ~np.isnan(model)
    obs_clean = obs[mask]
    model_clean = model[mask]

    # Fit regression line
    X = obs_clean.reshape(-1, 1)
    y = model_clean
    reg = LinearRegression().fit(X, y)
    y_fit = reg.predict(X)

    slope = reg.coef_[0]
    intercept = reg.intercept_

    # Plot log-scaled 2D histogram
    plt.figure(figsize=(5, 5))
    plt.hist2d(obs_clean, model_clean, bins=bins, norm=plt.matplotlib.colors.LogNorm(), cmap='turbo')

    # 1:1 reference line
    lims = [min(obs_clean.min(), model_clean.min()), max(obs_clean.max(), model_clean.max())]
    plt.plot(lims, lims, 'w-', label='1:1 line')

    # Regression line
    x_vals = np.linspace(lims[0], lims[1], 100)
    y_vals = slope * x_vals + intercept
    plt.plot(x_vals, y_vals, 'r-', label=f'Fit: y={slope:.2f}x+{intercept:.2f}')

    # Labels and title
    plt.xlabel(f"Drifter {component}")
    plt.ylabel(f"OSCAR {component}")
    plt.title(f"Log histogram SL={slope:.2f}")
    plt.axis('square')
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_binned_locations(ds):
    # Mask for valid velocity points
    valid_mask = np.isfinite(ds['ve']) & np.isfinite(ds['vn'])

    # Extract latitude and longitude values
    lats = ds['latitude'].values
    lons = ds['longitude'].values

    # Create 2D lat/lon mesh
    lon2d, lat2d = np.meshgrid(lons, lats)

    # Apply the mask
    valid_lats = lat2d[valid_mask]
    valid_lons = lon2d[valid_mask]

    # Normalize longitudes to [-180, 180]
    valid_lons = (valid_lons + 360) % 360
    valid_lons[valid_lons > 180] -= 360

    # Plot
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.set_global()

    ax.scatter(valid_lons, valid_lats, s=0.5, color='blue', transform=ccrs.PlateCarree())

    plt.title("Locations with Valid Velocity Data")
    plt.tight_layout()
    os.makedirs(os.path.dirname('/Users/stewarta/Desktop/oscarpy/diagnostics/drifters/binned_loc.png'), exist_ok=True)
    plt.savefig('/Users/stewarta/Desktop/oscarpy/diagnostics/drifters/binned_loc.png', dpi=300)
    plt.close()


#==========================================================================#
#def plot_daily_drifter_movement()

def get_unique_drifter_daily_avg(ds, output_folder):
    '''
    Input: month of 6 hourly drifters in one netcdf
    Output: individual drifter daily average  points for the year
    '''
    os.makedirs(output_folder, exist_ok=True)
    unique_ids = np.unique(ds['ID'].values)


    for drifter_id in unique_ids:
        
        mask = ds['ID'] == drifter_id
        drifter_ds = ds.where(mask, drop=True)
        drifter_ds = drifter_ds.sortby('time')
        daily_ds = drifter_ds.resample(time='1D').mean()
        daily_ds['ID'] = xr.DataArray(np.full(daily_ds.dims['time'], drifter_id), dims='time')
       
        outfile = os.path.join(output_folder, f"drifter_{drifter_id}_daily_avg.nc")
        daily_ds.to_netcdf(outfile)

        


def get_daily_avg_all_drifters(unique_drifters_dir, daily_avg_dir, selected_date):

    selected_date = np.datetime64(selected_date)

    # Initialize containers
    latitudes = []
    longitudes = []
    ve_list = []
    vn_list = []
    ve_oscar_list = []
    vn_oscar_list = []
    ids = []

    # Get all .nc files in folder
    file_list = sorted(glob(os.path.join(unique_drifters_dir, "*.nc")))


    for file in file_list:
        ds = xr.open_dataset(file)
    


        if selected_date in ds['time']:
            try:
                ve = ds['ve'].sel(time=selected_date).values.item()
                vn = ds['vn'].sel(time=selected_date).values.item()
                lat = ds['latitude'].sel(time=selected_date).values.item()
                lon = ds['longitude'].sel(time=selected_date).values.item()
                drifter_id = ds['ID'].sel(time=selected_date).values.item()

        

                ve_list.append(ve)
                vn_list.append(vn)
                latitudes.append(lat)
                longitudes.append(lon)
                ids.append(str(drifter_id))


            except Exception as e:
                print(f"⚠️ Skipping {file}: {e}")

    # Combine into dataset
    snapshot_ds = xr.Dataset({
        'latitude': xr.DataArray(latitudes, dims='drifter'),
        'longitude': xr.DataArray(longitudes, dims='drifter'),
        've': xr.DataArray(ve_list, dims='drifter'),
        'vn': xr.DataArray(vn_list, dims='drifter'),
        'ID': xr.DataArray(ids, dims='drifter'),
        'time': xr.DataArray([selected_date] * len(ids), dims='drifter')
    })
    # Convert to datetime object for formatting
    date_str = str(selected_date)[:10]  # 'YYYY-MM-DD'
    year, month, day = date_str.split('-')
    folder_path = os.path.join(daily_avg_dir)
    os.makedirs(folder_path, exist_ok=True)
    output_name = f"drifters_{year}_{month}_{day}.nc"
    output_path = os.path.join(folder_path, output_name)

    snapshot_ds.to_netcdf(output_path)
  



def interpolate_oscar_to_drifters(ds_oscar, ds_drifter, u_var = 'u', v_var = 'v', prefix = 'oscar', method='linear'):
    # Get OSCAR lat/lon and fields
    lats = ds_oscar['lat'].values
    lons = ds_oscar['lon'].values
    u = ds_oscar['u'].values
    v = ds_oscar['v'].values

    

    # If OSCAR has time dimension, drop it
    if u.ndim == 3:
        u = u[0, :, :]
        v = v[0, :, :]

    u = u.T  # shape now (719, 1440)
    v = v.T

            # Get drifter lat/lon
    drifter_lats = ds_drifter['latitude'].values
    drifter_lons = ds_drifter['longitude'].values



    # Handle lon wraparound if needed
    if np.any(lons > 180) and np.any(drifter_lons < 0):
        drifter_lons = (drifter_lons + 360) % 360
    # Set up interpolators
    u_interp = RegularGridInterpolator((lats, lons), u, method=method, bounds_error=False, fill_value=np.nan)
    v_interp = RegularGridInterpolator((lats, lons), v, method=method, bounds_error=False, fill_value=np.nan)


    # Stack into interpolation points
    interp_points = np.column_stack((drifter_lats, drifter_lons))

    # Interpolate
    u_interp_vals = u_interp(interp_points)
    v_interp_vals = v_interp(interp_points)

    ds_out = ds_drifter.copy()
    ds_out[f'{prefix}_ve'] = (('drifter',), u_interp_vals)
    ds_out[f'{prefix}_vn'] = (('drifter',), v_interp_vals)

    return ds_out



#==============DRIFTER vs OSCAR ============#

def summarize_time_jumps(df, col="time"):
    s = df[col].astype(str).str.replace(r"\s+", "", regex=True)
    d = pd.to_datetime(s, errors="coerce").dropna().dt.normalize()
    u = d.sort_values().drop_duplicates().reset_index(drop=True)
    if u.empty:
        return {"ranges": [], "breakpoints": [], "summary": "", "years": [], "year_months": [], "unique_days": 0}

    diff = u.diff().dt.days
    starts = [0] + [i for i, v in enumerate(diff.iloc[1:], start=1) if v != 1]
    ends   = [i - 1 for i in starts[1:]] + [len(u) - 1]

    ranges = [(u.iloc[s], u.iloc[e]) for s, e in zip(starts, ends)]
    bps = [ranges[0][0]] + [r[0] for r in ranges[1:]] + [ranges[-1][1]]
    summary = " to ".join(ts.strftime("%Y-%m-%d") for ts in bps)

    return summary
        
 
def plot_velocity_comparison_scatter(df, temporal_range, ssh):
    def compute_metrics(y_true, y_pred):
        # Flatten and mask NaNs
        y_true = np.array(y_true).flatten()
        y_pred = np.array(y_pred).flatten()
        mask = ~np.isnan(y_true) & ~np.isnan(y_pred)
        
        y_true_clean = y_true[mask]
        y_pred_clean = y_pred[mask]

        # Metrics
        difference = np.mean(y_pred_clean - y_true_clean)
        rmsd = np.sqrt(np.mean((y_pred_clean - y_true_clean) ** 2))
        r2 = 1 - (np.sum((y_true_clean - y_pred_clean) ** 2) / 
                  np.sum((y_true_clean - np.mean(y_true_clean)) ** 2))

        return r2, rmsd, difference
    
    fig = plt.figure(figsize=(14, 6))

    # ---- Zonal (U / ve) ----
    ax1 = plt.subplot(1, 2, 1)
    ax1.scatter(df['ve'], df[f'{ssh}_ve'], label= ssh.upper(), alpha=0.5, s=5)
    ax1.plot([-2, 2], [-2, 2], 'k--', linewidth=1)  # 1:1 line
    ax1.set_xlabel("Drifter Zonal Velocity(m/s)")
    ax1.set_ylabel(f"{ssh.upper()} Zonal Velocity (m/s)")
    #ax1.set_title(f"Zonal Velocity Comparison {temporal_range}")
    ax1.legend()
    ax1.grid(True)

    r2_u, rmsd_u, diff_u = compute_metrics(df['ve'], df[f'{ssh}_ve'])
    ax1.text(
        0.05, 0.95,
        f"R² = {r2_u:.3f}\nRMSD = {rmsd_u:.3f}\nBias = {diff_u:.3f}",
        transform=ax1.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7)
    )

    # ---- Meridional (V / vn) ----
    ax2 = plt.subplot(1, 2, 2)
    ax2.scatter(df['vn'], df[f'{ssh}_vn'], label=ssh.upper(), alpha=0.5, s=5)
    ax2.plot([-2, 2], [-2, 2], 'k--', linewidth=1)  # 1:1 line
    ax2.set_xlabel("Drifter Meridonal Velocity (m/s)")
    ax2.set_ylabel(f"{ssh.upper()} Meridonal Velocity (m/s)")
    #ax2.set_title(f"Meridional Velocity Comparison {temporal_range}")
    ax2.legend()
    ax2.grid(True)

    # Compute stats for meridional
    r2_v, rmsd_v, diff_v = compute_metrics(df['vn'], df[f'{ssh}_vn'])
    ax2.text(
        0.05, 0.95,
        f"R² = {r2_v:.3f}\nRMSD = {rmsd_v:.3f}\nBias = {diff_v:.3f}",
        transform=ax2.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7)
    )

    plt.suptitle(f"{ssh.upper()} vs Drifter Velocity Comparison\n{temporal_range}", fontsize=14, y=1.02)

    plt.tight_layout()
    #plt.show()

    return fig



def plot_validation_metrics(df, temporal_range, models, components):
    def compute_metrics(y_true, y_pred):
        # Flatten and mask NaNs
        y_true = np.array(y_true).flatten()
        y_pred = np.array(y_pred).flatten()
        mask = ~np.isnan(y_true) & ~np.isnan(y_pred)
        
        y_true_clean = y_true[mask]
        y_pred_clean = y_pred[mask]

        # Metrics
        difference = np.mean(y_pred_clean - y_true_clean)
        rmsd = np.sqrt(np.mean((y_pred_clean - y_true_clean) ** 2))
        r2 = 1 - (np.sum((y_true_clean - y_pred_clean) ** 2) / 
                  np.sum((y_true_clean - np.mean(y_true_clean)) ** 2))

        return r2, rmsd, difference

    # Compute metrics for each model & component
    metrics = {}
    for model in models:
        for component in components:
            r2, rmsd, difference = compute_metrics(
                df[component].values, 
                df[f"{model}_{component}"].values
            )
            metrics[f"{model}_{component}"] = {
                'R²': r2,
                'RMSD': rmsd,
                'Difference': difference
            }

    # Plot metrics
    metric_names = ['R²', 'RMSD', 'Difference']
    fig, axs = plt.subplots(1, 3, figsize=(16, 5))

    for i, metric_name in enumerate(metric_names):
        ax = axs[i]
        values = []

        for comp in components:
            for model in models:
                values.append(metrics[f"{model}_{comp}"][metric_name])

        x_labels = [f"{model.upper()} {comp}" for comp in components for model in models]
        colors = ['blue', 'red'] * len(components)

        ax.bar(x_labels, values, color=colors)
        ax.set_title(metric_name)
        ax.set_ylabel(metric_name)
        ax.grid(True)

    plt.suptitle(f"Validation Metrics: CMEMS & NEUROST vs Drifters {temporal_range}", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    return fig



def compute_binned_correlations(
    df,
    models=['cmems'],
    components=('ve', 'vn'),
    lat_range=(-60, 60),
    lon_range=(0, 360),
    bin_size=2,
    min_samples=5
):
    """
    Returns a dict of 2D correlation DataFrames keyed by (model, component).
    """
    def safe_corr(x, y, min_n=min_samples):
        m = ~np.isnan(x) & ~np.isnan(y)
        if np.count_nonzero(m) < min_n:
            return np.nan
        return pearsonr(x[m], y[m])[0]

    # Bins
    lat_bins = np.arange(lat_range[0], lat_range[1] + bin_size, bin_size)
    lon_bins = np.arange(lon_range[0], lon_range[1] + bin_size, bin_size)

    data = df.copy()
    data['lon_mod'] = (data['longitude'] % 360)
    data['lat_bin'] = pd.cut(data['latitude'], bins=lat_bins, labels=lat_bins[:-1])
    data['lon_bin'] = pd.cut(data['lon_mod'],  bins=lon_bins, labels=lon_bins[:-1])

    grouped = data.groupby(['lat_bin', 'lon_bin'], observed=True)

    corr_maps = {}
    for model in models:
        for comp in components:
            truth_col = comp
            pred_col = f"{model}_{comp}"

            corr_df = grouped.apply(lambda g: safe_corr(g[truth_col].values, g[pred_col].values)).unstack()
            corr_df.index = corr_df.index.astype(float)
            corr_df.columns = corr_df.columns.astype(float)
            corr_df = corr_df.sort_index().sort_index(axis=1)

            corr_maps[(model, comp)] = corr_df

    return corr_maps


def plot_binned_correlation_map(
    temporal_range,
    corr_df,
    model,
    component,
    zoom_extent=None,
    cmap='jet',
    vmin=0.0,
    vmax=1.0,
    title_prefix=None
):
  
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    lons = corr_df.columns.values
    lats = corr_df.index.values

    mesh = ax.pcolormesh(lons, lats, corr_df.values, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.coastlines()

    # Gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    # Optional zoom
    if zoom_extent is not None:
        lon_min, lon_max, lat_min, lat_max = zoom_extent
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # Title
    comp_name = 'Zonal' if component == 've' else 'Meridional'
    if title_prefix:
        title = f"{title_prefix} — {model.upper()} vs Drifters ({comp_name} on {temporal_range})"
    else:
        title = f"{model.upper()} vs Drifters ({comp_name} on {temporal_range})"
    ax.set_title(title, fontsize=13)

    cbar = plt.colorbar(mesh, ax=ax, orientation='vertical', fraction=0.03, pad=0.02)
    cbar.set_label('Correlation (R)')

    plt.tight_layout()
    return fig




def plot_density(
    df,
    temporal_range,
    model_prefix: str,            # 'cmems' or 'neurost'
    component: str,               # 've' or 'vn'
    model_label: str,             # 'CMEMS' or 'NEUROST'
    xlim=(-0.6, 0.6),
    ylim=(-0.6, 0.6),
    bins=150,
    log=False,
    cmap="turbo",
    title_suffix="1993 - Present",
    save_path: str | None = None, # e.g. os.path.join(plots_dir, "cmems_ve_lin_hist2d.png")
    show: bool = False,
    dpi: int = 300
):
   

    assert component in ("ve", "vn"), "component must be 've' or 'vn'"
    comp_label = "u" if component == "ve" else "v"

    # ground truth and model arrays
    x = np.asarray(df[component].values)
    y = np.asarray(df[f"{model_prefix}_{component}"].values)

    # clean + clip
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    m = (x >= xlim[0]) & (x <= xlim[1]) & (y >= ylim[0]) & (y <= ylim[1])
    x = x[m]; y = y[m]

    # figure/axes
    fig, ax = plt.subplots(figsize=(8, 6))
    norm = LogNorm() if log else None

    # density histogram
    h = ax.hist2d(x, y, bins=bins, range=[xlim, ylim], cmap=cmap, norm=norm)
    mappable = h[3]

    # y=x reference
    ax.plot(xlim, xlim, color="white", lw=1, zorder=3)

    # OLS fit (only if enough points)
    if x.size >= 3:
        slope, intercept, *_ = linregress(x, y)
        xx = np.array(xlim)
        ax.plot(xx, slope * xx + intercept, color="red", lw=2, zorder=4)
        ax.text(
            0.02, 0.95, f"SL={slope:.2f}",
            transform=ax.transAxes, ha="left", va="top",
            color="w", fontsize=10,
            bbox=dict(facecolor='black', alpha=0.5, edgecolor='none', pad=1)
        )

    # labels, limits, titles
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.set_title(f"{'Log ' if log else ''}histogram {model_label} {comp_label} {temporal_range}", fontsize=12)
    ax.set_xlabel(f"Drifter {comp_label} (m/s)")
    ax.set_ylabel(f"OSCAR {comp_label} (m/s)")

    # colorbar
    cb = fig.colorbar(mappable, ax=ax, orientation='vertical', fraction=0.035, pad=0.04)
    cb.set_label("Counts" + (" (log scale)" if log else ""))

    fig.suptitle(f"{model_label} vs Drifters ({title_suffix})", fontsize=13)
    fig.tight_layout()

    # optional save/show
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig.savefig(save_path, dpi=dpi)
    if show:
        plt.show()

    return fig


def make_slope_map(df_in,truth_col, model_col, title):
    """
    Compute per-grid-cell OLS slope between truth_col and model_col,
    plot the slope field, and return the figure object.
    """
    # ---- hard-coded config ----
    lat_bins = np.arange(-60, 61, 2)
    lon_bins = np.arange(0, 361, 2)
    MIN_SAMPLES = 5
    vmin, vmax = 0, 1
    cmap = 'jet'

    # ---- prep & bin ----
    df = df_in.copy()
    df['longitude'] = df['longitude'] % 360
    df['lat_bin'] = pd.cut(df['latitude'], bins=lat_bins, labels=lat_bins[:-1])
    df['lon_bin'] = pd.cut(df['longitude'], bins=lon_bins, labels=lon_bins[:-1])
    grouped = df.groupby(['lat_bin', 'lon_bin'])

    # ---- slope helper ----
    def _safe_slope(x, y, nmin=MIN_SAMPLES):
        x = np.asarray(x); y = np.asarray(y)
        m = np.isfinite(x) & np.isfinite(y)
        x = x[m]; y = y[m]
        if x.size < nmin or np.nanstd(x) == 0:
            return np.nan
        return linregress(x, y).slope

    # ---- compute slope grid ----
    s = grouped.apply(lambda g: _safe_slope(g[truth_col], g[model_col]))
    grd = s.unstack()
    grd.index = grd.index.astype(float)
    grd.columns = grd.columns.astype(float)
    grd = grd.sort_index().sort_index(axis=1)

    # ---- plot ----
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    lons = grd.columns.values
    lats = grd.index.values
    im = ax.pcolormesh(lons, lats, grd.values, vmin=vmin, vmax=vmax, cmap=cmap)

    ax.coastlines()
    ax.set_title(title, fontsize=14)

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    cb = plt.colorbar(im, ax=ax, orientation='vertical', fraction=0.01, pad=0.04)
    cb.set_label('Slope')

    plt.tight_layout()

    return fig

def make_rms_map(df_in, model_prefix, component, title, vmin=0.0, vmax=0.5, zoom=None):
    # --- hard-coded config ---
    LAT_BINS = np.arange(-60, 61, 2)
    LON_BINS = np.arange(0, 361, 2)
    MIN_SAMPLES = 5

    # --- prep & bin ---
    df = df_in.copy()
    df['longitude'] = df['longitude'] % 360
    df['lat_bin'] = pd.cut(df['latitude'],  bins=LAT_BINS, labels=LAT_BINS[:-1])
    df['lon_bin'] = pd.cut(df['longitude'], bins=LON_BINS, labels=LON_BINS[:-1])
    grouped = df.groupby(['lat_bin','lon_bin'])

    # --- RMSD helper ---
    def _safe_rmsd(truth, pred, nmin=MIN_SAMPLES):
        t = np.asarray(truth); p = np.asarray(pred)
        m = np.isfinite(t) & np.isfinite(p)
        t, p = t[m], p[m]
        if t.size < nmin:
            return np.nan
        return float(np.sqrt(np.mean((p - t)**2)))

    # --- compute RMSD grid ---
    truth_col = 've' if component == 've' else 'vn'
    model_col = f'{model_prefix}_{component}'
    s = grouped.apply(lambda g: _safe_rmsd(g[truth_col], g[model_col]))
    grid = s.unstack()
    grid.index = grid.index.astype(float)
    grid.columns = grid.columns.astype(float)
    grid = grid.sort_index().sort_index(axis=1)

    # --- plot ---
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    lons = grid.columns.values
    lats = grid.index.values
    im = ax.pcolormesh(lons, lats, grid.values, vmin=vmin, vmax=vmax,
                       cmap='jet', transform=ccrs.PlateCarree())

    ax.coastlines()
    ax.set_title(title, fontsize=14)

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    if zoom:
        lon_min, lon_max, lat_min, lat_max = zoom
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    cb = plt.colorbar(im, ax=ax, orientation='vertical', fraction=0.01, pad=0.04)
    cb.set_label('RMSD (m/s)', fontsize=12)

    plt.tight_layout()
    # no save, no show
    return fig


def make_residual_corr_map(df, model_prefix, component, title, vmin=-0.2, vmax=0.2, zoom=None):
 
    # --- hard-coded grid config ---
    LAT_BINS = np.arange(-60, 61, 2)
    LON_BINS = np.arange(0, 361, 2)
    MIN_SAMPLES = 5

    # --- prep & bin ---
    d = df.copy()
    d["longitude"] = d["longitude"] % 360
    d["lat_bin"] = pd.cut(d["latitude"],  bins=LAT_BINS, labels=LAT_BINS[:-1])
    d["lon_bin"] = pd.cut(d["longitude"], bins=LON_BINS, labels=LON_BINS[:-1])
    grouped = d.groupby(["lat_bin", "lon_bin"])

    # --- correlation helper (residual vs truth) ---
    def _safe_corr(a, b, nmin=MIN_SAMPLES):
        a = np.asarray(a); b = np.asarray(b)
        m = np.isfinite(a) & np.isfinite(b)
        if m.sum() < nmin:
            return np.nan
        a = a[m]; b = b[m]
        # guard against zero variance
        if np.nanstd(a) == 0 or np.nanstd(b) == 0:
            return np.nan
        return float(np.corrcoef(a, b)[0, 1])

    # --- compute residual correlation grid ---
    truth_col = "ve" if component == "ve" else "vn"
    model_col = f"{model_prefix}_{component}"
    s = grouped.apply(lambda g: _safe_corr(g[model_col] - g[truth_col], g[truth_col]))
    grid = s.unstack()
    grid.index = grid.index.astype(float)
    grid.columns = grid.columns.astype(float)
    grid = grid.sort_index().sort_index(axis=1)

    # --- plot ---
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    lons = grid.columns.values
    lats = grid.index.values
    im = ax.pcolormesh(
        lons, lats, grid.values,
        cmap="RdBu_r", vmin=vmin, vmax=vmax,
        transform=ccrs.PlateCarree()
    )

    ax.coastlines()
    ax.set_title(title, fontsize=14)

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color="gray", alpha=0.5, linestyle="--")
    gl.top_labels = False
    gl.right_labels = False

    if zoom:
        lon_min, lon_max, lat_min, lat_max = zoom
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    cb = plt.colorbar(im, ax=ax, orientation="vertical", fraction=0.035, pad=0.04)
    cb.set_label("Correlation(model−drifter, drifter)", fontsize=11)

    plt.tight_layout()
    # no save, no show
    return fig



def clean_validation_df(ds, ssh_mode):
    """
    Convert dataset to dataframe with core vars and filter unrealistic drifter velocities.
    Keeps rows where -2 <= ve <= 2 and -2 <= vn <= 2.
    """
    df = ds[['latitude', 'longitude', 've', 'vn',
             f'{ssh_mode}_ve', f'{ssh_mode}_vn',
             'time']].to_dataframe().dropna()

    df = df[(df['ve'].between(-2, 2)) & (df['vn'].between(-2, 2))]
    return df
