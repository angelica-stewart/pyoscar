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


def plot_drifter_map(ds, region='global'):
    """
    Plot drifter positions on a world map, with zoom options for specific regions.

    """
  

    # Flatten and clean lat/lon
    lat = ds['latitude'].values.flatten()
    lon = ds['longitude'].values.flatten()

    # Normalize longitudes to [-180, 180]
    lon = (lon + 360) % 360
    lon[lon > 180] -= 360

    # Setup plot
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Add map features
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')

    # Define zoom regions
    regions = {
        'global': [-180, 180, -90, 90],
        'gulf_stream': [-85, -30, 10, 50],
        'caribbean': [-90, -55, 5, 30],
        'indian_ocean': [20, 120, -40, 30]
    }

    # Set extent
    extent = regions.get(region.lower(), regions['global'])
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    # Plot points
    ax.scatter(lon, lat, s=0.5, color='blue', transform=ccrs.PlateCarree())

    plt.title(f"Drifter Positions - {region.replace('_', ' ').title()}")
    plt.tight_layout()
    #plt.show()

    os.makedirs(os.path.dirname("/Users/stewarta/Desktop/oscarpy/diagnostics/drifters/drifter_loc_map.png"), exist_ok=True)
    plt.savefig("/Users/stewarta/Desktop/oscarpy/diagnostics/drifters/drifter_loc_map.png", dpi=300)
    plt.close()
    
    return


def plot_drifter_currents(ds, region='global'):
    """
    Plot drifter positions as points colored by speed (no arrows).
    """
    lat = ds['latitude'].values
    lon = ds['longitude'].values
    u = ds['ve'].values
    v = ds['vn'].values

    # Normalize longitude to [-180, 180]
    lon = (lon + 360) % 360
    lon[lon > 180] -= 360

    # Compute speed magnitude
    speed = np.sqrt(u**2 + v**2)

    # Define regions
    regions = {
        'global': [-180, 180, -90, 90],
        'gulf_stream': [-85, -30, 10, 50],
        'caribbean': [-90, -55, 5, 30],
        'indian_ocean': [20, 120, -40, 30]
    }
    extent = regions.get(region.lower(), regions['global'])

    # Plot
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(ccrs.cartopy.feature.BORDERS, linestyle=':')
    ax.add_feature(ccrs.cartopy.feature.LAND, facecolor='lightgray')

    sc = ax.scatter(lon, lat, c=speed, cmap='turbo', s=10, transform=ccrs.PlateCarree())
    plt.colorbar(sc, ax=ax, label='Current Speed (m/s)')
    plt.title(f"Drifter Positions and Magnitude - {region.replace('_', ' ').title()}")
    plt.tight_layout()
    #plt.show()

    os.makedirs(os.path.dirname("/Users/stewarta/Desktop/oscarpy/diagnostics/drifters/currents.png"), exist_ok=True)
    plt.savefig("/Users/stewarta/Desktop/oscarpy/diagnostics/drifters/currents.png", dpi=300)
    plt.close()
   


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


def plot_drifter_component(binned_ds, component='ve',region='global', cmap='turbo'):
    """
    Plot a heatmap of the selected vector component from gridded drifter data with coastlines.
    
    Parameters:
        binned_ds (xarray.Dataset): Gridded dataset with 've' or 'vn'
        component (str): Component to plot ('ve' or 'vn')
        cmap (str): Colormap to use (default: 'turbo')
        region (str): Geographic region to plot (e.g. 'global', 'gulf_stream', 'caribbean')
    """

    if component not in binned_ds:
        raise ValueError(f"Component '{component}' not found in dataset.")

    lon = binned_ds['longitude']
    lat = binned_ds['latitude']
    data = binned_ds[component]

    # Predefined zoom regions
    regions = {
        'global': [-180, 180, -90, 90],
        'gulf_stream': [-85, -30, 10, 50],
        'caribbean': [-90, -55, 5, 30],
        'indian_ocean': [20, 120, -40, 30],
        'north_atlantic': [-80, 10, 0, 60],
        'custom': None  # use for future flexibility
    }

    fig = plt.figure(figsize=(14, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Add map features
    ax.coastlines()
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.BORDERS, linestyle=':')

    if region in regions and regions[region]:
        ax.set_extent(regions[region], crs=ccrs.PlateCarree())

    mesh = ax.pcolormesh(
        lon,
        lat,
        data,
        cmap=cmap,
        shading='auto',
        transform=ccrs.PlateCarree()
    )

    plt.colorbar(mesh, ax=ax, orientation='vertical', label=f"{component.upper()} (m/s)")
    ax.set_title(f"{component.upper()} Component of Drifter Velocity")
    plt.tight_layout()
    plt.show()

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
                'RMSE': np.nan,
                'Bias': np.nan,
                'Pearson': np.nan,
                'N_points': 0
            }
            continue

        # Compute metrics
        rmse = mean_squared_error(obs_clean, model_clean)
        r2 = r2_score(obs_clean, model_clean)
        bias = np.mean(model_clean - obs_clean)
        corr, _ = pearsonr(model_clean, obs_clean)

        results[label] = {
            'R2': r2,
            'RMSE': rmse,
            'Bias': bias,
            'Pearson': corr,
            'N_points': len(obs_clean)
        }

    return results

def plot_velocity_comparison(obs, model, component='u', limits=(-1, 1)):
    """
    Scatter plot comparing drifter (obs) vs model (OSCAR) velocity component.
    Computes and displays R², RMSE, Bias, and Pearson correlation.

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
    rmse = np.sqrt(mean_squared_error(obs_clean, model_clean))
    bias = np.mean(model_clean - obs_clean)
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
        f"RMSE: {rmse:.3f} m/s\n"
        f"Bias: {bias:.3f} m/s\n"
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
