import numpy as np
import xarray as xr
import glob
import os
import datetime as dt
import re
import yaml
from ..config import setup
from datetime import datetime, timedelta



def load_ds(date_str, var):

    date_obj = datetime.strptime(date_str, "%Y-%m-%d")
    next_day_obj = date_obj + timedelta(days=1)
    year = date_obj.strftime("%Y")
    month = date_obj.strftime("%m")
    day_str = date_obj.strftime("%Y%m%d")

    # Determine source path and file pattern
    if var == "ssh":

        search_dir = os.path.join(setup.ssh_src, year, month)
        pattern = setup.ssh_pattern.replace("*", f"{day_str}_*")
    elif var == "sst":
        search_dir = os.path.join(setup.sst_src, year, month)
        pattern = f"{day_str}*CMC-L4_GHRSST-SSTfnd*.nc"

    elif var == "wind":
        search_dir = os.path.join(setup.wind_src, year)
        pattern = setup.wind_pattern.replace("*", f"{year}{month}*")

    else:
        raise ValueError(f"Unsupported variable: {var}")
    
    search_path = os.path.join(search_dir, pattern)
    matches = glob.glob(search_path)

    if not matches:
        raise FileNotFoundError(f"No dataset found for {var} on {date_str} at {search_path}")

    if len(matches) > 1:
        raise RuntimeError(f"Multiple files found for {var} on {date_str}. Expected only one:\n{matches}")

    ds = xr.open_dataset(matches[0])


    if var == "ssh":
        ds = ds.rename({'adt': 'ssh'})

    elif var == "sst":
        ds = ds.rename({'analysed_sst': 'sst', 'lat': 'latitude', 'lon': 'longitude'})
        if 'mask' in ds:
            ocean_temp = np.where((ds.mask.values == 1), ds['sst'].values - 273.15, float("NaN"))
        else:
            ocean_temp = ds['sst'].values - 273.15
        ds['sst'].values = ocean_temp

    elif var == "wind":
        start_dt = np.datetime64(date_obj)
        end_dt = np.datetime64(next_day_obj)  - np.timedelta64(1, 'ns')
        ds = ds.sel(time=slice(start_dt, end_dt))
        ds = xr.merge([ds['u10'], ds['v10']])

    # === Normalize longitudes ===
    if 'longitude' in ds.coords:
        ds = ds.assign_coords(longitude=(ds.longitude % 360))
        ds = ds.sortby('longitude')

    if var in ["ssh", "sst"]:
        ds = ds[[var]]
    elif var == "wind":
        ds = ds[["u10", "v10"]]
        ds = ds.resample(time='1D').mean()

    if 'time' in ds.dims:
        if ds.sizes['time'] > 1:
            ds = ds.mean(dim='time', keep_attrs=True)

        ds['time'] = [np.datetime64(date_str)] 

    return ds

        
def write_netcdf_file(ds, path, filename):
    filename = filename + '.nc'
    print(filename)
    full_path = os.path.join(path, filename)
    print(full_path)
    os.makedirs(path, exist_ok=True)
    ds.to_netcdf(full_path)

def write_oscar_to_podaac(refDs,Ug,Uw,Ub,PODAACFILE,PODAACDIR,SSHLONGDESC,WINDLONGDESC,SSTLONGDESC,OSCARLONGDESC,OSCARSUMMARY,OSCARID,DOI, ssh_mode):
    for ii in range(len(refDs.time)):
        dsii=refDs.sel(time=refDs.time[ii])
        dt = str(dsii.time.to_numpy())
        dateDash = dt[0:10]
        year = dateDash[0:4]
        month = dateDash[5:7]
        day = dateDash[8:10]  
       
        date = year + month + day
        print("now processing output for " + dateDash)

        # Year/Month directory structure
        YRMONTHDIR = os.path.join(PODAACDIR, year, month)
        # YRDIR=PODAACDIR+'/'+year 

        # if not os.path.exists(YRDIR):
        #     print("Creating new directory: " + YRDIR)
        #     os.makedirs(YRDIR)
    
        # FILENM=YRDIR+'/'+PODAACFILE+date+'.nc'
        # if os.path.exists(FILENM):
        #     print("Removing existing file: " + FILENM)
        #     os.remove(FILENM)

        if not os.path.exists(YRMONTHDIR):
            print("Creating new directory: " + YRMONTHDIR)
            os.makedirs(YRMONTHDIR)

        FILENM = os.path.join(YRMONTHDIR, PODAACFILE + date + '.nc')

        if os.path.exists(FILENM):
            print("Removing existing file: " + FILENM)
            os.remove(FILENM)

        WRITE_METADATA(dsii,dateDash,Ug[ii,:,:],Uw[ii,:,:],Ub[ii,:,:],FILENM,SSHLONGDESC,WINDLONGDESC,SSTLONGDESC,OSCARLONGDESC,OSCARSUMMARY,OSCARID,DOI, ssh_mode)
    
        #display('REMOVING md5sum checksum from WRITE_OSCAR_PODAAC')
        # %eval(['! md5sum ',FILENM,' > ',FILENM,'.md5']) 
def WRITE_METADATA(refDs,dateDash,Ug,Uw,Ub,FILENM,SSHLONGDESC,WINDLONGDESC,SSTLONGDESC,OSCARLONGDESC,OSCARSUMMARY,OSCARID,DOI, ssh_mode):
    print("***********************************")
    print("writing output for " + FILENM)

    Ut = Ug+Uw+Ub

    longitude = refDs.coords['longitude'].to_numpy()
    latitude = refDs.coords['latitude'].to_numpy()
    time = np.atleast_1d(refDs.coords['time'].to_numpy())

    u = np.expand_dims(Ut.real.transpose(), 0)
    v = np.expand_dims(Ut.imag.transpose(), 0)
    ug = np.expand_dims(Ug.real.transpose(), 0)
    vg = np.expand_dims(Ug.imag.transpose(), 0)
    uw = np.expand_dims(Uw.real.transpose(), 0)
    vw = np.expand_dims(Uw.imag.transpose(), 0)
    ub = np.expand_dims(Ub.real.transpose(), 0)
    vb = np.expand_dims(Ub.imag.transpose(), 0)

    coords={'time': (['time'], time), 'lon': (['longitude'], longitude), 'lat': (['latitude'], latitude)}

    dsOut = xr.Dataset({'u': (['time', 'longitude', 'latitude'], u), 'v': (['time', 'longitude', 'latitude'], v), \
        'ug' : (['time', 'longitude', 'latitude'], ug), 'vg': (['time', 'longitude', 'latitude'], vg), \
        'uw' : (['time', 'longitude', 'latitude'], uw), 'vw': (['time', 'longitude', 'latitude'], vw), \
        'ub' : (['time', 'longitude', 'latitude'], ub), 'vb': (['time', 'longitude', 'latitude'], vb), \
        }, coords=coords)

    addMetadata(dsOut,dateDash,SSHLONGDESC,WINDLONGDESC,SSTLONGDESC,OSCARLONGDESC,OSCARSUMMARY,OSCARID,DOI, ssh_mode)

    epoch = {'units':'days since 1990-01-01'}
    fill = {'_FillValue':-999.0}
    encoding = {'time':epoch, 'u':fill, 'v':fill, 'ug':fill, 'vg':fill, 'uw':fill, 'vw':fill, 'ub':fill, 'vb':fill}

    dsOut.to_netcdf(FILENM, encoding=encoding)
def addMetadata(dsOut,dateDash,SSHLONGDESC,WINDLONGDESC,SSTLONGDESC,OSCARLONGDESC,OSCARSUMMARY,OSCARID,DOI, ssh_mode):
    if ssh_mode == 'cmems':
        lat_res = '0.25 degree'
        lon_res = '0.25 degree'
    else:
        lat_res = '0.1 degree'
        lon_res = '0.1 degree'
    # latitude
    dsOut.lat.attrs['long_name'] = 'latitude'
    dsOut.lat.attrs['standard_name'] = 'latitude'
    dsOut.lat.attrs['units'] = 'degrees_north'
    dsOut.lat.attrs['axis'] = 'Y'
    dsOut.lat.attrs['valid_min'] = -89.75
    dsOut.lat.attrs['valid_max'] =  89.75
    dsOut.lat.attrs['bounds'] = '[-89.75,89.75]'

    # longitude
    dsOut.lon.attrs['long_name'] = 'longitude'
    dsOut.lon.attrs['standard_name'] = 'longitude'
    dsOut.lon.attrs['units'] = 'degrees_east'
    dsOut.lon.attrs['axis'] =  'X'
    dsOut.lon.attrs['valid_min'] = 0.0
    dsOut.lon.attrs['valid_max'] = 359.75
    dsOut.lon.attrs['bounds'] = '[0,359.75]'

    dsOut.time.attrs['long_name'] = 'time centered on the day'
    dsOut.time.attrs['standard_name'] = 'time'
    # these are filled by epoch
    #ncwriteatt(FILENM, 'time', 'units', 'days since 1990-1-1');
    #ncwriteatt(FILENM, 'time', 'calendar', 'julian');
    dsOut.time.attrs['axis'] = 'T'
    dsOut.time.attrs['time_bounds'] = dateDash + ' 00:00:00 to ' + dateDash + ' 23:59:59'
    dsOut.time.attrs['comment'] = 'Data is averaged over the day'

    dsOut.u.attrs['long_name'] = 'zonal total surface current'
    dsOut.u.attrs['standard_name']  = 'eastward_sea_water_velocity'

    dsOut.v.attrs['long_name'] = 'meridional total surface current'
    dsOut.v.attrs['standard_name'] = 'northward_sea_water_velocity'

    for field in ['u', 'v']:
        dsOut[field].attrs['coverage_content_type'] = 'modelResult'
        dsOut[field].attrs['comment'] = 'Velocities are an average over the top 30m of the mixed layer'
        dsOut[field].attrs['source'] = 'SSH source: '+SSHLONGDESC+' ; WIND source: '+WINDLONGDESC+' ; SST source: '+SSTLONGDESC

    dsOut.ug.attrs['long_name'] = 'zonal geostrophic surface current'
    dsOut.ug.attrs['standard_name'] = 'geostrophic_eastward_sea_water_velocity'

    dsOut.vg.attrs['long_name'] = 'meridional geostrophic surface current'
    dsOut.vg.attrs['standard_name'] = 'geostrophic_northward_sea_water_velocity'

    for field in ['ug', 'vg']:
        dsOut[field].attrs['comment'] = 'Geostrophic velocities calculated from absolute dynamic topography'
        dsOut[field].attrs['source'] = 'SSH source: '+SSHLONGDESC

    dsOut.uw.attrs['long_name'] = 'NEED UPDATE'
    dsOut.uw.attrs['standard_name']  = 'NEED UPDATE'

    dsOut.vw.attrs['long_name'] = 'NEED UPDATE'
    dsOut.vw.attrs['standard_name'] = 'NEED UPDATE'

    for field in ['ug', 'vg']:
        dsOut[field].attrs['comment'] = 'NEED UPDATE'
        dsOut[field].attrs['source'] = 'NEED UPDATE'

    dsOut.ub.attrs['long_name'] = 'NEED UPDATE'
    dsOut.ub.attrs['standard_name']  = 'NEED UPDATE'

    dsOut.vb.attrs['long_name'] = 'NEED UPDATE'
    dsOut.vb.attrs['standard_name'] = 'NEED UPDATE'

    for field in ['ug', 'vg']:
        dsOut[field].attrs['comment'] = 'NEED UPDATE'
        dsOut[field].attrs['source'] = 'NEED  UPDATE'
    
    for field in ['u', 'v', 'ug', 'vg', 'uw', 'vw', 'ub', 'vb']:
        dsOut[field].attrs['units'] = 'm s-1'
        dsOut[field].attrs['valid_min'] = -3.0
        dsOut[field].attrs['valid_max'] = 3.0
        dsOut[field].attrs['depth'] = '15m'
    
    # GLOBAL ATTRIBUTES
    dsOut.attrs['title'] = OSCARLONGDESC
    dsOut.attrs['summary'] = 'Global, daily, 0.25 degree geostrophic and total mixed layer currents averaged over the top 30m. ' + OSCARSUMMARY;
    dsOut.attrs['keywords'] = 'ocean currents,ocean circulation,surface currents,ekman,geostrophic'
    #dsOut.attrs['keywords_vocabulary'] = 'CF: NetCDF COARDS Climate and Forecast Standard Names'
    dsOut.attrs['Conventions'] = 'CF-1.8 Standard Names v77, ACDD-1.3, netcdf 4.7.3, hdf5 1.8.12'
    dsOut.attrs['id'] = OSCARID
    dsOut.attrs['history'] = 'OSCAR 0.25 degree daily version 2.0 replaces OSCAR third degree 5 day'
    dsOut.attrs['source'] = 'OSCAR is based on simplified physics using satellite data; SSH source: '+SSHLONGDESC+' ; WIND source: '+WINDLONGDESC+' ; SST source: '+SSTLONGDESC
    dsOut.attrs['processing_level'] = 'L4'
    #dsOut.attrs['comment'] = 'RECOMMENDED - Provide useful additional information here'
    dsOut.attrs['standard_name_vocabulary'] = 'NetCDF Climate and Forecast (CF) Metadata Convention'
    dsOut.attrs['acknowledgment'] = 'OSCAR products are supported by NASA and may be freely distributed.'
    dsOut.attrs['product_version'] = 'v2.0'
    dsOut.attrs['creator_name'] = 'Kathleen Dohan'
    dsOut.attrs['creator_email'] = 'kdohan@esr.org'
    dsOut.attrs['creator_url'] = 'www.esr.org/research/oscar/'
    dsOut.attrs['creator_type'] = 'person'
    dsOut.attrs['creator_institution'] = 'ESR'
    dsOut.attrs['institution'] = 'Earth & Space Research'
    dsOut.attrs['references'] = 'www.esr.org/research/oscar/, PO.DAAC user guide, DOI: '+DOI
    dsOut.attrs['project'] = 'Ocean Surface Current Analyses Real-time (OSCAR)'
    dsOut.attrs['program'] = 'OSCAR'
    dsOut.attrs['publisher_name'] = 'NASA Physical Oceanography Distributed Active Archive Center (PO.DAAC)'
    dsOut.attrs['publisher_email'] = 'podaac@podaac.jpl.nasa.gov'
    dsOut.attrs['publisher_url'] = 'podaac.jpl.nasa.gov'
    dsOut.attrs['publisher_type'] = 'institution'
    dsOut.attrs['publisher_institution'] = 'PO.DAAC'
    dsOut.attrs['geospatial_lat_min'] = -89.75
    dsOut.attrs['geospatial_lat_max'] = 89.75
    dsOut.attrs['geospatial_lat_units'] = "degrees_north"
    dsOut.attrs['geospatial_lat_resolution'] = lat_res
    dsOut.attrs['geospatial_lon_min'] = 0.0
    dsOut.attrs['geospatial_lon_max'] = 359.75
    dsOut.attrs['geospatial_lon_units'] = "degrees_east"
    dsOut.attrs['geospatial_lon_resolution'] = lon_res
    dsOut.attrs['time_coverage_start'] = dateDash + 'T00:00:00'
    dsOut.attrs['time_coverage_end'] = dateDash + 'T23:59:59'
    


    dtNow = dt.datetime.utcnow()
    #Use ISO 8601:2004 for date and time
    dsOut.attrs['date_created'] = dtNow.strftime('%Y-%m-%d')
def get_last_processed_day(output_root, file_pattern="oscar_currents_nrt_*.nc"):
    """
    Recursively search the output directory for processed files and return the latest date found.
    """
    # Search all subdirectories
    search_path = os.path.join(output_root, "**", file_pattern)
    processed_files = sorted(glob.glob(search_path, recursive=True))

    if not processed_files:
        return None

    # Extract dates from filenames
    dates = []
    for file in processed_files:
        base_name = os.path.basename(file)
        try:
            # Extract the YYYYMMDD part
            date_str = base_name.split("_")[-1].replace(".nc", "")
            date = np.datetime64(f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:]}")
            dates.append(date)
        except Exception as e:
            print(f"Warning: could not parse date from {base_name}: {e}")

    if not dates:
        return None

    return max(dates)  # Return the latest date
def get_existing_dates(output_root, file_pattern="oscar_currents_nrt_*.nc"):
    """
    Recursively get a set of all dates for which output files already exist.
    """
    search_path = os.path.join(output_root, "**", file_pattern)
    processed_files = sorted(glob.glob(search_path, recursive=True))

    dates = set()
    for file in processed_files:
        base_name = os.path.basename(file)
        try:
            date_str = base_name.split("_")[-1].replace(".nc", "")
            date = np.datetime64(f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:]}")
            dates.add(date)
        except Exception as e:
            print(f"Warning: could not parse date from {base_name}: {e}")
    return dates

def extract_date(path):
    filename = os.path.basename(path)
    date_str = ''.join(filter(str.isdigit, filename))
    year = date_str[0:4]
    month = date_str[4:6]
    return year, month

def get_existing_dates(output_dir, start_date, end_date):
    existing = set()
    start = np.datetime64(start_date)
    end = np.datetime64(end_date)

    current = start
    while current <= end:
        year = str(current.astype('datetime64[Y]'))[:4]
        month = f"{int(current.astype('datetime64[M]').astype(int) % 12 + 1):02}"

        year_month_path = os.path.join(output_dir, year, month)
        if os.path.exists(year_month_path):
            for root, _, files in os.walk(year_month_path):
                for f in files:
                    match = re.search(r"(\d{8})", f)
                    if match:
                        try:
                            dt = np.datetime64(f"{match.group(1)[:4]}-{match.group(1)[4:6]}-{match.group(1)[6:]}")
                            existing.add(dt)
                        except Exception:
                            pass
        current = (current.astype('datetime64[M]') + 1).astype('datetime64[D]')


    return existing

def get_dates_to_process(start_date, end_date, output_dir, override=False):
    """
    Returns a list of np.datetime64 dates to process from start_date to end_date (inclusive).
    Skips already processed dates unless override=True.
    """
    all_dates = np.arange(np.datetime64(start_date), np.datetime64(end_date) + np.timedelta64(1, 'D'))
    existing_dates = get_existing_dates(output_dir, start_date, end_date)

    if override:
        return [str(d) for d in all_dates]

    else:
        return [str(d) for d in all_dates if d not in existing_dates]


def get_date_range(start_date: str, end_date: str):
    """
    Returns a list of date strings between start_date and end_date inclusive.
    
    Args:
        start_date (str): Start date in 'YYYY-MM-DD' format.
        end_date (str): End date in 'YYYY-MM-DD' format.
        
    Returns:
        List[str]: List of dates as strings in 'YYYY-MM-DD' format.
    """
    dates = np.arange(np.datetime64(start_date), np.datetime64(end_date) + np.timedelta64(1, 'D'))
    return [str(d) for d in dates]



import os
from datetime import datetime

import os
from datetime import datetime


def date_exists(date_str: str, ssh_mode) -> bool:
    """
    Check if a current file exists for the given date based on SSH mode.

    Args:
        date_str (str): Date string in format 'YYYY-MM-DD'.
        ssh_mode (str): Either 'cmems' or 'neurost'.

    Returns:
        bool: True if file exists, False otherwise.

    """

    with open("/Users/stewarta/Desktop/oscarpy/pyoscar/config/io_config.yaml", "r") as f:
        config = yaml.safe_load(f)

    if ssh_mode == "cmems":
        output_dir = config['podaac']['output_dir_cmems']
    else:
        output_dir = config['podaac']['output_dir_neurost']




    # Parse date string
    dt = datetime.strptime(date_str, "%Y-%m-%d")
    year = dt.strftime("%Y")
    month = dt.strftime("%m")
    yyyymmdd = dt.strftime("%Y%m%d")

    # Construct expected file path
    file_path = os.path.join(
        output_dir,
        year,
        month,
        f"oscar_currents_nrt{yyyymmdd}.nc"
    )

    return os.path.isfile(file_path)


# ds = load_ds("2020-01-01", var="wind")
# print(ds)

