
import subprocess
from datetime import datetime, timedelta
from copernicusmarine import subset
import cdsapi
from ..config.setup import *
import requests
from .data_utils import get_month_info
import calendar
from dateutil.relativedelta import relativedelta
import os, requests, calendar

import glob




# --- Helpers ----------------------------------------------------------

def _coerce_dt(d):
    """Return a datetime.datetime from either a datetime/date or 'YYYY-MM-DD' string."""
    if hasattr(d, "strftime"):
        return d
    return datetime.strptime(d, "%Y-%m-%d")


def _has_nc_for_date(folder: str, dt) -> bool:
    """
    True if there is at least one non-empty .nc file in `folder`
    whose name contains YYYYMMDD for `dt` (works for NEUROST/CMC patterns).
    """
    d8 = _coerce_dt(dt).strftime("%Y%m%d")
    pattern = os.path.join(folder, f"*{d8}*.nc")
    for p in glob.glob(pattern):
        try:
            if os.path.getsize(p) > 0:
                return True
        except OSError:
            pass
    return False


# --- CMEMS SSH --------------------------------------------------------

def download_ssh_cmems(dates):
    """
    Download daily CMEMS SSH files for each date.
    Respects skip_existing_cmems_download flag and writes 'ssh_YYYYMMDD.nc'.

    Uses:
      - oscar_mode to choose dataset + root dir (final vs interim)
      - ssh_src_cmems_final / ssh_src_cmems_interim
      - skip_existing_cmems_download
    """
    # pick dataset + root
    if oscar_mode == "final":
        dataset_id = "cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D"
        root_dir   = ssh_src_cmems_final
    elif oscar_mode == "interim":
        dataset_id = "cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D"
        root_dir   = ssh_src_cmems_interim
    else:
        raise ValueError(f"Unknown oscar_mode: {oscar_mode}")

    os.makedirs(root_dir, exist_ok=True)

    for d in dates:
        dt = _coerce_dt(d)
        y  = dt.strftime("%Y")
        m  = dt.strftime("%m")
        d8 = dt.strftime("%Y%m%d")

        out_dir  = os.path.join(root_dir, y, m)
        os.makedirs(out_dir, exist_ok=True)
        out_file = os.path.join(out_dir, f"ssh_{d8}.nc")

        # skip logic
        if skip_existing_cmems_download and os.path.exists(out_file) and os.path.getsize(out_file) > 0:
            print(f"⏭️  CMEMS skip existing: {out_file}")
            continue

        # perform subset
        subset(
            dataset_id=dataset_id,
            variables=[
                "adt","err_sla","err_ugosa","err_vgosa","flag_ice","sla",
                "tpa_correction","ugos","ugosa","vgos","vgosa"
            ],
            minimum_longitude=-179.9375,
            maximum_longitude= 179.9375,
            minimum_latitude= -89.9375,
            maximum_latitude=  89.9375,
            start_datetime=f"{d8}T00:00:00",
            end_datetime=f"{d8}T00:00:00",
            output_directory=out_dir,
            output_filename=f"ssh_{d8}.nc",
            file_format="netcdf",
        )
        print(f"✅ CMEMS saved: {out_file}")


# --- NEUROST SSH (PO.DAAC CLI) ---------------------------------------

def download_ssh_neurost(dates):
    """
    Download NEUROST SSH/SST L4 files for each date via podaac-data-downloader CLI.
    Respects skip_existing_neurost_download and keeps any provider filenames.
    """
    os.makedirs(ssh_src_neurost_final, exist_ok=True)

    for date in dates:
        start = _coerce_dt(date)
        end   = start + timedelta(days=1)
        start_str = start.strftime("%Y-%m-%dT00:00:00Z")
        end_str   = end.strftime("%Y-%m-%dT00:00:00Z")

        year  = start.strftime("%Y")
        month = start.strftime("%m")
        target_folder = os.path.join(ssh_src_neurost_final, year, month)

        # skip if any .nc for this YYYYMMDD is already present
        if skip_existing_neurost_download and _has_nc_for_date(target_folder, start):
            print(f"⏭️  NEUROST skip existing for {start.strftime('%Y-%m-%d')}: {target_folder}")
            continue

        os.makedirs(target_folder, exist_ok=True)

        cmd = [
            "podaac-data-downloader",
            "-c", "NEUROST_SSH-SST_L4_V2024.0",
            "-d", target_folder,
            "--start-date", start_str,
            "--end-date",   end_str,
            "-e", ".nc",
        ]
        subprocess.run(cmd, check=True)

        # clean up extra .txt (manifest/citation) files, if any
        for f in os.listdir(target_folder):
            if f.endswith(".txt"):
                try:
                    os.remove(os.path.join(target_folder, f))
                except OSError:
                    pass


# --- ERA5 WIND (CDS API) ---------------------------------------------

def download_wind_era5(dates):
    """
    Download daily ERA5 10m U/V wind files as single NetCDF per day.
    Respects skip_existing_era5_download and writes 'era5_YYYYMMDD.nc'.
    """
    client = cdsapi.Client()

    for date in dates:
        dt = _coerce_dt(date)
        year  = dt.strftime("%Y")
        month = dt.strftime("%m")
        day   = dt.strftime("%d")

        out_dir  = os.path.join(wind_src, year, month)
        os.makedirs(out_dir, exist_ok=True)
        out_file = os.path.join(out_dir, f"era5_{year}{month}{day}.nc")

        if skip_existing_era5_download and os.path.exists(out_file) and os.path.getsize(out_file) > 0:
            print(f"⏭️  ERA5 skip existing: {out_file}")
            continue

        request = {
            "product_type": "reanalysis",
            "variable": ["10m_u_component_of_wind", "10m_v_component_of_wind"],
            "year": year,
            "month": month,
            "day": day,
            "time": [f"{hour:02d}:00" for hour in range(24)],
            "format": "netcdf",
        }
        client.retrieve("reanalysis-era5-single-levels", request, out_file)
        print(f"✅ ERA5 saved: {out_file}")


# --- CMC SST (PO.DAAC CLI) -------------------------------------------

def download_sst_cmc(dates):
    """
    Download CMC daily SST via podaac-data-downloader CLI.
    Respects skip_existing_cmc_download. Uses collection based on download_mode + date.
    """
    os.makedirs(sst_src, exist_ok=True)

    for date in dates:
        dt = _coerce_dt(date)

        # choose collection
        if download_mode == "final" and dt < datetime(2016, 1, 1):
            collection = "CMC0.2deg-CMC-L4-GLOB-v2.0"
        else:
            collection = "CMC0.1deg-CMC-L4-GLOB-v3.0"

        start = dt
        end   = dt + timedelta(days=1)
        start_str = start.strftime("%Y-%m-%dT00:00:00Z")
        end_str   = end.strftime("%Y-%m-%dT00:00:00Z")

        year  = dt.strftime("%Y")
        month = dt.strftime("%m")
        target_folder = os.path.join(sst_src, year, month)

        if skip_existing_cmc_download and _has_nc_for_date(target_folder, start):
            print(f"⏭️  CMC SST skip existing for {dt.strftime('%Y-%m-%d')}: {target_folder}")
            continue

        os.makedirs(target_folder, exist_ok=True)

        cmd = [
            "podaac-data-downloader",
            "-c", collection,
            "-d", target_folder,
            "--start-date", start_str,
            "--end-date",   end_str,
            "-e", ".nc",
        ]
        subprocess.run(cmd, check=True)
        print(f"✅ CMC SST saved in: {target_folder}")


def get_drifter_data(year_month):
    year, month, start_day, end_day = get_month_info(year_month)
    # Construct date range
    print(month)
    start_str = f"{year}-{month}-{start_day}T00:00:00Z"
    end_str = f"{year}-{month}-{end_day}T23:59:59Z"

    # Save location
    
    target_dir = os.path.join(drifter_src_dir, year, month)
    os.makedirs(target_dir, exist_ok=True)

    filename = f"drifter_6hour_qc_{year}_{month}.nc"
    target_path = os.path.join(target_dir, filename)

    # Skip if file exists
    if os.path.exists(target_path):
        print(f"Already exists: {target_path}")
        return target_path

    # ERDDAP subset URL for that year
    base_url = "https://erddap.aoml.noaa.gov/gdp/erddap/tabledap/drifter_6hour_qc.nc"
    query = (
        "?ID,time,latitude,longitude,lon360,ve,vn,sst,err_lat,err_lon,err_sst,"
        "typedeath,deploy_date,deploy_lat,deploy_lon,start_date,start_lat,start_lon,"
        "end_date,end_lat,end_lon,drogue_lost_date"
        f"&time>={start_str}&time<={end_str}"
    )
    full_url = base_url + query

    try:
        response = requests.get(full_url)
        response.raise_for_status()
        with open(target_path, "wb") as f:
            f.write(response.content)
    except Exception as e:
        return None

    return target_path


def download_drifter_data(month_list):
    if not isinstance(month_list, (list, tuple)):
        raise ValueError("month_list must be a list or tuple of 'YYYY-MM' strings.")
    if len(month_list) not in (1, 2):
        raise ValueError("month_list must have length 1 or 2.")

    def norm(s):
        s = s.replace(" ", "")
        # Validate and normalize via datetime
        dt = datetime.strptime(s, "%Y-%m")
        return dt

    dts = [norm(s) for s in month_list]

    # Single month
    if len(dts) == 1:
        ym = dts[0].strftime("%Y-%m")
        return [get_drifter_data(ym)]

    # Range (inclusive)
    start, end = dts
    if end < start:
        raise ValueError(f"End month {end.strftime('%Y-%m')} is earlier than start month {start.strftime('%Y-%m')}.")

    results = []
    cur = start
    while cur <= end:
        ym = cur.strftime("%Y-%m")
        results.append(get_drifter_data(ym))
        cur += relativedelta(months=1)
    return results
