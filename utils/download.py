import os
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


def download_ssh_cmems(dates):
    if download_mode == "final":
        os.makedirs(ssh_src_cmems_final, exist_ok=True)
        for date in dates:
            #date_obj = datetime.strptime(date, "%Y-%m-%d")

            date_obj = date
            date_str = date_obj.strftime("%Y%m%d")

            subset(
                dataset_id="cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D",
                variables=["adt", "err_sla", "err_ugosa", "err_vgosa", "flag_ice", "sla", "tpa_correction", "ugos", "ugosa", "vgos", "vgosa"],
                minimum_longitude=-179.9375,
                maximum_longitude=179.9375,
                minimum_latitude=-89.9375,
                maximum_latitude=89.9375,
                start_datetime=f"{date_str}T00:00:00",
                end_datetime=f"{date_str}T00:00:00",
                output_directory=os.path.join(ssh_src_cmems_final, date_obj.strftime("%Y"), date_obj.strftime("%m")),
                output_filename=f"ssh_{date_str}.nc",
                file_format="netcdf"
                )
        
    if oscar_mode == "interim":
        os.makedirs(ssh_src_cmems_interim, exist_ok=True)
        for date in dates:
            date_obj = datetime.strptime(date, "%Y-%m-%d")
            date_str = date_obj.strftime("%Y%m%d")
            subset(
                dataset_id="cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D",
                variables=["adt", "err_sla", "err_ugosa", "err_vgosa", "flag_ice", "sla", "tpa_correction", "ugos", "ugosa", "vgos", "vgosa"],
                minimum_longitude=-179.9375,
                maximum_longitude=179.9375,
                minimum_latitude=-89.9375,
                maximum_latitude=89.9375,
                start_datetime=f"{date_str}T00:00:00",
                end_datetime=f"{date_str}T00:00:00",
                output_directory=os.path.join(ssh_src_cmems_interim, date_obj.strftime("%Y"), date_obj.strftime("%m")),
                output_filename=f"ssh_{date_str}.nc",
                file_format="netcdf"
                )
        return

def download_ssh_neurost(dates):
    
    os.makedirs(ssh_src_neurost_final, exist_ok=True)

    for date in dates:
        start = date
        end = date + timedelta(days=1)
        start_str = start.strftime("%Y-%m-%dT00:00:00Z")
        end_str = end.strftime("%Y-%m-%dT00:00:00Z")

        year = start.strftime("%Y")
        month = start.strftime("%m")
        target_folder = os.path.join(ssh_src_neurost_final, year, month)
        os.makedirs(target_folder, exist_ok=True)


        command = [
            "podaac-data-downloader",
            "-c", "NEUROST_SSH-SST_L4_V2024.0",
            "-d", target_folder,
            "--start-date", start_str,
            "--end-date", end_str,
            "-e", ".nc"
        ]

        
        subprocess.run(command, check=True)
        # After subprocess.run(...)
        for f in os.listdir(target_folder):
            if f.endswith(".txt"):
                os.remove(os.path.join(target_folder, f))

    return

#wind is same for both final and interim
def download_wind_era5(dates):
    client = cdsapi.Client()
    
    for date in dates:
        date_obj = date  # assume it's already a date object
        year = date_obj.strftime("%Y")
        month = date_obj.strftime("%m")
        day = date_obj.strftime("%d")
        
        target_file = os.path.join(wind_src, year, month, f"era5_{year}{month}{day}.nc")
        
        # ✅ This ensures the output folder exists before download
        os.makedirs(os.path.dirname(target_file), exist_ok=True)

        request = {
            "product_type": "reanalysis",
            "variable": ["10m_u_component_of_wind", "10m_v_component_of_wind"],
            "year": year,
            "month": month,
            "day": day,
            "time": [f"{hour:02d}:00" for hour in range(24)],
            "format": "netcdf",
        }

        
        client.retrieve("reanalysis-era5-single-levels", request, target_file)


#sst is same for final and interim only from 2016-present
def download_sst_cmc(dates):
    os.makedirs(sst_src, exist_ok = True)
    for date in dates:

        if download_mode == "final" and date < datetime(2016, 1, 1):
            c = "CMC0.2deg-CMC-L4-GLOB-v2.0"
        else:
            c = "CMC0.1deg-CMC-L4-GLOB-v3.0"
        
        start = date
        end = date + timedelta(days=1)
        start_str = start.strftime("%Y-%m-%dT00:00:00Z")
        end_str = end.strftime("%Y-%m-%dT00:00:00Z")

        date_obj = date
        year = date_obj.strftime("%Y")
        month = date_obj.strftime("%m")
        day = date_obj.strftime("%d")
        target_folder = os.path.join(sst_src, year, month)
            
        command = [
                "podaac-data-downloader",
                "-c", c,
                "-d", target_folder,
                "--start-date", start_str,
                "--end-date", end_str,
                "-e", ".nc"
            ]
        subprocess.run(command, check=True)

    return


def get_drifter_data(year_month):
    year, month, start_day, end_day = get_month_info(year_month)
    # Construct date range
    start_str = f"{year}-01-{start_day}T00:00:00Z"
    end_str = f"{year}-01-{end_day}T23:59:59Z"

    # Save location
    
    target_dir = os.path.join(drifter_src_dir, year, month)
    os.makedirs(target_dir, exist_ok=True)

    filename = f"drifter_6hour_qc_{year}_{month}.nc"
    target_path = os.path.join(target_dir, filename)

    # Skip if file exists
    if os.path.exists(target_path):
        print(f"✅ Already exists: {target_path}")
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
