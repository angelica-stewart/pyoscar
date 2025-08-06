import os
import subprocess
from datetime import datetime, timedelta
from copernicusmarine import subset
import cdsapi
from ..config.setup import *


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


