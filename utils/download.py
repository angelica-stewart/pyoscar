
import subprocess
from datetime import datetime, timedelta
from copernicusmarine import subset
import cdsapi
from ..config.setup import *
import requests
# from .data_utils import get_month_info
import calendar
from dateutil.relativedelta import relativedelta
import os, requests, calendar
import glob
import shutil
from .date_helpers import get_month_info

FINAL_ID   = "cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D"
INTERIM_ID = "cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D"
MIN_BYTES    = 10_000                   # guard against tiny placeholder files


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

def _nonempty(path: str) -> bool:
    try:
        return os.path.exists(path) and os.path.getsize(path) >= MIN_BYTES
    except OSError:
        return False
# --- CMEMS SSH --------------------------------------------------------


def download_ssh_cmems(dates):

    for d in dates:
        dt = _coerce_dt(d)
        y  = dt.strftime("%Y")
        m  = dt.strftime("%m")
        d8 = dt.strftime("%Y%m%d")


        final_dir   = os.path.join(ssh_src_cmems_final,   y, m)
        interim_dir = os.path.join(ssh_src_cmems_interim, y, m)
        os.makedirs(final_dir, exist_ok=True)
        os.makedirs(interim_dir, exist_ok=True)

        final_file   = os.path.join(final_dir,   f"ssh_{d8}.nc")
        interim_file = os.path.join(interim_dir, f"ssh_{d8}.nc")


        if  OVERRIDE_DOWNLOAD and _nonempty(final_file):
            print(f"CMEMS skip (final exists):   {final_file}")
            continue
        if  OVERRIDE_DOWNLOAD and _nonempty(interim_file):
            print(f"CMEMS skip (interim exists): {interim_file}")
            continue

         # ---- Try FINAL first
        tried_final = False
        try:
            tried_final = True
            subset(
                dataset_id=FINAL_ID,
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
                output_directory=final_dir,
                output_filename=f"ssh_{d8}.nc",
                file_format="netcdf",
            )
            if _nonempty(final_file):
                print(f"CMEMS FINAL saved: {final_file}")
                continue
            else:
                # clean up tiny/empty, fall through to interim
                try: os.remove(final_file)
                except OSError: pass
                print(f"FINAL returned empty/invalid file for {d8}, trying INTERIM…")
        except Exception as e:
            print(f"FINAL not available for {d8}: {e}. Trying INTERIM…")

        # ---- Fall back to INTERIM
        try:
            subset(
                dataset_id=INTERIM_ID,
                variables=[
                    "adt","err_sla","err_ugosa","err_vgosa","flag_ice","sla","ugos","ugosa","vgos","vgosa"
                ],
                minimum_longitude=-179.9375,
                maximum_longitude= 179.9375,
                minimum_latitude= -89.9375,
                maximum_latitude=  89.9375,
                start_datetime=f"{d8}T00:00:00",
                end_datetime=f"{d8}T00:00:00",
                output_directory=interim_dir,
                output_filename=f"ssh_{d8}.nc",
                file_format="netcdf",
            )
            if _nonempty(interim_file):
                tier = "INTERIM" if tried_final else "INTERIM (direct)"
                print(f"CMEMS {tier} saved: {interim_file}")
            else:
                try: os.remove(interim_file)
                except OSError: pass
                raise RuntimeError(f"INTERIM returned empty/invalid file for {d8}.")
        except Exception as e:
            print(f"Failed to download INTERIM/FINAL for {d8}: {e}")
        

        



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
        if OVERRIDE_DOWNLOAD and _has_nc_for_date(target_folder, start):
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

        if OVERRIDE_DOWNLOAD and os.path.exists(out_file) and os.path.getsize(out_file) > 0:
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
    os.makedirs(sst_src, exist_ok=True)

    for date in dates:
        dt = _coerce_dt(date)

        # FINAL-first choice based on date, with the other as INTERIM fallback
        final_collection   = "CMC0.2deg-CMC-L4-GLOB-v2.0" if dt < datetime(2016, 1, 1) else "CMC0.1deg-CMC-L4-GLOB-v3.0"
        interim_collection = "CMC0.1deg-CMC-L4-GLOB-v3.0" if final_collection == "CMC0.2deg-CMC-L4-GLOB-v2.0" else "CMC0.2deg-CMC-L4-GLOB-v2.0"

        start = dt
        end   = dt + timedelta(days=1)
        start_str = start.strftime("%Y-%m-%dT00:00:00Z")
        end_str   = end.strftime("%Y-%m-%dT00:00:00Z")

        year  = dt.strftime("%Y")
        month = dt.strftime("%m")
        target_folder = os.path.join(sst_src, year, month)

        # skip if we already have any .nc for the date in the target folder
        if OVERRIDE_DOWNLOAD and _has_nc_for_date(target_folder, start):
            print(f"CMC SST skip existing for {dt.strftime('%Y-%m-%d')}: {target_folder}")
            continue

        os.makedirs(target_folder, exist_ok=True)

        # --- Try FINAL first ---
        try:
            cmd = [
                "podaac-data-downloader",
                "-c", final_collection,
                "-d", target_folder,
                "--start-date", start_str,
                "--end-date",   end_str,
                "-e", ".nc",
            ]
            subprocess.run(cmd, check=True)
        except Exception as e:
            print(f"FINAL CMC SST attempt failed for {dt.strftime('%Y-%m-%d')} ({final_collection}): {e}")

        # Check if FINAL produced a file; if not, try INTERIM
        if not _has_nc_for_date(target_folder, start):
            try:
                cmd = [
                    "podaac-data-downloader",
                    "-c", interim_collection,
                    "-d", target_folder,
                    "--start-date", start_str,
                    "--end-date",   end_str,
                    "-e", ".nc",
                ]
                subprocess.run(cmd, check=True)
            except Exception as e:
                print(f"INTERIM CMC SST attempt failed for {dt.strftime('%Y-%m-%d')} ({interim_collection}): {e}")

            # Final check
            if not _has_nc_for_date(target_folder, start):
                print(
                    f"CMC SST not available for {dt.strftime('%Y-%m-%d')}"
                )
            else:

                print(f"CMC SST saved in: {target_folder}")


def get_drifter_data(year_month):
    year, month, start_day, end_day = get_month_info(year_month)
    # Construct date range
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

    # ERDDAP subset URL for that year/month
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
        print(f"⚠️ Failed to download {year_month}: {e}")
        return None

    return target_path


def download_drifter_data(dates):
    """
    Given a list of daily dates (YYYY-MM-DD), download drifter data 
    grouped by unique months.
    """
    if not isinstance(dates, (list, tuple)):
        raise ValueError("dates must be a list or tuple of 'YYYY-MM-DD' strings.")

    # Normalize and extract unique months
    months = set()
    for s in dates:
        dt = datetime.strptime(s.strip(), "%Y-%m-%d")
        months.add(dt.strftime("%Y-%m"))

    # Download each unique month
    results = []
    for ym in sorted(months):
        results.append(get_drifter_data(ym))

    return results