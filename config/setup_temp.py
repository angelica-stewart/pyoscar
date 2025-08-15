import yaml


with open("/Users/stewarta/Desktop/oscarpy/pyoscar/config/io_config.yaml", "r") as f:
    config = yaml.safe_load(f)

# ==== DOWNLOAD TOGGLES ==== #

DOWNLOAD = config['download']
DOWNLOAD_MODE = DOWNLOAD['download_mode']

#SSH - CMEMS SETUP
DOWNLOAD_SSH_CMEMS = DOWNLOAD['ssh']['cmems']
DOWNLOAD_CMEMS_FILES = DOWNLOAD_SSH_CMEMS['do_download']
OVERRIDE_CMEMS_FILES_DOWNLOAD = DOWNLOAD_SSH_CMEMS['override_download']

#SSH - NEUROST SETUP
DOWNLOAD_SSH_NEUROST = DOWNLOAD['ssh']['neurost']
DOWNLOAD_NEUROST_FILES = DOWNLOAD_SSH_NEUROST['do_download']
OVERRIDE_NEUROST_FILES_DOWNLOAD = DOWNLOAD_SSH_NEUROST['override_download']

#SST - CMC SETUP
DOWNLOAD_SST_CMC = DOWNLOAD['sst']['cmc']
DOWNLOAD_CMC_FILES = DOWNLOAD_SST_CMC['do_download']
OVERRIDE_CMC_FILES_DOWNLOAD = DOWNLOAD_SST_CMC['override_download']

#Wind - ERA5 Setup
DOWNLOAD_WIND_ERA5 = DOWNLOAD['wind']['era5']
DOWNLOAD_ERA5_FILES = DOWNLOAD_WIND_ERA5['do_download']
OVERRIDE_ERA5_FILES_DOWNLOAD = DOWNLOAD_WIND_ERA5['override_download']

#Drifter Set Up
DOWNLOAD_DRIFTER = DOWNLOAD['drifter']
DOWNLOAD_DRIFTER_FILES = DOWNLOAD_DRIFTER['do_download']
OVERRIDE_DRIFTER_FILES_DOWNLOAD = DOWNLOAD_DRIFTER['override_download']






