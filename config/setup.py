import os, yaml

here = os.path.dirname(__file__)   # /Users/stewarta/Desktop/oscarpy/pyoscar/config
config_path = os.path.join(here, "io_config.yaml")
with open(config_path, "r") as f:
    config = yaml.safe_load(f)


PLOT_CURRENTS = config['general']['plot_currents'] #change to do_check_plots
DO_VALIDATION = config['general']['do_validation']


START_DATE = config["general"]["start_date"]
END_DATE = config["general"]["end_date"]

DATADIR = os.path.abspath(os.path.join(here, "../../datasets"))
OVERRIDE_DOWNLOAD = config['general']['override_download']
OVERRIDE_CURRENT = config['general']['override_current']
filename_interp = "interpolated_{var}_{day}"
filename_gradient = "grad_{var}_{day}"
override = config['general']['override_current']


SSH_MODE = config['general']['ssh_mode']
if SSH_MODE == 'cmems':
    SSH_PATTERN = "ssh_*.nc"
elif SSH_MODE == 'neurost':
    SSH_PATTERN = "NeurOST_SSH-SST_*.nc"
else:
    SSH_PATTERN = None



ssh_src_cmems_interim = DATADIR + "/SSH/CMEMS/INTERIM/SRC/"
ssh_src_cmems_final = DATADIR + "/SSH/CMEMS/FINAL/SRC/"
ssh_src_neurost_interim = DATADIR + "/SSH/NEUROST/INTERIM/SRC/"
ssh_src_neurost_final = DATADIR + "/SSH/NEUROST/FINAL/SRC"


sst_src = DATADIR + "/SST/CMC/SRC/"



wind_src = DATADIR + "/WIND/ERA5/SRC/"
WIND_PATTERN = "ERA5_*.nc"



do_eq = config['general']['do_eq']


podaacfile = "oscar_currents_" 

ssh_long_desc = "CMEMS SSALTO/DUACS SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS_008_046 DOI: 10.48670/moi-00149"
wind_long_desc = "ECMWF ERA5 10m wind DOI: 10.24381/cds.adbb2d47"
sst_long_desc =  "CMC 0.1 deg SST V3.0 DOI: 10.5067/GHCMC-4FM03"
oscar_long_desc = "Ocean Surface Current Analyses Real-time (OSCAR) Surface Currents - Interim 0.25 Degree (Version 2.0)"
oscar_summary = "Higher quality than NRT currents, but lesser quality than final currents."
oscar_id = "OSCAR_L4_OC_INTERIM_V2.0"
doi =  "10.5067/OSCAR-25I20"

if SSH_MODE == "cmems":
    output_dir = DATADIR + "/CURRENTS/FINAL/CMEMS"
    podaacdir = DATADIR + "/CURRENTS/FINAL/CMEMS"
else:
    output_dir = DATADIR + "/CURRENTS/FINAL/NEUROST"
    podaacdir = DATADIR + "/CURRENTS/FINAL/NEUROST"


REGION = config['plot_currents']['region']
FIG_ROOT = os.path.abspath(os.path.join(here, "../../diagnostics/figures"))
VALIDATION_DIR = os.path.abspath(os.path.join(here, "../../diagnostics/validation"))


if SSH_MODE == "cmems":
    search_path_plots = DATADIR + "/CURRENTS/FINAL/CMEMS"
else:
    search_path_plots =  DATADIR + "/CURRENTS/FINAL//NEUROST"

search_path_cmems = DATADIR + "/CURRENTS/FINAL/CMEMS"
search_path_neurost = DATADIR + "/CURRENTS/FINAL/NEUROST"

DRIFTER_SRC_DIR = DATADIR + "/DRIFTERS"


currents_cmems_final = DATADIR + "/CURRENTS/FINAL/CMEMS"
currents_neurost_final = DATADIR + "/CURRENTS/FINAL/NEUROST"
currents_cmems_interim = DATADIR + "/CURRENTS/INTERIM/CMEMS"
currents_neurost_interim = DATADIR + "/CURRENTS/INTERIM/NEUROST"

explanations = [
    "explanation1",
    "explanation2",
    "explanation3",
    "explanation4",
    "explanation5",
    "explanation6",
    "explanation7",
    "explanation8",
    "explanation9",
    "explanation10",
    "explanation11"
]
