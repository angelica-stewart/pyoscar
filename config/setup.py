import yaml


with open("/Users/stewarta/Desktop/oscarpy/pyoscar/config/io_config.yaml", "r") as f:
    config = yaml.safe_load(f)

do_currents = config['general']['do_currents']
do_plots = config['general']['do_plots']
do_validation = config['general']['do_validation']
do_download = config['general']['do_download']

start_date = config["general"]["start_date"]
end_date = config["general"]["end_date"]
oscar_mode = config['general']['oscar_mode']
download_mode = config['general']['download_mode']
diagnostics = config['diagnostics']
save_dir_fig = diagnostics['nrt']['interp_grad_fig'] 
filename_interp = "interpolated_{var}_{day}"
filename_gradient = "grad_{var}_{day}"
override = config['general']['override']


ssh_mode = config['general']['ssh_mode']
ssh = config['ssh'][ssh_mode]['input']
ssh_var = ssh['var']
ssh_pattern = ssh['file_pattern']
ssh_src_cmems_interim = config['ssh']['cmems']['input']['src_dir_cmems_interim']
ssh_src_cmems_final = config['ssh']['cmems']['input']['src_dir_cmems_final']
ssh_src_neurost_interim = config['ssh']['neurost']['input']['src_dir_neurost_interim']
ssh_src_neurost_final = config['ssh']['neurost']['input']['src_dir_neurost_final']
download_dates_cmems = config['ssh']['cmems']['download']['download_dates_cmems']
download_dates_neurost = config['ssh']['neurost']['download']['download_dates_neurost']
download_cmems = config['ssh']['cmems']['download']['download_cmems']
download_neurost = config['ssh']['neurost']['download']['download_neurost']

sst_mode = config['general']['sst_mode']
sst = config['sst'][sst_mode]['input']
sst_var = sst['var']
sst_pattern = sst['file_pattern']
sst_src = sst['src_dir']
download_dates_sst = config['sst'][sst_mode]['download']['download_dates_sst']
download_sst = config['sst'][sst_mode]['download']['download_sst']

wind_mode = config['general']['wind_mode']
wind = config['wind'][wind_mode]['input']
wind_var_u = wind['var_u']
wind_var_v = wind['var_v']
wind_pattern = wind['file_pattern']
wind_src = wind['src_dir']
download_dates_wind = config['wind'][wind_mode]['download']['download_dates_wind']
download_wind = config['wind'][wind_mode]['download']['download_wind']

do_eq = config['general']['do_eq']


podaacfile = config['podaac']['output_file_prefix'] + oscar_mode

ssh_long_desc = config['podaac']['description']['ssh']
wind_long_desc = config['podaac']['description']['wind']
sst_long_desc = config['podaac']['description']['sst']
oscar_long_desc = config['podaac']['description']['oscar']
oscar_summary = config['podaac']['summary']
oscar_id = config['podaac']['id']
ssh_long_desc = config['podaac']['description']['ssh']
doi = config['podaac']['doi']

if ssh_mode == "cmems":
    output_dir = config['podaac']['output_dir_cmems']
    podaacdir = config['podaac']['output_dir_cmems']
else:
    output_dir = config['podaac']['output_dir_neurost']
    podaacdir = config['podaac']['output_dir_neurost']


fig_root = config['diagnostics']['currents']['fig_root']


# plot_currents = config['plots']['currents']['do_plot']
# plot_validation = config['plots']['validation']['do_validate']
region = config['plots']['currents']['region']
dates_to_plot = config['plots']['currents']['dates_to_plot']
dates_to_validate = config['plots']['validation']['dates_to_validate']

ssh_mode_plots = config['plots']['ssh_mode']
if ssh_mode_plots == "cmems":
    search_path_plots = config['podaac']['output_dir_cmems']
else:
    search_path_plots = config['podaac']['output_dir_neurost']




search_path_cmems = config['podaac']['output_dir_cmems']
search_path_neurost = config['podaac']['output_dir_neurost']
