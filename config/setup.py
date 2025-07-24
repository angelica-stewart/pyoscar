import yaml


with open("/Users/stewarta/Desktop/oscarpy/pyoscar/config/io_config.yaml", "r") as f:
    config = yaml.safe_load(f)

start_date = config["general"]["start_date"]
end_date = config["general"]["end_date"]
oscar_mode = config['general']['oscar_mode']
diagnostics = config['diagnostics']
save_dir_fig = diagnostics[oscar_mode]['interp_grad_fig'] 
filename_interp = "interpolated_{var}_{day}"
filename_gradient = "grad_{var}_{day}"
override = config['general']['override']


ssh_mode = config['general']['ssh_mode']
ssh = config['ssh'][ssh_mode]['input']
ssh_var = ssh['var']
ssh_pattern = ssh['file_pattern']
ssh_src = ssh['src_dir']



sst_mode = config['general']['sst_mode']
sst = config['sst'][sst_mode]['input']
sst_var = sst['var']
sst_pattern = sst['file_pattern']
sst_src = sst['src_dir']

wind_mode = config['general']['wind_mode']
wind = config['wind'][wind_mode]['input']
wind_var_u = wind['var_u']
wind_var_v = wind['var_v']
wind_pattern = wind['file_pattern']
wind_src = wind['src_dir']

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


region = config['diagnostics']['currents']['region']
fig_root = config['diagnostics']['currents']['fig_root']