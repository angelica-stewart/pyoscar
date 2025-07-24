from .plot_currents import extract_data, plot_currents
import xarray as xr
import numpy as np


cmems = xr.open_dataset('/Users/stewarta/Desktop/oscarpy/datasets/CURRENTS/INTERIM/PODAAC/CMEMS/2020/01/oscar_currents_nrt20200105.nc')
cmems_path ='/Users/stewarta/Desktop/oscarpy/datasets/CURRENTS/INTERIM/PODAAC/CMEMS/2020/01/oscar_currents_nrt20200105.nc'
neurost = xr.open_dataset('/Users/stewarta/Desktop/oscarpy/datasets/CURRENTS/INTERIM/PODAAC/NEUROST/2020/01/oscar_currents_nrt20200105.nc')
neurost_path = '/Users/stewarta/Desktop/oscarpy/datasets/CURRENTS/INTERIM/PODAAC/NEUROST/2020/01/oscar_currents_nrt20200105.nc'
test = xr.open_dataset('/Users/stewarta/Desktop/NASA_OSCAR/datasets_package/CURRENTS/INTERIM/PODAAC/2020/oscar_currents_interim_20200101.nc')
test_path = '/Users/stewarta/Desktop/NASA_OSCAR/datasets_package/CURRENTS/INTERIM/PODAAC/2020/oscar_currents_interim_20200101.nc'

datasets = {
    'cmems': cmems,
    'neurost': neurost,
    'test': test
}

extracted = {}


u_cmems, v_cmems, lon_cmems, lat_cmems, speed_cmems, lon2d_cmems, lat2d_cmems = extract_data(cmems)


# Extract and create meshgrid for NEUROST
u_neurost, v_neurost, lon_neurost, lat_neurost, speed_neurost, lon2d_neurost, lat2d_neurost = extract_data(neurost)

plot_currents(cmems_path, lon2d_cmems, lat2d_cmems, u_cmems, v_cmems, speed_cmems, 'CMEMS OSCAR Surface Currents 05-Jan-2020', 'cmems')
plot_currents(neurost_path, lon2d_neurost, lat2d_neurost, u_neurost, v_neurost, speed_neurost, 'NEUROST Surface Currents 05-Jan-2020', 'neurost')
#plot_currents(test_path, lon2d_test, lat2d_test, u_test, v_test, speed_test, 'NEUROST TEST Surface Currents 01-Jan-2020', 'test')
