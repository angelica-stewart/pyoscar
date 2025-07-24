import xarray as xr
import numpy as np

def compare_netcdf_files(file1: str, file2: str, var_map: dict = None, rtol=1e-5, atol=1e-8, check_attrs=False):
    """
    Compare two NetCDF files for equality in data, structure, and optionally metadata.

    Args:
        file1 (str): Path to the first NetCDF file.
        file2 (str): Path to the second NetCDF file.
        var_map (dict): Dictionary to rename variables in file2 to match file1 (e.g., {"sshx": "adtx"}).
        rtol (float): Relative tolerance for numerical comparison.
        atol (float): Absolute tolerance for numerical comparison.
        check_attrs (bool): If True, checks global and variable attributes as well.

    Returns:
        bool: True if files are equivalent (within tolerance), False otherwise.
    """
    try:
        ds1 = xr.open_dataset(file1)
        ds2 = xr.open_dataset(file2)

        # Rename variables in ds2 if mapping provided
        if var_map:
            ds2 = ds2.rename(var_map)

        # Check dimensions
        if ds1.dims != ds2.dims:
            print("Mismatch in dimensions.")
            return False

        # Check data variables
        if set(ds1.data_vars) != set(ds2.data_vars):
            print("Mismatch in data variables.")
            return False

        # Check coordinates
        if set(ds1.coords) != set(ds2.coords):
            print("Mismatch in coordinates.")
            return False

        # Check variable data
        for var in ds1.data_vars:
            if not np.allclose(ds1[var].values, ds2[var].values, rtol=rtol, atol=atol, equal_nan=True):
                print(f"Mismatch in values for variable: {var}")
                return False

        # Optionally check attributes
        if check_attrs:
            if ds1.attrs != ds2.attrs:
                print("Mismatch in global attributes.")
                return False
            for var in ds1.data_vars:
                if ds1[var].attrs != ds2[var].attrs:
                    print(f"Mismatch in attributes for variable: {var}")
                    return False

        print("✅ Files are equivalent.")
        return True

    except Exception as e:
        print(f"❌ Error comparing NetCDF files: {e}")
        return False

cmems_new = '/Users/stewarta/Desktop/oscarpy/neurost_b.nc'      # contains sshx/sshy
cmems_orig = '/Users/stewarta/Desktop/NASA_OSCAR/neurost_a.nc'  # contains adtx/adty

# Rename sshx → adtx, sshy → adty to match
var_mapping = {
    "sshx": "adtx",
    "sshy": "adty"
}

compare_netcdf_files(
    file1=cmems_orig,
    file2=cmems_new,
    var_map=var_mapping,
    check_attrs=True
)

