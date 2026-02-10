"""
Parse a Helios .exo output file and save to HDF5 format.

Author: Marcus Borscz
HB11 Energy
"""

from netCDF4 import Dataset
import h5py
from pathlib import Path
import numpy as np

PATH2DATA = Path(__file__).parent.parent / 'data' / 'helios-runs'


def main(fname, outname=None):
    helios_vars = {
        "time_whole": "t",
        "zone_boundaries": "r",
        "ion_density": "ni",
        "elec_density": "ne",
        "ion_temperature": "Ti",
        "elec_temperature": "Te",
        "fluid_velocity": "u",
        "mass_density": "rho",
    }

    helios_vars_info = {
        "t": "Time [sec]",
        "r": "Zone center [cm]",
        "ni": "Ion density [#/cm^3]",
        "ne": "Electron density [#/cm^3]",
        "Ti": "Ion temperature [eV]",
        "Te": "Electron temperature [eV]",
        "u": "Fluid velocity [cm/s]",
        "rho": "Mass density [g/cm^3]",
    }

    data = {}

    fpath = PATH2DATA / f"{fname}.exo"
    if not fpath.exists():
        raise FileNotFoundError(f"File {fpath} does not exist.")
    
    with Dataset(fpath, 'r') as nc_file:
        # Extracting variables
        for k, v in helios_vars.items():
            vardata = nc_file.variables[k][:]

            if np.ma.is_masked(vardata):
                print(f"Warning: Variable '{k}' contains masked values. Filling masked values with NaN.")
                vardata = vardata.filled(np.nan)
            else:
                vardata = np.ma.asarray(vardata)

            if v in ("r", "u"):
                data[v] = 0.5 * (vardata[:, :-1] + vardata[:, 1:]) # Zone centers
            else:
                data[v] = vardata

    if outname is None: outname = fname
    outpath = PATH2DATA / f"{outname}.h5"

    # Save to HDF5
    with h5py.File(outpath, 'w') as h5_file:
        for k, v in data.items():
            dset = h5_file.create_dataset(k, data=v)
            dset.attrs['info'] = helios_vars_info[k]


if __name__ == '__main__':
    fname = 'BDT_rho800Ti25a1H2.15M50'
    outname = None

    main(fname, outname)