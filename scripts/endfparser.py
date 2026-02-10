from endf_parserpy import EndfParserFactory
import numpy as np
from pathlib import Path
from plasmapy.particles import particle_mass
from scipy.constants import physical_constants
import h5py

PATH2DATA = Path(__file__).parent.parent / 'data'
KG2EV = physical_constants['kilogram-electron volt relationship'][0]


def endf_neutron_parser(target):
    incident = 'n'

    mtar = particle_mass(target).value * KG2EV
    minc = particle_mass(incident).value * KG2EV
    mu = 1/(1/mtar + 1/minc)
    lab2cm = mu / minc

    fdir = 'endf-neutron'
    fname = f'{target.upper()}({incident.upper()},EL)'
    fpath = PATH2DATA / fdir / f'{fname}.txt'
    outpath = PATH2DATA / fdir / f'{fname}.h5'

    endf = EndfParserFactory.create().parsefile(fpath)
    F3T2_data = endf[3][2]['xstable']

    with h5py.File(outpath, 'w') as h5_file:
        metadata = h5_file.create_group('metadata')
        metadata.attrs['target'] = target
        metadata.attrs['incident'] = incident
        metadata.attrs['mass_target_eV'] = mtar
        metadata.attrs['mass_incident_eV'] = minc
        metadata.attrs['reduced_mass_eV'] = mu

        F3T2 = h5_file.create_group('F3T2')
        F3T2.attrs['info'] = 'Reaction cross section, elastic scattering'

        F3T2_E = F3T2.create_dataset('E', data=np.array(F3T2_data['E']) * lab2cm)
        F3T2_E.attrs['info'] = 'CoM kinetic energy [eV]'
        F3T2_xs = F3T2.create_dataset('xs', data=np.array(F3T2_data['xs']))
        F3T2_xs.attrs['info'] = 'Cross section [barns]'
        F3T2_int = F3T2.create_dataset('int', data=np.array(F3T2_data['INT']))
        F3T2_int.attrs['info'] = 'Cross section interpolation scheme (1: constant, 2: lin-lin, 3: lin-log, 4: log-lin, 5: log-log)'
        F3T2_nbt = F3T2.create_dataset('NBT', data=np.array(F3T2_data['NBT']))
        F3T2_nbt.attrs['info'] = 'Index separating the interpolation regions (if multiple interpolation schemes are used)' 


def main():
    target = 'H-2' # Change this to parse a different target
    endf_neutron_parser(target)


if __name__ == '__main__':
    main()