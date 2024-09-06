"""
Plot scalars.

Usage:
    plot_spheroid.py <file1> <file2> <file3> <file4>

"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import pathlib
from docopt import docopt
from scipy.stats import linregress

# Parameters
Ekin = 'E_kin'
Edis = 'E_dis'
figsize = (8, 5)
log_scale = False

# Plot
fig = plt.figure(figsize=figsize)
args = docopt(__doc__)

with h5py.File(args['<file1>'], mode='r') as file:
    t = np.array(file['scales/sim_time'])
    dset_kin = file['tasks'][Ekin]
    dset_dis = file['tasks'][Edis]
    avg_Ekin1 = np.average(dset_kin[:].ravel())
    avg_Edis1 = np.average(dset_dis[:].ravel())
    plt.scatter(np.log(avg_Ekin1), np.log(avg_Edis1), label='Ek=1.6e-5')

with h5py.File(args['<file2>'], mode='r') as file:
    t = np.array(file['scales/sim_time'])
    dset_kin = file['tasks'][Ekin]
    dset_dis = file['tasks'][Edis]
    avg_Ekin2 = np.average(dset_kin[:].ravel())
    avg_Edis2 = np.average(dset_dis[:].ravel())
    plt.scatter(np.log(avg_Ekin2), np.log(avg_Edis2), label='Ek=3.2e-6')
   
with h5py.File(args['<file3>'], mode='r') as file:
    t = np.array(file['scales/sim_time'])
    dset_kin = file['tasks'][Ekin]
    dset_dis = file['tasks'][Edis]
    avg_Ekin3 = np.average(dset_kin[:].ravel())
    avg_Edis3 = np.average(dset_dis[:].ravel())
    plt.scatter(np.log(avg_Ekin3), np.log(avg_Edis3), label='Ek=5e-6')

with h5py.File(args['<file4>'], mode='r') as file:
    t = np.array(file['scales/sim_time'])
    dset_kin = file['tasks'][Ekin]
    dset_dis = file['tasks'][Edis]
    avg_Ekin4 = np.average(dset_kin[:].ravel())
    avg_Edis4 = np.average(dset_dis[:].ravel())
    plt.scatter(np.log(avg_Ekin4), np.log(avg_Edis4), label='Ek=8.7e-6')
'''
with h5py.File(args['<file5>'], mode='r') as file:
    t = np.array(file['scales/sim_time'])
    dset_kin = file['tasks'][Ekin]
    dset_dis = file['tasks'][Edis]
    avg_Ekin5 = np.average(dset_kin[:].ravel())
    avg_Edis5 = np.average(dset_dis[:].ravel())
    plt.scatter(np.log(avg_Ekin5), np.log(avg_Edis5), label='Ek=3.2e-6')    

with h5py.File(args['<file6>'], mode='r') as file:
    t = np.array(file['scales/sim_time'])
    dset_kin = file['tasks'][Ekin]
    dset_dis = file['tasks'][Edis]
    avg_Ekin6 = np.average(dset_kin[:].ravel())
    avg_Edis6 = np.average(dset_dis[:].ravel())
    plt.scatter(np.log(avg_Ekin6), np.log(avg_Edis6), label='Ek=5e-6')

with h5py.File(args['<file7>'], mode='r') as file:
    t = np.array(file['scales/sim_time'])
    dset_kin = file['tasks'][Ekin]
    dset_dis = file['tasks'][Edis]
    avg_Ekin7 = np.average(dset_kin[:].ravel())
    avg_Edis7 = np.average(dset_dis[:].ravel())
    plt.scatter(np.log(avg_Ekin7), np.log(avg_Edis7), label='Ek=8.7e-4')
    
with h5py.File(args['<file8>'], mode='r') as file:
    t = np.array(file['scales/sim_time'])
    dset_kin = file['tasks'][Ekin]
    dset_dis = file['tasks'][Edis]
    avg_Ekin8 = np.average(dset_kin[:].ravel())
    avg_Edis8 = np.average(dset_dis[:].ravel())
    plt.scatter(np.log(avg_Ekin8), np.log(avg_Edis8), label='Ek=8.7e-6')

with h5py.File(args['<file9>'], mode='r') as file:
    t = np.array(file['scales/sim_time'])
    dset_kin = file['tasks'][Ekin]
    dset_dis = file['tasks'][Edis]
    avg_Ekin9 = np.average(dset_kin[:].ravel())
    avg_Edis9 = np.average(dset_dis[:].ravel())
    plt.scatter(np.log(avg_Ekin9), np.log(avg_Edis9), label='Ek=9.7e-5')
'''
Edis = [avg_Edis1, avg_Edis2, avg_Edis3, avg_Edis4] 
Ekin = [avg_Ekin1, avg_Ekin2, avg_Ekin3, avg_Ekin4] 
slope, intercept, _, _, _ = linregress(np.log(Ekin), np.log(Edis))
fit_line = slope*np.log(Ekin) + intercept
constant = np.exp(intercept)
plt.plot(np.log(Ekin), fit_line, color = 'red', label=f'Slope: {slope:.3f}, Intercept: {constant:.3f}')
plt.xlabel('log(avg_Ekin)')
plt.ylabel('log(avg_Edis)')
if log_scale:
    plt.yscale('log')
    plt.xscale('log')
plt.legend()
plt.tight_layout()
plt.savefig('Ekin_vs_Edis_fit_4.pdf')
plt.show()

