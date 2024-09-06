"""
Plot scalars.

Usage:
    plot_spheroid.py <file1> <file2>

"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import pathlib
from docopt import docopt

# Parameters
tasks = ['E_kin']
figsize = (6, 4)
log_scale = True

# Plot
fig = plt.figure(figsize=figsize)
args = docopt(__doc__)
with h5py.File(args['<file1>'], mode='r') as file:
    t = np.array(file['scales/sim_time'])
    for task in tasks:
        dset = file['tasks'][task]
        plt.plot(t, dset[:].ravel(), label=task, color='C0')
with h5py.File(args['<file2>'], mode='r') as file:
    t = np.array(file['scales/sim_time'])
    for task in tasks:
        dset = file['tasks'][task]
        plt.plot(t, dset[:].ravel(), label=task, color='C0')
    
plt.xlabel('t')
if log_scale:
    plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.savefig('Ekin.pdf')
plt.show()

