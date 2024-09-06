"""
Plot disk outputs.

Usage:
    plot_disk.py <files>... [--output=<dir>]

Options:
    --output=<dir>  Output directory [default: ./frames_u]

"""

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dedalus.extras import plot_tools
import dedalus.public as d3
import matplotlib.colors as colors
#from dedalus.extras.plot_tools import FieldWrapper
#import plot_tools


def main(filename, start, count, output):
    """Save plot of specified tasks for given range of analysis writes."""

    # Plot settings
    tasks = ['u_3d']
    cmap = plt.cm.RdBu_r
    projection='polar'
    levels=50
    figsize = (8,8)
    savename_func = lambda write: 'write_{:06}.png'.format(write)
    title_func = lambda sim_time: 't = {:.3f}'.format(sim_time)
    dpi = 200
    func = lambda phi, r, data: (r*np.cos(phi), r*np.sin(phi), data)

    # Plotting loop
    with h5py.File(filename, mode='r') as file:
        for i in range(start, start+count):
            for n, task in enumerate(tasks): 
            	dset = file['tasks'][task]
            	phi = dset.dims[2]
            	theta = dset.dims[3]
            	r = dset.dims[4]
            	u_x = np.mean(dset[i:,0], axis=(0,3)) # averaged over z-axis and time
            	u_y = np.mean(dset[i:,1], axis=(0,3)) # axis(t, x, y, z)
            	u_z = np.mean(dset[i:,2], axis=(0,3))
            	y,x = np.meshgrid(theta[0],phi[0])
            	
            	fig, ax3 = plt.subplots(figsize=figsize, subplot_kw=dict(projection=projection))
            	#c=ax3.contourf(x, y, u_z, levels=levels, norm=colors.CenteredNorm(), cmap=cmap)
            	c=ax3.pcolormesh(x,y,u_z,cmap=cmap, rasterized=1)
            	fig.colorbar(c, ax=ax3 ,format='%.2e', shrink=0.45)
            	ax3.set_title(r'$u_{phi}$')
            	ax3.axis('off')
            # Add time title
            t1=dset.dims[0]['sim_time'][i]
            fig.suptitle(f'time from t = {t1}')
            # Save figure
            savename = savename_func(file['scales/write_number'][i])
            savepath = output.joinpath(savename)
            fig.savefig(str(savepath), dpi=dpi)
            fig.clear()
        plt.close(fig)
        	
        	
if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync

    args = docopt(__doc__)

    output_path = pathlib.Path(args['--output']).absolute()
    # Create output directory if needed
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_path.exists():
                output_path.mkdir()
    post.visit_writes(args['<files>'], main, output=output_path)

