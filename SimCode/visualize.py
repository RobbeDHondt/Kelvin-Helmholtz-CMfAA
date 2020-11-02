# General imports
import sys
import glob # regex
import numpy as np
# Data extraction and visualization
import amrvac_pytools as apt # Installation: http://amrvac.org/md_doc_python_setup.html
import matplotlib.pyplot as plt
import matplotlib.animation as animation

basename = "kh_2d_"
data_dir = "SimData/"
outname  = "../Animations/pyvistest.mp4"

def animate():

    fig, ax = plt.subplots()
   
    def update(iframe):
        ax.clear()

        # Extract the current frame
        filename = data_dir + basename + f"{iframe:04d}.dat"
        ds = apt.load_datfile(filename)
        ad = ds.load_all_data(regriddir=data_dir+"regridded_data/")
        xx = ds.get_coordinate_arrays()[0]
        yy = ds.get_coordinate_arrays()[1]
        rho = ad['rho']
        # omega = ad['omega'] # Do auxiliary variables not make it to the .dat file?
        time = ds.get_time()

        # Make a contourplot
        ax.contourf(xx, yy, rho) #vmin=minval, vmax=maxval, extend='both', cmap='Oranges', levels=200)
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        ax.set_xlim([0.0, 1.0]) # Unit square
        ax.set_ylim([0.0, 1.0]) # Unit square
        ax.set_title(f"Density at time {time:.2f}")
        plt.draw()

    nframes = len(glob.glob1(data_dir, basename+"*.dat"))
    anim = animation.FuncAnimation(fig, update, frames=nframes, repeat=False)
    anim.save(outname)

def yt_test():
    """Alternative way of plotting via yt. Not sure how to make animations though."""
    import yt # Installation: http://amrvac.org/md_doc_yt_usage.html
    ds = yt.load("SimData/kh_2d_0015.dat")
    pt = yt.plot_2d(ds, 'rho')
    pt.save("test")

animate()
# yt_test()