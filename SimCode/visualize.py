# General imports
import sys
import glob # regex
import numpy as np
# Data extraction and visualization
import amrvac_pytools as apt # Installation: http://amrvac.org/md_doc_python_setup.html
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable

data_dir = "./SimData/"
out_dir  = "../Animations/"

def animate(basename="kh_2d_", outname=None):
    if outname is None:
        outname = basename + ".mp4"

    fig, ax = plt.subplots()
    div = make_axes_locatable(ax)
    cax = div.append_axes('right', '5%', '5%')
   
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
        contourplot = ax.contourf(xx, yy, rho) #vmin=minval, vmax=maxval, extend='both', cmap='Oranges', levels=200)
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        ax.set_xlim([0.0, 1.0]) # Unit square
        ax.set_ylim([0.0, 1.0]) # Unit square
        ax.set_title(f"Density at time {time:.2f}")
        # Colorbar
        cax.cla()
        fig.colorbar(contourplot, cax=cax)
        plt.draw()

    nframes = len(glob.glob1(data_dir, basename+"*.dat"))
    anim = animation.FuncAnimation(fig, update, frames=nframes, repeat=False)
    anim.save(out_dir+outname)

def yt_test():
    """Alternative way of plotting via yt. Not sure how to make animations though."""
    import yt # Installation: http://amrvac.org/md_doc_yt_usage.html
    ds = yt.load("SimData/kh_2d_0015.dat")
    pt = yt.plot_2d(ds, 'rho')
    pt.save("test")

def logfile_reader():
    """Still a work in progress."""
    with open(data_dir + "kh_2d_fd_.log") as f:
        head = f.readline()
        data = []
        for line in f:
            data.append(line.split())
        data = np.array(data)
    return head, data

# animate("kh_2d_tvdmu_", "kh_2d_tvdmu.mp4")
# animate("kh_2d_fd_", "kh_2d_fd.mp4")
# animate("kh_2d_tvd_","test.mp4")
# yt_test()

if len(sys.argv) > 1:
    animate(sys.argv[1])