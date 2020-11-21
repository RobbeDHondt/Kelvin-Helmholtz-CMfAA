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

def animate(basename, outname=None, ext=".mp4", quantity="omega"):
    """Produces an animated contourplot of the simulation with given name.
    
    @param basename: Base filename, given in the .par file by &filelist, base_filename.
    @param outname : Name of the movie to be produced.
    @param ext     : Output format, depends on support of your VideoWriter (usually mp4 is OK).
    @param quantity: The parameter to visualize (any of the primitive or conserved variables,
        omega, grad_omega).
    """
    # No name given ==> produce one yourself
    if outname is None:
        outname = basename + ext
    # Make sure extension is added to the filename
    if outname[-len(ext):] != ext:
        outname += ext

    # Make figure axes and colorbar
    fig, ax = plt.subplots()
    div = make_axes_locatable(ax)
    cax = div.append_axes('right', '5%', '5%')
   
    def update(iframe):
        """Updating function for the MATLAB AnimationWriter."""
        ax.clear()

        # Extract the data for the current frame
        filename = data_dir + basename + f"{iframe:04d}.dat"
        ds = apt.load_datfile(filename)
        ad = ds.load_all_data(regriddir=data_dir+"regridded_data/")
        xx = ds.get_coordinate_arrays()[0]
        yy = ds.get_coordinate_arrays()[1]
        zz = ad[quantity].transpose()
        time = ds.get_time()

        # Make a contourplot #vmin=minval, vmax=maxval, extend='both', cmap='Oranges', levels=200)
        contourplot = ax.contourf(xx, yy, zz) 
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        ax.set_xlim([0.0, 1.0]) # Unit square
        ax.set_ylim([0.0, 1.0]) # Unit square
        ax.set_title(f"{quantity} at time {time:.2f}")

        # Update the colorbar
        cax.cla()
        fig.colorbar(contourplot, cax=cax)

    nframes = len(glob.glob1(data_dir, basename+"*.dat"))
    anim = animation.FuncAnimation(fig, update, frames=nframes, repeat=False)
    anim.save(out_dir+outname)

def animate_amr(basename, outname=None, ext=".mp4", quantity="omega"):
    """Produces an animated AMR contourplot of the simulation with given name.
    
    @param basename: Base filename, given in the .par file by &filelist, base_filename.
    @param outname : Name of the movie to be produced.
    @param ext     : Output format, depends on support of your VideoWriter (usually mp4 is OK).
    @param quantity: The parameter to visualize (any of the primitive or conserved variables,
        omega, grad_omega).
    """
    # No name given ==> produce one yourself
    if outname is None:
        outname = basename + ext
    # Make sure extension is added to the filename
    if outname[-len(ext):] != ext:
        outname += ext

    # Make figure axes and colorbar
    fig, ax = plt.subplots()
   
    def update(iframe):
        """Updating function for the MATLAB AnimationWriter."""
        ax.clear()

        # Extract the data for the current frame
        filename = data_dir + basename + f"{iframe:04d}.dat"
        ds = apt.load_datfile(filename)
        time = ds.get_time()

        # Plot and set figure
        ds.amrplot(
            quantity, fig=fig, ax=ax, draw_mesh=True, mesh_linewidth=0.5, 
            mesh_color='white', mesh_linestyle='solid', mesh_opacity=0.8
        )
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        ax.set_xlim([0.0, 1.0]) # Unit square
        ax.set_ylim([0.0, 1.0]) # Unit square
        ax.set_title(f"{quantity} at time {time:.2f}")

    nframes = len(glob.glob1(data_dir, basename+"*.dat"))
    anim = animation.FuncAnimation(fig, update, frames=nframes, repeat=False)
    anim.save(out_dir+outname)

def quantities(basename):
    # =========
    # Computing
    nframes = len(glob.glob1(data_dir, basename+"*.dat"))
    time      = np.zeros(nframes)
    enstrophy = np.zeros(nframes)
    for iframe in range(nframes):
        # Load the data
        filename = data_dir + basename + f"{iframe:04d}.dat"
        ds = apt.load_datfile(filename)
        ad = ds.load_all_data(regriddir=data_dir+"regridded_data/")
        omega = ad['omega']

        # Compute quantities
        time[iframe] = ds.get_time()
        enstrophy[iframe] = 1/2 * np.linalg.norm(omega)**2

    print(enstrophy)

    # ======
    # Saving (doesn't seem to be necessary, loads very fast when regridded)
    # np.save(data_dir + basename + "time.npy", time)
    # np.save(data_dir + basename + "enstr.npy", enstrophy)

    # ========
    # Plotting
    plt.semilogy(time, enstrophy)
    plt.xlabel("Time")
    plt.ylabel("Enstrophy")
    plt.show()
    return time, enstrophy

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

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Read datfile basename from commandline
        animate(sys.argv[1])
    else:
        # animate("kh_2d_tvdmu_", "kh_2d_tvdmu.mp4")
        # animate("kh_2d_fd_", "kh_2d_fd.mp4")
        # animate("kh_2d_tvd_","test.mp4")
        # yt_test()
        # animate("test_dat_", "temp")

        # animate_amr("kh2d_bad-amr_") # 700 seconds
        # animate_amr("kh2d_good-amr_", outname="kh2d_good-amr_testlessblocks") # 1000 seconds
        # animate_amr("kh2d_good-amr_", outname="kh2d_good-amr_specialrefine") # 1370 seconds (4 levels)

        # animate_amr("kh2d_SETUP-NAME_", "kh2d_compressible_amr") # 1564 seconds (4 levels)
        # animate("kh2d_SETUP-NAME_", "kh2d_compressible")
        
        # quantities("kh2d_SETUP-NAME_")
        pass