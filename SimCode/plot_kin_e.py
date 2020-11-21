# import necessary packages
import yt
from yt.funcs import mylog as ytlog
import matplotlib.pyplot as plt
import numpy as np
import glob # regex
import amrvac_pytools as apt

data_dir = "../Test FD/"
basename = "kh_2d_fd(2)"
out_dir  = "./"

# second attempt
def plot_kin_e():
    # set ytlog.setLevel(40) to suppress logging up to a certain level
    # set ytlog.setLevel(20) to enable more logging
    ytlog.setLevel(40)

    # load data as time series of individual data sets
    # would yt recognise our data as 2d? Since override_geometry does not seem to work?
    ts = yt.load(data_dir+basename+'????.dat')

    # initialise arrays
    times = []
    E_kin = []

    grid_level = 2 # maximal refinement level used

    # add field: kinetic energy per cell
    def _kinetic_energy_per_cell(field, data):
        return data["kinetic_energy_density"] * data["cell_volume"]
    yt.add_field(("gas", "kinetic_energy_per_cell"), function=_kinetic_energy_per_cell, units='g*cm**2/s**2')

    # compute kinetic energy for each time step
    for ds in ts:
        #print(ds)
        dd = ds.all_data()
        #TODO use smoothed_covering_grid() for top down regridding of data
        # or find out whether dd.max_level uses interpolation
        dd.max_level = grid_level
        times.append(ds.current_time.in_units("s"))
        E_kin.append(dd.quantities.total_quantity([("gas", "kinetic_energy_per_cell")]))
    E_kin = np.array(E_kin)

    # plot data
    plt.semilogy(times, E_kin)
    plt.ylabel('Kinetic energy ($g\cdot cm^2/s^2$)')
    plt.xlabel('Time (s)')
    plt.title('Kinetic energy')
    plt.show()
    plt.savefig(out_dir+'kinetic_energy_'+basename+'.png')


if __name__ == '__main__':
    # test plot_kin_e()
    #plot_kin_e()
    print('main accessed')

