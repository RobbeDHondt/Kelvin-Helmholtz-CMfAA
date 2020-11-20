# import necessary packages
import yt
from yt.funcs import mylog as ytlog
import matplotlib.pyplot as plt
import numpy as np
import math
import glob # regex
#import amrvac_pytools as apt

data_dir = "C:/cygwin64/home/eigenaar/Kelvin-Helmholtz-CMfAA/SimCode/SimData/"
basename = "kh2d_test2011_"
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
    Palinstrophy = []
    Enstrophy = []

    grid_level = 2 # maximal refinement level used

    # add field: kinetic energy per cell
    # the cell-volume is important as weight, s.t. data with different mesh sizes can be compared
    # TODO yt does not recognise our data as 2d
    def _kinetic_energy_per_cell(field, data):
        return data["kinetic_energy_density"] * data["cell_volume"]
    yt.add_field(("gas", "_kinetic_energy_per_cell"), function=_kinetic_energy_per_cell, units='g*cm**2/s**2')

    # add field: kinetic energy per cell
    # the cell-volume is important as weight, s.t. data with different mesh sizes can be compared
    # TODO yt does not recognise our data as 2d
    # TODO dimensionless?
    def _quad_omega(field, data):
        return data["omega"] ** 2
    yt.add_field(("amrvac", "_quad_omega"), function=_quad_omega, units='dimensionless')

    # add field: kinetic energy per cell
    # the cell-volume is important as weight, s.t. data with different mesh sizes can be compared
    # TODO yt does not recognise our data as 2d
    # TODO dimensionless?
    def _norm_grad_omega(field, data):
        return data["grad_omega_x_"]**2 + data["grad_omega_y_"]**2
    yt.add_field(("gas", "_norm_grad_omega"), function=_norm_grad_omega, units='dimensionless')

    # compute kinetic energy for each time step
    for ds in ts:
        #print(ds.shape)
        dd = ds.all_data()
        #print(dd)
        #TODO use smoothed_covering_grid() for top down regridding of data
        # or find out whether dd.max_level uses interpolation
        dd.max_level = grid_level
        times.append(ds.current_time.in_units("s"))
        E_kin.append(dd.quantities.total_quantity([("gas", "_kinetic_energy_per_cell")]))
        Enstrophy.append(0.5*dd.quantities.total_quantity([("amrvac","_quad_omega")]))
        Palinstrophy.append(0.5*dd.quantities.total_quantity([("gas", "_norm_grad_omega")]))
    E_kin = np.array(E_kin)
    Palinstrophy = np.array(Palinstrophy)

    # plot data
    fig, axs = plt.subplots(3)
    #fig.suptitle('Vertically stacked subplots')
    #axs[0].plot(x, y)
    #axs[1].plot(x, -y)
    axs[0].plot(times, E_kin)
    axs[0].set(xlabel='Time (s)', ylabel='Kinetic energy ($g\cdot cm^2/s^2$)')

    axs[1].plot(times, Enstrophy)
    axs[1].set(xlabel='Time (s)', ylabel='Enstrophy')

    axs[2].plot(times, Palinstrophy)
    axs[2].set(xlabel='Time (s)', ylabel='Palinstrophy')

    fig.show()
    fig.savefig(out_dir+'kinetic_energy_'+basename+'.pdf')


if __name__ == '__main__':
    # test plot_kin_e()
    plot_kin_e()
