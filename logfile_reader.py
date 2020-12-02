import numpy as np
import matplotlib.pyplot as plt


def file_len(fname):
    """Returns the number of lines in the given file."""
    i = -1 # For empty files
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def logfile_reader(fname):
    """Returns the column names and relevant data present in the logfile."""
    # (i.e. skips everything after the "|")
    if fname[-3:] != ".log":
        fname = fname + ".log"
    numLines = file_len(fname)
    with open(fname) as f:
        head = f.readline().split("|")[0].split()
        data = np.zeros((numLines-1,len(head)), dtype=float)
        for i,line in enumerate(f):
            entries = line.split("|")[0].split()
            entries = np.array(entries, dtype=float)
            data[i,:] = entries
    return head, data

def main():
    head, data = logfile_reader("SimCode/SimData/kh2d_robbe_tvd_roe_")
    data = data[:-1] # The last dt is 0, annoying for many calculations

    # Column indices of some interesting quantities (index "data" with it)
    # print(head)
    it          = 0
    global_time = 1
    dt          = 2

    # Boxplot of timesteps
    plt.boxplot(np.log10(data[:,dt])) # plt.yscale("log") hides ylabel :c
    plt.xticks([])
    plt.title("Comparison of timestep sizes")
    plt.ylabel(r"$\log_{10}(\Delta t$)")
    plt.show()

    # Analyse wavespeeds
    dx = 1./256
    courantpar = 0.8
    Smax = courantpar * dx / data[1:,dt] # Forward Euler: dt = courantpar*dx / Smax
    plt.plot(data[1:,global_time], Smax, "kx-")
    plt.title("Maximal wavespeeds estimated by AMRVAC")
    plt.xlabel(r"Simulation time $t$")
    plt.ylabel(r"$S_{max}$")
    plt.show()

if __name__ == "__main__":
    main()