#!/Users/ponet/tools/anaconda3/envs/PhD/bin/python

import DFTools.plotter as plotter

if __name__ == "__main__":
    import sys

    sysarg = sys.argv
    plot = plotter.Plotter()
    erange = []
    erange.append(float(sysarg[2]))
    erange.append(float(sysarg[3]))
    plot.plot_bands(str(sysarg[1]),energy_range=erange)
    plot.save(sysarg[4])