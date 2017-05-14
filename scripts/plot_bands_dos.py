#!/Users/ponet/tools/anaconda3/envs/PhD/bin/python
import DFTools.plotter as plotter

if __name__ == "__main__":
    import sys

    sysarg = sys.argv
    plot = plotter.Plotter()
    erange = []
    tmp = sysarg[2].split(',')
    erange.append(float(tmp[0]))
    erange.append(float(tmp[1]))
    plot.plot_bands(str(sysarg[1]),energy_range=erange)
    plot.overlay_pdos("{}".format(sysarg[3]),sysarg[4])
    plot.save(sysarg[5])
