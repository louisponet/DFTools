import plotter
import numpy as np

if __name__ == "__main__":
    import sys

    sysarg = sys.argv
    plot = plotter.Plotter()
    erange = []
    tmp = sysarg[2].split(',')
    erange.append(float(tmp[0]))
    erange.append(float(tmp[1]))
    plot.plot_bands(str(sysarg[1]),energy_range=erange)
    tmp_1 = [sysarg[4].split(";")[1].split(",")[0],sysarg[4].split(";")[1].split(",")[1]]
    tmp_2 = [sysarg[4].split(";")[1].split(",")[0],sysarg[4].split(";")[1].split(",")[2]]
    kys_list = np.arange(sysarg[4].split(";")[1].split(",")[0],sysarg[4].split(";")[1].split(",")[1],sysarg[4].split(";")[1].split(",")[2])
    kzs_list = np.arange(sysarg[4].split(";")[2].split(",")[0],sysarg[4].split(";")[2].split(",")[1],sysarg[4].split(";")[2].split(",")[2])
    plot.save(sysarg[5])

