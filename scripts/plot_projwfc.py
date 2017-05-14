#!/Users/ponet/tools/anaconda3/envs/PhD/bin/python
import DFTools.plotter as plotter

if __name__ == "__main__":
    import sys

    sysarg = sys.argv
    plot = plotter.Plotter()
    if len(sysarg)>4:
        energy_range=[sysarg[3],sysarg[4]]
        save_file = sysarg[5]
    else:
        energy_range=None
        save_file=sysarg[3]
    plot.plot_k_pdos(str(sysarg[1]),int(sysarg[2]),energy_range=None)
    plot.save(save_file)