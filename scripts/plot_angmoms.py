from DFTools import plotter, fileprc, wrap
import numpy as np

centres = [np.asarray([0.0, 0.0, 0.1114785], np.float64), np.asarray([0.0, 0.0, 0.1114785], np.float64), np.asarray([0.0, 0.0, 0.1114785], np.float64), np.asarray([0.0, 0.0, 0.1114785], np.float64), np.asarray([0.0, 0.0, 5.4186777], np.float64), np.asarray([0.0, 0.0, 5.4186777], np.float64), np.asarray([0.0, 0.0, 5.4186777], np.float64), np.asarray([0.0, 0.0, 5.4186777], np.float64), np.asarray([0.0, 0.0, 0.1114785], np.float64), np.asarray([0.0, 0.0, 5.4186777], np.float64)]
k_points = fileprc.extract_ks_from_file("/Users/ponet/Documents/Fysica/PhD/GeTe/fullrel/GeTe_bands_notrot.out")
kxs = np.hstack((np.linspace(k_points[95, 0], 0.5, 500), np.linspace(0.5, k_points[105, 0], 500)))
kys = np.hstack((np.linspace(k_points[95, 1], 0.5, 500), np.linspace(0.5, k_points[105, 1], 500)))
kzs = np.hstack((np.linspace(k_points[95, 2], 0.5, 500), np.linspace(0.5, k_points[105, 2], 500)))
k_points = [kxs, kys, kzs]
ticks = [0,499, 999]
ticklabels = ['<- A', 'Z', 'U ->']

plot = plotter.Plotter(w=15)
plot.plot_wan_angmom("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest/", 4, centres, l='Lx', k_points=k_points, label="First valence band", linewidth=2, c='r')
plot.plot_wan_angmom("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest/", 3, centres, l='Lx', k_points=k_points, label="Second valence band", linewidth=2, c='b')
plot.plot_wan_angmom("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest/", 2, centres, l='Lx', k_points=k_points, label="Third valence band", linewidth=2, c='g')

plot.add_plot((1, 2))

plot.plot_wan_angmom("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest/", 4, centres, l='Ly', k_points=k_points, axes_index=(1, 2), linewidth=2, c='r')
plot.plot_wan_angmom("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest/", 3, centres, l='Ly', k_points=k_points, axes_index=(1, 2), linewidth=2, c='b')
plot.plot_wan_angmom("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest/", 2, centres, l='Ly', k_points=k_points, axes_index=(1, 2), linewidth=2, c='g')

plot.show_legend(fontsize=10)
plot.title("Orbital Angular Momentum without SOC", size=30)
plot.stretch_plot(-0.05, "Top", axes_index=(1, 1))
plot.stretch_plot(-0.05, "Top", axes_index=(1, 2))
plot.set_xaxis_margin(0)
plot.set_xaxis_margin(0, axes_index=(1, 2))
plot.set_title("Lx", axes_index=(1, 1), size=20)
plot.set_title("Ly", axes_index=(1, 2), size=20)
plot.set_xaxis_ticks(ticks, ticklabels, fontsize=15)
plot.set_xaxis_ticks(ticks, ticklabels, fontsize=15, axes_index=(1, 2))
plot.set_yaxis_tick_params(labelsize=15)
plot.set_yaxis_tick_params(labelsize=15, axes_index=(1, 2))
plot.set_xaxis_grid((1, 1), c='#444444',linewidth=2,linestyle = '-')
plot.set_xaxis_grid((1, 2), c='#444444',linewidth=2,linestyle = '-')
plot.stretch_plot(-0.04, "t", axes_index=(1, 1))
plot.stretch_plot(-0.04, "t", axes_index=(1, 2))
plot.save("/Users/ponet/Documents/Fysica/PhD/GeTe/eigenLyLx.png")
