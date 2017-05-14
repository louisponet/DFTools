#%%
from DFTools import plotter,fileprc,wrap
import numpy as np

# centres = [np.asarray([0.0,0.0,0.1114785],np.float64),np.asarray([0.0,0.0,0.1114785],np.float64),np.asarray([0.0,0.0,0.1114785],np.float64),np.asarray([0.0,0.0,0.1114785],np.float64),np.asarray([0.0,0.0,5.4186777],np.float64),np.asarray([0.0,0.0,5.4186777],np.float64),np.asarray([0.0,0.0,5.4186777],np.float64),np.asarray([0.0,0.0,5.4186777],np.float64),np.asarray([0.0,0.0,0.1114785],np.float64),np.asarray([0.0,0.0,5.4186777],np.float64)]
# k_points = fileprc.extract_ks_from_file("/Users/ponet/Documents/Fysica/PhD/GeTe/fullrel/GeTe_bands_notrot.out")
# kxs = np.hstack((np.linspace(k_points[95,0],0.5,500),np.linspace(0.5,k_points[105,0],500)))
# kys = np.hstack((np.linspace(k_points[95,1],0.5,500),np.linspace(0.5,k_points[105,1],500)))
# kzs = np.hstack((np.linspace(k_points[95,2],0.5,500),np.linspace(0.5,k_points[105,2],500)))
# k_points = [kxs,kys,kzs]

# ls = [-0.1,-0.1,-0.1,-0.1,-0.32,-0.32,-0.32,-0.32]
# ticks = [0,499,999]
# tick_labels=["<- A","Z","U ->"]

#---------------3rd valence band----------------------------#

# plot = plotter.Plotter(w=15)
# plot.settings["h_spacing"] = 0.09

# plot.plot_wan_band_dip_SOC("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest",4,ls,centres,k_points=k_points,label="",c='r',linewidth=2)
# plot.plot_wan_band_dip_SOC("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest",5,ls,centres,k_points=k_points,label="",c='b',linewidth=2)

# plot.plot_wan_band_eigval_SOC("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest",4,ls,centres,k_points=k_points,axes_index=(1,2),c=(1,0,0,0.4),linewidth=2)
# plot.plot_wan_band_eigval_SOC("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest",5,ls,centres,k_points=k_points,axes_index=(1,2),c=(0,0,1,0.4),linewidth=2)

# dipoles4_ = []
# dipoles5_ = []
# for dip in plot.data_storage.dipoles:
#     dipoles4_.append(dip[4][2].real)
#     dipoles5_.append(dip[5][2].real)

# mid = dipoles4_[499]

# dipoles4 = []
# dipoles5 = []

# for dip4,dip5 in zip(dipoles4_,dipoles5_):
#     dipoles4.append(dip4-mid)
#     dipoles5.append(dip5 - mid)

# ct = 15.058788*0.7/15

# energy_mid = plot.data_storage.eigenvalues[499][4]
# energy4 = []
# energy5 = []

# for dip4,dip5 in zip(dipoles4,dipoles5):
#     energy4.append(energy_mid+dip4*ct)
#     energy5.append(energy_mid + dip5 * ct)

# plot.plot_data(energy4,axes_index=(1,2),linewidth=2,c='r')
# plot.plot_data(energy5,axes_index=(1,2),linewidth=2,c='b')

# plot.set_data_label("Interpolated Eigenvalues",0,axes_index=(1,2))
# plot.set_data_label("Electrostatic energy",2,axes_index=(1,2))

# plot.show_legend(axes_index=(1,2),indices=[0,2],fontsize=15)

# plot.set_xaxis_margin(0)
# plot.set_xaxis_ticks(ticks,tick_labels)
# plot.set_xaxis_grid(c='#444444')
# plot.set_yaxis_label("Rz (Angstrom)",size=20)
# plot.set_title("Band center",size=20)

# plot.set_xaxis_margin(0,axes_index=(1,2))
# plot.set_xaxis_ticks(ticks,tick_labels,axes_index=(1,2))
# plot.set_xaxis_grid(axes_index=(1,2),c='#444444')
# plot.set_yaxis_label("energy (eV)",axes_index=(1,2),size=20)
# plot.set_title("Band dispersion",axes_index=(1,2),size=20)

# plot.title("Band center vs dispersion (3rd v-band)",size=30)

# plot.stretch_plot(-0.07,"top")
# plot.stretch_plot(-0.07,"top",axes_index=(1,2))
# plot.move_plot(0.01,"right")
# plot.move_plot(0.01,"right",axes_index=(1,2))
# plot.save("/Users/ponet/Documents/Fysica/PhD/GeTe/dip_tests/dipoles45_2ndneighbour1.png")

#---------------2nd valence band----------------------------#

# plot = plotter.Plotter(w=15)
# plot.settings["h_spacing"] = 0.09
# plot.plot_wannier_dipoles("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest","","",6,2,ls,centres,k_points=k_points,label="",c='r',linewidth=2)
# plot.plot_wannier_dipoles("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest","","",7,2,ls,centres,k_points=k_points,label="",c='b',linewidth=2)
#
# plot.plot_wannier_eigenvalues("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest","",6,ls,centres,k_points=k_points,axes_index=(1,2),c=(1,0,0,0.4),linewidth=2)
# plot.plot_wannier_eigenvalues("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest","",7,ls,centres,k_points=k_points,axes_index=(1,2),c=(0,0,1,0.4),linewidth=2)
#
# dipoles6_ = []
# dipoles7_ = []
# for dip in plot.data_storage.dipoles:
#     dipoles6_.append(dip[6][2].real)
#     dipoles7_.append(dip[7][2].real)
#
# mid = dipoles6_[499]
#
# dipoles6 = []
# dipoles7 = []
#
# for dip6,dip7 in zip(dipoles6_,dipoles7_):
#     dipoles6.append(dip6-mid)
#     dipoles7.append(dip7 - mid)
#
# ct = 15.058788*0.7/29
#
# energy_mid = plot.data_storage.eigenvalues[499][6]
# energy6 = []
# energy7 = []
#
# for dip6,dip7 in zip(dipoles6,dipoles7):
#     energy6.append(energy_mid+dip6*ct)
#     energy7.append(energy_mid + dip7 * ct)
#
# plot.plot_set(energy6,axes_index=(1,2),linewidth=2,c='r')
# plot.plot_set(energy7,axes_index=(1,2),linewidth=2,c='b')
#
# plot.set_bands_label("Interpolated Eigenvalues",0,axes_index=(1,2))
# plot.set_bands_label("Electrostatic energy",2,axes_index=(1,2))
#
# plot.show_legend(axes_index=(1,2),indices=[0,2],fontsize=15)
#
# plot.set_xaxis_margin(0)
# plot.set_xaxis_ticks(ticks,tick_labels)
# plot.set_xaxis_grid(c='#444444')
# plot.set_yaxis_label("Rz (Angstrom)",size=20)
# plot.set_title("Band center",size=20)
#
# plot.set_xaxis_margin(0,axes_index=(1,2))
# plot.set_xaxis_ticks(ticks,tick_labels,axes_index=(1,2))
# plot.set_xaxis_grid(axes_index=(1,2),c='#444444')
# plot.set_yaxis_label("energy (eV)",axes_index=(1,2),size=20)
# plot.set_title("Band dispersion",axes_index=(1,2),size=20)
#
# plot.title("Band center vs dispersion (2nd v-band)",size=30)
#
# plot.stretch_plot(-0.07,"top")
# plot.stretch_plot(-0.07,"top",axes_index=(1,2))
# plot.move_plot(0.01,"right")
# plot.move_plot(0.01,"right",axes_index=(1,2))
# plot.save("/Users/ponet/Documents/Fysica/PhD/GeTe/dipoles67.png")

#---------------1st valence band----------------------------#

# plot = plotter.Plotter(w=15)
# plot.settings["h_spacing"] = 0.09
# plot.plot_wan_band_dip_SOC("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest",8,ls,centres,k_points=k_points,label="",c='r',linewidth=2)
# plot.plot_wan_band_dip_SOC("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest",9,ls,centres,k_points=k_points,label="",c='b',linewidth=2)

# plot.plot_wan_band_eigval_SOC("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest",8,ls,centres,k_points=k_points,axes_index=(1,2),c=(1,0,0,0.4),linewidth=2)
# plot.plot_wan_band_eigval_SOC("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest",9,ls,centres,k_points=k_points,axes_index=(1,2),c=(0,0,1,0.4),linewidth=2)

# dipoles8_ = []
# dipoles9_ = []
# for dip in plot.data_storage.dipoles:
#     dipoles8_.append(dip[8][2].real)
#     dipoles9_.append(dip[9][2].real)

# mid = dipoles8_[499]

# dipoles8 = []
# dipoles9 = []

# for dip8,dip9 in zip(dipoles8_,dipoles9_):
#     dipoles8.append(dip8-mid)
#     dipoles9.append(dip9 - mid)

# ct = 15.058788*0.7/10

# energy_mid = plot.data_storage.eigenvalues[499][8]
# energy8 = []
# energy9 = []

# for dip8,dip9 in zip(dipoles8,dipoles9):
#     energy8.append(energy_mid+dip8*ct)
#     energy9.append(energy_mid + dip9 * ct)

# plot.plot_data(energy8,axes_index=(1,2),linewidth=2,c='r')
# plot.plot_data(energy9,axes_index=(1,2),linewidth=2,c='b')

# plot.set_data_label("Interpolated Eigenvalues",0,axes_index=(1,2))
# plot.set_data_label("Electrostatic energy",2,axes_index=(1,2))

# plot.show_legend(axes_index=(1,2),indices=[0,2],fontsize=15)

# plot.set_xaxis_margin(0)
# plot.set_xaxis_ticks(ticks,tick_labels)
# plot.set_xaxis_grid(c='#444444')
# plot.set_yaxis_label("Rz (Angstrom)",size=20)
# plot.set_title("Band center",size=20)

# plot.set_xaxis_margin(0,axes_index=(1,2))
# plot.set_xaxis_ticks(ticks,tick_labels,axes_index=(1,2))
# plot.set_xaxis_grid(axes_index=(1,2),c='#444444')
# plot.set_yaxis_label("energy (eV)",axes_index=(1,2),size=20)
# plot.set_title("Band dispersion",axes_index=(1,2),size=20)

# plot.title("Band center vs dispersion (1st v-band)",size=30)

# plot.stretch_plot(-0.07,"top")
# plot.stretch_plot(-0.07,"top",axes_index=(1,2))
# plot.move_plot(0.01,"right")
# plot.move_plot(0.01,"right",axes_index=(1,2))
# plot.save("/Users/ponet/Documents/Fysica/PhD/GeTe/dip_tests/dipoles89_test.png")

# mesh1 = plot.data_storage.wannier_meshes[0]
# mesh2 = plot.data_storage.wannier_meshes[1]
# #%%
# import time
# t = time.time()
# print(plot.data_storage.dipoles[19])
# print("time:{}".format(time.time()-t))
#---------------1st valence band relative----------------------------#

# ls0 = [0.0,0,0,0,0,0,0,0]
# data_nosoc = plotter.DataStorage(wan_dir="/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest/",ls = ls0,k_points=k_points,centres=centres,computation='dipoles_eigenvalues')

# plot = plotter.Plotter(w=15)
# plot.settings["h_spacing"] = 0.09
# plot.plot_wan_band_dip_SOC("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest",8,ls,centres,k_points=k_points,label="",c='r',linewidth=2)
# plot.plot_wan_band_dip_SOC("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest",9,ls,centres,k_points=k_points,label="",c='b',linewidth=2)

# energy9 = []
# energy8 = []

# for i in range(0,len(data_nosoc.eigenvalues)):
#     energy9.append(plot.data_storage.eigenvalues[i][9]-data_nosoc.eigenvalues[i][9])
#     energy8.append(plot.data_storage.eigenvalues[i][8] - data_nosoc.eigenvalues[i][8])

# dipoles8 = []
# dipoles9 = []
# for i in range(0,len(data_nosoc.eigenvalues)):
#     dipoles8.append(plot.data_storage.dipoles[i][8][2].real-data_nosoc.dipoles[i][8][2].real)
#     dipoles9.append(plot.data_storage.dipoles[i][9][2].real-data_nosoc.dipoles[i][9][2].real)

# mid = dipoles8[499]

# ct = 15.058788*0.7/25

# energy_mid = plot.data_storage.eigenvalues[499][8]
# dip_mid = dipoles8[499]
# dip_energy8 = []
# dip_energy9 = []

# for dip8,dip9 in zip(dipoles8,dipoles9):
#     dip_energy8.append((dip8-dip_mid)*ct+energy8[499])
#     dip_energy9.append((dip9-dip_mid) * ct+energy9[499])

# plot.plot_data(energy8,axes_index=(1,2),linewidth=2,c=(1,0,0,0.4))
# plot.plot_data(energy9,axes_index=(1,2),linewidth=2,c=(0,0,1,0.4))
# plot.plot_data(dip_energy8,axes_index=(1,2),linewidth=2,c='r')
# plot.plot_data(dip_energy9,axes_index=(1,2),linewidth=2,c='b')

# plot.set_data_label("Relative eigenvalues",0,axes_index=(1,2))
# plot.set_data_label("Relative electrostatic energy",2,axes_index=(1,2))

# plot.show_legend(axes_index=(1,2),indices=[0,2],fontsize=15)

# plot.set_xaxis_margin(0)
# plot.set_xaxis_ticks(ticks,tick_labels)
# plot.set_xaxis_grid(c='#444444')
# plot.set_yaxis_label("Rz (Angstrom)",size=20)
# plot.set_title("Band center",size=20)

# plot.set_xaxis_margin(0,axes_index=(1,2))
# plot.set_xaxis_ticks(ticks,tick_labels,axes_index=(1,2))
# plot.set_xaxis_grid(axes_index=(1,2),c='#444444')
# plot.set_yaxis_label("energy (eV)",axes_index=(1,2),size=20)
# plot.set_title("Relative band dispersion",axes_index=(1,2),size=20)

# plot.title("Band center vs dispersion (1st v-band)",size=30)

# plot.stretch_plot(-0.07,"top")
# plot.stretch_plot(-0.07,"top",axes_index=(1,2))
# plot.move_plot(0.01,"right")
# plot.move_plot(0.01,"right",axes_index=(1,2))
# plot.save("/Users/ponet/Documents/Fysica/PhD/GeTe/dip_tests/dipoles89_rel.png")
ls0 = [0.0,0,0,0,0,0,0,0]
centers = 
plot = plotter.Plotter(w=15)
plot.plot_bands("/Users/ponet/Documents/Fysica/PhD/CsPbF3/colin/",energy_range=[0,15],c="b",linewidth=2)
plot.plot_wan_band_eigval_SOC()