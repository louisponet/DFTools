#%%
import DFTools

# This tutorial will give example calculations for the Wannier portion of the plotting tool.
# This is intended as an interactive tutorial, but one can put a plot.save at any point to do it non-interactive.
plot = DFTools.Plotter(w=15)
#%%
# First let us check how well the interpolated Wannier band_structure looks compared to the original one.
plot.plot_bands("../tests/colin/GeTe_bands.out", energy_range=[0, 15], c='b', linewidth=2)
# Now we overlay the Wannier interpolation --> read the tight-binding Hamiltonian from Wannier interpolation output, compute the band structure and overlay
plot.overlay_wannier("../tests/colin/xsftest/03GeTe_up_hr.dat", c='r', s=50)
plot.set_title("Collinear Bandstructure", size=20)
plot.set_data_label("DFT", 0, fontsize=15)
plot.set_data_label("Wannier interpolation", 2, fontsize=15)
# The first input is the tight binding Hamiltonian that we got out of the Wannier90 program.
#%%
# We see that the Wannier interpolation worked quite well. Now let's see what happens with SOC.
plot.plot_bands("../tests/fullrel/GeTe_bands.out", energy_range=[0, 15], axes_index=(1, 2), c='b', linewidth=2)
plot.set_title("Noncollinear Bandstructure", axes_index=(1, 2), size=20)
plot.set_data_label("DFT", 0, axes_index=(1, 2), fontsize=15)
plot.set_data_label("Wannier interpolation", 2, axes_index=(1, 2), fontsize=15)
# As you may have noticed, we didn't explicitely add another plot.
# In theory each function is set up to create a new plot at the specified grid position if it doesn't exist yet.

# Since we need to go through extra steps now, we can't just call the overlay_Wannier method.
# The way SOC is calculated is by defining the directory where the program can find the Wannier meshes,
# the centres around which each Wannier function is localized, the lambdas that should be used for each Wannier function,
# and the k_points for which the interpolation should be carried out.
import numpy as np

Ge_center = np.asarray([0.0, 0.0, 0.1114785], np.float64)
Te_center = np.asarray([0.0, 0.0, 5.4186777], np.float64)

centers = [Ge_center, Ge_center, Ge_center, Ge_center, Te_center, Te_center, Te_center, Te_center]
lambdas = [-0.1, -0.1, -0.1, -0.1, -0.32, -0.32, -0.32, -0.32]
plot.plot_wan_eigvals_SOC("../tests/colin/xsftest/", lambdas, centers, k_points='plot', axes_index=(1, 2), c='r', linestyle='None', marker='.', markersize=4)
plot.set_data_label("DFT", 0, fontsize=15)
plot.set_data_label("Wannier interpolation", 2, fontsize=15)
# k_points = 'plot' is a predefined switch that will cause the function to look to the first bands_set in the ax
# and take the k-points that are used there. This could also be defined as a list [kxs,kys,kzs],
# with kxs a list of relative k-values with regards to the first spanning vector,kys to the second,etc.
# We could also define k_points_file = .... with the filename of a bands.out file where the k_points should be taken from.
# We see that things agree quite nicely.

# Now I will showcase some other things that can be done with the Wannier code namely: dipoles and angular momenta
ks = plot.get_k_points(plot.axes[1])  # this is just to say that we want to get the k-points from the first bandstructure on the first plot.
plot.plot_wan_band_dip_eigval_SOC("../tests/colin/xsftest/", 8, lambdas, centers, k_points=ks, axes_indices=[(2, 1), (2, 2)], direction=2, c='b', linewidth=2)
# The 4 here is to indicate that we want to plot things for the 4th band. There is not functionality to plotting everything since it's very messy.
# Direction = 2 means that we want to have the component of the dipoles along the z-axis, this is the default value. 0 = x , 1 = y
plot.plot_wan_band_dip_eigval_SOC("../tests/colin/xsftest/", 9, lambdas, centers, k_points=ks, axes_indices=[(2, 1), (2, 2)], direction=2, c='r', linewidth=2)

plot.set_data_label("band 8", 0, axes_index=(2, 1), fontsize=15)
plot.set_data_label("band 9", 1, axes_index=(2, 1), fontsize=15)
plot.set_data_label("band 8", 0, axes_index=(2, 2), fontsize=15)
plot.set_data_label("band 9", 1, axes_index=(2, 2), fontsize=15)
plot.set_yaxis_label("Rz (Angstrom)", axes_index=(2, 1), size=15)
plot.set_yaxis_label("energy (eV)", axes_index=(2, 2), size=15)
plot.set_title("Band center", axes_index=(2, 1), size=20)
plot.set_title("Band energy", axes_index=(2, 2), size=20)
# Let's clear the bottom two plots to display the angular momenta
plot.clear(axes_index=(1, 1))
plot.clear(axes_index=(1, 2))
plot.plot_wan_angmom("../tests/colin/xsftest/", 4, centers, l='Lx', k_file="../tests/colin/GeTe_bands.out", c='b', linewidth=2)
plot.plot_wan_angmom("../tests/colin/xsftest/", 4, centers, l='Ly', k_file="../tests/colin/GeTe_bands.out", axes_index=(1, 2), c='b', linewidth=2)
# as you can see, we used this time instead of k_points = ... k_file = ... to use the k-points specified in the file.
# Since for now I didn't implement a function that calculates the angular momentum for SOC included calculation, there is no need to define Lambda_soc.
# There are also half the amount of bands, which means that band index 4 corresponds to the band that is spin-split in the upper two plots
# However we still need to define around which centers the Lx's will be calculated.
plot.set_title("Lx", size=20)
plot.set_data_label("Band 8-9", 0, fontsize=15)
plot.set_title("Ly", axes_index=(1, 2), size=20)
plot.set_data_label("Band 8-9", 0, axes_index=(1, 2), fontsize=15)

# Interesting thing to note: Ly doesn't vary when kx or kz don't vary; indices 0->100 only vary ky
plot.save("wannier_example.png")
