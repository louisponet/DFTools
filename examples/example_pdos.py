#%%
import DFTools
##-------------------------------Band Structure and PDOS overlay---------------------------##

# Create a plotting window.
plot = DFTools.Plotter(w=15)

# First we want to again get a nice band structure plot
plot.plot_bands("../tests/colin/GeTe_bands.out", energy_range=[6, 15], c='r', linewidth=2)

#Now we add a second plot that displays the total DOS. Energy will be on the y-axes and the DOS will be on the x-axis.
plot.add_plot((1, 2))
plot.plot_pdos("../tests/fullrel/pwo.pdos_tot", 0, energy_range=[6, 15], axes_index=(1, 2))
# The 0 signifies that we want to plot the zeroth density column = total density, in the pdos file.

# It's always nice to have the zero-energy match up with the Fermi-level, so we shift our graphs by using set_Fermi_level again.
plot.set_fermi_level(filename="../tests/fullrel/GeTe_scf.out", axes_index=(1, 1))
# Remark: there is a index (e.g, index=(1)) variable that you can use to selectively shift a plot to the fermi-level.
plot.set_fermi_level(filename="../tests/fullrel/GeTe_scf.out", axes_index=(1, 2))

# finalize the picture a bit
plot.set_title("Bandstructure GeTe", size=20)
plot.set_title("Total DOS GeTe", axes_index=(1, 2), size=20)

plot.set_xaxis_label("DOS (Arbitrary units)",axes_index=(1,2),size=15)
plot.set_yaxis_label("energy (eV)",axes_index=(1,2),size=15)
# Since I know the future, I know we will have to stretch the second plot to allow for space for the x-axis label.
# We also stretch the first plot to allow for better comparison.
plot.stretch_plot(-0.03,"bottom",axes_index=(1,1))
plot.stretch_plot(-0.03,"bottom",axes_index=(1,2))
# Now we save the finalized product, see that everything is good, and proceed with the rest of our daily business.
plot.save("example_pdos.png")