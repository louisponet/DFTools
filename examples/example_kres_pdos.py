#%%
import DFTools

##-------------------------------Band Structure and DOS ---------------------------##

# Create a plotting window.
plot = DFTools.Plotter(w=15)

# First we want to again get a nice band structure plot
plot.plot_bands("../tests/colin/GeTe_bands.out", energy_range=[6, 15], c='r', linewidth=2)

# To analyse what orbitals contribute the most we can add a Partial density of states that gives the component of
# the p-orbitals of each atom
plot.add_plot((1, 2))
plot.add_plot((1, 3))
plot.plot_k_pdos("../tests/colin/pwo.pdos_atm#2(Ge)_wfc#2(p)", 0, energy_range=[6, 15], axes_index=(1, 2), vmin=0, vmax=10)
# The 0 signifies that we want to plot the zeroth density column = total density, in the pdos file.
# vmin,vmax gives relative limits on the colourscale used in the 2d histogram.
plot.plot_k_pdos("../tests/colin/pwo.pdos_atm#1(Te)_wfc#2(p)", 0, energy_range=[6, 15], axes_index=(1, 3), vmin=0, vmax=10)

# It's always nice to have the zero-energy match up with the Fermi-level,
# so we shift our graphs by using set_Fermi_level again.
plot.set_fermi_level(filename="../tests/fullrel/GeTe_scf.out", axes_index=(1, 1))
# Remark: there is a index (e.g, index=(1)) variable that you can use to
# selectively shift a plot to the fermi-level.
plot.set_fermi_level(filename="../tests/fullrel/GeTe_scf.out", axes_index=(1, 2))
plot.set_fermi_level(filename="../tests/fullrel/GeTe_scf.out", axes_index=(1, 3))

# We see that there is redundant information, namely there is an extra
# colourbar.
# ass_plot denotes that we are removing an assisting plot to the main axes.
plot.remove_ass_plot(0, axes_index=(1, 2))
# Now we can stretch the pdos plots to fill the freed space.
plot.stretch_plot(0.03, 'r', axes_index=(1, 2))
plot.stretch_plot(0.03, 'l', axes_index=(1, 3))

# finalize the picture a bit
plot.set_title("Bandstructure GeTe", size=20)
plot.set_title("Ge(p)", axes_index=(1, 2), size=20)
plot.set_title("Te(p)", axes_index=(1, 3), size=20)

# Now we save the finalized product, see that everything is good, and
# proceed with the rest of our daily business.
plot.save("example_kres_pdos.png")
