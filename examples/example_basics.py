#%%
# First import the basic plotter package
import DFTools

#-----------------------------Basic Band Structure Plot-------------------------##
# Create a Plotter object, w=15 results in a window that has width 15.
# This creates a window with an empty "Axes" object, which represents a graph.

plot = DFTools.Plotter(w=15)
# Now we plot the output of a Quantum-Espresso "bands" run.
plot.plot_bands("../tests/fullrel/GeTe_bands.out")

# The resulting graph is not really displaying what we want, so let's add a zoomed-in plot around the 6-15 eV range,
# while we're at it, we also make the bands a bit thicker and give them all the same colour.

plot.add_plot((1, 2))  # This adds a plot at grid index 1,2 and automatically resizes the plots to fit the window.
plot.plot_bands("../tests/fullrel/GeTe_bands.out", energy_range=[6, 15], axes_index=(1, 2), c="r", linewidth=2)

# Now let's take the Fermi level from an earlier scf output file.
plot.set_fermi_level(filename="../tests/fullrel/GeTe_scf.out", axes_index=(1, 2))
# the axes_index variable denotes the subplot to apply the functions to. This is present in most functions.

# Now let's overlay a bandstructure that does not include SOC.
plot.plot_bands("../tests/colin/GeTe_bands.out", energy_range=[6, 15], axes_index=(1, 2), c="b", linewidth=2)

# To make the plot have a little story we can add various titles to it, and also produce a legend.
plot.set_title("Full Bandstructure", axes_index=(1, 1), size=20)  # axes_index=(1,1) is always implied so we could omit this.
plot.set_title("SOC vs No SOC", axes_index=(1, 2), size=20)  # size is a kwarg that denotes the fontsize of the title.
plot.title("Bandstructure tests", size=30)  # this sets the title of the complete figure.
plot.set_data_label("SOC", 0, axes_index=(1, 2), fontsize=15)  # this assigns labels to the various bands_sets and shows the legend.
plot.set_data_label("No SOC", 1, axes_index=(1, 2), fontsize=15)
# Because I can predict the future, we will have to rescale our graphs a bit to allow some room for our big title.
plot.stretch_plot(-0.07, 'Top', axes_index=(1, 1))
plot.stretch_plot(-0.07, 'Top', axes_index=(1, 2))

# Now that we have everything set up like we want it, let's save the figure.
plot.save("example_basics.png")
