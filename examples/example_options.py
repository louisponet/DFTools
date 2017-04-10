import DFTools

# This example shows how to make your plot possibly more beautiful

plot = DFTools.Plotter(w=15)

# Let's first get our information plotted, we'll plot a bandstructure without SOC and one with SOC,
# add ticks and ticklabels, customize the grid layout, and get a nice legend going.

plot.plot_bands("../tests/colin/GeTe_bands.out", energy_range=[6, 15], c='r', linewidth=2)
plot.plot_bands("../tests/fullrel/GeTe_bands.out", energy_range=[6, 15], c='b', linewidth=2)

plot.set_fermi_level(filename="../tests/fullrel/GeTe_scf.out")  # axes_index=(1,1) is always implied

# Let's start by adding some ticklabels and axis labels
x_ticks = [0, 100, 200]
x_ticklabels = ['A', 'Z', 'U']
plot.set_xaxis_ticks(x_ticks, x_ticklabels, size=20)
plot.set_yaxis_tick_params(labelsize=20)
plot.set_yaxis_label("energy (eV)", size=20)
# Now there is some glaring discrepancy between "size" and "labelsize",
# this is because this functions are basically forwarding a load of **kwargs to the underlying matplotlib functions
# and for some reason this is how matplotlib likes to keep things extra "clear". If you are getting errors with regards
# to fontsizes there are 3 kwargs to try: labelsize,size,fontsize.

# now let's get a grid going for the xaxis
plot.set_xaxis_grid(c="#444444")
# c stands for colour in general, this can be passed through in a ton of different formats,
# hex,rgba,'b' for blue, etc

plot.title("Bandstructure GeTe", size=30)
plot.stretch_plot(-0.05, 't')

# We can play around with the details of the plotted bands_sets, these take kwargs that get called when the plot is constructed.
plot.add_set_kwargs(0, linestyle='', marker='x')
# Since we had a colinear band_plot there are 2 sets, one for each spin, so we have to change each of them
plot.add_set_kwargs(1, linestyle='', marker='x')
plot.add_set_kwargs(2, marker='.')
plot.refresh()  # we need to refresh to redraw the elements with the new details

# If you ever add a kwarg that matplotlib doesn't like, you can always remove them again
# plot.add_set_kwargs(0,line=None)
# plot.refresh()
# plot.remove_set_kwargs(0,axes_index=(1,1),line=None)
# plot.refresh()

plot.save("example_options.png")
