from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import matplotlib.image as mpimg
from . import fileprc
from . import wrap


def normalize_data(data_set):
    energy = data_set.energy_zero
    if data_set.type != "Image":
        data = data_set.data
        if np.ndim(data) > 2:
            for band in data:
                i = 0
                while i < len(band):
                    band[i, 1] -= energy
                    i += 1
        else:
            i = 0
            while i < len(data):
                data[i, 1] -= energy
                i += 1
        data_set.data = data
    else:
        return


class DataSet:
    def __init__(self, data_type, data_set, label="", energy_zero=0, **kwargs):
        self.type = data_type
        self.data = data_set
        self.energy_zero = energy_zero
        self.label = label
        self.handle = None
        self.kwargs = kwargs


class Axes(plt.Axes):
    """Main graphing object. It represents a subplot in the plotting entity. Extends the functionality of the basic
    matplotlib axes object. Init param are the parent figure and position rectangle."""

    def __init__(self, fig, rect):
        # TODO: Test everything!!!!!!
        # TODO: add change settings
        # TODO: think about how to handle colorbars
        super().__init__(fig, rect)
        self.settings = {}
        self.settings["linewidth"] = 1
        self.settings["scatter_size"] = 40
        self.plot_title = ""
        self.ass_axes = []
        self.index = None
        self.legend_kwargs = {}
        # this includes normal bands (up/down), dos overlays and wannier bands
        self.data_sets = []
        self.fig = fig
        self.energy_zero = 0
        self.x_range = None
        self.x_ticks = None
        self.x_tick_labels = None
        self.x_tick_params = {}
        self.x_label = ""
        self.x_label_params = {}
        self.x_grid_params = {}
        self.y_ticks = None
        self.y_tick_labels = None
        self.y_tick_params = {}
        self.y_label = ""
        self.y_label_params = {}
        self.y_grid_params = {}
        self.tick_params(labelsize=15)

    def set_fermi_level(self, energy, data_index=None):
        """Sets the fermi level (energy zero point). Set index parameter selects which subset gets shifted with respect
        to the zero point."""
        self.energy_zero = energy
        if data_index is None:
            for data_set in self.data_sets:
                if data_set.energy_zero != energy:
                    data_set.energy_zero = energy
                    normalize_data(data_set)
        else:
            if abs(data_index) < len(self.data_sets):
                if self.data_sets[data_index].energy_zero != energy:
                    self.data_sets[data_index].energy_zero = energy
                    normalize_data(self.data_sets[data_index])
            else:
                print("Error: Data index exceeds amount of plotted data sets.")
 
        self.energy_range[0]-=self.energy_zero
        self.energy_range[1]-=self.energy_zero         
        self.display()

    def display(self):
        """Main plotting function, gets called after each addition of a data_set or any change to the graphed objects."""

        super().clear()

        labels = []

        for data_set in self.data_sets:
            d_type = data_set.type
            if "KPDOS" == d_type:
                data_set_tmp = data_set.data
                k_max = np.max(data_set_tmp[:, 0])
                e_bin = int(len(data_set_tmp) / k_max)
                data_set.handle = self.hist2d(data_set_tmp[:, 0], data_set_tmp[:, 1], weights=data_set_tmp[:, 2], bins=[
                    k_max, e_bin], cmap=plt.get_cmap('plasma'), **data_set.kwargs)[-1]

        for data_set in self.data_sets:
            kwargs = data_set.kwargs
            d_type = data_set.type
            labels.append(data_set.label)
            if "Bandstructure" in d_type and "PDOS" not in d_type:
                for band in data_set.data:
                    if "linewidth" in data_set.kwargs:
                        data_set.handle, = self.plot(
                            band[:, 0], band[:, 1], **kwargs)
                    else:
                        data_set.handle, = self.plot(
                            band[:, 0], band[:, 1], linewidth=self.settings["linewidth"], **kwargs)

            elif "PDOS" in d_type and d_type != "KPDOS":
                band = data_set.data
                if 'linewidth' in kwargs:
                    data_set.handle, = self.plot(
                        band[:, 1], band[:, 0], **kwargs)
                else:
                    data_set.handle, = self.plot(band[:, 1], band[:, 0], linewidth=self.settings["linewidth"], **kwargs)

            elif "Wannier" in d_type:
                for band in data_set.data:
                    if "s" in kwargs:
                        data_set.handle = self.scatter(
                            band[:, 0], band[:, 1], marker='.', edgecolors='none', **kwargs)
                    else:
                        data_set.handle = self.scatter(
                            band[:, 0], band[:, 1], marker='.', edgecolors='none', s=self.settings["scatter_size"], **kwargs)
            elif "Image" in d_type:
                self.imshow(data_set.data)
                self.spines['left'].set_color('none')
                self.spines['right'].set_color('none')
                self.spines['top'].set_color('none')
                self.spines['bottom'].set_color('none')
                self.xaxis.set_ticks([])
                self.yaxis.set_ticks([])
            elif "Fatbands" in d_type:
                for band in data_set.data:
                    if "s" in kwargs:
                        data_set.handle = self.scatter(band[:, 0], band[:, 1], cmap=plt.get_cmap('plasma'), c=band[:, 2],
                                                       marker='.', edgecolors='none',
                                                       **kwargs)
                    else:
                        data_set.handle = self.scatter(band[:, 0], band[:, 1], cmap=plt.get_cmap(
                            'plasma'), c=band[:, 2], marker='.', edgecolors='none', s=self.settings["scatter_size"], **kwargs)
            elif "Plot" in d_type or "plot" in d_type:
                if "linewidth" in data_set.kwargs:
                    data_set.handle, = self.plot(
                        data_set.data[:, 0], data_set.data[:, 1], **kwargs)
                else:
                    data_set.handle, = self.plot(
                        data_set.data[:, 0], data_set.data[:, 1], linewidth=self.settings["linewidth"], **kwargs)

            elif "Scatter" in d_type or "scatter" in d_type:
                data_set.handle = self.scatter(data_set.data[:, 0], data_set.data[:, 1], **kwargs)

        if self.energy_range is not None:
            self.set_ylim([self.energy_range[0],
                           self.energy_range[1]])
        if self.x_range is not None:
            self.set_xlim(self.x_range)
        self.set_title(self.plot_title)
        self.config_axis()
        self.fig.show()

    def add_data_set(self, data_set):
        """Adds data_set and normalizes it if there was a zero_energy set beforehand."""
        data_set.energy_zero = self.energy_zero
        normalize_data(data_set)
        self.data_sets.append(data_set)

    def change_band_colour(self, band, colour):
        """Change one of the data_set's colour."""
        if type(band) == int and band < len(self.data_sets):
            self.data_sets[band].colour = colour
        elif type(band) == str:
            for band_set in self.data_sets:
                if band in band_set.type:
                    band_set.kwargs['c'] = colour

    def remove_data_set(self, index):
        """Remove a data_set from the graph."""
        if abs(index) < len(self.data_sets):
            self.data_sets.pop(index)
        else:
            print("Error: Index exceeds amount of data sets.")

    def set_data_label(self, set_index=0, label=''):
        """Sets the legend label of a displayed data_set."""
        if abs(set_index) < len(self.data_sets):
            self.data_sets[set_index].label = label

    def show_legend(self, indices=None, **kwargs):
        """Displays/changes the legend of the graph. The indices parameter allows for selective display in the legend."""
        labels = []
        handles = []
        for key, val in kwargs.items():
            self.legend_kwargs[key] = val
        if indices is None:
            for data_set in self.data_sets:
                if data_set.type == "KPDOS":
                    continue
                if data_set.label == "":
                    continue
                labels.append(data_set.label)
                handles.append(data_set.handle)
        else:
            for i in indices:
                labels.append(self.data_sets[i].label)
                handles.append(self.data_sets[i].handle)
        if len(handles) > 0:
            self.legend(handles=handles, labels=labels, **self.legend_kwargs)
        else:
            pass

    def create_colourbar(self):
        """Semi-private function that creates a colourbar when a DOS/PDOS styled graph is displayed."""
        for asax in self.ass_axes:
            if type(asax) == matplotlib.colorbar.Colorbar:
                return
        data_set = None
        for d_set in self.data_sets:
            if d_set.type == "KPDOS" or d_set.type == "Fatbands":
                data_set = d_set
        self.ass_axes.append(self.fig.colorbar(data_set.handle, ax=self))

    def set_xlabel(self, label, **kwargs):
        super().set_xlabel(label, **kwargs)
        self.x_label = label
        self.x_label_params = kwargs

    def set_ylabel(self, label, **kwargs):
        super().set_ylabel(label, **kwargs)
        self.y_label = label
        self.y_label_params = kwargs

    def set_xticks(self, ticks):
        super().set_xticks(ticks)
        self.x_ticks = ticks

    def set_yticks(self, ticks):
        super().set_yticks(ticks)
        self.y_ticks = ticks

    def set_xticklabels(self, labels, **kwargs):
        super().set_xticklabels(labels, **kwargs)
        self.x_tick_labels = labels
        self.x_tick_params = kwargs

    def set_yticklabels(self, labels, **kwargs):
        super().set_yticklabels(labels, **kwargs)
        self.y_tick_labels = labels
        self.y_tick_params = kwargs

    def set_xtick_params(self, **kwargs):
        self.tick_params(axis="x", **kwargs)
        self.x_tick_params = kwargs

    def set_ytick_params(self, **kwargs):
        self.tick_params(axis="y", **kwargs)
        self.y_tick_params = kwargs

    def set_xgrid(self, **kwargs):
        self.xaxis.grid(**kwargs)
        self.x_grid_params = kwargs

    def set_ygrid(self, **kwargs):
        self.yaxis.grid(**kwargs)
        self.y_grid_params = kwargs

    def config_axis(self):
        try:
            super().set_xlabel(self.x_label, **self.x_label_params)
        except:
            pass
        try:
            super().set_ylabel(self.y_label, **self.y_label_params)
        except:
            pass
        try:
            super().set_xticks(self.x_ticks)
        except:
            pass
        try:
            super().set_yticks(self.y_ticks)
        except:
            pass
        try:
            super().set_xticklabels(self.x_tick_labels, **self.x_tick_params)
        except:
            pass
        try:
            super().set_yticklabels(self.y_tick_labels, **self.y_tick_params)
        except:
            pass
        try:
            super().xaxis.grid(**self.x_grid_params)
        except:
            pass
        try:
            super().yaxis.grid(**self.y_grid_params)
        except:
            pass
        self.show_legend(**self.legend_kwargs)

    def clear(self):
        """Clears the entire axes object, and also clears the list of data_sets."""
        super().clear()
        self.data_sets.clear()
        self.energy_zero = 0
        self.energy_range = None
        self.x_range = None
        self.x_ticks = None
        self.x_tick_labels = None
        self.x_tick_params = {}
        self.x_label = ""
        self.x_label_params = {}
        self.x_grid_params = {}
        self.y_ticks = None
        self.y_tick_labels = None
        self.y_tick_params = {}
        self.y_label = ""
        self.y_label_params = {}
        self.y_grid_params = {}


class DataStorage:
    def __init__(self, **kwargs):
        self.wan_dir = None
        self.wannier_meshes = None
        self.hami_raw = None
        self.grid = None
        self.k_points = None
        self.dipoles = None
        self.eigenvalues = None
        self.eigenvectors = None
        self.centres = None
        self.ls = None
        self.angmoms = None

        computation = False
        tmp_band_index = None
        for key, value in kwargs.items():
            if key == "wan_dir":
                self.wan_dir = value
                self.load_wannier_info(self.wan_dir)
            elif key == "k_points":
                self.load_kpoints(k_points=value)
            elif key == "k_file":
                self.load_kpoints(k_file=value)
            elif key == "ls":
                self.ls = value
            elif key == "centres":
                self.centres = value
            elif key == "band_index":
                tmp_band_index = value
            elif key == "computation":
                computation = value
        if computation is not None:
            if computation == "dipoles_eigenvalues":
                self.calculate_dipoles_eigenvalues(self.ls, self.centres)
            elif computation == "eigenvectors_SOC":
                self.calculate_eigenvectors_SOC(
                    tmp_band_index, self.ls, self.centres)
            elif computation == "angmoms":
                self.calc_angmoms(self.centres)

    def load_wannier_info(self, wan_dir):
        import glob
        wan_files = None
        hami_file = None
        if wan_dir[-1] == "/":
            wan_files = glob.glob(wan_dir + "*.xsf")
            hami_file = glob.glob(wan_dir + "*hr.dat")
        else:
            wan_files = glob.glob(wan_dir + "/*.xsf")
            hami_file = glob.glob(wan_dir + "/*hr.dat")

        if len(wan_files) == 0:
            print("Error: no xsf files found in dir {}!".format(wan_dir))
            return
        if len(hami_file) == 0:
            print("Error: no hamiltonian file found in dir {}!".format(wan_dir))
            return
        elif len(hami_file) > 1:
            print("Error: ambiguity concerning which hamiltonian file to use.\nSupplied files:{}.".format(
                hami_file))
            return
        self.wannier_meshes = []

        prim_cells = []
        origins = []
        a_span_vecs = []
        b_span_vecs = []
        c_span_vecs = []
        a_span_arrays = []
        b_span_arrays = []
        c_span_arrays = []

        for file in wan_files:
            mesh, prim_cell, origin, a_vec, b_vec, c_vec, a_array, b_array, c_array = fileprc.read_wannier_function(
                file)
            self.wannier_meshes.append(mesh)
            prim_cells.append(prim_cell)
            origins.append(origin)
            a_span_vecs.append(a_vec)
            b_span_vecs.append(b_vec)
            c_span_vecs.append(c_vec)
            a_span_arrays.append([a * a_vec for a in a_array])
            b_span_arrays.append([b * b_vec for b in b_array])
            c_span_arrays.append([c * c_vec for c in c_array])
        origin = origins[0]
        a_span_array = a_span_arrays[0]
        b_span_array = b_span_arrays[0]
        c_span_array = c_span_arrays[0]

        self.grid = wrap.generate_grid(
            origin, a_span_array, b_span_array, c_span_array)
        self.hami_raw = fileprc.read_wan_hami(hami_file[0])
        self.dir = wan_dir

    def load_kpoints(self, k_file=None, k_points=None):
        if k_points is None:
            if k_file is None:
                print(
                    "Error: specify either a bands.out file or a preconfigured set of k_points!")
                return
            bands = fileprc.create_bands_from_file(k_file)
            kxs = bands[0][0][:, 5]
            kys = bands[0][0][:, 6]
            kzs = bands[0][0][:, 7]
        else:
            kxs = k_points[0]
            kys = k_points[1]
            kzs = k_points[2]

        self.k_points = np.empty((len(kxs), 3))
        self.k_points[:, 0] = kxs
        self.k_points[:, 1] = kys
        self.k_points[:, 2] = kzs

    def calculate_dipoles_eigenvalues(self, lambda_socs, centres):
        if self.k_points is None:
            print("Error: load k-points first!")
        if self.hami_raw is None:
            print("Error: load wannier info first!")
        results = wrap.calculate_tb_dipoles(
            self.hami_raw, self.k_points[:, 0], self.k_points[:, 1], self.k_points[:, 2], self.wannier_meshes, self.grid, centres, lambda_socs, lambda_socs, lambda_socs)
        self.dipoles = []
        self.eigenvalues = []
        for result in results:
            self.eigenvalues.append(result[-1].real)
            self.dipoles.append(result[0])

    def calculate_eigenvectors_SOC(self, band_index, lambda_socs, centres):
        self.eigenvectors = wrap.generate_wannier_bands_SOC(
            self.hami_raw, self.k_points[:, 0], self.k_points[:, 1], self.k_points[:, 2], self.wannier_meshes, self.grid,
            centres, lambda_socs, lambda_socs, lambda_socs)[band_index][1]

    def calc_angmoms(self, centres):
        self.angmoms = wrap.calculate_tb_angmoms(
            self.hami_raw, self.k_points[:, 0], self.k_points[:, 1], self.k_points[:, 2], self.wannier_meshes, self.grid, centres)


class Plotter:
    # TODO make functions specifically for multiple actions at once!
    # TODO sort out how to apply settings and the like!
    # TODO Allow for selection of multiple subplot indices to plot something on
    """The plotting class, supplies all plotting functions and allows control over the various subplots which are represented by the Axes class.
    Initialization variables w and h denote weight and height. Further kwargs will be passed through to the matplotlib fig() function, which creates
    the matplotlib figure."""

    def __init__(self, w=8, h=8, **kwargs):
        self.fig = plt.figure(figsize=(w, h), frameon=False, **kwargs)
        self.axes = []
        self.settings = {}
        self.settings["h_spacing"] = 0.07
        self.settings["v_spacing"] = 0.07
        self.settings["h_border"] = 0.05
        self.settings["v_border"] = 0.05
        # plt.ion()
        ax = Axes(self.fig, [0.05, 0.05, 0.90, 0.90])
        ax.index = [1, 1]
        self.fig.add_axes(ax)
        self.axes.append(ax)

        self.subplot_dim = [1, 1]
        self.fig.show()
        self.data_storage = DataStorage()
##-------------------------Start of basic plot manipulation section ---------------------------##
    def add_plot(self, axes_index):
        """Add a subplot in the subplot grid. The grid is set up so that the rows rise in number from bottom to top, while the columns from left to right."""

        if len(axes_index) != 2:
            print("Error: please provide valid indices in the correct format: (x,y).")
            return
        ax = self.find_axes(axes_index)
        if ax != 0:
            print("Error: that subplot already exists!")
            return

        if axes_index[0] > self.subplot_dim[0]:
            self.subplot_dim[0] = axes_index[0]

        if axes_index[1] > self.subplot_dim[1]:
            self.subplot_dim[1] = axes_index[1]

        new_height = (1 - 2 * self.settings["v_border"] - self.settings["v_spacing"] * (
            self.subplot_dim[0] - 1)) / self.subplot_dim[0]
        new_width = (1 - 2 * self.settings["h_border"] - self.settings["h_spacing"] * (
            self.subplot_dim[1] - 1)) / self.subplot_dim[1]
        new_ax = Axes(self.fig,
                      [self.settings["h_border"] + self.settings["h_spacing"] * (axes_index[1] - 1) + new_width * (axes_index[1] - 1),
                       self.settings["h_border"] + self.settings["v_spacing"] * (axes_index[0] - 1) + new_height * (axes_index[0] - 1),
                       new_width, new_height])
        new_ax.index = [axes_index[0], axes_index[1]]
        self.fig.add_axes(new_ax)
        self.axes.append(new_ax)
        self._resize_axes(self.subplot_dim[0], self.subplot_dim[1])

    def remove_plot(self, axes_index):
        """Removes a subplot and all related quantities from the figure and plotter object. Same indexing as add_plot."""
        if len(axes_index) != 2:
            print("Error: please provide valid indices in the correct format: (x,y).")
            return

        i = 0
        while i < len(self.axes):
            if self.axes[i].index[0] == axes_index[0] and self.axes[i].index[1] == axes_index[1]:
                self.fig.delaxes(self.axes[i])
                self.axes.pop(i)
                break

            i += 1
            if i == len(self.axes):
                print("Error: no such subplot index in use!")
                return

        max_x_index = 1
        max_y_index = 1
        for ax in self.axes:
            if ax.index[0] > max_x_index:
                max_x_index = ax.index[0]

            if ax.index[1] > max_y_index:
                max_y_index = ax.index[1]
        self._resize_axes(max_x_index, max_y_index)

    def _resize_axes(self, x_index, y_index):
        self.subplot_dim = [x_index, y_index]
        new_height = (1 - 2 * self.settings["v_border"] - self.settings["v_spacing"] * (
            self.subplot_dim[0] - 1)) / self.subplot_dim[0]
        new_width = (1 - 2 * self.settings["h_border"] - self.settings["h_spacing"] * (
            self.subplot_dim[1] - 1)) / self.subplot_dim[1]
        for ax in self.axes:
            for ass_ax in ax.ass_axes:
                ass_ax.remove()
            ax.ass_axes.clear()

            ax.set_position([self.settings["h_border"] + self.settings["h_spacing"] * (ax.index[1] - 1) + new_width * (ax.index[1] - 1),
                             self.settings["v_border"] + self.settings["v_spacing"] * (ax.index[0] - 1) + new_height * (ax.index[0] - 1), new_width, new_height])
            for data_set in ax.data_sets:
                if data_set.type == "PDOS" or data_set.type == "Fatbands":
                    ax.create_colourbar()

    def stretch_plot(self, amount, direction, axes_index=(1, 1)):
        ax = self.find_axes(axes_index)
        if ax != 0:
            origpos = ax.get_position()
            if direction == "left" or direction == "Left" or direction == "l":
                ax.set_position([origpos.x0 - amount, origpos.y0,
                                 origpos.width + amount, origpos.height])
            elif direction == "right" or direction == "Right" or direction == "r":
                ax.set_position(
                    [origpos.x0, origpos.y0, origpos.width + amount, origpos.height])
            elif direction == "top" or direction == "Top" or direction == "t":
                ax.set_position(
                    [origpos.x0, origpos.y0, origpos.width, origpos.height + amount])
            elif direction == "bottom" or direction == "Bottom" or direction == "b":
                ax.set_position([origpos.x0, origpos.y0 - amount,
                                 origpos.width, origpos.height + amount])
        else:
            print("Error: axes doesn't exist!")

    def move_plot(self, amount, direction, axes_index=(1, 1)):
        ax = self.find_axes(axes_index)
        if ax != 0:
            origpos = ax.get_position()
            if direction == "left" or direction == "Left" or direction == "l":
                ax.set_position(
                    [origpos.x0 - amount, origpos.y0, origpos.width, origpos.height])
            elif direction == "right" or direction == "Right" or direction == "r":
                ax.set_position(
                    [origpos.x0 + amount, origpos.y0, origpos.width, origpos.height])
            elif direction == "top" or direction == "Top" or direction == "t":
                ax.set_position(
                    [origpos.x0, origpos.y0 + amount, origpos.width, origpos.height])
            elif direction == "bottom" or direction == "Bottom" or direction == "b":
                ax.set_position(
                    [origpos.x0, origpos.y0 - amount, origpos.width, origpos.height])
        else:
            print("Error: axes doesn't exist!")

    def copy_plot(self, start_index, end_index, overwrite=False):
        ax_to_copy = None
        ax_to_overwrite = None
        for ax in self.axes:
            if ax.index[0] == start_index[0] and ax.index[1] == start_index[1]:
                ax_to_copy = ax
            if ax.index[0] == end_index[0] and ax.index[1] == end_index[1]:
                if not overwrite:
                    print(
                        "Error: That index is already used by another subplot!\n To overwrite already existing plots set 'overwrite = True'.")
                    return
                else:
                    ax_to_overwrite = ax

        if ax_to_copy is None:
            print("Error: No plot exists with that index!")
            return

        if ax_to_overwrite is None:
            self.add_plot(end_index)
            self.axes[-1].plot_title = ax_to_copy.plot_title
            self.axes[-1].data_sets = ax_to_copy.data_sets
            self.axes[-1].x_range = ax_to_copy.x_range
            self.axes[-1].energy_zero = ax_to_copy.energy_zero
            self.axes[-1].plotted_indices = ax_to_copy.plotted_indices
            self.axes[-1].energy_range = ax_to_copy.energy_range
            self.axes[-1].index = [end_index[0], end_index[1]]
            self.axes[-1].display()

        else:
            self.axes[self.axes.index(ax_to_overwrite)
                      ].plot_title = ax_to_copy.plot_title
            self.axes[self.axes.index(ax_to_overwrite)
                      ].data_sets = ax_to_copy.data_sets
            self.axes[self.axes.index(ax_to_overwrite)
                      ].x_range = ax_to_copy.x_range
            self.axes[self.axes.index(ax_to_overwrite)
                      ].energy_zero = ax_to_copy.energy_zero
            self.axes[self.axes.index(
                ax_to_overwrite)].plotted_indices = ax_to_copy.plotted_indices
            self.axes[self.axes.index(ax_to_overwrite)
                      ].energy_range = ax_to_copy.energy_range
            self.axes[self.axes.index(ax_to_overwrite)].index = [
                end_index[0], end_index[1]]
            self.axes[self.axes.index(ax_to_overwrite)].display()

    def remove_ass_plot(self, index, axes_index=(1, 1)):
        ax = self.find_axes(axes_index)
        if ax != 0:
            ax.ass_axes[index].remove()
        else:
            print("Error: That ax doesn't exist.")
            return

    def find_axes(self, index):
        for ax in self.axes:
            if hasattr(ax, "index"):
                if index[0] == ax.index[0] and index[1] == ax.index[1]:
                    return ax
            else:
                continue
        return 0

    def set_size(self, w, h):
        self.fig.set_size_inches(w, h, forward=True)

    def save(self, filename, **kwargs):
        plt.savefig('{}'.format(filename), **kwargs)

    def title(self, title, **kwargs):
        self.fig.suptitle(title, **kwargs)

    def clear(self, axes_index=(1, 1)):
        ax = self.find_axes(axes_index)
        ax.clear()
        ax.data_sets.clear()
        for ass_ax in ax.ass_axes:
            ass_ax.remove()
        ax.ass_axes.clear()
        self._resize_axes(self.subplot_dim[0], self.subplot_dim[1])

    def show(self):
        self.fig.show()

    def refresh(self):
        for ax in self.axes:
            ax.display()

    def resize(self, w, h):
        self.fig.set_size_inches(w, h)
        self._resize_axes(self.subplot_dim[0], self.subplot_dim[1])

 ##-------------------------Start of plot options section ---------------------------##
    def set_data_label(self, label, set_index, axes_index=(1, 1), **kwargs):
        ax = self.find_axes(axes_index)
        if ax != 0:
            ax.set_data_label(set_index, label)
            ax.show_legend(**kwargs)
        else:
            print("Error: there is no plot with that index.")

    def remove_data_label(self, set_index, axes_index=(1, 1)):
        ax = self.find_axes(axes_index)
        if ax != 0:
            ax.data_sets[set_index].label = ''
            ax.show_legend()
        else:
            print("Error: there is no plot with that index.")

    def set_xaxis_label(self, label, axes_index=(1, 1), **kwargs):
        ax = self.find_axes(axes_index)
        if ax == 0:
            print("Error:That ax doesn't exist.")
            return
        ax.set_xlabel(label, **kwargs)

    def set_yaxis_label(self, label, axes_index=(1, 1), **kwargs):
        ax = self.find_axes(axes_index)
        if ax == 0:
            print("Error:That ax doesn't exist.")
            return
        ax.set_ylabel(label, **kwargs)

    def set_xaxis_ticks(self, ticks, labels, axes_index=(1, 1), **kwargs):
        ax = self.find_axes(axes_index)
        if ax == 0:
            print("Error:That ax doesn't exist.")
            return
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, **kwargs)

    def set_yaxis_ticks(self, ticks, labels, axes_index=(1, 1), **kwargs):
        ax = self.find_axes(axes_index)
        if ax == 0:
            print("Error:That ax doesn't exist.")
            return
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels, **kwargs)

    def set_xaxis_tick_params(self, axes_index=(1, 1), **kwargs):
        ax = self.find_axes(axes_index)
        if ax == 0:
            print("Error:That ax doesn't exist.")
            return
        ax.set_xtick_params(**kwargs)

    def set_yaxis_tick_params(self, axes_index=(1, 1), **kwargs):
        ax = self.find_axes(axes_index)
        if ax == 0:
            print("Error:That ax doesn't exist.")
            return
        ax.set_ytick_params(**kwargs)

    def set_xaxis_grid(self, axes_index=(1, 1), **kwargs):
        ax = self.find_axes(axes_index)
        if ax != 0:
            ax.set_xgrid(**kwargs)
        else:
            print("Error:That ax doesn't exist.")

    def set_yaxis_grid(self, axes_index=(1, 1), **kwargs):
        ax = self.find_axes(axes_index)
        if ax != 0:
            ax.set_ygrid(**kwargs)
        else:
            print("Error:That ax doesn't exist.")

    def set_xaxis_margin(self, margin, axes_index=(1, 1), **kwargs):
        ax = self.find_axes(axes_index)
        if ax != 0:
            ax.set_xmargin(margin)
            ax.autoscale()
        else:
            print("Error: axes doesn't exist!")

    def set_yaxis_margin(self, margin, axes_index=(1, 1), **kwargs):
        ax = self.find_axes(axes_index)
        if ax != 0:
            ax.set_ymargin(margin)
            ax.autoscale()
        else:
            print("Error: axes doesn't exist!")

    def show_legend(self, axes_index=(1, 1), **kwargs):
        ax = self.find_axes(axes_index)
        if ax == 0:
            print("Error: There is no plot with that index.")
            return
        else:
            ax.show_legend(**kwargs)

    def add_set_kwargs(self, set_index, axes_index=(1, 1), **kwargs):
        ax = self.find_axes(axes_index)
        if ax == 0:
            print("Error: There is no plot with that index.")
            return
        else:
            bset = ax.data_sets[set_index]
            for key, value in kwargs.items():
                bset.kwargs[key] = value

    def remove_set_kwargs(self, set_index, axes_index=(1, 1), **kwargs):
        ax = self.find_axes(axes_index)
        if ax == 0:
            print("Error: There is no plot with that index.")
            return
        else:
            bset = ax.data_sets[set_index]
            for key, value in kwargs.items():
                bset.kwargs.pop(key)

    def remove_legend_kwargs(self, axes_index=(1, 1), **kwargs):
        ax = self.find_axes(axes_index)
        if ax == 0:
            print("Error: There is no plot with that index.")
            return
        else:
            for key, value in kwargs.items():
                ax.legend_kwargs.pop(key)

    def set_title(self, title, axes_index=(1, 1), size=None):
        for ax in self.axes:
            if ax.index[0] == axes_index[0] and ax.index[1] == axes_index[1]:
                if size is not None:
                    ax.set_title(title, size=size)
                else:
                    ax.set_title(title)

    def remove_data_set(self, band_index, axes_index=(1, 1)):
        ax = self.find_axes(axes_index)
        if ax != 0:
            ax.remove_data_set(band_index)
            ax.display()
        else:
            print("Error: there is no plot with that index.")

    def set_energy_range(self, energy_range, axes_index=(1, 1)):
        ax = self.find_axes(axes_index)
        if ax != 0:
            ax.energy_range = energy_range
            ax.display()
        else:
            print("Error: there is no plot with that index.")

 ##-------------------------Start of plotting section ---------------------------##
 # TODO how to manage spin polarized stuff???? ask for input?
    def plot_data(self, val_list, d_type='Plot', axes_index=(1, 1), label="", **kwargs):
        indices = np.arange(0, len(val_list))
        set_ = np.transpose(np.vstack((indices, val_list)))
        d_set = DataSet(d_type, set_, label, **kwargs)
        ax = self.find_axes(axes_index)
        if ax != 0:
            ax.add_data_set(d_set)
            ax.display()
        else:
            self.add_plot(axes_index)
            self.axes[-1].add_data_set(d_set)
            self.axes[-1].display()

    def plot_bands(self, filename, axes_index=(1, 1), energy_range=[-15,15], label=None, **kwargs):
        bands_array = fileprc.create_bands_from_file(filename)
        ax = self.find_axes(axes_index)
        if ax == 0:
            self.add_plot(axes_index)
            ax = self.axes[-1]

        legend_index = 1
        for data_set in ax.data_sets:
            if "Bandstructure" in data_set.type:
                legend_index += 1
        ax.energy_range = energy_range
        for i in range(0, len(bands_array)):
            tmp_set = DataSet("Bandstructure", bands_array[i], **kwargs)
            if label is not None:
                tmp_set.label = label
            ax.add_data_set(tmp_set)
        ax.x_range = [0, len(bands_array[0][0])]
        ax.display()

    def change_colour(self, subplot, band, colour):
        for ax in self.axes:
            if ax.index[0] == subplot[0] and ax.index[1] == subplot[1]:
                ax.change_band_colour(band, colour)
                ax.display()

    def overlay_dos(self, proj_out_filename, atoms=None, axes_indices=None):
        indices_to_plot = []

        if axes_indices is None:
            indices_to_plot = [(1, 1)]
        else:
            indices_to_plot = axes_indices

        axes_to_plot = []
        for ax in self.axes:
            if tuple(ax.index) in indices_to_plot and len(["Bandstructure" in x for x in iter(ax.data_sets.keys())]) > 0:
                axes_to_plot.append(ax)

        if len(axes_to_plot) == 0:
            print("Error: Plot Bandstructure first!")
            return

        wfc, k_point_arrays = fileprc.read_pdos_colin(proj_out_filename, atoms)

        for ax in axes_to_plot:
            indices = np.where(["Bandstructure" in x for x in ax.data_sets])
            for index in indices[0]:
                bands_array = ax.data_sets[index].data_set
                i = 0
                while i < len(bands_array):
                    i2 = 0
                    while i2 < len(bands_array[i]):
                        bands_array[i][i2, 2:5] = k_point_arrays[i2][i][1:]
                        i2 += 1
                    i += 1
                    ax.data_sets[index].type = ax.data_sets[index].type + "_DOS"

            ax.plot_title = ax.plot_title + " with DOS"
            if atoms is not None:
                ax.plot_title += " (" + atoms + ")"
            ax.display()

    def overlay_wannier(self, hami_file, axes_indices=[(1, 1)], **kwargs):
        hami = fileprc.read_wan_hami(hami_file)
        for axes_index in axes_indices:
            for ax in self.axes:
                if ax.index[0] == axes_index[0] and ax.index[1] == axes_index[1]:
                    if len(ax.data_sets) == 0:
                        print("Error: Plot Bandstructure first!")
                        return
                    else:
                        band_set = None
                        for data_set in ax.data_sets:
                            if "Bandstructure" in data_set.type:
                                band_set = np.asarray(
                                    data_set.data[0][:, 5:], dtype=np.double)
                        bands_tmp = wrap.generate_wannier_bands(
                            hami, band_set[:, 0], band_set[:, 1], band_set[:, 2])
                        bands_t = []
                        for bands in bands_tmp:
                            bands_t.append(np.transpose(
                                np.vstack((np.arange(0, len(bands[0])), bands[0]))))
                        wannier_bands = DataSet("Wannier", bands_t, **kwargs)
                        ax.add_data_set(wannier_bands)
                        ax.display()

    def show_png(self, filename, axes_index=(1, 1), overwrite=False, **kwargs):

        image = mpimg.imread(filename)

        drawn = False
        for ax in self.axes:
            if ax.index[0] == axes_index[0] and ax.index[1] == axes_index[1]:
                if len(ax.data_sets) != 0 and not overwrite:
                    print(
                        "Error: This plot is already in use.\nIf you want to overwrite set 'overwrite=False'.")
                    return
                else:
                    drawn = True
                    ax.clear()
                    ax.add_data_set("Image", (image, None), **kwargs)
                    ax.plot_title = filename.split("/")[-1]
                    ax.display()

        if not drawn:
            self.add_plot(axes_index)
            self.axes[-1].clear()
            self.axes[-1].plot_title = filename.split("/")[-1]
            self.axes[-1].add_data_set("Image", (image, None))
            self.axes[-1].display()

    def plot_pdos(self, filename, column, energy_range=None, label=None, axes_index=(1, 1), **kwargs):
        data_set = DataSet("PDOS", fileprc.read_pdos(filename, column=column + 1), **kwargs)
        if label is None:
            data_set.label = "PDOS(Column {})".format(column)
        else:
            data_set.label = label
        ax = self.find_axes(axes_index)
        if ax != 0:
            ax.data_sets.clear()
            ax.add_data_set(data_set)
            ax.energy_range = energy_range
            ax.display()
        else:
            self.add_plot(axes_index)
            self.axes[-1].add_data_set(data_set)
            self.axes[-1].energy_range = energy_range
            self.axes[-1].display()

    def plot_k_pdos(self, filename, column, energy_range=None, label=None, axes_index=(1, 1), **kwargs):
        import time
        t = time.time()
        data_set = DataSet("KPDOS", fileprc.read_k_pdos(
            filename, column + 2), **kwargs)
        print("time was:{}".format(time.time()-t))
        if label is None:
            data_set.label = "PDOS(Column {})".format(column)
        else:
            data_set.label = label
        ax = self.find_axes(axes_index)
        if ax != 0:
            ax.data_sets.clear()
            ax.add_data_set(data_set)
            ax.energy_range = energy_range
            ax.display()
            ax.create_colourbar()
        else:
            self.add_plot(axes_index)
            self.axes[-1].add_data_set(data_set)
            self.axes[-1].energy_range = energy_range
            self.axes[-1].display()
            self.axes[-1].create_colourbar()

    def plot_fatbands(self, filename, energy_range=None, axes_index=(1, 1), **kwargs):
        data_set = DataSet(
            "Fatbands", fileprc.read_fatbands(filename), **kwargs)
        data_set.label = "Fatbands"
        ax = self.find_axes(axes_index)
        if ax != 0:
            ax.add_data_set(data_set)
            ax.x_range = [0, len(data_set.data[0])]
            ax.energy_range = energy_range
            ax.display()
            ax.create_colourbar()
        else:
            self.add_plot(axes_index)
            self.axes[-1].add_data_set(data_set)
            self.axes[-1].x_range = [0, len(data_set.data[0])]
            self.axes[-1].energy_range = energy_range
            self.axes[-1].display()
            self.axes[-1].create_colourbar()

    def plot_wan_band_dip_SOC(self, wan_dir, band_index, lambda_socs, centres, direction=2, k_points=None, k_file=None, axes_index=(1, 1), label="", **kwargs):
        ax = self.find_axes(axes_index)
        if self.data_storage.eigenvalues is None:
            self.data_storage.load_wannier_info(wan_dir)
            if k_points == 'plot' and ax != 0:
                k = self.get_k_points(ax)
                self.data_storage.load_kpoints(k_file, k)
            else:
                self.data_storage.load_kpoints(k_file, k_points)
            self.data_storage.calculate_dipoles_eigenvalues(
                lambda_socs, centres)
        indices = np.arange(0, len(self.data_storage.dipoles))
        dipoles = []
        for dipole in self.data_storage.dipoles:
            dipoles.append(dipole[band_index][direction].real)
        set_ = np.transpose(np.vstack((indices, dipoles)))
        data_set = DataSet("Plot", set_, **kwargs)
        data_set.label = label
        if ax != 0:
            ax.add_data_set(data_set)
            ax.display()
        else:
            self.add_plot(axes_index)
            self.axes[-1].add_data_set(data_set)
            self.axes[-1].display()

    def plot_wan_band_eigval_SOC(self, wan_dir, band_index, lambda_socs, centres, k_points=None, k_file=None, axes_index=(1, 1), label="", **kwargs):
        ax = self.find_axes(axes_index)
        if self.data_storage.eigenvalues is None:
            self.data_storage.load_wannier_info(wan_dir)
            if k_points == 'plot' and ax != 0:
                k = self.get_k_points(ax)
                self.data_storage.load_kpoints(k_file, k)
            else:
                self.data_storage.load_kpoints(k_file, k_points)
            self.data_storage.calculate_dipoles_eigenvalues(
                lambda_socs, centres)
        indices = np.arange(0, len(self.data_storage.dipoles))
        eigenvalues = []
        for eigenval in self.data_storage.eigenvalues:
            eigenvalues.append(eigenval[band_index])
        set_ = np.transpose(np.vstack((indices, eigenvalues)))
        data_set = DataSet("Plot", set_, **kwargs)
        data_set.label = label
        if ax != 0:
            ax.add_data_set(data_set)
            ax.display()
        else:
            self.add_plot(axes_index)
            self.axes[-1].add_data_set(data_set)
            self.axes[-1].display()

    def plot_wan_band_dip_eigval_SOC(self, wan_dir, band_index, lambda_socs, centres, direction=2, k_points=None, k_file=None, axes_indices=[(1, 1), (1, 2)], label="", **kwargs):
        self.plot_wan_band_dip_SOC(wan_dir, band_index,
                                   lambda_socs, centres, k_points, k_file, direction, axes_indices[0], label, **kwargs)
        self.plot_wan_band_eigval_SOC(wan_dir, band_index,
                                      lambda_socs, centres, k_points, k_file, axes_indices[1], label, **kwargs)

    def plot_wan_band_eigvec_SOC(self, wan_dir, band_index, c_index, lambda_socs, centres, k_points=None, k_file=None, axes_index=(1, 1), label="", **kwargs):

        ax = self.find_axes(axes_index)
        if self.data_storage.eigenvalues is None:
            self.data_storage.load_wannier_info(wan_dir)
            if k_points == 'plot' and ax != 0:
                k = self.get_k_points(ax)
                self.data_storage.load_kpoints(k_file, k)
            else:
                self.data_storage.load_kpoints(k_file, k_points)
            self.data_storage.calculate_dipoles_eigenvalues(
                lambda_socs, centres)
        indices = np.arange(0, len(self.data_storage.eigenvectors))
        tmp_vec = []
        for vec in self.data_storage.eigenvectors:
            tmp_vec.append(
                np.sqrt(vec[c_index].real**2 + vec[c_index].imag**2))
        set_ = np.transpose(np.vstack((indices, tmp_vec)))
        data_set = DataSet("Plot", set_, **kwargs)
        data_set.label = label
        if ax != 0:
            ax.add_data_set(data_set)
            ax.display()
        else:
            self.add_plot(axes_index)
            self.axes[-1].add_data_set(data_set)
            self.axes[-1].display()

    def plot_wan_eigvals_SOC(self, wan_dir, lambda_socs, centres, k_points=None, k_file=None, axes_index=(1, 1), label="", **kwargs):
        ax = self.find_axes(axes_index)
        if self.data_storage.eigenvalues is None:
            self.data_storage.load_wannier_info(wan_dir)
            if k_points == 'plot' and ax != 0:
                k = self.get_k_points(ax)
                self.data_storage.load_kpoints(k_file, k)
            else:
                self.data_storage.load_kpoints(k_file, k_points)
            self.data_storage.calculate_dipoles_eigenvalues(
                lambda_socs, centres)
        indices = np.arange(0, len(self.data_storage.eigenvalues))
        if ax == 0:
            self.add_plot(axes_index)
        ax = self.axes[-1]
        for i in range(0, len(self.data_storage.eigenvalues[0])):
            tmp_set = []
            for val in self.data_storage.eigenvalues:
                tmp_set.append(val[i])
            tmp_set1 = np.transpose(np.vstack((indices, tmp_set)))
            data_set = DataSet("Plot", tmp_set1, **kwargs)
            data_set.label = label
            ax.add_data_set(data_set)
        ax.display()

    def plot_wan_angmom(self, wan_dir, band_index, centres, l="Lx", k_points=None, k_file=None, axes_index=(1, 1), **kwargs):
        if self.data_storage.angmoms is None:
            self.data_storage.load_wannier_info(wan_dir)
            self.data_storage.load_kpoints(k_file, k_points)
            self.data_storage.calc_angmoms(centres)
        tmp_moms = []
        for mom in self.data_storage.angmoms:
            if l == "Lx":
                tmp_moms.append(mom[band_index][0].real)
            elif l == "Ly":
                tmp_moms.append(mom[band_index][1].real)
            elif l == "Lz":
                tmp_moms.append(mom[band_index][2].real)

        self.plot_data(tmp_moms, d_type='Plot', axes_index=axes_index, **kwargs)

    def get_k_points(self, ax):
        return [ax.data_sets[0].data[0][:, 5], ax.data_sets[0].data[0][:, 6], ax.data_sets[0].data[0][:, 7]]

    def set_fermi_level(self, index=None, filename=None, energy=None, axes_index=(1, 1)):
        if filename is not None:
            found = False
            with open(filename, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if "Fermi" in line:
                        found = True
                        fermi_level = float(line.split()[4])
                        ax = self.find_axes(axes_index)
                        if ax == 0:
                            print("Error: There is no plot with that index.")
                        else:
                            ax.set_fermi_level(fermi_level, index)
                            ax.display()
                if not found:
                    print("Error: could not find Fermi level in specified file!")
                    return
        elif energy is not None:
            ax = self.find_axes(axes_index) 
            if ax == 0:
                print("Error: There is no plot with that index.")
                return
            else:
                ax.set_fermi_level(energy, index)
                ax.display()
        else:
            print("Error: specify either an scf output file or the fermi level!")
            return
