import numpy as np
import glob
import DFTools.plotter as plotter
import DFTools.fileprc as fileprc

kxs = np.hstack(([0.5 for i in range(0, 100)], np.linspace(0.5, 0.4875, 100)))
kys = np.hstack((np.linspace(0.525, 0.5, 100), np.linspace(0.5, 0.525, 100)))
kzs = np.hstack((np.linspace(0.475, 0.5, 100), np.linspace(0.5, 0.4875, 100)))
k_points = [kxs, kys, kzs]

meshes = []
meshes1 = []
prim_cells = []
origins = []
a_span_vecs = []
b_span_vecs = []
c_span_vecs = []
a_span_arrays = []
b_span_arrays = []
c_span_arrays = []


for i in range(1, 201):
    mesh, prim_cell, origin, a_vec, b_vec, c_vec, a_array, b_array, c_array = fileprc.read_wannier_function("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/cut3D/b19_20/wfcnorm_k{}_b19_spinor1".format(i))

    meshes.append(mesh)
    prim_cells.append(prim_cell)
    origins.append(origin)
    a_span_vecs.append(a_vec)
    b_span_vecs.append(b_vec)
    c_span_vecs.append(c_vec)
    a_span_arrays.append([a * a_vec for a in a_array])
    b_span_arrays.append([b * b_vec for b in b_array])
    c_span_arrays.append([c * c_vec for c in c_array])
    mesh_, prim_cell_, origin_, a_vec_, b_vec_, c_vec_, a_array_, b_array_, c_array_ = fileprc.read_wannier_function(
        "/Users/ponet/Documents/Fysica/PhD/GeTe/colin/cut3D/b19_20/wfcnorm_k{}_b20_spinor1".format(i))
    meshes1.append(mesh_)

origin = origins[0]
a_span_array = a_span_arrays[0]
b_span_array = b_span_arrays[0]
c_span_array = c_span_arrays[0]

grid = plotter.wrap.generate_grid(origin, a_span_array, b_span_array, c_span_array)

dipoles = []
dipoles1 = []
for mesh in meshes:
    dipoles.append(plotter.wrap.calculate_intracell_dipole(mesh, mesh, grid))
for mesh in meshes1:
    dipoles1.append(plotter.wrap.calculate_intracell_dipole(mesh, mesh, grid))
set_ = []
set_1 = []
for dip in dipoles:
    set_.append(dip[2].real)
for dip in dipoles1:
    set_1.append(dip[2].real)
plot = plotter.Plotter(w=15)
plot.plot_set(set_, c='b')
plot.plot_set(set_1, c='r')
# plot.plot_set(set_01,colour='m')
# plot.plot_set(set_11,colour='m')
plot.save("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/cut3D/b19_20/dipoles_s1.png")
