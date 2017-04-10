import DFTools.fileprc as fileprc
import DFTools.wrap as wrap
import glob
import numpy as np

bands = fileprc.create_bands_from_file("/Users/ponet/Documents/Fysica/PhD/GeTe/fullrel/GeTe_bands.out")
kxs = bands[0][0][:,5]
kys = bands[0][0][:,6]
kzs = bands[0][0][:,7]

files = glob.glob("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/bestxsf/*.xsf")
meshes = []
prim_cells = []
origins = []
a_span_vecs = []
b_span_vecs = []
c_span_vecs = []
a_span_arrays = []
b_span_arrays = []
c_span_arrays = []
for file in files:
   mesh,prim_cell,origin,a_vec,b_vec,c_vec,a_array,b_array,c_array = fileprc.read_wannier_function(file)
   meshes.append(mesh)
   prim_cells.append(prim_cell)
   origins.append(origin)
   a_span_vecs.append(a_vec)
   b_span_vecs.append(b_vec)
   c_span_vecs.append(c_vec)
   a_span_arrays.append([a*a_vec for a in a_array])
   b_span_arrays.append([b*b_vec for b in b_array])
   c_span_arrays.append([c*c_vec for c in c_array])


origin = origins[0]
a_span_array = a_span_arrays[0]
b_span_array = b_span_arrays[0]
c_span_array = c_span_arrays[0]
grid = wrap.generate_grid(origin,a_span_array,b_span_array,c_span_array)

hami_raw = fileprc.read_wan_hami("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/bestxsf/03GeTe_up_hr.dat")

hamiltonian = wrap.hami_from_k(hami_raw, kxs[100],kys[100],kzs[100])
x = -0.09
l = [x for i in range(0,10)]
full_hamis = wrap.construct_SO_hamiltonians(hami_raw,kxs[99:100],kys[99:100],kzs[99:100],meshes,grid,l,l,l)

eigenvalues,eigenvectors = np.linalg.eig(full_hamis[0])
print(eigenvalues)
for i in [-5,-6,-9,-10,-11,-12]:
    eigenvector = eigenvectors[i]
    eigenmesh = np.empty_like(meshes[0])
    for i1 in range(0,len(meshes)):
        eigenmesh += eigenvector[i1]*meshes[i1]+eigenvector[i1+len(meshes)]*meshes[i1]

    fileprc.save_wannier_xsf("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/eigenmeshesSO/eigenmesh_v{}".format(i),eigenmesh,grid)
