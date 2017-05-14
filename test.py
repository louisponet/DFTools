#%%
from DFTools import wrap,fileprc
import time
#%%
mesh, prim_cell, origin, a_vec, b_vec, c_vec, a_array, b_array, c_array = fileprc.read_wannier_function("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest/03GeTe_up_00001.xsf")
a_span_array=[a * a_vec for a in a_array]
b_span_array=[b * b_vec for b in b_array]
c_span_array=[c * c_vec for c in c_array]

grid = wrap.generate_grid(
            origin, a_span_array, b_span_array, c_span_array)
mesh1,prim_cell, origin, a_vec, b_vec, c_vec, a_array, b_array, c_array = fileprc.read_wannier_function("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/xsftest/03GeTe_up_00004.xsf")
t = time.time()
print(wrap.calculate_angular_momentum(mesh,mesh1,grid,[0.0,0.0,0.1114785],wrap.AngMomentum.Lx))
print(time.time()-t)

#%%


