#%%
import os
os.chdir('./')
import numpy as np
import DFTools.wrap as wrap
import DFTools.fileprc as fileprc

#%%
mesh1 = wrap.normalize(fileprc.read_wannier_function("tests/colin/xsftest/03GeTe_up_00001.xsf")[0])
mesh2 = wrap.normalize(fileprc.read_wannier_function("tests/colin/xsftest/03GeTe_up_00002.xsf")[0])
print(wrap.calculate_overlap_pot(mesh1,mesh2,wrap.Direction.END))