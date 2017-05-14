import numpy as np
from subprocess import call
import os
import cmath

def create_bands_from_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        bands = False
        spin_ud = False
        i = 0
        k_points = np.empty((0, 3), np.double)
        arrays = []
        while i < len(lines):
            if "cryst." in lines[i].split() and len(lines[i].split()) == 2:
                i1 = 1
                while lines[i + i1] != '\n':
                    line = lines[i + i1].split()
                    kx = np.double(line[4])
                    ky = np.double(line[5])
                    kz = np.double(line[6][:-2])
                    k_points = np.vstack((k_points, np.asarray([kx, ky, kz], np.double)))
                    i1 += 1
                i += i1 - 1
            if "k =" in lines[i] and "bands (ev)" in lines[i]:
                tmp_array = []
                i1 = 2
                while lines[i + i1] != '\n':
                    tmp_array.append(lines[i + i1].split())
                    i1 += 1
                i += i1 - 1
                tmp_array2 = []
                for x in tmp_array[:]:
                    for i2 in x[:]:
                        tmp_array2.append(np.double(i2))

                arrays.append(tmp_array2)

            if " ------ SPIN DOWN ----------" in lines[i]:
                spin_ud = True
            if "End of band structure calculation" in lines[i]:
                bands = True
            i += 1
        if not bands:
            print("Error: please select a valid bands output file!")
            return

        if not spin_ud:
            bands_tmp = []
            i = 0
            while i < len(arrays[0]):
                i2 = 0
                band = []
                while i2 < len(arrays):
                    band.append(arrays[i2][i])
                    i2 += 1
                np_band = np.arange(0, len(band), 1, dtype=np.double)
                np_band = np.vstack((np_band, np.asarray(band, dtype=np.double)))
                np_band = np.transpose(np.vstack((np_band, np.zeros((3, len(band)), dtype=np.double))))
                np_band = np.hstack((np_band, k_points))
                bands_tmp.append(np_band)
                i += 1
            bands = [bands_tmp]
        else:
            bands1 = []
            bands2 = []
            i = 0
            while i < len(arrays[0]):
                i2 = 0
                band = []
                while i2 < len(arrays) / 2:
                    band.append(arrays[i2][i])
                    i2 += 1
                np_band = np.arange(0, len(band), 1, dtype=np.double)
                np_band = np.vstack((np_band, np.asarray(band, dtype=np.double)))
                np_band = np.transpose(np.vstack((np_band, np.zeros((3, len(band)), dtype=np.double))))
                np_band = np.hstack((np_band, k_points))
                bands1.append(np_band)
                i += 1
            i = 0
            while i < len(arrays[0]):
                i2 = int(len(arrays) / 2)
                band = []
                while i2 < len(arrays):
                    band.append(arrays[i2][i])
                    i2 += 1
                np_band = np.arange(0, len(band), 1, dtype=np.double)
                np_band = np.vstack((np_band, np.asarray(band, dtype=np.double)))
                np_band = np.transpose(np.vstack((np_band, np.zeros((3, len(band)), dtype=np.double))))
                np_band = np.hstack((np_band, k_points))
                bands2.append(np_band)
                i += 1
            bands = [bands1, bands2]
        return bands

def extract_ks_from_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        k_points = np.empty((0, 3), np.double)
        bands = False
        i = 0
        while i < len(lines):
            if "cryst." in lines[i].split() and len(lines[i].split()) == 2:
                i1 = 1
                while lines[i + i1] != '\n':
                    line = lines[i + i1].split()
                    kx = np.double(line[4])
                    ky = np.double(line[5])
                    kz = np.double(line[6][:-2])
                    k_points = np.vstack((k_points, np.asarray([kx, ky, kz], np.double)))
                    i1 += 1
                i += i1 - 1
            if "End of band structure calculation" in lines[i]:
                bands = True
            i += 1
        if bands:
            return k_points
        else:
            print("Error:Please enter a valid bandscalculation output file!")


def read_pdos_colin(filename, atoms=None):
    wfc = {}
    k_point_arrays = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            if "state #" in lines[i]:
                i1 = 0
                while "state #" in lines[i + i1]:
                    wfc[i1 + 1] = (lines[i + i1].split()[5][1:], int(lines[i + i1].split()[8]),
                                   int(lines[i + i1].split()[-1][:-1]))
                    i1 += 1
                i += i1

            if " k =" in lines[i]:
                if len(wfc) > 18:
                    k_point_array = np.empty((0, 5))
                else:
                    k_point_array = np.empty((0, 4))
                i1 = 1
                tmp_array = np.empty((1, 4))
                while " k =" not in lines[i + i1] and lines[i + i1] != "\n":
                    if " e = " in lines[i + i1]:
                        tmp_array[0, 0] = lines[i + i1].split()[-2]
                    elif " psi = " in lines[i + i1]:
                        line = lines[i + i1].split("]+")[:-1]
                        s_value = 0
                        p_value = 0
                        d_value = 0
                        for part in line:
                            segs = part.split()
                            contr = np.double(segs[-2][:-3])
                            if atoms is None:
                                if wfc[int(segs[-1])][1] == 1:
                                    s_value += np.double(contr)
                                elif wfc[int(segs[-1])][1] == 2:
                                    p_value += np.double(contr)
                                elif wfc[int(segs[-1])][1] == 3:
                                    d_value += np.double(contr)
                            else:
                                if wfc[int(segs[-1])][1] == 1 and wfc[int(segs[-1])][0] in atoms:
                                    s_value += np.double(contr)
                                elif wfc[int(segs[-1])][1] == 2 and wfc[int(segs[-1])][0] in atoms:
                                    p_value += np.double(contr)
                                elif wfc[int(segs[-1])][1] == 3 and wfc[int(segs[-1])][0] in atoms:
                                    d_value += np.double(contr)
                        if s_value > 1:
                            s_value -= 0.0011
                        if p_value > 1:
                            p_value -= 0.0011
                        if d_value > 1:
                            d_value -= 0.0011
                        tmp_array[0, 1] = s_value
                        tmp_array[0, 2] = p_value
                        tmp_array[0, 3] = d_value

                        k_point_array = np.vstack((k_point_array, tmp_array))
                        tmp_array = np.empty((1, 4))
                    i1 += 1

                k_point_arrays.append(k_point_array)
                i += i1 - 1
            i += 1
    return wfc, k_point_arrays


def read_pdos_SOC(filename, atoms=None):
    os.system('/Users/louisponet/Documents/Fysica/PhD/Python/Python/projwfc_processing {0} {1}'.format(filename, filename + '.pp'))
    wfc_array = None
    atom_array = []
    full_array = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            if "state #" in lines[i]:
                if len(lines[i]) > 54:
                    wfc_array = np.empty((0, 5))
                    i1 = 0
                    while "state #" in lines[i + i1]:
                        line = lines[i + i1].split()
                        if len(line) == 12:
                            array = [line[4], line[8], line[9][-3:], line[10][-1], line[11][-5:-1]]
                        else:
                            array = [line[4], line[8], line[9][-3:], line[10][-1], line[12][:-1]]
                        wfc_array = np.vstack((wfc_array, np.asarray(array, dtype=np.np.double)))

                        i1 += 1
                    i += i1
                else:
                    wfc_array = np.empty((0, 4))
                    i1 = 0
                    while "state #" in lines[i + i1]:
                        line = lines[i + i1].split()
                        array = [line[4], line[8], line[9][-1], line[10][-2]]
                        wfc_array = np.vstack((wfc_array, np.asarray(array, dtype=np.int)))
                        i1 += 1
                    i += i1
            i += 1

    with open(filename + '.pp', 'r') as f:
        lines = f.readlines()
        k_array = []
        i = 0
        while i < len(lines):
            line = lines[i]
            if line == '\n':
                band_array = []
                i1 = 1
                while i + i1 < len(lines) and lines[i + i1] != '\n':
                    band = int(np.double(lines[i + i1].split()[0]))
                    line = lines[i + i1].split()[1:]
                    tmp = []
                    for (x, y) in zip(line[::2], line[1::2]):
                        tmp.append((np.double(x), np.double(y)))
                    band_array.append((band, tmp))
                    i1 += 1
                i += i1 - 1
                k_array.append(band_array)
            i += 1
    return wfc_array, k_array


def read_wan_hami(filename):
    out_list = []
    with open(filename, 'r') as f:
        for line in f.readlines():
            segs = line.split()
            if len(segs) == 7:
                tmp_array = (int(segs[0]), int(segs[1]), int(segs[2]), int(segs[3]), int(segs[4]), np.double(segs[5]),
                             np.double(segs[6]))
                out_list.append(tmp_array)
    return out_list


def read_wannier_function(filename):
    lines = None
    with open(filename, 'r') as f:
        lines = f.readlines()
    prim_cell = np.empty((0, 3), dtype=np.double)
    numbers = []
    nx, ny, nz = 0, 0, 0
    origin = None
    a_vec = None
    b_vec = None
    c_vec = None
    i = 0

    while i < len(lines):
        if lines[i] == "PRIMVEC\n":
            i += 1
            prim_cell = np.vstack((prim_cell, np.asarray(lines[i].split(), dtype=np.double)))
            prim_cell = np.vstack((prim_cell, np.asarray(lines[i + 1].split(), dtype=np.double)))
            prim_cell = np.vstack((prim_cell, np.asarray(lines[i + 2].split(), dtype=np.double)))
            i += 1
        if lines[i] == " DATAGRID_3D_DENSITY\n":
            i += 1
            nx, ny, nz = int(lines[i].split()[0]), int(lines[i].split()[1]), int(lines[i].split()[2])
            origin = np.asarray(lines[i + 1].split(), dtype=np.double)
            a_vec = np.asarray(lines[i + 2].split(), dtype=np.double)
            b_vec = np.asarray(lines[i + 3].split(), dtype=np.double)
            c_vec = np.asarray(lines[i + 4].split(), dtype=np.double)
            i += 5
            i1 = i
            while lines[i1] != "\n":
                numbers.append(np.complex128(float(lines[i1])))
                i1 += 1
            i += i1 - i
        if lines[i] == "BEGIN_DATAGRID_3D_UNKNOWN\n":
            i += 1
            nx, ny, nz = int(lines[i].split()[0]), int(lines[i].split()[1]), int(lines[i].split()[2])
            origin = np.asarray(lines[i + 1].split(), dtype=np.double)
            a_vec = np.asarray(lines[i + 2].split(), dtype=np.double)
            b_vec = np.asarray(lines[i + 3].split(), dtype=np.double)
            c_vec = np.asarray(lines[i + 4].split(), dtype=np.double)
            i += 5
            i1 = i
            while lines[i1] != "END_DATAGRID_3D\n":
                i2 = 0
                while i2 < len(lines[i1].split()):
                    numbers.append(np.complex128(float(lines[i1].split()[i2])))
                    i2 += 1
                i1 += 1
            i += i1 - i
        i += 1
    mesh = np.empty((nz, ny, nx), dtype=np.complex128)
    i = 0
    while i < nz:
        i1 = 0
        while i1 < ny:
            i2 = 0
            while i2 < nx:
                mesh[i, i1, i2] = numbers[nx * ny * i + nx * i1 + i2]
                i2 += 1
            i1 += 1
        i += 1
    a_array = np.linspace(0, 1, nx, dtype=np.double)
    b_array = np.linspace(0, 1, ny, dtype=np.double)
    c_array = np.linspace(0, 1, nz, dtype=np.double)

    return mesh, prim_cell, origin, a_vec, b_vec, c_vec, a_array, b_array, c_array

def read_pdos(filename,column=0):
    array = np.genfromtxt(filename, skip_header=True)
    return np.asarray([[x,y] for (x, y) in zip(array[:,0],array[:,column])],dtype=np.double)

def read_k_pdos(filename, column=0):
    array = np.genfromtxt(filename, skip_header=True)
    array1 = np.asarray([[x, y, z] for (x, y, z) in zip(array[:, 0], array[:, 1], array[:, column])], dtype=np.double)
    return array1

def read_fatbands(filename):
    bands_array = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            line = lines[i]
            if "BAND number" in line:
                array = []
                i1 = 1
                while lines[i + i1] != '\n' and "&" not in lines[i + i1]:
                    line_s = lines[i + i1].split()
                    array.append([i1 - 1, line_s[1], line_s[2]])
                    i1 += 1
                bands_array.append(np.asarray(array, dtype=np.double))
                i += i1 - 1
            i += 1
    return bands_array


def save_wannier_xsf(filename, mesh, grid):
    shape = np.shape(grid)
    with open(filename + '.xsf', 'w') as f:
        f.writelines(["# Generated from PhD calculations\n", "", "BEGIN_BLOCK_DATAGRID_3D\n", "3D_FIELD\n",
                      "BEGIN_DATAGRID_3D_UNKNOWN\n"])
        f.write("{0}    {1}    {2}\n".format(shape[0], shape[1], shape[2]))
        f.write("{0}    {1}    {2}\n".format(grid[0, 0, 0, 0], grid[0, 0, 0, 1], grid[0, 0, 0, 2]))
        f.write("{0}    {1}    {2}\n".format(grid[0, 0, -1, 0] - grid[0, 0, 0, 0], grid[0, 0, -1, 1] - grid[0, 0, 0, 1], grid[0, 0, -1, 2] - grid[0, 0, 0, 2]))
        f.write("{0}    {1}    {2}\n".format(grid[0, -1, 0, 0] - grid[0, 0, 0, 0], grid[0, -1, 0, 1] - grid[0, 0, 0, 1], grid[0, -1, 0, 2] - grid[0, 0, 0, 2]))
        f.write("{0}    {1}    {2}\n".format(grid[-1, 0, 0, 0] - grid[0, 0, 0, 0], grid[-1, 0, 0, 1] - grid[0, 0, 0, 1], grid[-1, 0, 0, 2] - grid[0, 0, 0, 2]))
        for i in range(0, shape[0]):
            for i1 in range(0, shape[1]):
                for i2 in range(0, shape[2]):
                    if isinstance(mesh[0, 0, 0], np.complex128):
                        f.write("{} ".format(np.linalg.norm(mesh[i, i1, i2])))
                    else:
                        f.write("{} ".format(mesh[i, i1, i2]))
                f.write('\n')
            f.write('\n')
        f.writelines(["END_DATAGRID_3D\n", "END_BLOCK_DATAGRID_3D\n"])
    if not isinstance(mesh[0, 0, 0], np.complex128):
        return
    with open(filename + '.ph.xsf', 'w') as f:
        f.writelines(["# Generated from PhD calculations\n", "", "BEGIN_BLOCK_DATAGRID_3D\n", "3D_FIELD\n",
                      "BEGIN_DATAGRID_3D_UNKNOWN\n"])
        f.write("{0}    {1}    {2}\n".format(shape[0], shape[1], shape[2]))
        f.write("{0}    {1}    {2}\n".format(grid[0, 0, 0, 0], grid[0, 0, 0, 1], grid[0, 0, 0, 2]))
        f.write("{0}    {1}    {2}\n".format(grid[0, 0, -1, 0] - grid[0, 0, 0, 0], grid[0, 0, -1, 1] - grid[0, 0, 0, 1], grid[0, 0, -1, 2] - grid[0, 0, 0, 2]))
        f.write("{0}    {1}    {2}\n".format(grid[0, -1, 0, 0] - grid[0, 0, 0, 0], grid[0, -1, 0, 1] - grid[0, 0, 0, 1], grid[0, -1, 0, 2] - grid[0, 0, 0, 2]))
        f.write("{0}    {1}    {2}\n".format(grid[-1, 0, 0, 0] - grid[0, 0, 0, 0], grid[-1, 0, 0, 1] - grid[0, 0, 0, 1], grid[-1, 0, 0, 2] - grid[0, 0, 0, 2]))
        for i in range(0, shape[0]):
            for i1 in range(0, shape[1]):

                for i2 in range(0, shape[2]):
                    f.write("{} ".format(cmath.phase(mesh[i, i1, i2])))

                f.write('\n')
            f.write('\n')
        f.writelines(["END_DATAGRID_3D\n", "END_BLOCK_DATAGRID_3D\n"])
