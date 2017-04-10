//
//  typedefs.h
//  Python
//
//  Created by Louis Ponet on 22/01/2017.
//  Copyright © 2017 Louis Ponet. All rights reserved.
//

#ifndef typedefs_h
#define typedefs_h

#include <vector>
#include "pybind11/pybind11.h"
#include "pybind11/eigen.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"



namespace py = pybind11;

typedef unsigned long Ul;

typedef std::complex<float> Cf;
typedef std::complex<double> Cd;

typedef std::vector<float> vec_f;
typedef std::vector<double> vec_d;
typedef std::vector<Cf> vec_Cf;
typedef std::vector<Cd> vec_Cd;

typedef py::array_t<float> array_f;
typedef py::array_t<double> array_d;
typedef py::array_t<Cd,py::array::c_style> array_Cd;

enum Direction{
    Ex,
    Ex_,
    Ey,
    Ey_,
    Ez,
    Ez_,
    ExEx,
    Ex_Ex_,
    ExEy,
    Ex_Ey_,
    Ex_Ey,
    ExEy_,
    ExEz,
    Ex_Ez_,
    ExEz_,
    Ex_Ez,
    EyEz,
    Ey_Ez_,
    Ey_Ez,
    EyEz_,
    ExEyEz,
    Ex_Ey_Ez_,
    END
};
enum AngMomentum{
    Lx,
    Ly,
    Lz
};

enum Momentum{
    Px,
    Py,
    Pz
};
#endif /* typedefs_h */
