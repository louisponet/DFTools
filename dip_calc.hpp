//
//  dipole_calculations.hpp
//  Python
//
//  Created by Louis Ponet on 14/02/2017.
//  Copyright © 2017 Louis Ponet. All rights reserved.
//

#ifndef dipole_calculations_hpp
#define dipole_calculations_hpp

#include <stdio.h>
#include "typedefs.hpp"


vec_Cd calculate_intracell_dipole(array_Cd mesh1, array_Cd mesh2, array_d grid);
vec_Cd calculate_intracell_dipole_(Cd* meshPtr1, Cd* meshPtr2, double* grid,std::vector<Ul> shape);
vec_Cd calculate_intercell_dipole(array_Cd mesh1, array_Cd mesh2, array_d grid, Direction d);
vec_Cd calculate_intercell_dipole_(Cd* meshPtr1, Cd* meshPtr2, double* gridPtr,std::vector<Ul> shape, Direction d);

std::tuple<std::vector<vec_Cd>,Eigen::VectorXcd> calculate_eigen_dipoles(Eigen::Matrix<Cd,Eigen::Dynamic,Eigen::Dynamic> hami, std::vector<vec_Cd> dipoles,std::vector<vec_Cd> intra_dip, float kx, float ky, float kz);
void init_dip_calc(py::module& m);
#endif /* dipole_calculations_hpp */
