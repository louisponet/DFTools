//
//  sup_calc.hpp
//  Python
//
//  Created by Louis Ponet on 14/02/2017.
//  Copyright © 2017 Louis Ponet. All rights reserved.
//

#ifndef sup_calc_hpp
#define sup_calc_hpp

#include <stdio.h>
#include "typedefs.hpp"

std::vector<Cd *> create_meshPtrs(std::vector<array_Cd> meshes_,std::vector<Ul> shape);
std::vector<double *> create_centrePtrs(std::vector<array_d> centres_);

array_Cd construct_px(array_d grid);
array_Cd construct_py(array_d grid);
array_Cd construct_pz(array_d grid);


array_d generate_grid(array_d origin,array_d a_span_array,array_d b_span_array,array_d c_span_array);

Cd inner_product_3D(array_Cd mesh1,array_Cd mesh2);
array_Cd normalize(array_Cd mesh);
array_Cd construct_bloch_sum_center(array_Cd mesh, double kx, double ky, double kz);
array_Cd construct_bloch_sum_origin(array_Cd mesh, double kx, double ky, double kz);

Cd calculate_angular_momentum(array_Cd mesh1,array_Cd mesh2, array_d grid, array_d centre, AngMomentum angMom);
Cd calculate_angular_momentum_(Cd* mesh1Ptr,Cd* mesh2Ptr, double* gridPtr,double* centrePtr, std::vector<Ul> shape, AngMomentum angMom);
array_Cd project_on_angular_momentum(array_Cd mesh, array_d grid,array_d centre,AngMomentum angMom);
array_Cd project_on_momentum(array_Cd mesh1, array_d grid,Momentum mom);
Cd calculate_overlap_pot(array_Cd mesh1, array_Cd mesh2,Direction d);

vec_Cd calculate_center(array_Cd mesh,array_d grid);

void init_sup_calc(py::module& m);
#endif /* sup_calc_hpp */
