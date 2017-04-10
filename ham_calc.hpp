//
//  ham_calc.hpp
//  Python
//
//  Created by Louis Ponet on 14/02/2017.
//  Copyright © 2017 Louis Ponet. All rights reserved.
//

#ifndef ham_calc_hpp
#define ham_calc_hpp

#include <stdio.h>
#include "typedefs.hpp"
Eigen::Matrix<Cf,Eigen::Dynamic,Eigen::Dynamic> hami_from_k(std::vector<std::tuple<int,int,int,int,int,float,float> > hami_raw,float kx,float ky, float kz);

std::vector<Eigen::Matrix<Cd,Eigen::Dynamic,Eigen::Dynamic> > calculate_L_hamiltonians(std::vector<array_Cd> meshes_, array_d grid,std::vector<array_d>centres,vec_d lam_x,vec_d lam_y, vec_d lam_z);
std::vector<Eigen::Matrix<Cd, Eigen::Dynamic, Eigen::Dynamic> > construct_SO_hamiltonians(std::vector<std::tuple<int,int,int,int,int,float,float> > hami_raw,vec_f kxs,vec_f kys, vec_f kzs,std::vector<array_Cd> meshes, array_d grid,std::vector<array_d>centres, vec_d lam_x,vec_d  lam_y,vec_d lam_z);

void init_ham_calc(py::module& m);
#endif /* ham_calc_hpp */
