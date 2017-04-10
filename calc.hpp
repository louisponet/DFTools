//
//  gen_calc.hpp
//  Python
//
//  Created by Louis Ponet on 14/02/2017.
//  Copyright © 2017 Louis Ponet. All rights reserved.
//

#ifndef gen_calc_hpp
#define gen_calc_hpp

#include <stdio.h>
#include "typedefs.hpp"

std::vector<std::tuple<Eigen::VectorXd,std::vector<Eigen::VectorXcd> > > generate_wannier_bands(std::vector<std::tuple<int,int,int,int,int,float,float> > hami_raw, vec_f kxs,vec_f kys, vec_f kzs);
std::vector<std::tuple<Eigen::VectorXd,std::vector<Eigen::VectorXcd> > >  generate_wannier_bands_SOC(std::vector<std::tuple<int,int,int,int,int,float,float> > hami_raw,vec_f kxs,vec_f kys, vec_f kzs,std::vector<array_Cd> meshes_, array_d grid, std::vector<array_d> centres, vec_d lam_x,vec_d  lam_y,vec_d lam_z);
std::vector<std::tuple<std::vector<vec_Cd>,Eigen::VectorXcd > > calculate_tb_dipoles(std::vector<std::tuple<int,int,int,int,int,float,float> > hami_raw,vec_f kxs,vec_f kys, vec_f kzs,std::vector<array_Cd> meshes_, array_d grid, std::vector<array_d> centres,vec_d lam_x,vec_d  lam_y,vec_d lam_z);
std::vector<std::vector<vec_Cd> > calculate_tb_angmoms(std::vector<std::tuple<int,int,int,int,int,float,float> > hami_raw,vec_f kxs,vec_f kys, vec_f kzs,std::vector<array_Cd> meshes_, array_d grid, std::vector<array_d> centres_);
std::vector<array_Cd> construct_eigenmesh_SOC(std::vector<std::tuple<int,int,int,int,int,float,float> > hami_raw, float kx, float ky, float kz, std::vector<array_Cd> meshes_, array_d grid, vec_d ls, std::vector<array_d> centres_);
void init_calc(py::module& m);
#endif /* gen_calc_hpp */
