//
//  ham_calc.cpp
//  Python
//
//  Created by Louis Ponet on 14/02/2017.
//  Copyright © 2017 Louis Ponet. All rights reserved.
//
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <omp.h>
#include "ham_calc.hpp"
#include "sup_calc.hpp"

Eigen::Matrix< Cf, Eigen::Dynamic, Eigen::Dynamic > hami_from_k(std::vector< std::tuple< int, int, int, int, int, float, float > > hami_raw, float kx, float ky, float kz) {

    int dim = 0;
    for (unsigned long i = 0; i < hami_raw.size(); i++) {
        std::tuple< int, int, int, int, int, float, float > h = hami_raw[i];
        int dim_temp = std::get< 3 >(h);
        if (dim_temp > dim)
            dim = dim_temp;
        else
            break;
    }
    Eigen::Matrix< Cf, Eigen::Dynamic, Eigen::Dynamic > out = Eigen::Matrix< Cf, Eigen::Dynamic, Eigen::Dynamic >::Zero(dim, dim);
    for (int i = 0; i < (int)hami_raw.size(); i++) {
        std::tuple< int, int, int, int, int, float, float > h = hami_raw[i];
        float complex_part = 2 * M_PI * (kx * std::get< 0 >(h) + ky * std::get< 1 >(h) + kz * std::get< 2 >(h));
        out(std::get< 3 >(h) - 1, std::get< 4 >(h) - 1) += Cf(std::get< 5 >(h), std::get< 6 >(h)) * std::exp(Cf(0, complex_part));
    }
    return out;
}

std::vector< Eigen::Matrix< Cd, Eigen::Dynamic, Eigen::Dynamic > > calculate_L_hamiltonians(std::vector< array_Cd > meshes_, array_d grid, std::vector< array_d > centres_, vec_d lam_x, vec_d lam_y, vec_d lam_z) {

    py::buffer_info gridInfo = grid.request();
    std::vector< Ul > shape = gridInfo.shape;
    double* gridPtr = static_cast< double* >(gridInfo.ptr);

    std::vector< Cd* > meshes = create_meshPtrs(meshes_, shape);
    std::vector< double* > centres = create_centrePtrs(centres_);
    Ul dim = meshes.size();
    Eigen::MatrixXcd hami_Lx(dim, dim);
    Eigen::MatrixXcd hami_Ly(dim, dim);
    Eigen::MatrixXcd hami_Lz(dim, dim);
#pragma omp parallel for shared(hami_Lx, hami_Ly, hami_Lz)
    for (Ul i = 0; i < dim; i++) {
        for (Ul i1 = 0; i1 < dim; i1++) {
            hami_Lx(i, i1) = lam_x[i1] * calculate_angular_momentum_(meshes[i], meshes[i1], gridPtr, centres[i1], shape, AngMomentum::Lx);
            hami_Ly(i, i1) = lam_y[i1] * calculate_angular_momentum_(meshes[i], meshes[i1], gridPtr, centres[i1], shape, AngMomentum::Ly);
            hami_Lz(i, i1) = lam_z[i1] * calculate_angular_momentum_(meshes[i], meshes[i1], gridPtr, centres[i1], shape, AngMomentum::Lz);
        }
    }
    for (Ul i = 0; i < meshes.size(); i++) {
        delete meshes[i];
        delete centres[i];
    }
    std::vector< Eigen::MatrixXcd > out({hami_Lx, hami_Ly, hami_Lz});
    return out;
}

std::vector< Eigen::Matrix< Cd, Eigen::Dynamic, Eigen::Dynamic > > construct_SO_hamiltonians(std::vector< std::tuple< int, int, int, int, int, float, float > > hami_raw, vec_f kxs, vec_f kys, vec_f kzs, std::vector< array_Cd > meshes, array_d grid, std::vector< array_d > centres, vec_d lam_x, vec_d lam_y, vec_d lam_z) {

    //    py::print("calculating SO hamiltonians -> started");
    Ul dim_0 = meshes.size();
    Eigen::MatrixXcf hami_0(dim_0, dim_0);

    std::vector< Eigen::MatrixXcd > ang_hamis = calculate_L_hamiltonians(meshes, grid, centres, lam_x, lam_y, lam_z);
    Eigen::MatrixXcd hami_Lx = ang_hamis[0];
    Eigen::MatrixXcd hami_Ly = ang_hamis[1];
    Eigen::MatrixXcd hami_Lz = ang_hamis[2];

    //    py::print(hami_Lx);
    //    py::print(hami_Ly);
    //    py::print(hami_Lz);

    Eigen::MatrixXcd tmp(2 * dim_0, 2 * dim_0);

    std::vector< Eigen::MatrixXcd > tot_hamis(kxs.size(), tmp);
    for (Ul i = 0; i < kxs.size(); i++) {
        hami_0 = hami_from_k(hami_raw, kxs[i], kys[i], kzs[i]);
        Eigen::MatrixXcd tot_hami(2 * dim_0, 2 * dim_0);
        Eigen::MatrixXcd tmp_hami = hami_0.cast< Cd >();

        tot_hami.block(0, 0, dim_0, dim_0) = tmp_hami + hami_Lz;
        tot_hami.block(0, dim_0, dim_0, dim_0) = hami_Lx + hami_Ly * Cd(0, 1);
        tot_hami.block(dim_0, 0, dim_0, dim_0) = hami_Lx - hami_Ly * Cd(0, 1);
        tot_hami.block(dim_0, dim_0, dim_0, dim_0) = tmp_hami - hami_Lz;
        tot_hamis[i] = tot_hami;
    }
    //    py::print("calculating SO hamiltonians -> done");
    return tot_hamis;
}

void init_ham_calc(py::module& m) {
    m.def("hami_from_k", &hami_from_k, "Get the bloch hamiltonian from k values", "hami_raw", "kx", "ky", "kz");
    m.def("calculate_L_hamiltonians", &calculate_L_hamiltonians, "Calculate the angular momentum matrix elements.", "meshes", "grid", "centres", "lam_x", "lam_y", "lam_z");
    m.def("construct_SO_hamiltonians", &construct_SO_hamiltonians, "hami_raw", "kxs", "kys", "kzs", "meshes", "grid", "lam_x", "lam_y", "lam_z");
}
