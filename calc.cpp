//
//  gen_calc.cpp
//  Python
//
//  Created by Louis Ponet on 14/02/2017.
//  Copyright © 2017 Louis Ponet. All rights reserved.
//

#include <Eigen/Dense>
#include "calc.hpp"
#include "ham_calc.hpp"
#include "dip_calc.hpp"
#include "sup_calc.hpp"

std::vector< std::tuple< Eigen::VectorXd, std::vector< Eigen::VectorXcd > > > generate_wannier_bands(std::vector< std::tuple< int, int, int, int, int, float, float > > hami_raw, vec_f kxs, vec_f kys, vec_f kzs) {
    Ul length = kxs.size();
    std::vector< Eigen::VectorXcd > tmp_eig(length);
    std::vector< Eigen::MatrixXcd > tmp_vec(length);
    //#pragma omp parallel for //watchout if this doesn't result in weirdness test
    for (Ul i = 0; i < length; i++) {
        Eigen::MatrixXcd hamiltonian = hami_from_k(hami_raw, kxs[i], kys[i], kzs[i]).cast< Cd >();

        Eigen::ComplexEigenSolver< Eigen::MatrixXcd > eigensolver(hamiltonian, Eigen::ComputeEigenvectors);
        eigensolver.compute(hamiltonian);
        tmp_eig[i] = eigensolver.eigenvalues();
        tmp_vec[i] = eigensolver.eigenvectors();
    }
    std::vector< std::tuple< Eigen::VectorXd, std::vector< Eigen::VectorXcd > > > out;
    Ul nr_bands = tmp_eig[0].size();

    for (Ul i = 0; i < nr_bands; i++) {
        Eigen::VectorXd tmp_band(length);
        std::vector< Eigen::VectorXcd > tmp_eigvs(length);
        for (Ul i1 = 0; i1 < length; i1++) {
            tmp_band(i1) = tmp_eig[i1](i).real();
            tmp_eigvs[i1] = tmp_vec[i1].col(i);
        }

        out.push_back(std::make_tuple(tmp_band, tmp_eigvs));
    }
    return out;
}

std::vector< std::tuple< Eigen::VectorXd, std::vector< Eigen::VectorXcd > > > generate_wannier_bands_SOC(std::vector< std::tuple< int, int, int, int, int, float, float > > hami_raw, vec_f kxs, vec_f kys, vec_f kzs, std::vector< array_Cd > meshes_, array_d grid, std::vector< array_d > centres, vec_d lam_x, vec_d lam_y, vec_d lam_z) {
    std::vector< Eigen::MatrixXcd > socHamis = construct_SO_hamiltonians(hami_raw, kxs, kys, kzs, meshes_, grid, centres, lam_x, lam_y, lam_z);

    Ul length = kxs.size();
    std::vector< Eigen::VectorXcd > tmp_eig(length);
    std::vector< Eigen::MatrixXcd > tmp_vec(length);
    for (Ul i = 0; i < length; i++) {

        Eigen::ComplexEigenSolver< Eigen::MatrixXcd > eigensolver(socHamis[i], Eigen::ComputeEigenvectors);
        eigensolver.compute(socHamis[i]);
        tmp_eig[i] = eigensolver.eigenvalues();
        tmp_vec[i] = eigensolver.eigenvectors();
    }
    std::vector< std::tuple< Eigen::VectorXd, std::vector< Eigen::VectorXcd > > > out;
    Ul nr_bands = tmp_eig[0].size();

    for (Ul i = 0; i < nr_bands; i++) {
        Eigen::VectorXd tmp_band(length);
        std::vector< Eigen::VectorXcd > tmp_eigvs(length);
        for (Ul i1 = 0; i1 < length; i1++) {
            tmp_band(i1) = tmp_eig[i1](i).real();
            tmp_eigvs[i1] = tmp_vec[i1].col(i);
        }

        out.push_back(std::make_tuple(tmp_band, tmp_eigvs));
    }
    return out;
}

std::vector< std::tuple< std::vector< vec_Cd >, Eigen::VectorXcd > > calculate_tb_dipoles(std::vector< std::tuple< int, int, int, int, int, float, float > > hami_raw, vec_f kxs, vec_f kys, vec_f kzs, std::vector< array_Cd > meshes_, array_d grid, std::vector< array_d > centres, vec_d lam_x, vec_d lam_y, vec_d lam_z) {
    std::vector< Eigen::Matrix< Cd, Eigen::Dynamic, Eigen::Dynamic > > socHamis = construct_SO_hamiltonians(hami_raw, kxs, kys, kzs, meshes_, grid, centres, lam_x, lam_y, lam_z);
    std::vector< array_Cd > socMeshes;
    for (int i = 0; i < 2; i++) {
        for (Ul i1 = 0; i1 < meshes_.size(); i1++) {
            socMeshes.push_back(meshes_[i1]);
        }
    }
    py::buffer_info gridInfo = grid.request();
    double* gridPtr = static_cast< double* >(gridInfo.ptr);

    std::vector< Ul > shape = gridInfo.shape;
    std::vector< Cd* > meshes = create_meshPtrs(socMeshes, shape);

    vec_Cd tmp({Cd(0, 0), Cd(0, 0), Cd(0, 0)});
    Ul dimEig = socMeshes.size();
    std::vector< vec_Cd > dipoles(dimEig * dimEig * Direction::END, tmp);
//    py::print("Calculating intercell dipoles -> started");
#pragma omp parallel for
    for (Ul orb1 = 0; orb1 < dimEig; orb1++) {
        Cd* mesh1 = meshes[orb1];
        for (Ul orb2 = 0; orb2 < dimEig; orb2++) {
            Cd* mesh2 = meshes[orb2];
            for (int i = Direction::Ex; i < Direction::END; i++) {
                if (orb1 < dimEig / 2 && orb2 >= dimEig / 2) {
                    dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + i] = vec_Cd({Cd(0, 0), Cd(0, 0), Cd(0, 0)});
                }
                else if (orb1 >= dimEig / 2 && orb2 < dimEig / 2) {
                    dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + i] = vec_Cd({Cd(0, 0), Cd(0, 0), Cd(0, 0)});
                }
                else {
                    vec_Cd tmp_ = calculate_intercell_dipole_(mesh1, mesh2, gridPtr, shape, (Direction)i);
                    dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + i] = tmp_;
                }
            }
        }
    }
    //    py::print("Calculating intercell dipoles -> done");
    //    py::print("Calculating intracell dipoles -> started");
    std::vector< vec_Cd > intra_dip(dimEig * dimEig, tmp);
#pragma omp parallel for
    for (Ul orb1 = 0; orb1 < dimEig; orb1++) {
        Cd* mesh1 = meshes[orb1];
        for (Ul orb2 = 0; orb2 < dimEig; orb2++) {
            Cd* mesh2 = meshes[orb2];
            if (orb1 < dimEig / 2 && orb2 >= dimEig / 2) {
                intra_dip[orb1 * dimEig + orb2] = vec_Cd({Cd(0, 0), Cd(0, 0), Cd(0, 0)});
            }
            else if (orb1 >= dimEig / 2 && orb2 < dimEig / 2) {
                intra_dip[orb1 * dimEig + orb2] = vec_Cd({Cd(0, 0), Cd(0, 0), Cd(0, 0)});
            }
            else {
                vec_Cd tmp_ = calculate_intracell_dipole_(mesh1, mesh2, gridPtr, shape);
                intra_dip[orb1 * dimEig + orb2] = tmp_;
            }
        }
    }

    //    py::print("Calculating intracell dipoles -> done");
    std::vector< vec_Cd > tmp2(dimEig * dimEig, tmp);
    Eigen::VectorXcd tmp3 = Eigen::VectorXcd::Random(dimEig);
    std::vector< std::tuple< std::vector< vec_Cd >, Eigen::VectorXcd > > out(kxs.size(), std::make_tuple(tmp2, tmp3));
//    py::print("Calculating total K-point dipoles -> started");
#pragma omp parallel for
    for (Ul i = 0; i < socHamis.size(); i++) {
        out[i] = calculate_eigen_dipoles(socHamis[i], dipoles, intra_dip, kxs[i], kys[i], kzs[i]);
    }
    //    py::print("Calculating total K-point dipoles -> done");
    for (Ul i = 0; i < meshes.size(); i++) {
        delete meshes[i];
    }
    return out;
}

std::vector< std::vector< vec_Cd > > calculate_tb_angmoms(std::vector< std::tuple< int, int, int, int, int, float, float > > hami_raw, vec_f kxs, vec_f kys, vec_f kzs, std::vector< array_Cd > meshes_, array_d grid, std::vector< array_d > centres_) {

    py::buffer_info gridInfo = grid.request();
    double* gridPtr = static_cast< double* >(gridInfo.ptr);

    std::vector< Ul > shape = gridInfo.shape;
    std::vector< Cd* > meshes = create_meshPtrs(meshes_, shape);

    std::vector< double* > centres = create_centrePtrs(centres_);

    vec_Cd tmp({Cd(0, 0), Cd(0, 0), Cd(0, 0)});
    std::vector< vec_Cd > angMoms(meshes_.size() * meshes_.size(), tmp);
    std::vector< vec_Cd > outTmp(meshes_.size());
    std::vector< std::vector< vec_Cd > > outMoms(kxs.size(), outTmp);

#pragma omp parallel for
    for (Ul i = 0; i < meshes_.size(); i++) {
        for (Ul i1 = 0; i1 < meshes_.size(); i1++) {
            vec_Cd tmpAng(3);
            tmpAng[0] = calculate_angular_momentum_(meshes[i], meshes[i1], gridPtr, centres[i1], shape, AngMomentum::Lx);
            tmpAng[1] = calculate_angular_momentum_(meshes[i], meshes[i1], gridPtr, centres[i1], shape, AngMomentum::Ly);
            tmpAng[2] = calculate_angular_momentum_(meshes[i], meshes[i1], gridPtr, centres[i1], shape, AngMomentum::Lz);
            angMoms[i * meshes_.size() + i1] = tmpAng;
        }
    }
#pragma omp parallel for
    for (Ul i = 0; i < kxs.size(); i++) {
        Eigen::MatrixXcd hami = hami_from_k(hami_raw, kxs[i], kys[i], kzs[i]).cast< Cd >();
        Eigen::ComplexEigenSolver< Eigen::MatrixXcd > eigensolver(hami);
        Eigen::MatrixXcd eigenvectors = eigensolver.eigenvectors();

        Ul dimEig = eigenvectors.cols();
        vec_Cd tmpAng({Cd(0, 0), Cd(0, 0), Cd(0, 0)});
        for (Ul band = 0; band < dimEig; band++) {
            Eigen::VectorXcd eigVec = eigenvectors.col(band);
            for (Ul orb1 = 0; orb1 < dimEig; orb1++) {
                for (Ul orb2 = 0; orb2 < dimEig; orb2++) {
                    tmpAng[0] += std::conj(eigVec[orb1]) * eigVec[orb2] * angMoms[orb1 * dimEig + orb2][0];
                    tmpAng[1] += std::conj(eigVec[orb1]) * eigVec[orb2] * angMoms[orb1 * dimEig + orb2][1];
                    tmpAng[2] += std::conj(eigVec[orb1]) * eigVec[orb2] * angMoms[orb1 * dimEig + orb2][2];
                }
            }
            outMoms[i][band] = tmpAng;
        }
    }
    for (Ul i = 0; i < centres.size(); i++) {
        delete centres[i];
    }
    for (Ul i = 0; i < meshes.size(); i++) {
        delete meshes[i];
    }
    return outMoms;
}

std::vector< array_Cd > construct_eigenmesh_SOC(std::vector< std::tuple< int, int, int, int, int, float, float > > hami_raw, float kx, float ky, float kz, std::vector< array_Cd > meshes_, array_d grid, vec_d ls, std::vector< array_d > centres_) {
    Eigen::MatrixXcd hami = construct_SO_hamiltonians(hami_raw, vec_f({kx}), vec_f({ky}), vec_f({kz}), meshes_, grid, centres_, ls, ls, ls)[0];
    Eigen::ComplexEigenSolver< Eigen::MatrixXcd > eigensolver(hami);
    Eigen::MatrixXcd eigenvectors = eigensolver.eigenvectors();

    Ul dimEig = eigenvectors.cols();
    py::buffer_info gridInfo = grid.request();

    Ul dim0 = gridInfo.shape[0];
    Ul dim1 = gridInfo.shape[1];
    Ul dim2 = gridInfo.shape[2];

    std::vector< array_Cd > socMeshes;
    for (int i = 0; i < 2; i++) {
        for (Ul i1 = 0; i1 < meshes_.size(); i1++) {
            socMeshes.push_back(meshes_[i1]);
        }
    }

    std::vector< Cd* > meshes = create_meshPtrs(socMeshes, gridInfo.shape);

    std::vector< Cd* > out_tmp(dimEig);
    for (Ul i = 0; i < dimEig; i++) {
        out_tmp[i] = new Cd[dim0 * dim1 * dim2];
    }
    std::vector< array_Cd > out;
#pragma omp parallel for
    for (Ul i = 0; i < dimEig; i++) {
        Cd* tmpMeshPtr = out_tmp[i];
        Eigen::VectorXcd eigenvector = eigenvectors.col(i);
        for (Ul orb1 = 0; orb1 < dimEig; orb1++) {
            Cd c1 = eigenvector[orb1];
            Cd* meshPtr1 = meshes[orb1];
            for (Ul i1 = 0; i1 < dim0; i1++) {
                for (Ul i2 = 0; i2 < dim1; i2++) {
                    for (Ul i3 = 0; i3 < dim2; i3++) {
                        if (orb1 == 0) {
                            tmpMeshPtr[i1 * dim1 * dim2 + i2 * dim2 + i3] = c1 * meshPtr1[i1 * dim1 * dim2 + i2 * dim2 + i3];
                        }
                        else {
                            tmpMeshPtr[i1 * dim1 * dim2 + i2 * dim2 + i3] += c1 * meshPtr1[i1 * dim1 * dim2 + i2 * dim2 + i3];
                        }
                    }
                }
            }
        }
    }
    for (Ul i = 0; i < out_tmp.size(); i++) {
        out.push_back(array_Cd(std::vector< size_t >({dim0, dim1, dim2}), out_tmp[i]));
        delete out_tmp[i];
    }
    for (Ul i = 0; i < meshes.size(); i++) {
        delete meshes[i];
    }
    return out;
}

void init_calc(py::module& m) {
    m.def("generate_wannier_bands", &generate_wannier_bands, "Calculate the wannier bands for a sequence of k-values", "hami_raw", "kxs", "kys", "kzs");
    m.def("generate_wannier_bands_SOC", &generate_wannier_bands_SOC, "Calculate the wannier bands with SOC", "hami_raw", "kxs", "kys", "kzs", "meshes_", "grid", "centres", "lam_x", "lam_y", "lam_z");
    m.def("calculate_tb_dipoles", &calculate_tb_dipoles, "hami_raw", "kxs", "kys", "kzs", "meshes", "grid", "centres", "lam_x", "lam_y", "lam_z");
    m.def("calculate_tb_angmoms", &calculate_tb_angmoms, "hami_raw", "kxs", "kys", "kzs", "meshes", "grid", "centres");
    m.def("construct_eigenmesh_SOC", &construct_eigenmesh_SOC, "hami_raw", "kx", "ky", "kz", "meshes", "grid", "ls", "centres");
}
