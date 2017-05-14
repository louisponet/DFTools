//
//  dipole_calculations.cpp
//  Python
//
//  Created by Louis Ponet on 14/02/2017.
//  Copyright © 2017 Louis Ponet. All rights reserved.
//
#include <Eigen/Dense>
#include <omp.h>
#include "dip_calc.hpp"
#include "sup_calc.hpp"

vec_Cd calculate_intracell_dipole(array_Cd mesh1, array_Cd mesh2, array_d grid) {

    py::buffer_info gridInfo = grid.request();
    py::buffer_info mesh1Info = mesh1.request();
    py::buffer_info mesh2Info = mesh2.request();

    Ul dim0 = gridInfo.shape[0];
    Ul dim1 = gridInfo.shape[1];
    Ul dim2 = gridInfo.shape[2];
    Ul dim3 = gridInfo.shape[3];

    auto gridPtr = static_cast< double* >(gridInfo.ptr);
    auto meshPtr1 = static_cast< Cd* >(mesh1Info.ptr);
    auto meshPtr2 = static_cast< Cd* >(mesh2Info.ptr);

    Cd dx(0, 0);
    Cd dy(0, 0);
    Cd dz(0, 0);
    double n1 = 0;
    double n2 = 0;

    for (Ul i = 0; i < dim0; i++) {
        for (Ul i1 = 0; i1 < dim1; i1++) {
            for (Ul i2 = 0; i2 < dim2; i2++) {
                dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                n2 += std::norm(meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2]);
            }
        }
    }
    double n = std::sqrt(n1 * n2);
    dx /= n;
    dy /= n;
    dz /= n;
    return vec_Cd({dx, dy, dz});
}

vec_Cd calculate_intracell_dipole_(Cd* meshPtr1, Cd* meshPtr2, double* gridPtr, std::vector< Ul > shape) {


    Ul dim0 = shape[0];
    Ul dim1 = shape[1];
    Ul dim2 = shape[2];
    Ul dim3 = shape[3];

    Cd dx(0, 0);
    Cd dy(0, 0);
    Cd dz(0, 0);
    double n1 = 0;
    double n2 = 0;

    for (Ul i = 0; i < dim0; i++) {
        for (Ul i1 = 0; i1 < dim1; i1++) {
            for (Ul i2 = 0; i2 < dim2; i2++) {
                dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                n2 += std::norm(meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2]);
            }
        }
    }
    double n = std::sqrt(n1 * n2);
    dx /= n;
    dy /= n;
    dz /= n;
    return vec_Cd({dx, dy, dz});
}


vec_Cd calculate_intercell_dipole(array_Cd mesh1, array_Cd mesh2, array_d grid, Direction d) {
    py::buffer_info gridInfo = grid.request();
    py::buffer_info mesh1Info = mesh1.request();
    py::buffer_info mesh2Info = mesh2.request();

    Ul dim0 = gridInfo.shape[0];
    Ul dim1 = gridInfo.shape[1];
    Ul dim2 = gridInfo.shape[2];
    Ul dim3 = gridInfo.shape[3];

    auto gridPtr = static_cast< double* >(gridInfo.ptr);
    auto meshPtr1 = static_cast< Cd* >(mesh1Info.ptr);
    auto meshPtr2 = static_cast< Cd* >(mesh2Info.ptr);

    Cd dx(0, 0);
    Cd dy(0, 0);
    Cd dz(0, 0);
    double n1 = 0;
    double n2 = 0;
    double n = 0;
    switch (d) {
        case Ex:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = dim2 / 3 + 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3]);
                    }
                }
            }
            break;
        case Ey:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = dim1 / 3 + 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < dim0; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2]);
                    }
                }
            }
            break;
        case Ez:
            for (Ul i = dim0 / 3 + 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2]));
                    }
                }
            }
            break;
        case Ex_:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < 2 * dim2 / 3 + 2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2]);
                    }
                }
            }
            break;
        case Ey_:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < 2 * dim1 / 3 + 2; i1++) {
                    for (Ul i2 = 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2]));
                    }
                }
            }
            break;
        case Ez_:
            for (Ul i = 1; i < 2 * dim0 / 3 + 2; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2]));
                    }
                }
            }
            break;
        case ExEx:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = 2 * dim2 / 3; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - 2 * dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - 2 * dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - 2 * dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - 2 * dim2 / 3]));
                    }
                }
            }
            break;
        case Ex_Ex_:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < dim2 / 3 + 1; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + 2 * dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + 2 * dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + 2 * dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + 2 * dim2 / 3]));
                    }
                }
            }
            break;
        case ExEy:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = dim1 / 3 + 1; i1 < dim1; i1++) {
                    for (Ul i2 = dim2 / 3 + 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3]));
                    }
                }
            }
            break;
        case Ex_Ey:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = dim1 / 3 + 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < 2 * dim2 / 3 + 2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 + dim2 / 3 - 2]));
                    }
                }
            }
            break;
        case ExEy_:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < 2 * dim1 / 3 + 2; i1++) {
                    for (Ul i2 = dim2 / 3 + 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 - dim2 / 3]));
                    }
                }
            }
            break;
        case Ex_Ey_:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < 2 * dim1 / 3 + 2; i1++) {
                    for (Ul i2 = 1; i2 < 2 * dim2 / 3 + 2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2]));
                    }
                }
            }
            break;
        case ExEz:
            for (Ul i = dim0 / 3 + 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = dim2 / 3 + 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3]));
                    }
                }
            }
            break;
        case Ex_Ez:
            for (Ul i = dim0 / 3 + 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < 2 * dim2 / 3 + 2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2]));
                    }
                }
            }
            break;
        case ExEz_:
            for (Ul i = 1; i < 2 * dim0 / 3 + 2; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = dim2 / 3 + 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3]));
                    }
                }
            }
            break;
        case Ex_Ez_:
            for (Ul i = 1; i < 2 * dim0 / 3 + 2; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < 2 * dim2 / 3 + 2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2]));
                    }
                }
            }
            break;
        case EyEz:
            for (Ul i = dim0 / 3 + 1; i < dim0; i++) {
                for (Ul i1 = dim1 / 3 + 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < dim0; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2]));
                    }
                }
            }
            break;
        case Ey_Ez:
            for (Ul i = dim0 / 3 + 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < 2 * dim1 / 3 + 2; i1++) {
                    for (Ul i2 = 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2]));
                    }
                }
            }
            break;
        case EyEz_:
            for (Ul i = 1; i < 2 * dim0 / 3 + 2; i++) {
                for (Ul i1 = dim1 / 3 + 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < dim0; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2]));
                    }
                }
            }
            break;
        case Ey_Ez_:
            for (Ul i = 1; i < 2 * dim0 / 3 + 2; i++) {
                for (Ul i1 = 1; i1 < 2 * dim1 / 3 + 2; i1++) {
                    for (Ul i2 = 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2]));
                    }
                }
            }
            break;
        case ExEyEz:
            for (Ul i = dim0 / 3 + 1; i < dim0; i++) {
                for (Ul i1 = dim1 / 3 + 1; i1 < dim1; i1++) {
                    for (Ul i2 = dim2 / 3 + 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3]));
                    }
                }
            }
            break;
        case Ex_Ey_Ez_:
            for (Ul i = 1; i < 2 * dim0 / 3 + 2; i++) {
                for (Ul i1 = 1; i1 < 2 * dim1 / 3 + 2; i1++) {
                    for (Ul i2 = 1; i2 < 2 * dim2 / 3 + 2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n += std::sqrt(std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * std::norm(meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2]));
                    }
                }
            }
            break;

        default:
            break;
    }

    n = std::sqrt(n1 * n2);
    dx /= n;
    dy /= n;
    dz /= n;

    return vec_Cd({dx, dy, dz});
}

vec_Cd calculate_intercell_dipole_(Cd* meshPtr1, Cd* meshPtr2, double* gridPtr, std::vector< Ul > shape, Direction d) {
    //we take the integral to be only over the unit cell -> only terms within reach of r+2 count
    // indices for center wannier mesh : 26->53, indices off centered wannier 53->80, 0->26 in each direction (a,b,c)
    //Issue 1: cell is not given by these mesh coordinates! origin is at 28,28,28 cell boundary is at 55,55,55
    //      2: Something is going on with z?

    Ul dim0 = shape[0];
    Ul dim1 = shape[1];
    Ul dim2 = shape[2];
    Ul dim3 = shape[3];


    Cd dx(0, 0);
    Cd dy(0, 0);
    Cd dz(0, 0);
    double n1 = 0;
    double n2 = 0;
    switch (d) {
        case Ex:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = dim2 / 3 + 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3]);
                    }
                }
            }
            break;
        case Ey:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = dim1 / 3 + 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < dim0; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2]);
                    }
                }
            }
            break;
        case Ez:
            for (Ul i = dim0 / 3 + 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2]);
                    }
                }
            }
            break;
        case Ex_:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < 2 * dim2 / 3 + 2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2]);
                    }
                }
            }
            break;
        case Ey_:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < 2 * dim1 / 3 + 2; i1++) {
                    for (Ul i2 = 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2]);
                    }
                }
            }
            break;
        case Ez_:
            for (Ul i = 1; i < 2 * dim0 / 3 + 2; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2]);
                    }
                }
            }
            break;
        case ExEx:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = 2 * dim2 / 3; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - 2 * dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - 2 * dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - 2 * dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 - 2 * dim2 / 3]);
                    }
                }
            }
            break;
        case Ex_Ex_:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < dim2 / 3 + 1; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + 2 * dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + 2 * dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + 2 * dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[i * dim1 * dim2 + i1 * dim2 + i2 + 2 * dim2 / 3]);
                    }
                }
            }
            break;
        case ExEy:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = dim1 / 3 + 1; i1 < dim1; i1++) {
                    for (Ul i2 = dim2 / 3 + 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3]);
                    }
                }
            }
            break;
        case Ex_Ey:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = dim1 / 3 + 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < 2 * dim2 / 3 + 2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[i * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 + dim2 / 3 - 2]);
                    }
                }
            }
            break;
        case ExEy_:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < 2 * dim1 / 3 + 2; i1++) {
                    for (Ul i2 = dim2 / 3 + 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 - dim2 / 3]);
                    }
                }
            }
            break;
        case Ex_Ey_:
            for (Ul i = 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < 2 * dim1 / 3 + 2; i1++) {
                    for (Ul i2 = 1; i2 < 2 * dim2 / 3 + 2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[i * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2]);
                    }
                }
            }
            break;
        case ExEz:
            for (Ul i = dim0 / 3 + 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = dim2 / 3 + 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3]);
                    }
                }
            }
            break;
        case Ex_Ez:
            for (Ul i = dim0 / 3 + 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < 2 * dim2 / 3 + 2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[(i - dim0 / 3) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2]);
                    }
                }
            }
            break;
        case ExEz_:
            for (Ul i = 1; i < 2 * dim0 / 3 + 2; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = dim2 / 3 + 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 - dim2 / 3]);
                    }
                }
            }
            break;
        case Ex_Ez_:
            for (Ul i = 1; i < 2 * dim0 / 3 + 2; i++) {
                for (Ul i1 = 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < 2 * dim2 / 3 + 2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + i1 * dim2 + i2 + dim2 / 3 - 2]);
                    }
                }
            }
            break;
        case EyEz:
            for (Ul i = dim0 / 3 + 1; i < dim0; i++) {
                for (Ul i1 = dim1 / 3 + 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < dim0; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2]);
                    }
                }
            }
            break;
        case Ey_Ez:
            for (Ul i = dim0 / 3 + 1; i < dim0; i++) {
                for (Ul i1 = 1; i1 < 2 * dim1 / 3 + 2; i1++) {
                    for (Ul i2 = 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2]);
                    }
                }
            }
            break;
        case EyEz_:
            for (Ul i = 1; i < 2 * dim0 / 3 + 2; i++) {
                for (Ul i1 = dim1 / 3 + 1; i1 < dim1; i1++) {
                    for (Ul i2 = 1; i2 < dim0; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2]);
                    }
                }
            }
            break;
        case Ey_Ez_:
            for (Ul i = 1; i < 2 * dim0 / 3 + 2; i++) {
                for (Ul i1 = 1; i1 < 2 * dim1 / 3 + 2; i1++) {
                    for (Ul i2 = 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2]);
                    }
                }
            }
            break;
        case ExEyEz:
            for (Ul i = dim0 / 3 + 1; i < dim0; i++) {
                for (Ul i1 = dim1 / 3 + 1; i1 < dim1; i1++) {
                    for (Ul i2 = dim2 / 3 + 1; i2 < dim2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[(i - dim0 / 3) * dim1 * dim2 + (i1 - dim1 / 3) * dim2 + i2 - dim2 / 3]);
                    }
                }
            }
            break;
        case Ex_Ey_Ez_:
            for (Ul i = 1; i < 2 * dim0 / 3 + 2; i++) {
                for (Ul i1 = 1; i1 < 2 * dim1 / 3 + 2; i1++) {
                    for (Ul i2 = 1; i2 < 2 * dim2 / 3 + 2; i2++) {
                        dx += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3];
                        dy += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 1];
                        dz += std::conj(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]) * meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2] * gridPtr[i * dim1 * dim2 * dim3 + i1 * dim2 * dim3 + i2 * dim3 + 2];
                        n1 += std::norm(meshPtr1[i * dim1 * dim2 + i1 * dim2 + i2]);
                        n2 += std::norm(meshPtr2[(i + dim0 / 3 - 2) * dim1 * dim2 + (i1 + dim1 / 3 - 2) * dim2 + i2 + dim2 / 3 - 2]);
                    }
                }
            }
            break;

        default:
            break;
    }
    double n = std::sqrt(n1 * n2);
    dx /= n;
    dy /= n;
    dz /= n;
    return vec_Cd({dx, dy, dz});
}


std::tuple< std::vector< vec_Cd >, Eigen::VectorXcd > calculate_eigen_dipoles(Eigen::Matrix< Cd, Eigen::Dynamic, Eigen::Dynamic > hami, std::vector< vec_Cd > dipoles, std::vector< vec_Cd > intra_dip, float kx, float ky, float kz) {
    //P=\frac{1}{V_uc}\int_{\rm{uc}}d^3r\quad r|\psi(r)|^2=\frac{1}{V_uc}\int_{\rm{uc}}d^3r\quad \sum_{nR} c_n^* w_n^*(r+R)e^{-ikR} r \sum_{mR'} c_m w_m(r+R')e^{ikR'}


    //not sure about this but it seems to be needed since we also did this to construct the hamiltonian from k
    kx = 2 * M_PI * kx;
    ky = 2 * M_PI * ky;
    kz = 2 * M_PI * kz;

    //decompose the hamiltonian
    Eigen::ComplexEigenSolver< Eigen::MatrixXcd > eigensolver(hami);
    Eigen::MatrixXcd eigenvectors = eigensolver.eigenvectors();
    Ul dimEig = eigenvectors.cols();
    //create outEig for also returning the eigenvalues of the hamiltonian
    Eigen::VectorXcd outEig = eigensolver.eigenvalues();
    //initialize list of dipoles so they stay in order during parallel processing
    vec_Cd tmp({Cd(0, 0), Cd(0, 0), Cd(0, 0)});
    std::vector< vec_Cd > outDipoles(dimEig, tmp);

    //precacalculate the bloch phases
    std::vector< Cd > factors;
    factors.push_back(std::exp(Cd(0, (double)kx)));
    factors.push_back(std::exp(Cd(0, (double)-kx)));
    factors.push_back(std::exp(Cd(0, (double)ky)));
    factors.push_back(std::exp(Cd(0, (double)-ky)));
    factors.push_back(std::exp(Cd(0, (double)kz)));
    factors.push_back(std::exp(Cd(0, (double)-kz)));
    factors.push_back(std::exp(Cd(0, (double)2 * kx)));
    factors.push_back(std::exp(Cd(0, (double)-2 * kx)));
    factors.push_back(std::exp(Cd(0, (double)kx + ky)));
    factors.push_back(std::exp(Cd(0, (double)-kx - ky)));
    factors.push_back(std::exp(Cd(0, (double)-kx + ky)));
    factors.push_back(std::exp(Cd(0, (double)kx - ky)));
    factors.push_back(std::exp(Cd(0, (double)kx + kz)));
    factors.push_back(std::exp(Cd(0, (double)-kx - kz)));
    factors.push_back(std::exp(Cd(0, (double)-kx + kz)));
    factors.push_back(std::exp(Cd(0, (double)kx - kz)));
    factors.push_back(std::exp(Cd(0, (double)ky + kz)));
    factors.push_back(std::exp(Cd(0, (double)-ky - kz)));
    factors.push_back(std::exp(Cd(0, (double)-ky + kz)));
    factors.push_back(std::exp(Cd(0, (double)ky - kz)));
    factors.push_back(std::exp(Cd(0, (double)kx + ky + kz)));
    factors.push_back(std::exp(Cd(0, (double)-kx - ky - kz)));

    //#pragma omp parallel for
    for (Ul i = 0; i < dimEig; i++) {
        Cd dip_x(0, 0);
        Cd dip_y(0, 0);
        Cd dip_z(0, 0);
        Eigen::VectorXcd eigenvector = eigenvectors.col(i);
        for (Ul orb1 = 0; orb1 < dimEig; orb1++) {
            Cd c_orb1 = eigenvector(orb1);
            for (Ul orb2 = 0; orb2 < dimEig; orb2++) {
                Cd c_orb2 = eigenvector(orb2);
                //                 exp(i*k_1*1)*[(conj(w_n(r-R_1))*w_m(r))*r+conj(w_n(r))*w_m(r-R_1)*r]+exp(-i*k_1*1)*[conj(w_n(r+R_1)*w_m(r))*r+conj(w_n(r))*w_m(r-R_1)]
                dip_x += std::conj(c_orb1) * c_orb2 * (factors[0] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 1][0]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END][0]) + factors[1] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 1][0] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END][0])));
                dip_x += std::conj(c_orb1) * c_orb2 * (factors[2] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 3][0]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 2][0]) + factors[3] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 3][0] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 2][0])));
                dip_x += std::conj(c_orb1) * c_orb2 * (factors[4] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 5][0]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 4][0]) + factors[5] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 5][0] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 4][0])));
                //
                //                dip_x += std::conj(c_orb1)*c_orb2*(factors[6]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+7][0])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+6][0])+factors[7]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+7][0]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+6][0])));
                //                dip_x += std::conj(c_orb1)*c_orb2*(factors[8]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+9][0])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+8][0])+factors[9]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+9][0]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+8][0])));
                //                dip_x += std::conj(c_orb1)*c_orb2*(factors[10]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+11][0])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+10][0])+factors[11]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+11][0]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+10][0])));
                //                dip_x += std::conj(c_orb1)*c_orb2*(factors[12]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+13][0])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+12][0])+factors[13]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+13][0]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+12][0])));
                //                dip_x += std::conj(c_orb1)*c_orb2*(factors[14]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+15][0])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+14][0])+factors[15]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+15][0]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+14][0])));
                //                dip_x += std::conj(c_orb1)*c_orb2*(factors[16]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+17][0])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+16][0])+factors[17]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+17][0]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+16][0])));
                //                dip_x += std::conj(c_orb1)*c_orb2*(factors[18]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+19][0])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+18][0])+factors[19]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+19][0]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+18][0])));
                //                dip_x += std::conj(c_orb1)*c_orb2*(factors[20]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+21][0])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+20][0])+factors[21]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+21][0]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+20][0])));
                //
                dip_y += std::conj(c_orb1) * c_orb2 * (factors[0] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 1][1]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END][1]) + factors[1] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 1][1] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END][1])));
                dip_y += std::conj(c_orb1) * c_orb2 * (factors[2] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 3][1]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 2][1]) + factors[3] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 3][1] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 2][1])));
                dip_y += std::conj(c_orb1) * c_orb2 * (factors[4] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 5][1]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 4][1]) + factors[5] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 5][1] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 4][1])));
                //
                //                dip_y += std::conj(c_orb1)*c_orb2*(factors[6]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+7][1])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+6][1])+factors[7]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+7][1]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+6][1])));
                //                dip_y += std::conj(c_orb1)*c_orb2*(factors[8]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+9][1])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+8][1])+factors[9]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+9][1]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+8][1])));
                //                dip_y += std::conj(c_orb1)*c_orb2*(factors[10]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+11][1])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+10][1])+factors[11]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+11][1]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+10][1])));
                //                dip_y += std::conj(c_orb1)*c_orb2*(factors[12]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+13][1])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+12][1])+factors[13]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+13][1]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+12][1])));
                //                dip_y += std::conj(c_orb1)*c_orb2*(factors[14]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+15][1])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+14][1])+factors[15]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+15][1]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+14][1])));
                //                dip_y += std::conj(c_orb1)*c_orb2*(factors[16]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+17][1])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+16][1])+factors[17]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+17][1]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+16][1])));
                //                dip_y += std::conj(c_orb1)*c_orb2*(factors[18]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+19][1])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+18][1])+factors[19]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+19][1]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+18][1])));
                //                dip_y += std::conj(c_orb1)*c_orb2*(factors[20]*(std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+21][1])+dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+20][1])+factors[21]*(dipoles[orb1*dimEig*Direction::END+orb2*Direction::END+21][1]+std::conj(dipoles[orb2*dimEig*Direction::END+orb1*Direction::END+20][1])));
                // //
                dip_z += std::conj(c_orb1) * c_orb2 * (factors[0] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 1][2]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END][2]) + factors[1] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 1][2] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END][2])));
                dip_z += std::conj(c_orb1) * c_orb2 * (factors[2] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 3][2]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 2][2]) + factors[3] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 3][2] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 2][2])));
                dip_z += std::conj(c_orb1) * c_orb2 * (factors[4] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 5][2]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 4][2]) + factors[5] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 5][2] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 4][2])));
                ////
                // dip_z += std::conj(c_orb1) * c_orb2 * (factors[6] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 7][2]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 6][2]) + factors[7] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 7][2] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 6][2])));
                // dip_z += std::conj(c_orb1) * c_orb2 * (factors[8] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 9][2]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 8][2]) + factors[9] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 9][2] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 8][2])));
                // dip_z += std::conj(c_orb1) * c_orb2 * (factors[10] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 11][2]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 10][2]) + factors[11] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 11][2] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 10][2])));
                // dip_z += std::conj(c_orb1) * c_orb2 * (factors[12] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 13][2]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 12][2]) + factors[13] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 13][2] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 12][2])));
                // dip_z += std::conj(c_orb1) * c_orb2 * (factors[14] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 15][2]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 14][2]) + factors[15] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 15][2] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 14][2])));
                // dip_z += std::conj(c_orb1) * c_orb2 * (factors[16] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 17][2]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 16][2]) + factors[17] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 17][2] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 16][2])));
                // dip_z += std::conj(c_orb1) * c_orb2 * (factors[18] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 19][2]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 18][2]) + factors[19] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 19][2] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 18][2])));
                // dip_z += std::conj(c_orb1) * c_orb2 * (factors[20] * (std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 21][2]) + dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 20][2]) + factors[21] * (dipoles[orb1 * dimEig * Direction::END + orb2 * Direction::END + 21][2] + std::conj(dipoles[orb2 * dimEig * Direction::END + orb1 * Direction::END + 20][2])));

                dip_x += std::conj(c_orb1) * c_orb2 * intra_dip[orb1 * dimEig + orb2][0];
                dip_y += std::conj(c_orb1) * c_orb2 * intra_dip[orb1 * dimEig + orb2][1];
                dip_z += std::conj(c_orb1) * c_orb2 * intra_dip[orb1 * dimEig + orb2][2];
            }
        }
        outDipoles[i] = vec_Cd({dip_x, dip_y, dip_z});
    }
    return std::make_tuple(outDipoles, outEig);
}

void init_dip_calc(py::module& m) {
    m.def("calculate_intracell_dipole", &calculate_intracell_dipole, "Calculate dipole between two meshes within the same unit cell", "mesh1", "mesh2", "grid");
    m.def("calculate_intercell_dipole", &calculate_intercell_dipole, "Calculate dipole between two meshes 1 cell apart.", "mesh1", "mesh2", "grid", "direction");
}
