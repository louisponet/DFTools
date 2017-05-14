//
//  sup_calc.cpp
//  Python
//
//  Created by Louis Ponet on 14/02/2017.
//  Copyright © 2017 Louis Ponet. All rights reserved.
//
#include <Eigen/Geometry>
#include "sup_calc.hpp"

namespace py = pybind11;

double h_bar = 1.0 ;

std::vector<Cd *> create_meshPtrs(std::vector<array_Cd> meshes_,std::vector<Ul> shape){
    Ul dim0 = shape[0];
    Ul dim1 = shape[1];
    Ul dim2 = shape[2];
    std::vector<Cd*> meshes;
    for (Ul i =0; i<meshes_.size(); i++) {
        Cd* tmp = new Cd[dim0*dim1*dim2];
        Cd* meshPtr = static_cast<Cd *>(meshes_[i].request().ptr);
        memcpy(tmp,meshPtr,sizeof(Cd)*dim1*dim2*dim0);
        meshes.push_back(tmp);
    }
    return meshes;
}

std::vector<double *> create_centrePtrs(std::vector<array_d> centres_){
    std::vector<double*> centres;
    for (Ul i = 0; i<centres_.size(); i++) {
        double* tmp = new double[3];
        double* centrePtr = static_cast<double *>(centres_[i].request().ptr);
        memcpy(tmp,centrePtr,sizeof(double)*3);
        centres.push_back(tmp);
    }
    return centres;
}

array_Cd construct_px(array_d grid){
    py::buffer_info gridInfo = grid.request();
    Ul dim0 = gridInfo.shape[0];
    Ul dim1 = gridInfo.shape[1];
    Ul dim2 = gridInfo.shape[2];
    Ul dim3 = gridInfo.shape[3];
    auto gridPtr = static_cast<double *>(gridInfo.ptr);
    array_Cd out(std::vector<size_t>({dim0,dim1,dim2}));
    auto outPtr = static_cast<Cd *>(out.request().ptr);
    
    for (Ul i = 0; i<dim0; i++) {
        for (Ul i1 = 0; i1<dim1; i1++) {
            for (Ul i2 = 0; i2<dim2 ;i2++) {
                double x = gridPtr[i*dim1*dim2*dim3+i1*dim2*dim3+i2*dim3];
                double y = gridPtr[i*dim1*dim2*dim3+i1*dim2*dim3+i2*dim3+1];
                double z = gridPtr[i*dim1*dim2*dim3+i1*dim2*dim3+i2*dim3+2];
                outPtr[i*dim1*dim2+i1*dim2+i2] = Cd(x*exp(-sqrt(x*x+y*y+z*z)),0);
            }
        }
    }
    return out;
    
}

array_Cd construct_py(array_d grid){
    py::buffer_info gridInfo = grid.request();
    Ul dim0 = gridInfo.shape[0];
    Ul dim1 = gridInfo.shape[1];
    Ul dim2 = gridInfo.shape[2];
    Ul dim3 = gridInfo.shape[3];
    auto gridPtr = static_cast<double *>(gridInfo.ptr);
    array_Cd out(std::vector<size_t>({dim0,dim1,dim2}));
    auto outPtr = static_cast<Cd *>(out.request().ptr);
    
    for (Ul i = 0; i<dim0; i++) {
        for (Ul i1 = 0; i1<dim1; i1++) {
            for (Ul i2 = 0; i2<dim2 ;i2++) {
                double x = gridPtr[i*dim1*dim2*dim3+i1*dim2*dim3+i2*dim3];
                double y = gridPtr[i*dim1*dim2*dim3+i1*dim2*dim3+i2*dim3+1];
                double z = gridPtr[i*dim1*dim2*dim3+i1*dim2*dim3+i2*dim3+2];
                outPtr[i*dim1*dim2+i1*dim2+i2] = Cd(y*exp(-sqrt(x*x+y*y+z*z)),0);
            }
        }
    }
    return out;
    
}

array_Cd construct_pz(array_d grid){
    py::buffer_info gridInfo = grid.request();
    Ul dim0 = gridInfo.shape[0];
    Ul dim1 = gridInfo.shape[1];
    Ul dim2 = gridInfo.shape[2];
    Ul dim3 = gridInfo.shape[3];
    auto gridPtr = static_cast<double *>(gridInfo.ptr);
    array_Cd out(std::vector<size_t>({dim0,dim1,dim2}));
    auto outPtr = static_cast<Cd *>(out.request().ptr);
    
    for (Ul i = 0; i<dim0; i++) {
        for (Ul i1 = 0; i1<dim1; i1++) {
            for (Ul i2 = 0; i2<dim2 ;i2++) {
                double x = gridPtr[i*dim1*dim2*dim3+i1*dim2*dim3+i2*dim3];
                double y = gridPtr[i*dim1*dim2*dim3+i1*dim2*dim3+i2*dim3+1];
                double z = gridPtr[i*dim1*dim2*dim3+i1*dim2*dim3+i2*dim3+2];
                outPtr[i*dim1*dim2+i1*dim2+i2] = Cd(z*exp(-sqrt(x*x+y*y+z*z)),0);
            }
        }
    }
    return out;
    
}


array_d generate_grid(array_d origin,array_d a_span_array,array_d b_span_array,array_d c_span_array){
    py::buffer_info span_inf_0 = a_span_array.request();
    py::buffer_info span_inf_1 = b_span_array.request();
    py::buffer_info span_inf_2 = c_span_array.request();
    py::buffer_info origin_inf = origin.request();
    
    auto span_ptr_0 = static_cast<double *>(span_inf_0.ptr);
    auto span_ptr_1 = static_cast<double *>(span_inf_1.ptr);
    auto span_ptr_2 = static_cast<double *>(span_inf_2.ptr);
    auto origin_ptr = static_cast<double *>(origin_inf.ptr);
    
    array_d out(std::vector<size_t>({span_inf_0.shape[0],span_inf_1.shape[0],span_inf_2.shape[0],3}));
    auto out_ptr = static_cast<double *>(out.request().ptr);
    
    for (int i = 0 ; i<(int)span_inf_0.shape[0]*3;i+=3){
        for(int i1 = 0 ; i1<(int)span_inf_1.shape[0]*3;i1+=3){
            for(int i2= 0 ; i2<(int)span_inf_2.shape[0]*3;i2+=3){
                out_ptr[i*span_inf_1.shape[0]*span_inf_2.shape[0]+i1*span_inf_2.shape[0]+i2] = origin_ptr[0]+span_ptr_2[i]+span_ptr_1[i1]+span_ptr_0[i2];
                out_ptr[i*span_inf_1.shape[0]*span_inf_2.shape[0]+i1*span_inf_2.shape[0]+i2+1] = origin_ptr[1]+span_ptr_2[i+1]+span_ptr_1[i1+1]+span_ptr_0[i2+1];
                out_ptr[i*span_inf_1.shape[0]*span_inf_2.shape[0]+i1*span_inf_2.shape[0]+i2+2] = origin_ptr[2]+span_ptr_2[i+2]+span_ptr_1[i1+2]+span_ptr_0[i2+2];
            }
        }
    }
    
    return out ;
}

Cd inner_product_3D(array_Cd mesh1,array_Cd mesh2){
    py::buffer_info mesh1Info = mesh1.request();
    py::buffer_info mesh2Info = mesh2.request();
    
    assert(mesh1Info.shape[0]==mesh2Info.shape[0]);
    assert(mesh1Info.shape[1]==mesh2Info.shape[1]);
    assert(mesh1Info.shape[2]==mesh2Info.shape[2]);
    
    auto mesh1Ptr = static_cast<Cd *>(mesh1Info.ptr);
    auto mesh2Ptr = static_cast<Cd *>(mesh2Info.ptr);
    
    Ul dim0 = mesh1Info.shape[0];
    Ul dim1 = mesh1Info.shape[1];
    Ul dim2 = mesh1Info.shape[2];
    Cd inner_product = 0;
    double n1 = 0;
    double n2 = 0;
    for (Ul i = 0; i<dim0; i++) {
        for(Ul i1 = 0; i1<dim1;i1++){
            for (Ul i2 = 28; i2<dim2; i2++) {
                inner_product += std::conj(mesh1Ptr[dim1*dim2*i+dim2*i1+i2]) * mesh2Ptr[dim2*dim1*i+dim2*i1+i2-28];
                n1 += std::norm(mesh1Ptr[i*dim1*dim2+i1*dim2+i2]);
                n2 += std::norm(mesh2Ptr[dim2*dim1*i+dim2*i1+i2-28]);
            }
        }
    }
    double n = std::sqrt(n1*n2);
    return inner_product/n;
}

array_Cd normalize(array_Cd mesh){
    //Normalizes to total probability of 1 over the grid. No volume considerations are taken into account.
    py::buffer_info meshInfo = mesh.request();
    
    auto meshPtr = static_cast<Cd  *>(meshInfo.ptr);
    array_Cd out({meshInfo.shape[0],meshInfo.shape[1],meshInfo.shape[2]});
    py::buffer_info outInfo = out.request();
    auto outPtr = static_cast<Cd  *>(outInfo.ptr);
    double norm_ = 0;
    for(Ul i = 0; i < meshInfo.shape[0];i++){
        for(Ul i1 = 0; i1 < meshInfo.shape[1];i1++){
            for(Ul i2 = 0; i2 < meshInfo.shape[2];i2++){
                norm_+= std::norm(meshPtr[meshInfo.shape[1]*meshInfo.shape[2]*i+meshInfo.shape[2]*i1+i2]);
            }
        }
    }
    norm_ = sqrt(norm_);
    for(Ul i = 0; i < meshInfo.shape[0];i++){
        for(Ul i1 = 0; i1 < meshInfo.shape[1];i1++){
            for(Ul i2 = 0; i2 < meshInfo.shape[2];i2++){
                outPtr[meshInfo.shape[1]*meshInfo.shape[2]*i+meshInfo.shape[2]*i1+i2] = meshPtr[meshInfo.shape[1]*meshInfo.shape[2]*i+meshInfo.shape[2]*i1+i2]/norm_;
            }
        }
    }
    return out;
}

array_Cd construct_bloch_sum_center(array_Cd mesh, double kx, double ky, double kz){
    py::buffer_info meshInfo = mesh.request();
    auto meshPtr = static_cast<Cd *>(meshInfo.ptr);
    
    Ul dim0 = meshInfo.shape[0];
    Ul dim1 = meshInfo.shape[1];
    Ul dim2 = meshInfo.shape[2];
    
    Cd fac_x1 = std::exp(Cd(0,2*M_PI*kx));
    Cd fac_y1 = std::exp(Cd(0,2*M_PI*ky));
    Cd fac_z1 = std::exp(Cd(0,2*M_PI*kz));
    Cd fac_x1_ = std::exp(Cd(0,-2*M_PI*kx));
    Cd fac_y1_ = std::exp(Cd(0,-2*M_PI*ky));
    Cd fac_z1_ = std::exp(Cd(0,-2*M_PI*kz));

    array_Cd out({dim0,dim1,dim2});
    auto outPtr = static_cast<Cd *>(out.request().ptr);
    for(Ul i = 0; i < dim0;i++){
        for(Ul i1 = 0; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < dim2;i2++){
                outPtr[i*dim1*dim2+i1*dim2+i2] = meshPtr[i*dim1*dim2+i1*dim2+i2];
            }
        }
    }
// exp(+ikR)w(r+R)
    for(Ul i = 0; i < dim0;i++){
        for(Ul i1 = 0; i1 < dim1;i1++){
            for(Ul i2 = dim2/3+1; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1*meshPtr[i*dim1*dim2+i1*dim2+i2-dim2/3];
            }
        }
    }
    for(Ul i = 0; i < dim0;i++){
        for(Ul i1 = dim1/3+1; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_y1*meshPtr[i*dim1*dim2+(i1-dim1/3)*dim2+i2];
            }
        }
    }
    for(Ul i = dim0/3+1; i < dim0;i++){
        for(Ul i1 = 0; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_z1*meshPtr[(i-dim0/3)*dim1*dim2+i1*dim2+i2];
            }
        }
    }
    for(Ul i = 0; i < dim0;i++){
        for(Ul i1 = dim1/3+1; i1 < dim1;i1++){
            for(Ul i2 = dim2/3+1; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1*fac_y1*meshPtr[i*dim1*dim2+(i1-dim1/3)*dim2+i2-dim2/3];
            }
        }
    }
    for(Ul i = dim0/3+1; i < dim0;i++){
        for(Ul i1 = 0; i1 < dim1;i1++){
            for(Ul i2 = dim2/3+1; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1*fac_z1*meshPtr[(i-dim0/3)*dim1*dim2+i1*dim2+i2-dim2/3];
            }
        }
    }
    for(Ul i = dim0/3+1; i < dim0;i++){
        for(Ul i1 = dim1/3+1; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_y1*fac_z1*meshPtr[(i-dim0/3)*dim1*dim2+(i1-dim1/3)*dim2+i2];
            }
        }
    }
    for(Ul i = dim0/3+1; i < dim0;i++){
        for(Ul i1 = dim1/3+1; i1 < dim1;i1++){
            for(Ul i2 = dim2/3+1; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1*fac_y1*fac_z1*meshPtr[(i-dim0/3)*dim1*dim2+(i1-dim1/3)*dim2+(i2-dim2/3)];
            }
        }
    }

// exp(-ikR)w(r-R)
    for(Ul i = 0; i < dim0;i++){
        for(Ul i1 = 0; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < 2*dim2/3-1;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1_*meshPtr[i*dim1*dim2+i1*dim2+i2+dim2/3+1];
            }
        }
    }
    for(Ul i = 0; i < dim0;i++){
        for(Ul i1 = 0; i1 < 2*dim1/3-1;i1++){
            for(Ul i2 = 0; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_y1_*meshPtr[i*dim1*dim2+(i1+dim1/3+1)*dim2+i2];
            }
        }
    }
    for(Ul i = 0; i < 2*dim0/3-1;i++){
        for(Ul i1 = 0; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_z1_*meshPtr[(i+dim0/3+1)*dim1*dim2+i1*dim2+i2];
            }
        }
    }
    for(Ul i = 0; i < dim0;i++){
        for(Ul i1 = 0; i1 < 2*dim1/3-1;i1++){
            for(Ul i2 = 0; i2 < 2*dim2/3-1;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1_*fac_y1_*meshPtr[i*dim1*dim2+(i1+dim1/3+1)*dim2+i2+dim2/3+1];
            }
        }
    }
    for(Ul i = 0; i < 2*dim0/3-1;i++){
        for(Ul i1 = 0; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < 2*dim2/3-1;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1_*fac_z1_*meshPtr[(i+dim0/3+1)*dim1*dim2+i1*dim2+i2+dim2/3+1];
            }
        }
    }
    for(Ul i = 0; i < 2*dim0/3-1;i++){
        for(Ul i1 = 0; i1 < 2*dim1/3-1;i1++){
            for(Ul i2 = 0; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_y1_*fac_z1_*meshPtr[(i+dim0/3+1)*dim1*dim2+(i1+dim1/3+1)*dim2+i2];
            }
        }
    }
    for(Ul i = 0; i < 2*dim0/3-1;i++){
        for(Ul i1 = 0; i1 < 2*dim1/3-1;i1++){
            for(Ul i2 = 0; i2 < 2*dim2/3-1;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1_*fac_y1_*fac_z1_*meshPtr[(i+dim0/3+1)*dim1*dim2+(i1+dim1/3+1)*dim2+(i2+dim2/3+1)];
            }
        }
    }
    
    // exp(-ikR1+ikR2)w(r-R1+R2)
    for(Ul i = 0; i < dim0;i++){
        for(Ul i1 = 0; i1 < 2*dim1/3-1;i1++){
            for(Ul i2 = dim2/3+1; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1*fac_y1_*meshPtr[i*dim1*dim2+(i1+dim1/3+1)*dim2+i2-dim2/3];
            }
        }
    }
    for(Ul i = 0; i < 2*dim0/3-1;i++){
        for(Ul i1 = 0; i1 < dim1;i1++){
            for(Ul i2 = dim2/3+1; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1*fac_z1_*meshPtr[(i+dim0/3+1)*dim1*dim2+i1*dim2+i2-dim2/3];
            }
        }
    }
    for(Ul i = 0; i < dim0;i++){
        for(Ul i1 = dim1/3+1; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < 2*dim2/3-1;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1_*fac_y1*meshPtr[i*dim1*dim2+(i1-dim1/3)*dim2+(i2+dim2/3+1)];
            }
        }
    }
    for(Ul i = dim0/3+1; i < dim0;i++){
        for(Ul i1 = 0; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < 2*dim2/3-1;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1_*fac_z1*meshPtr[(i-dim0/3)*dim1*dim2+i1*dim2+(i2+dim2/3+1)];
            }
        }
    }
    for(Ul i = 0; i < 2*dim0/3-1;i++){
        for(Ul i1 = dim1/3+1; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_y1*fac_z1_*meshPtr[(i+dim0/3+1)*dim1*dim2+(i1-dim1/3)*dim1+i2];
            }
        }
    }
    for(Ul i = dim0/3+1; i < dim0;i++){
        for(Ul i1 = 0; i1 < 2*dim1/3-1;i1++){
            for(Ul i2 = 0; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_y1_*fac_z1*meshPtr[(i-dim0/3)*dim1*dim2+(i1+dim1/3+1)*dim1+i2];
            }
        }
    }
    
    // exp(-ikR1-ikR2+ikR3)w(r-R1-R2+R3)
    for(Ul i = dim0/3+1; i < dim0;i++){
        for(Ul i1 = 0; i1 < 2*dim1/3-1;i1++){
            for(Ul i2 = dim2/3+1; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1*fac_y1_*fac_z1*meshPtr[(i-dim0/3)*dim1*dim2+(i1+dim1/3+1)*dim2+i2-dim2/3];
            }
        }
    }
    for(Ul i = 0; i < 2*dim0/3-1;i++){
        for(Ul i1 = 0; i1 < 2*dim1/3-1;i1++){
            for(Ul i2 = dim2/3+1; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1*fac_y1_*fac_z1_*meshPtr[(i+dim0/3+1)*dim1*dim2+(i1+dim1/3+1)*dim2+i2-dim2/3];
            }
        }
    }
    for(Ul i = dim0/3+1; i < dim0;i++){
        for(Ul i1 = dim1/3+1; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < 2*dim2/3-1;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1_*fac_y1*fac_z1*meshPtr[(i-dim0/3)*dim1*dim2+(i1-dim1/3)*dim2+(i2+dim2/3+1)];
            }
        }
    }
    for(Ul i = 0; i < 2*dim0/3-1;i++){
        for(Ul i1 = dim1/3+1; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < 2*dim2/3-1;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1_*fac_y1*fac_z1_*meshPtr[(i+dim0/3+1)*dim1*dim2+(i1-dim1/3)*dim2+(i2+dim2/3+1)];
            }
        }
    }
    for(Ul i = dim0/3+1; i < dim0;i++){
        for(Ul i1 = 0; i1 < 2*dim1/3-1;i1++){
            for(Ul i2 = 0; i2 < 2*dim2/3-1;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1_*fac_y1_*fac_z1*meshPtr[(i-dim0/3)*dim1*dim2+(i1+dim1/3+1)*dim2+i2+dim2/3+1];
            }
        }
    }
    for(Ul i = 0; i < 2*dim0/3-1;i++){
        for(Ul i1 = dim1/3+1; i1 < dim1;i1++){
            for(Ul i2 = dim2/3+1; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1*fac_y1*fac_z1_*meshPtr[(i+dim0/3+1)*dim1*dim2+(i1-dim1/3)*dim2+i2-dim2/3];
            }
        }
    }
    return out;
}

array_Cd construct_bloch_sum_origin(array_Cd mesh, double kx, double ky, double kz){
    py::buffer_info meshInfo = mesh.request();
    auto meshPtr = static_cast<Cd *>(meshInfo.ptr);
    
    Ul dim0 = meshInfo.shape[0];
    Ul dim1 = meshInfo.shape[1];
    Ul dim2 = meshInfo.shape[2];
    
    Cd fac_x1 = std::exp(Cd(0,2*M_PI*kx));
    Cd fac_y1 = std::exp(Cd(0,2*M_PI*ky));
    Cd fac_z1 = std::exp(Cd(0,2*M_PI*kz));
    Cd fac_x2 = std::exp(Cd(0,4*M_PI*kx));
    Cd fac_y2 = std::exp(Cd(0,4*M_PI*ky));
    Cd fac_z2 = std::exp(Cd(0,4*M_PI*kz));
    Cd fac_x1_ = std::exp(Cd(0,-2*M_PI*kx));
    Cd fac_y1_ = std::exp(Cd(0,-2*M_PI*ky));
    Cd fac_z1_ = std::exp(Cd(0,-2*M_PI*kz));
    Cd fac_x2_ = std::exp(Cd(0,-4*M_PI*kx));
    Cd fac_y2_ = std::exp(Cd(0,-4*M_PI*ky));
    Cd fac_z2_ = std::exp(Cd(0,-4*M_PI*kz));
    
    
    array_Cd out({dim0,dim1,dim2});
    auto outPtr = static_cast<Cd *>(out.request().ptr);
    for(Ul i = 0; i < dim0;i++){
        for(Ul i1 = 0; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < dim2;i2++){
                outPtr[i*dim1*dim2+i1*dim2+i2] = meshPtr[i*dim1*dim2+i1*dim2+i2];
            }
        }
    }
    // exp(+ikR)w(r+R)
    for(Ul i = 0; i < dim0;i++){
        for(Ul i1 = 0; i1 < dim1;i1++){
            for(Ul i2 = dim2/3+1; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1*meshPtr[i*dim1*dim2+i1*dim2+i2-dim2/3];
            }
        }
    }
    for(Ul i = 0; i < dim0;i++){
        for(Ul i1 = dim1/3+1; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_y1*meshPtr[i*dim1*dim2+(i1-dim1/3)*dim2+i2];
            }
        }
    }
    for(Ul i = dim0/3+1; i < dim0;i++){
        for(Ul i1 = 0; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_z1*meshPtr[(i-dim0/3)*dim1*dim2+i1*dim2+i2];
            }
        }
    }
    for(Ul i = 0; i < dim0;i++){
        for(Ul i1 = dim1/3+1; i1 < dim1;i1++){
            for(Ul i2 = dim2/3+1; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1*fac_y1*meshPtr[i*dim1*dim2+(i1-dim1/3)*dim2+i2-dim2/3];
            }
        }
    }
    for(Ul i = dim0/3+1; i < dim0;i++){
        for(Ul i1 = 0; i1 < dim1;i1++){
            for(Ul i2 = dim2/3+1; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1*fac_z1*meshPtr[(i-dim0/3)*dim1*dim2+i1*dim2+i2-dim2/3];
            }
        }
    }
    for(Ul i = dim0/3+1; i < dim0;i++){
        for(Ul i1 = dim1/3+1; i1 < dim1;i1++){
            for(Ul i2 = 0; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_y1*fac_z1*meshPtr[(i-dim0/3)*dim1*dim2+(i1-dim1/3)*dim2+i2];
            }
        }
    }
    for(Ul i = dim0/3+1; i < dim0;i++){
        for(Ul i1 = dim1/3+1; i1 < dim1;i1++){
            for(Ul i2 = dim2/3+1; i2 < dim2;i2++){
                outPtr[i*dim0*dim1+i1*dim0+i2] += fac_x1*fac_y1*fac_z1*meshPtr[(i-dim0/3)*dim1*dim2+(i1-dim1/3)*dim2+(i2-dim2/3)];
            }
        }
    }
    
    return out;
}

Cd calculate_angular_momentum(array_Cd mesh1,array_Cd mesh2, array_d grid,array_d centre,AngMomentum angMom){
    
    py::buffer_info mesh1Info = mesh1.request();
    py::buffer_info mesh2Info = mesh2.request();
    py::buffer_info gridInfo = grid.request();
    py::buffer_info centreInfo = centre.request();
    
    unsigned long dim0 = gridInfo.shape[0];
    unsigned long dim1 = gridInfo.shape[1];
    unsigned long dim2 = gridInfo.shape[2];
    unsigned long dim3 = gridInfo.shape[3];
    
    auto mesh1Ptr = static_cast<Cd *>(mesh1Info.ptr);
    auto mesh2Ptr = static_cast<Cd *>(mesh2Info.ptr);
    auto gridPtr = static_cast<double *>(gridInfo.ptr);
    auto centrePtr = static_cast<double *>(centreInfo.ptr);
    
    Cd angMomOut(0,0);
    double n1 = 0;
    double n2 = 0;
    
    Eigen::Vector3d origin = Eigen::Vector3d({gridPtr[0],gridPtr[1],gridPtr[2]});
    Eigen::Vector3d a_1 = Eigen::Vector3d({gridPtr[dim3],gridPtr[dim3+1],gridPtr[dim3+2]})-origin;
    Eigen::Vector3d b_1 = Eigen::Vector3d({gridPtr[dim2*dim3],gridPtr[dim2*dim3+1],gridPtr[dim2*dim3+2]})-origin;
    Eigen::Vector3d c_1 = Eigen::Vector3d({gridPtr[dim1*dim2*dim3],gridPtr[dim1*dim2*dim3+1],gridPtr[dim1*dim2*dim3+2]})-origin;
    Eigen::Matrix3d V ;
    V.row(0)= a_1;
    V.row(1)=b_1;
    V.row(2)=c_1;
    
    Eigen::Matrix3d W = (V.transpose()).inverse();
    
    double dadx = W(0,0);
    double dbdx = W(1,0);
    double dcdx = W(2,0);
    
    double dady =W(0,1);
    double dbdy =W(1,1);
    double dcdy =W(2,1);
    
    double dadz = W(0,2);
    double dbdz = W(1,2);
    double dcdz = W(2,2);

    Cd ddax;
    Cd ddbx;
    Cd ddcx;
    Cd ddx ;
    
    Cd dday;
    Cd ddby;
    Cd ddcy;
    Cd ddy ;
    
    Cd ddaz;
    Cd ddbz;
    Cd ddcz;
    Cd ddz ;
    
    Ul i;
    Ul i1 ;
    Ul i2;
    switch (angMom) {
        case Lx:
            //#pragma omp parallel for shared(angMomOut)
            for (i =0; i<dim0-1; i++) /* span_vec_c */ {
                //                printf("%d\n",omp_get_num_threads());
                for (i1 = 0 ; i1<dim1-1; i1++) /* span_vec_b */ {
                    for (i2 = 0 ; i2<dim2-1; i2++) /* span_vec_a */{
                        dday = (mesh2Ptr[dim1*dim2*i+dim2*i1+i2+1]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dady;
                        ddby = (mesh2Ptr[dim1*dim2*i+dim2*(i1+1)+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dbdy;
                        ddcy = (mesh2Ptr[dim1*dim2*(i+1)+dim2*i1+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dcdy;
                        ddy =  dday+ddby+ddcy;
                        
                        ddaz = (mesh2Ptr[dim1*dim2*i+dim2*i1+i2+1]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dadz;
                        ddbz = (mesh2Ptr[dim1*dim2*i+dim2*(i1+1)+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dbdz;
                        ddcz = (mesh2Ptr[dim1*dim2*(i+1)+dim2*i1+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dcdz;
                        ddz = ddaz + ddbz + ddcz;
                        angMomOut += std::conj(mesh1Ptr[dim1*dim2*i+dim2*i1+i2]) * Cd(0,-h_bar)*(/* z */-(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3+2]-centrePtr[2])*ddy + /* y */(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3+1]-centrePtr[1])*ddz );
                        n1 += std::norm(mesh1Ptr[dim1*dim2*i+dim2*i1+i2]);
                        n2 += std::norm(mesh2Ptr[dim1*dim2*i+dim2*i1+i2]);
                    }
                }
            }
            break;
        case Ly:
            //#pragma omp parallel for private(ddax,ddbx,ddcx,ddx,ddaz,ddbz,ddcz,ddz) shared(angMomOut)
            for (i =0; i<dim0-1; i++) /* span_vec_c */ {
                for (i1 = 0 ; i1<dim1-1; i1++) /* span_vec_b */ {
                    for (i2 = 0; i2<dim2-1; i2++) /* span_vec_a */{
                        
                        ddax = (mesh2Ptr[dim1*dim2*i+dim2*i1+i2+1]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dadx;
                        ddbx = (mesh2Ptr[dim1*dim2*i+dim2*(i1+1)+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dbdx;
                        ddcx = (mesh2Ptr[dim1*dim2*(i+1)+dim2*i1+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dcdx;
                        ddx =  ddax+ddbx+ddcx;
                        
                        ddaz = (mesh2Ptr[dim1*dim2*i+dim2*i1+i2+1]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dadz;
                        ddbz = (mesh2Ptr[dim1*dim2*i+dim2*(i1+1)+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dbdz;
                        ddcz = (mesh2Ptr[dim1*dim2*(i+1)+dim2*i1+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dcdz;
                        ddz = ddaz+ddbz + ddcz;
                        
                        angMomOut += std::conj(mesh1Ptr[dim1*dim2*i+dim2*i1+i2]) * Cd(0,-h_bar)*(/* z */(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3+2]-centrePtr[2])*ddx - /* x */(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3]-centrePtr[0])*ddz );
                        n1 += std::norm(mesh1Ptr[dim1*dim2*i+dim2*i1+i2]);
                        n2 += std::norm(mesh2Ptr[dim1*dim2*i+dim2*i1+i2]);
                    }
                }
            }
            break;
        case Lz:
            //#pragma omp parallel for private(ddax,ddbx,ddcx,ddx,dday,ddby,ddcy,ddy) shared(angMomOut)
            for (i =0; i<dim0-1; i++) /* span_vec_c */ {
                for (i1 = 0 ; i1<dim1-1; i1++) /* span_vec_b */ {
                    for (i2 = 0; i2<dim2-1; i2++) /* span_vec_a */{
                        
                        ddax = (mesh2Ptr[dim1*dim2*i+dim2*i1+i2+1]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dadx;
                        ddbx = (mesh2Ptr[dim1*dim2*i+dim2*(i1+1)+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dbdx;
                        ddcx = (mesh2Ptr[dim1*dim2*(i+1)+dim2*i1+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dcdx;
                        ddx =  ddax+ddbx+ddcx;
                        
                        dday = (mesh2Ptr[dim1*dim2*i+dim2*i1+i2+1]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dady;
                        ddby = (mesh2Ptr[dim1*dim2*i+dim2*(i1+1)+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dbdy;
                        ddcy = (mesh2Ptr[dim1*dim2*(i+1)+dim2*i1+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dcdy;
                        ddy =  dday+ddby+ddcy;
                        
                        angMomOut += std::conj(mesh1Ptr[dim1*dim2*i+dim2*i1+i2]) * Cd(0,-h_bar)*(/* x */(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3]-centrePtr[0])*ddy -/* y */(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3+1]-centrePtr[1])*ddx );
                        n1 += std::norm(mesh1Ptr[dim1*dim2*i+dim2*i1+i2]);
                        n2 += std::norm(mesh2Ptr[dim1*dim2*i+dim2*i1+i2]);
                    }
                }
            }
            
    }
    double n = std::sqrt(n1*n2);
    return angMomOut/n;
}

Cd calculate_angular_momentum_(Cd* mesh1Ptr,Cd* mesh2Ptr, double* gridPtr, double* centrePtr, std::vector<Ul> shape, AngMomentum angMom){
    
    Ul dim0 = shape[0];
    Ul dim1 = shape[1];
    Ul dim2 = shape[2];
    Ul dim3 = shape[3];
    
    Cd angMomOut(0,0);
    double n1 = 0;
    double n2 = 0;
    
    Eigen::Vector3d origin = Eigen::Vector3d({gridPtr[0],gridPtr[1],gridPtr[2]});
    Eigen::Vector3d a_1 = Eigen::Vector3d({gridPtr[dim3],gridPtr[dim3+1],gridPtr[dim3+2]})-origin;
    Eigen::Vector3d b_1 = Eigen::Vector3d({gridPtr[dim2*dim3],gridPtr[dim2*dim3+1],gridPtr[dim2*dim3+2]})-origin;
    Eigen::Vector3d c_1 = Eigen::Vector3d({gridPtr[dim1*dim2*dim3],gridPtr[dim1*dim2*dim3+1],gridPtr[dim1*dim2*dim3+2]})-origin;
    Eigen::Matrix3d V ;
    V.row(0)= a_1;
    V.row(1)= b_1;
    V.row(2)= c_1;
    
    Eigen::Matrix3d W = (V.transpose()).inverse();
    
    double dadx = W(0,0);
    double dbdx = W(1,0);
    double dcdx = W(2,0);
    
    double dady =W(0,1);
    double dbdy =W(1,1);
    double dcdy =W(2,1);
    
    double dadz = W(0,2);
    double dbdz = W(1,2);
    double dcdz = W(2,2);
    
    Cd ddax;
    Cd ddbx;
    Cd ddcx;
    Cd ddx ;
    
    Cd dday;
    Cd ddby;
    Cd ddcy;
    Cd ddy ;
    
    Cd ddaz;
    Cd ddbz;
    Cd ddcz;
    Cd ddz ;
    
    Ul i;
    Ul i1 ;
    Ul i2;
    switch (angMom) {
        case Lx:
            for (i =0; i<dim0-1; i++) /* span_vec_c */ {
                for (i1 = 0 ; i1<dim1-1; i1++) /* span_vec_b */ {
                    for (i2 = 0 ; i2<dim2-1; i2++) /* span_vec_a */{
                        dday = (mesh2Ptr[dim1*dim2*i+dim2*i1+i2+1]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dady;
                        ddby = (mesh2Ptr[dim1*dim2*i+dim2*(i1+1)+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dbdy;
                        ddcy = (mesh2Ptr[dim1*dim2*(i+1)+dim2*i1+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dcdy;
                        ddy =  dday+ddby+ddcy;
                        
                        ddaz = (mesh2Ptr[dim1*dim2*i+dim2*i1+i2+1]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dadz;
                        ddbz = (mesh2Ptr[dim1*dim2*i+dim2*(i1+1)+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dbdz;
                        ddcz = (mesh2Ptr[dim1*dim2*(i+1)+dim2*i1+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dcdz;
                        ddz = ddaz + ddbz + ddcz;
                        
                        angMomOut += std::conj(mesh1Ptr[dim1*dim2*i+dim2*i1+i2]) * Cd(0,-h_bar)*(/* z */-(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3+2]-centrePtr[2])*ddy + /* y */(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3+1]-centrePtr[1])*ddz );
//                        n+=std::sqrt(std::abs(mesh1Ptr[i*dim1*dim2+i1*dim2+i2])*std::abs(mesh2Ptr[i*dim1*dim2+i1*dim2+i2]));
                        n1 += std::norm(mesh1Ptr[i*dim1*dim2+i1*dim2+i2]);
                        n2 += std::norm(mesh2Ptr[i*dim1*dim2+i1*dim2+i2]);

                    }
                }
            }
            break;
        case Ly:
            for (i =0; i<dim0-1; i++) /* span_vec_c */ {
                for (i1 = 0 ; i1<dim1-1; i1++) /* span_vec_b */ {
                    for (i2 = 0; i2<dim2-1; i2++) /* span_vec_a */{
                        
                        ddax = (mesh2Ptr[dim1*dim2*i+dim2*i1+i2+1]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dadx;
                        ddbx = (mesh2Ptr[dim1*dim2*i+dim2*(i1+1)+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dbdx;
                        ddcx = (mesh2Ptr[dim1*dim2*(i+1)+dim2*i1+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dcdx;
                        ddx =  ddax+ddbx+ddcx;
                        
                        ddaz = (mesh2Ptr[dim1*dim2*i+dim2*i1+i2+1]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dadz;
                        ddbz = (mesh2Ptr[dim1*dim2*i+dim2*(i1+1)+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dbdz;
                        ddcz = (mesh2Ptr[dim1*dim2*(i+1)+dim2*i1+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dcdz;
                        ddz = ddaz+ddbz + ddcz;
                        
                        angMomOut += std::conj(mesh1Ptr[dim1*dim2*i+dim2*i1+i2]) * Cd(0,-h_bar)*(/* z */(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3+2]-centrePtr[2])*ddx - /* x */(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3]-centrePtr[0])*ddz );
//                        n+=std::sqrt(std::abs(mesh1Ptr[i*dim1*dim2+i1*dim2+i2])*std::abs(mesh2Ptr[i*dim1*dim2+i1*dim2+i2]));
                        n1 += std::norm(mesh1Ptr[i*dim1*dim2+i1*dim2+i2]);
                        n2 += std::norm(mesh2Ptr[i*dim1*dim2+i1*dim2+i2]);
                    }
                }
            }
            break;
        case Lz:
            for (i =0; i<dim0-1; i++) /* span_vec_c */ {
                for (i1 = 0 ; i1<dim1-1; i1++) /* span_vec_b */ {
                    for (i2 = 0; i2<dim2-1; i2++) /* span_vec_a */{
                        
                        ddax = (mesh2Ptr[dim1*dim2*i+dim2*i1+i2+1]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dadx;
                        ddbx = (mesh2Ptr[dim1*dim2*i+dim2*(i1+1)+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dbdx;
                        ddcx = (mesh2Ptr[dim1*dim2*(i+1)+dim2*i1+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dcdx;
                        ddx =  ddax+ddbx+ddcx;
                        
                        dday = (mesh2Ptr[dim1*dim2*i+dim2*i1+i2+1]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dady;
                        ddby = (mesh2Ptr[dim1*dim2*i+dim2*(i1+1)+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dbdy;
                        ddcy = (mesh2Ptr[dim1*dim2*(i+1)+dim2*i1+i2]-mesh2Ptr[dim1*dim2*i+dim2*i1+i2])*dcdy;
                        ddy =  dday+ddby+ddcy;
                        
                        angMomOut += std::conj(mesh1Ptr[dim1*dim2*i+dim2*i1+i2]) * Cd(0,-h_bar)*(/* x */(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3]-centrePtr[0])*ddy -/* y */(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3+1]-centrePtr[1])*ddx );
//                        n+=std::sqrt(std::abs(mesh1Ptr[i*dim1*dim2+i1*dim2+i2])*std::abs(mesh2Ptr[i*dim1*dim2+i1*dim2+i2]));
                        n1 += std::norm(mesh1Ptr[i*dim1*dim2+i1*dim2+i2]);
                        n2 += std::norm(mesh2Ptr[i*dim1*dim2+i1*dim2+i2]);
                        
                    }
                }
            }
            
    }
    double n = std::sqrt(n1*n2);
    return angMomOut/n;
}

array_Cd project_on_angular_momentum(array_Cd mesh, array_d grid,array_d centre, AngMomentum angMom){
    
    py::buffer_info meshInfo = mesh.request();
    py::buffer_info gridInfo = grid.request();
    py::buffer_info centreInfo = centre.request();
    
    unsigned long dim0 = gridInfo.shape[0];
    unsigned long dim1 = gridInfo.shape[1];
    unsigned long dim2 = gridInfo.shape[2];
    unsigned long dim3 = gridInfo.shape[3];
    
    auto meshPtr = static_cast<Cd *>(meshInfo.ptr);
    auto gridPtr = static_cast<double *>(gridInfo.ptr);
    auto centrePtr = static_cast<double *>(centreInfo.ptr);
    
    array_Cd out(std::vector<size_t>({dim0,dim1,dim2}));
    py::buffer_info outInfo = out.request();
    auto outPtr = static_cast<Cd *>(outInfo.ptr);
    
    Eigen::Vector3d origin = Eigen::Vector3d({gridPtr[0],gridPtr[1],gridPtr[2]});
    Eigen::Vector3d a_1 = Eigen::Vector3d({gridPtr[dim3],gridPtr[dim3+1],gridPtr[dim3+2]})-origin;
    Eigen::Vector3d b_1 = Eigen::Vector3d({gridPtr[dim2*dim3],gridPtr[dim2*dim3+1],gridPtr[dim2*dim3+2]})-origin;
    Eigen::Vector3d c_1 = Eigen::Vector3d({gridPtr[dim1*dim2*dim3],gridPtr[dim1*dim2*dim3+1],gridPtr[dim1*dim2*dim3+2]})-origin;
    Eigen::Matrix3d V ;
    V.row(0)= a_1;
    V.row(1)=b_1;
    V.row(2)=c_1;
    
    
    Eigen::Matrix3d W = (V.transpose()).inverse();
    
    double dadx = W(0,0);
    double dbdx = W(1,0);
    double dcdx = W(2,0);
    
    double dady = W(0,1);
    double dbdy = W(1,1);
    double dcdy = W(2,1);
    
    double dadz =W(0,2);
    double dbdz =W(1,2);
    double dcdz =W(2,2);
    
    switch (angMom) {
        case Lx:
            for (unsigned long i =0; i<dim0-1; i++) /* span_vec_c */ {
                for (unsigned long i1 = 0 ; i1<dim1-1; i1++) /* span_vec_b */ {
                    for (unsigned long i2 = 0; i2<dim2-1; i2++) /* span_vec_a */{
                        
                        Cd dday = (meshPtr[dim1*dim2*i+dim2*i1+i2+1]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dady;
                        Cd ddby = (meshPtr[dim1*dim2*i+dim2*(i1+1)+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dbdy;
                        Cd ddcy = (meshPtr[dim1*dim2*(i+1)+dim2*i1+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dcdy;
                        Cd ddy =  dday+ddby+ddcy;
                        
                        Cd ddaz = (meshPtr[dim1*dim2*i+dim2*i1+i2+1]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dadz;
                        Cd ddbz = (meshPtr[dim1*dim2*i+dim2*(i1+1)+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dbdz;
                        Cd ddcz = (meshPtr[dim1*dim2*(i+1)+dim2*i1+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dcdz;
                        Cd ddz = ddaz+ddbz + ddcz;
                        
                        outPtr[dim1*dim2*i+dim2*i1+i2] =  Cd(0,-h_bar)*(-(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3+2]-centrePtr[2])*ddy + /* y */(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3+1]-centrePtr[1])*ddz );
                    }
                }
            }
            break;
        case Ly:
            for (unsigned long i =0; i<dim0-1; i++) /* span_vec_c */ {
                for (unsigned long i1 = 0 ; i1<dim1-1; i1++) /* span_vec_b */ {
                    for (unsigned long i2 = 0; i2<dim2-1; i2++) /* span_vec_a */{
                        
                        Cd ddax = (meshPtr[dim1*dim2*i+dim2*i1+i2+1]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dadx;
                        Cd ddbx = (meshPtr[dim1*dim2*i+dim2*(i1+1)+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dbdx;
                        Cd ddcx = (meshPtr[dim1*dim2*(i+1)+dim2*i1+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dcdx;
                        Cd ddx =  ddax+ddbx+ddcx;
                        
                        Cd ddaz = (meshPtr[dim1*dim2*i+dim2*i1+i2+1]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dadz;
                        Cd ddbz = (meshPtr[dim1*dim2*i+dim2*(i1+1)+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dbdz;
                        Cd ddcz = (meshPtr[dim1*dim2*(i+1)+dim2*i1+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dcdz;
                        Cd ddz = ddaz+ddbz + ddcz;
                        
                        outPtr[dim1*dim2*i+dim2*i1+i2] = Cd(0,-h_bar)*(/* z */(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3+2]-centrePtr[2])*ddx - /* x */(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3]-centrePtr[0])*ddz );
                    }
                }
            }
            break;
        case Lz:
            for (unsigned long i =0; i<dim0-1; i++) /* span_vec_c */ {
                for (unsigned long i1 = 0 ; i1<dim1-1; i1++) /* span_vec_b */ {
                    for (unsigned long i2 = 0; i2<dim2-1; i2++) /* span_vec_a */{
                        
                        Cd ddax = (meshPtr[dim1*dim2*i+dim2*i1+i2+1]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dadx;
                        Cd ddbx = (meshPtr[dim1*dim2*i+dim2*(i1+1)+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dbdx;
                        Cd ddcx = (meshPtr[dim1*dim2*(i+1)+dim2*i1+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dcdx;
                        Cd ddx =  ddax+ddbx+ddcx;
                        
                        Cd dday = (meshPtr[dim1*dim2*i+dim2*i1+i2+1]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dady;
                        Cd ddby = (meshPtr[dim1*dim2*i+dim2*(i1+1)+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dbdy;
                        Cd ddcy = (meshPtr[dim1*dim2*(i+1)+dim2*i1+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dcdy;
                        Cd ddy =  dday+ddby+ddcy;
                        
                        outPtr[dim1*dim2*i+dim2*i1+i2] = Cd(0,-h_bar)*(/* x */(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3]-centrePtr[0])*ddy - /* y */(gridPtr[dim1*dim2*dim3*i+dim2*dim3*i1+i2*dim3+1]-centrePtr[1])*ddx );
                        
                    }
                }
            }
            
    }
    
    return out;
}

array_Cd project_on_momentum(array_Cd mesh, array_d grid,Momentum mom){
    py::buffer_info meshInfo = mesh.request();
    py::buffer_info gridInfo = grid.request();
    
    assert(meshInfo.shape[0]==gridInfo.shape[0]);
    assert(meshInfo.shape[1]==gridInfo.shape[1]);
    assert(meshInfo.shape[2]==gridInfo.shape[2]);
    
    auto meshPtr = static_cast<Cd *>(meshInfo.ptr);
    auto gridPtr = static_cast<double *>(gridInfo.ptr);
    
    unsigned long dim0 = gridInfo.shape[0];
    unsigned long dim1 = gridInfo.shape[1];
    unsigned long dim2 = gridInfo.shape[2];
    unsigned long dim3 = gridInfo.shape[3];
    
    Eigen::Vector3d origin = Eigen::Vector3d({gridPtr[0],gridPtr[1],gridPtr[2]});
    Eigen::Vector3d a_1 = Eigen::Vector3d({gridPtr[dim3],gridPtr[dim3+1],gridPtr[dim3+2]})-origin;
    Eigen::Vector3d b_1 = Eigen::Vector3d({gridPtr[dim2*dim3],gridPtr[dim2*dim3+1],gridPtr[dim2*dim3+2]})-origin;
    Eigen::Vector3d c_1 = Eigen::Vector3d({gridPtr[dim1*dim2*dim3],gridPtr[dim1*dim2*dim3+1],gridPtr[dim1*dim2*dim3+2]})-origin;
    Eigen::Matrix3d V ;
    V.row(0)= a_1;
    V.row(1)=b_1;
    V.row(2)=c_1;
    
    
    Eigen::Matrix3d W = (V.transpose()).inverse();
    
    double dadx = W(0,0);
    double dbdx = W(1,0);
    double dcdx = W(2,0);
    
    double dady = W(0,1);
    double dbdy = W(1,1);
    double dcdy = W(2,1);
    
    double dadz =W(0,2);
    double dbdz =W(1,2);
    double dcdz =W(2,2);
    
    array_Cd out(std::vector<size_t>({dim0,dim1,dim2}));
    py::buffer_info outInfo = out.request();
    auto outPtr = static_cast<Cd *>(outInfo.ptr);
    
    switch (mom) {
        case Px:
            for (unsigned long i =0; i<dim0-1; i++) /* span_vec_c */ {
                for (unsigned long i1 = 0 ; i1<dim1-1; i1++) /* span_vec_b */ {
                    for (unsigned long i2 = 0; i2<dim2-1; i2++) /* span_vec_a */{
                        
                        Cd ddax = (meshPtr[dim1*dim2*i+dim2*i1+i2+1]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dadx;
                        Cd ddbx = (meshPtr[dim1*dim2*i+dim2*(i1+1)+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dbdx;
                        Cd ddcx = (meshPtr[dim1*dim2*(i+1)+dim2*i1+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dcdx;
                        outPtr[dim1*dim2*i+dim2*i1+i2] = Cd(0,-h_bar)*(ddax+ddbx+ddcx);
                        
                    }
                }
            }
        case Py:
            for (unsigned long i =0; i<dim0-1; i++) /* span_vec_c */ {
                for (unsigned long i1 = 0 ; i1<dim1-1; i1++) /* span_vec_b */ {
                    for (unsigned long i2 = 0; i2<dim2-1; i2++) /* span_vec_a */{
                        Cd dday = (meshPtr[dim1*dim2*i+dim2*i1+i2+1]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dady;
                        Cd ddby = (meshPtr[dim1*dim2*i+dim2*(i1+1)+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dbdy;
                        Cd ddcy = (meshPtr[dim1*dim2*(i+1)+dim2*i1+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dcdy;
                        outPtr[dim1*dim2*i+dim2*i1+i2] = Cd(0,-h_bar)*(dday + ddby + ddcy);
                    }
                }
            }
        case Pz:
            for (unsigned long i =0; i<dim0-1; i++) /* span_vec_c */ {
                for (unsigned long i1 = 0 ; i1<dim1-1; i1++) /* span_vec_b */ {
                    for (unsigned long i2 = 0; i2<dim2-1; i2++) /* span_vec_a */{
                        
                        Cd ddaz = (meshPtr[dim1*dim2*i+dim2*i1+i2+1]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dadz;
                        Cd ddbz = (meshPtr[dim1*dim2*i+dim2*(i1+1)+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dbdz;
                        Cd ddcz = (meshPtr[dim1*dim2*(i+1)+dim2*i1+i2]-meshPtr[dim1*dim2*i+dim2*i1+i2])*dcdz;
                        outPtr[dim1*dim2*i+dim2*i1+i2] = Cd(0,-h_bar)*(ddaz+ddbz + ddcz);
                        
                    }
                }
            }
    }
    return out;
}

Cd calculate_overlap_pot(array_Cd mesh1, array_Cd mesh2,Direction d){
    auto shape = mesh1.request().shape;
    auto meshPtrs = create_meshPtrs(std::vector<array_Cd>({mesh1,mesh2}),shape);
    auto mesh1Ptr = meshPtrs[0];
    auto mesh2Ptr = meshPtrs[1];
    Cd out(0.0,0.0);
    switch(d){
        case END:
            for(Ul i = 0; i<shape[0];i++){
                for(Ul i1 = 0; i1<shape[1];i1++){
                    for(Ul i2 = 0; i2<shape[2];i2++){
                        out += std::conj(mesh1Ptr[i*shape[1]*shape[2]+i1*shape[2]+i2])*mesh2Ptr[i*shape[1]*shape[2]+i1*shape[2]+i2];
                    }
                }
            }
    }
    return out;
}

vec_Cd calculate_center(array_Cd mesh,array_d grid){
    py::buffer_info gridInfo = grid.request();
    py::buffer_info meshInfo = mesh.request();
    
    Ul dim0 = gridInfo.shape[0];
    Ul dim1 = gridInfo.shape[1];
    Ul dim2 = gridInfo.shape[2];
    Ul dim3 = gridInfo.shape[3];
    
    auto gridPtr = static_cast<double *>(gridInfo.ptr);
    auto meshPtr = static_cast<Cd *>(meshInfo.ptr);
    
    Cd x(0,0);
    Cd y(0,0);
    Cd z(0,0);
    
    for (Ul i = 0; i<dim0; i++) {
        for (Ul i1 = 0; i1<dim1; i1++) {
            for (Ul i2 = 0; i2<dim2 ;i2++) {
                x += meshPtr[i*dim1*dim2+i1*dim2+i2]*meshPtr[i*dim1*dim2+i1*dim2+i2]*gridPtr[i*dim1*dim2*dim3+i1*dim2*dim3+i2*dim3];
                y += meshPtr[i*dim1*dim2+i1*dim2+i2]*meshPtr[i*dim1*dim2+i1*dim2+i2]*gridPtr[i*dim1*dim2*dim3+i1*dim2*dim3+i2*dim3+1];
                z += meshPtr[i*dim1*dim2+i1*dim2+i2]*meshPtr[i*dim1*dim2+i1*dim2+i2]*gridPtr[i*dim1*dim2*dim3+i1*dim2*dim3+i2*dim3+2];
            }
        }
    }
    
    return vec_Cd({x,y,z});
}

void init_sup_calc(py::module& m ){
    m.def("construct_px",&construct_px,"construct a test px mesh","grid");
    m.def("construct_py",&construct_py,"construct a test py mesh","grid");
    m.def("construct_pz",&construct_pz,"construct a test pz mesh","grid");
    m.def("generate_grid",&generate_grid,"Generate grid given an origin and spanning arrays","origin","a_span_array","b_span_array","c_span_array");
    m.def("inner_product_3D",&inner_product_3D,"Calculate the inner product of two wavefunctions in 3D with the same dimensions","mesh1","mesh2");
    m.def("normalize",&normalize,"Normalize a wavefunction in 3D.","mesh");
    m.def("calculate_angular_momentum",&calculate_angular_momentum,"calculate angular momentum.","mesh1","mesh2","grid","centre","AngMom");
    m.def("project_on_angular_momentum",&project_on_angular_momentum,"Apply angular momentum operator to a wfc.","mesh","grid","centre","AngMom");
    m.def("project_on_momentum",&project_on_momentum,"Apply momentum operator to a wfc.","mesh","grid","mom");
    m.def("construct_bloch_sum_center",&construct_bloch_sum_center,"construct the bloch sum for a wannier function centered in the unit cell.","mesh","kx","ky","kz");
    m.def("construct_bloch_sum_origin",&construct_bloch_sum_origin,"construct the bloch sum for a wannier function centered in the unit cell.","mesh","kx","ky","kz");
    m.def("calculate_overlap_pot",&calculate_overlap_pot,"Calculate the overlap potential between two meshes","mesh1","mesh2","direction");
}

