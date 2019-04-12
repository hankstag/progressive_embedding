#include "local_smooth.h"
#include "auto_grad.hpp"
#include <igl/grad.h>
#include <igl/adjacency_list.h>
#include <igl/local_basis.h>
#include <igl/cat.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/jet.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <iostream>
#include <igl/Timer.h>
#include <stdlib.h>     
#include <time.h>  
#include "coloring_mesh.h"
#include <tbb/parallel_for.h>
#define DEBUG

void compute_dirichlet_energy(
  const double target_area,
  const Eigen::MatrixXd& uv,
  const Eigen::MatrixXi& F,
  const std::vector<int>& L,
  const double eps,
  Eigen::VectorXd& E
){
  E.resize(L.size());
  // calc gradient operator
  double h = std::sqrt(target_area/sin(M_PI / 3.0));
  Eigen::Matrix<double,3,3> Tx;
  Tx<<0,0,0,h,0,0,h/2.,(std::sqrt(3)/2.)*h,0;
  Eigen::Matrix<double,3,3> gx;
  grad_operator(Tx,gx);
  Eigen::Matrix<double,2,3> G_t = gx.topRows(2);
  for(int i=0;i<L.size();i++){
    Eigen::Matrix<double,3,2> Tuv;
    int f = L[i];
    Tuv<<uv.row(F(f,0)),uv.row(F(f,1)),uv.row(F(f,2));
    double e = autogen::sd_energy(Tuv,G_t);
    if(e < eps)
      e = 4;
    E(i) = e;
  }
}

// given a face in 3d, return the gradient operator for any double function defined on it
template <typename Scalar>
void grad_operator(
    const Eigen::Matrix<Scalar,3,3>& T, // reference shape
    Eigen::Matrix<Scalar,3,3>& g
){
    Eigen::Matrix<Scalar,1,3> v01 = T.row(1) - T.row(0);
    Eigen::Matrix<Scalar,1,3> v12 = T.row(2) - T.row(1);
    Eigen::Matrix<Scalar,1,3> v20 = T.row(0) - T.row(2);
    Eigen::Matrix<Scalar, 1, 3> n = v01.cross(v12);
    // area of parallelogram is twice area of triangle
    // area of parallelogram is || v1 x v2 ||
    // This does correct l2 norm of rows, so that it contains #F list of twice
    // triangle areas
    Scalar dblA = n.norm();
    Eigen::Matrix<Scalar, 1, 3> u(0,0,1);
    // now normalize normals to get unit normals
    // u = n / dblA;

    Eigen::Matrix<Scalar,1,3> eperp20,eperp01;
    // rotate each vector 90 degrees around normal
    eperp20 = u.cross(v20);
    eperp20.normalize();
    
    eperp20 *= v20.norm()/dblA;
    
    eperp01 = u.cross(v01);
    eperp01.normalize();
    eperp01 *= v01.norm()/dblA;
    g.col(0) = -(eperp20+eperp01);
    g.col(1) = eperp20;
    g.col(2) = eperp01;
}

void local_smoothing(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::VectorXi& M, // constraints marked as 1
    Eigen::MatrixXd& uv,
    int loop,
    double eps
){
    std::vector<std::vector<int>> VF,VFi;
    igl::vertex_triangle_adjacency(V,F,VF,VFi);
    // initialize smallest area
    double target_area = -1e5;
    int itr_g = 0;
    Eigen::VectorXi C;
    coloring_mesh(F,C);
    std::cout<<"# color "<<C.maxCoeff()<<std::endl;
    Eigen::Matrix<double,Eigen::Dynamic,1> Energy;
    //auto cmp = [](std::pair<double,int> a, std::pair<double,int> b) { return a.first > b.first || (a.first == b.first && a.second > b.second);};
    //std::set<std::pair<double,int>,decltype(cmp)> L(cmp);
    int multiplier = 1;
    Eigen::Matrix<double,Eigen::Dynamic,1> A;
    A.setZero(F.rows());
    Eigen::Matrix<double,2,3> G_t;
    std::ofstream mf;
    #ifdef DEBUG
    std::cout<<"export..."<<std::endl;
    //mf.open("stat_"+model_name+".csv");
    //mf<<"itr,min_area,target_area,Energy,flips\n";
    #endif
    // std::cout<<"# v:"<<V.rows()<<std::endl;
    // std::cout<<"# f:"<<F.rows()<<std::endl;
    // extern igl::Timer tom;
    // extern std::ofstream timelog;
    auto step = [&]()->int{
        Energy.setZero(F.rows());
        int c = 0;
        #ifdef DEBUG
        std::cout<<"itr_g "<<itr_g<<std::endl;
        #endif
        // calc average area
        Eigen::VectorXd area;
        target_area = 0;
        igl::doublearea(uv,F,area);
        for(int k=0;k<area.rows();k++){
            if(!std::isnan(area(k)) && !std::isinf(area(k))){
                target_area += area(k);
            }
        }
        //std::cout<<target_area<<" divide into "<<F.rows()<<std::endl;
        target_area /= F.rows();
        if(itr_g == 0){ // update min_area every 5 iterations
            double da = target_area;
            //std::cout<<"target area "<<da<<std::endl;
            // use equilateral with area dlA as target
            double h = std::sqrt( double(da)/sin(M_PI / 3.0));
            Eigen::Matrix<double,3,3> Tx;
            Tx<<0,0,0,h,0,0,h/2.,(std::sqrt(3)/2.)*h,0;
            Eigen::Matrix<double,3,3> gx;
            grad_operator(Tx,gx);
            G_t = gx.topRows(2);
        }
        #ifdef DEBUG
        std::cout<<"current min area "<<area.minCoeff()<<std::endl;
        double my_min = 1e5;
        int min_i = 0;
        for(int i=0;i<area.rows();i++){
            if(my_min > area(i)){
                my_min = area(i);
                min_i = i;
            }
        }
        std::cout<<"confirm min area "<<my_min<<" at "<<min_i<<" ("<<F(min_i,0)<<","<<F(min_i,1)<<","<<F(min_i,2)<<")"<<std::endl;
        std::cout<<"current target area "<<target_area<<std::endl;
        #endif
        std::vector<int> lt;
        //double t2 = timer.getElapsedTime();
        #ifdef DEBUG
        //std::cout<<"laps1 "<<t2-t1<<std::endl;
        #endif
        // mark the 2% of faces that has worst energy
        // sort vertices according to the energy of neighbor
        double max_init_energy;
        Eigen::VectorXi updated_indices;
        updated_indices.setZero(V.rows());
        for(int i=0;i<V.rows();i++){
            //double e = 0;
            double max_ring_energy = 0;
            for(int k=0;k<VF[i].size();k++){
                int f = VF[i][k];
                // Eigen::Matrix<double,3,2> Tuv;
                // for(int t=0;t<3;t++){
                //     // make sure vertex v is always the first in T
                //     Tuv.row(t) << uv.row(F(f,t));
                // }
                // double energy = autogen::sd_energy(Tuv,G_t); 
                Eigen::VectorXd SE;
                compute_dirichlet_energy(target_area,uv,F,{f},eps,SE);
                Energy(f) = SE(0);

                max_ring_energy = std::max(Energy(f),max_ring_energy);
                // if(energy < eps){
                //     Energy(f) = 4;
                //     e = 4;
                // }else{
                //     Energy(f) = energy;
                //     if(std::isinf(energy)){
                //         std::cout<<f<<std::endl;
                //         std::cout<<std::setprecision(17)<<Tuv<<std::endl;
                //         std::cout<<"inf\n";
                //         exit(0);
                //     }
                //     e = std::max(energy,e);
                // } 
            }
            if(max_ring_energy>4.0){
                //L.insert(std::make_pair(max_ring_energy,i));
                lt.push_back(i);
            }
            if(max_init_energy < max_ring_energy)
                max_init_energy = max_ring_energy;
        }
        
        int dt = 0;
        #ifdef DEBUG
        std::cout<<"size of L "<<lt.size()<<std::endl;
        std::cout<<"max init energy "<<max_init_energy<<std::endl;
        #endif
        // genearte a random sequence [0:L.size()]
        // srand (0);
        // std::vector<int> RD;
        // for(int r=0;r<L.size();r++){
        //     RD.push_back(rand() % L.size());
        // }
        // std::set<std::pair<int,std::pair<double,int>>> rand_L;
        // int r = 0;
        // for(auto q: L){
        //     rand_L.insert(std::make_pair(RD[r++],q));
        // }
        std::vector<std::vector<int>> groups;
        groups.resize(C.maxCoeff()+1);
        for(auto id: lt){
            groups[C[id]].push_back(id);
        } 
        for(auto group : groups){
            //int num = 0;
//        for(auto x: lt){
            //if(num++ > V.rows()/4) break;
            tbb::parallel_for(tbb::blocked_range<size_t>( 0, group.size()), 
            [&](const tbb::blocked_range<size_t>& r){
                for(size_t id=r.begin();id!=r.end();id++){
                    int i = group[id];
                    //int i = x;
                    if(M(i) == 1) continue;
                    Eigen::Matrix<double,1,2> J;
                    Eigen::Matrix<double,2,2> H;
                    J.setZero();
                    H.setZero();
                    double old_energy = 0.0;
                    c++;
                    double e_tt = 0.0;
                    for(int k=0;k<VF[i].size();k++){
                        int f = VF[i][k];
                        Eigen::Matrix<double,3,2> Tuv;
                        Eigen::Matrix<double,3,3> T;
                        for(int t=0;t<3;t++){
                            int s = (t+VFi[i][k])%3;
                            Tuv.row(t) << uv.row(F(f,s));
                        }
                        Eigen::Matrix<double,2,3> G_e;
                        if(Energy(f) > eps){
                            G_e = G_t;
                        }else{
                            // use uv itself as target -> meaning try to keep the same
                            T.block(0,0,3,2) = Tuv;
                            T.col(2).setZero();
                            Eigen::Matrix<double,3,3> g;
                            grad_operator(T,g);
                            G_e = g.topRows(2);
                        }
                        
                        Eigen::Matrix<double,1,2> Jl;
                        Eigen::Matrix<double,2,2> Hl;
                        autogen::sd_grad(Tuv,G_e,Jl);
                        autogen::sd_hess(Tuv,G_e,Hl);
                        Energy(f) = autogen::sd_energy(Tuv,G_e);
                        old_energy = std::max(Energy(f),old_energy);
                        J += Jl;
                        H += Hl;
                    }
                    // gradient descent/Newton's method on vertex i
                    double alpha = 1.0;
                    Eigen::Matrix<double,1,2> uv_n;
                    int MAX_IT = 100;
                    int it = 0;
                    for (; it < MAX_IT; it++) {
                        #define NEWTON
                        #ifdef NEWTON
                        Eigen::Matrix<double,Eigen::Dynamic,1> x0 = uv.row(i);
                        double step_size = (alpha*H.inverse()*J.transpose()).norm();
                        uv_n = uv.row(i) - (alpha*H.inverse()*J.transpose()).transpose();
                        #else
                        double step_size = (alpha*J).norm();
                        uv_n = uv.row(i) - alpha * J;
                        #endif
                        // save it for rollback
                        double posx = uv(i,0);
                        double posy = uv(i,1);
                        double max_ring_energy = 0;
                        bool flipped = false;
                        for(int k=0;k<VF[i].size();k++){
                            int f = VF[i][k];

                            uv.row(i) = uv_n;
                            Eigen::VectorXd SE;
                            compute_dirichlet_energy(target_area,uv,F,{f},eps,SE);
                            max_ring_energy = std::max(max_ring_energy,SE(0));

                            Eigen::Matrix<double,3,2> Tuv;
                            for(int t=0;t<3;t++){
                                // make sure vertex v is always the first in T
                                Tuv.row(t) << uv.row(F(f,t));
                            }
                            // trial update position of vertex i
                            Tuv.row(VFi[i][k]) << uv_n;
                            //double energy_trial = autogen::sd_energy(Tuv,G_t); 
                            Eigen::MatrixXi Fl(1,3);
                            Fl<<0,1,2;              
                            if(test_flip(Tuv,Fl) > 0)
                                flipped = true;
                            //std::cout<<new_energy<<std::endl;
                            //new_energy = std::max(energy_trial,new_energy);
                        }
                        if(std::isinf(max_ring_energy) || std::isnan(max_ring_energy) || max_ring_energy > old_energy || flipped == true){
                            alpha /= 2;
                            uv.row(i) << posx,posy;
                            continue;
                        }
                        updated_indices(i) = 1;
                        break;
                    }
                }
            });
            #ifdef DEBUG
            std::cout<<"sum of updated indices: "<<updated_indices.sum()<<std::endl;
            #endif
        }
        return 0; // test_flip(uv,F);
    };
    itr_g = 0;
    int loop_count = 0;
    for(int i=0;i<loop;i++){
        //if(loop_count++ > 500) break;
        int flip = step();
        // if(Energy.maxCoeff()>eps){
        //     i--;
        // }
        #ifdef DEBUG
        std::cout<<"max energy "<<Energy.maxCoeff()<<","<<eps<<std::endl;
        double my_energy = 0;
        int my_ei = 0;
        for(int i=0;i<Energy.rows();i++){
            if(my_energy <Energy(i)){
                my_energy = Energy(i);
                my_ei = i;
            }
        }
        std::cout<<"max energy is at face "<<my_ei<<" ("<<F(my_ei,0)<<","<<F(my_ei,1)<<","<<F(my_ei,2)<<")"<<std::endl;
        #endif
        // extern igl::Timer tom;
        // double tm = tom.getElapsedTime();
        // extern std::ofstream timelog;
        // timelog << tm <<","<<Energy.maxCoeff()<<std::endl;
        itr_g++;
    }
}
template void grad_operator<double>(Eigen::Matrix<double, 3, 3, 0, 3, 3> const&, Eigen::Matrix<double, 3, 3, 0, 3, 3>&);