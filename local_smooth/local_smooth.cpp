#include "local_smooth.h"
#include "../validity_check.h"
#include "auto_grad.hpp"
#include <igl/vertex_triangle_adjacency.h>
#include <iostream>
#include <igl/Timer.h>
#include "coloring_mesh.h"
#include <tbb/parallel_for.h>

// #define DEBUG

void local_smoothing(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::VectorXi& U, // constraints marked as 1
  Eigen::MatrixXd& uv,
  int loop,
  double eps,
  double avg
){
  
  std::vector<std::vector<int>> VF,VFi;
  igl::vertex_triangle_adjacency(V,F,VF,VFi);

  Eigen::Matrix<double,2,3> G_t;
  grad_to_eqtri(avg,G_t);

  Eigen::VectorXd Energy(F.rows());
  Eigen::VectorXd VE(V.rows());
  auto step = [&]()->int{
    // collect list of vertices to be updated
    // e.g. have max 1-ring energy > eps
    std::vector<int> lt;
    for(int f=0;f<F.rows();f++){
      Eigen::Matrix<double,3,2> Tuv;
      Tuv<<uv.row(F(f,0)),uv.row(F(f,1)),uv.row(F(f,2));
      Energy(f) = autogen::sd_energy(Tuv,G_t);
    }
    for(int i=0;i<uv.rows();i++){
      auto nb = VF[i];
      double max_energy = std::numeric_limits<double>::min();
      for(int f: nb){
        max_energy = std::max(max_energy,Energy(f));
        VE(i) = max_energy;
      }
      if(max_energy > eps && U(i) != 1)
        lt.push_back(i);
    }
        
    #ifdef DEBUG
    std::cout<<"size of L "<<lt.size()<<std::endl;
    #endif
  
    // put vertices of the same color into same group
    Eigen::VectorXi C;
    coloring_mesh(F,C);
    std::vector<std::vector<int>> groups;
    groups.resize(C.maxCoeff()+1);
    for(auto id: lt){
      groups[C[id]].push_back(id);
    } 
      
    // optmization begins
    Eigen::VectorXi updated;
    updated.setZero(V.rows());
    for(auto group : groups){
      tbb::parallel_for(tbb::blocked_range<size_t>(0, group.size()), 
      [&](const tbb::blocked_range<size_t>& r){
        for(size_t id=r.begin();id!=r.end();id++){
          int i = group[id];
          
          // compute Jacobian and Hessian for every vertex i
          // which is the sum of the J and H of adjacent faces
          Eigen::Matrix<double,1,2> J;
          Eigen::Matrix<double,2,2> H;
          J.setZero();
          H.setZero();
          double old_energy = 0.0;
          
          // compute J and H for vertex i
          for(int k=0;k<VF[i].size();k++){
            int f = VF[i][k];
            
            // prapare the target shape
            Eigen::Matrix<double,3,2> Tuv;
            for(int t=0;t<3;t++){
              int s = (t+VFi[i][k])%3;
              Tuv.row(t) << uv.row(F(f,s));
            }
            // prapare grad operator for energy above/below eps
            Eigen::Matrix<double,2,3> G_e;
            if(Energy(f) < eps){
              Eigen::Matrix<double,3,3> T0;
              Eigen::Matrix<double,3,3> g;
              T0 << Tuv, Eigen::Vector3d::Zero();
              grad_operator(T0,g);
              G_e = g.topRows(2);
            }else
              grad_to_eqtri(avg, G_e);
            
            Eigen::Matrix<double,1,2> Jl;
            Eigen::Matrix<double,2,2> Hl;
            autogen::sd_grad(Tuv,G_e,Jl);
            autogen::sd_hess(Tuv,G_e,Hl);
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
            Eigen::RowVector2d x0 = uv.row(i);
            double step_size = (alpha*H.inverse()*J.transpose()).norm();
            uv_n = uv.row(i) - (alpha*H.inverse()*J.transpose()).transpose();
            #else
            double step_size = (alpha*J).norm();
            uv_n = uv.row(i) - alpha * J;
            #endif
            double max_ring_energy = 0;
            bool invalid = false;
            for(int k=0;k<VF[i].size();k++){
              int f = VF[i][k];
              uv.row(i) = uv_n;
              Eigen::Matrix<double,3,2> T0;
              T0 << uv.row(F(f,0)),uv.row(F(f,1)),uv.row(F(f,2));
              max_ring_energy = std::max(max_ring_energy, 
                                         autogen::sd_energy(T0,G_t));
              if(!std::isfinite(max_ring_energy) || is_face_flipped(T0)){
                invalid = true;
                break;
              }
            }
            if(invalid || max_ring_energy > VE(i)){
              alpha /= 2;
              uv.row(i) = x0;
              continue;
            }
            
            updated(i) = 1;
            break;
            
          }
        }
      });
      #ifdef DEBUG
        std::cout<<"#updated: "<<updated.sum()<<std::endl;
      #endif
    }
    return 0;
  };
  int itr_g = 0;
  for(int i=0;i<loop;i++){
    #ifdef DEBUG
    std::cout<<"itr_g "<<itr_g<<std::endl;
    #endif
    int flip = step();
    itr_g++;
  }
}