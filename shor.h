#ifndef SHOR_H
#define SHOR_H
#include "embed_points.h"
#include <Eigen/Core>
#include <iostream>
#include <igl/triangle/triangulate.h>
#include <igl/slice.h>
#include <igl/list_to_matrix.h>

short orientation(const Eigen::Matrix<double,3,2>& P);
struct Angle
{
  Eigen::RowVector2d w1,v,w2;
  int r;
  Angle(){r=0;}
  Angle(const Eigen::RowVector2d& x,
        const Eigen::RowVector2d& y,
        const Eigen::RowVector2d& z): w1(x),v(y),w2(z),r(0){}
  Angle operator+(Angle K){
    assert(K.w1 == w2);
    auto w3 = K.w2;
    Angle S(w1,v,w3);

    Eigen::Matrix<double,3,2> tri1,tri2,tri3;
    tri1<<w1,v,w2;  //w1,v,w2
    tri2<<w2,v,w3;  //w2,v,w3
    tri3<<w1,v,w3;  //w1,v,w3

    short ori1 = orientation(tri1);
    short ori2 = orientation(tri2);
    
    bool reflex_or_180_1 = (ori1 < 0);
    bool reflex_or_180_2 = (ori2 < 0);
    bool is_zero_1 = false;
    bool is_zero_2 = false;
    if(ori1 == 0) // w1,v,w2 are collinear
    {
      Eigen::Matrix<double,3,2> T1;
      Eigen::Matrix<double,3,2> T2;
      Eigen::RowVector2d ref;
      ref<<0,0;
      T1<<ref,v,w1;
      if(orientation(T1) == 0){
        T1.row(0)<<1,0;
        ref<<1,0;
        if(orientation(T1) == 0){
          T1.row(0)<<0,1;
          ref<<0,1;
        }
      }
      T2<<ref,v,w2;
      if((orientation(T1) < 0) != (orientation(T2) < 0)) // w1 and w2 on different sides of v
        reflex_or_180_1 = true;
            else 
                is_zero_1 = true;
    }
    if(ori2 == 0) // w3,v,w2 are collinear
    {
      Eigen::Matrix<double,3,2> T1;
      Eigen::Matrix<double,3,2> T2;
      Eigen::RowVector2d ref;
      ref<<0,0;
      T1<<ref,v,w3;
      if(orientation(T1) == 0){
        T1.row(0)<<1,0;
        ref<<1,0;
        if(orientation(T1) == 0){
          T1.row(0)<<0,1;
          ref<<0,1;
        }
      }
      T2<<ref,v,w2;
      if((orientation(T1) < 0) != (orientation(T2) < 0)) // w1 and w2 on different sides of v
        reflex_or_180_2 = true;
            else
                is_zero_2 = true;
    }
        // COUNTER-EXAMPLE: PI + 0
    bool add_one = (is_zero_1 || is_zero_2) ? false : (((reflex_or_180_1 != reflex_or_180_2) && orientation(tri3)>=0 ) || (reflex_or_180_1 && reflex_or_180_2));

    if(add_one)
      S.r = r + K.r + 1;
    else
      S.r = r + K.r;
    return S;
  };
};

void simplify_triangulation(
  const Eigen::MatrixXd& V_i,
  const Eigen::MatrixXd& C,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F
);


// test for weakly-self-overlapping
bool weakly_self_overlapping(
  const Eigen::MatrixXd& P,
  const Eigen::VectorXi& R,
  Eigen::MatrixXi& F
);

void subdivide_polygon(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  std::vector<std::vector<int>>& L
);

void drop_colinear(
  const Eigen::MatrixXd& P,
  const Eigen::VectorXi& R,
  Eigen::VectorXi& B,
  Eigen::MatrixXd& mP,
  Eigen::VectorXi& mR
);

void add_colinear(
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXi& nF,
  const Eigen::VectorXi& B,
  Eigen::MatrixXi& F
);

void add_triangle(
  Eigen::MatrixXi& F, 
  int i, int j, 
  std::vector<std::vector<int>>& K,
  std::vector<std::vector<int>>& Q
);

void set_rotation_index(
  const Eigen::MatrixXd& uv,
  const Eigen::MatrixXi& F, 
  Eigen::VectorXi& R,
  int offset=0
);

bool Shor_van_wyck(
  const Eigen::MatrixXd& P,
  const Eigen::VectorXi& R,
  const std::string flags,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  bool do_refine=true // false meaning no internal vertices
);

bool Shor_van_wyck_v2(
  const Eigen::MatrixXd& P,
  const Eigen::VectorXi& R,
  const Eigen::VectorXi& mark,
  const std::string flags,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  bool do_refine
);

#endif