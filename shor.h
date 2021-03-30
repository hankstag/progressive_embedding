#ifndef SHOR2_H
#define SHOR2_H
#include "embed_points.h"
#include <Eigen/Core>
#include <iostream>
#include <igl/triangle/triangulate.h>
#include <igl/slice.h>
#include <unsupported/Eigen/MPRealSupport>
#include <igl/list_to_matrix.h>

template <typename Scalar>
short orientation(const Eigen::Matrix<Scalar,3,2>& P);

template <typename Scalar>
struct Angle
{
  Eigen::Matrix<Scalar, 1, 2> w1,v,w2;
  int r;
  Angle(){r=0;}
  Angle(const Eigen::Matrix<Scalar, 1, 2>& x,
        const Eigen::Matrix<Scalar, 1, 2>& y,
        const Eigen::Matrix<Scalar, 1, 2>& z): w1(x),v(y),w2(z),r(0){}
  Angle operator+(Angle K){
    assert(K.w1 == w2);
    auto w3 = K.w2;
    Angle<Scalar> S(w1,v,w3);

    Eigen::Matrix<Scalar,3,2> tri1,tri2,tri3;
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
      Eigen::Matrix<Scalar,3,2> T1;
      Eigen::Matrix<Scalar,3,2> T2;
      Eigen::Matrix<Scalar,1,2> ref;
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
      Eigen::Matrix<Scalar,3,2> T1;
      Eigen::Matrix<Scalar,3,2> T2;
      Eigen::Matrix<Scalar,1,2> ref;
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

void str_to_num(const std::vector<std::vector<std::string>> &uv_str, Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> &uv, int n_digits);
void num_to_str(Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> &pos, std::vector<std::vector<std::string>> &pos_str);

Eigen::VectorXi set_ri(
    const Eigen::MatrixXd& uv,
    const Eigen::MatrixXi& F
);

void simplify_triangulation(
  const Eigen::MatrixXd& V_i,
  const Eigen::MatrixXd& C,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F
);


// test for weakly-self-overlapping
template <typename Scalar>
bool weakly_self_overlapping(
  const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& P,
  const Eigen::VectorXi& R,
  Eigen::MatrixXi& F
);

// wrapper for weakly selfoverlappnig test - passing mpf as string
bool weakly_self_overlapping_str(
  const std::vector<std::vector<std::string>>& P_str,
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

template <typename Scalar>
void set_rotation_index(
   const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& uv,
   const Eigen::MatrixXi& F,
   Eigen::VectorXi& R,
   int offset = 0
);

template <typename Scalar>
Eigen::VectorXi set_all_ri(
    const Eigen::Matrix<Scalar, -1, -1> & uv,
    const Eigen::MatrixXi& F
);

Eigen::VectorXi set_all_ri_str(
    const std::vector<std::vector<std::string>> &uv_str,
    const Eigen::MatrixXi &F);

bool Shor_van_wyck(
  const Eigen::MatrixXd& P,
  const Eigen::VectorXi& R,
  const std::string flags,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  bool do_refine=true // false meaning no internal vertices
);

template <typename Scalar>
bool Shor_van_wyck_v2(
  const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& P,
  const Eigen::VectorXi& R,
  const Eigen::VectorXi& mark,
  const std::string flags,
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& V,
  Eigen::MatrixXi& F,
  Eigen::MatrixXi& Fn,
  bool do_refine
);

#endif