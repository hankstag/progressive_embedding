#include <Eigen/Dense>
#include "validity_check.h"
#include <igl/copyleft/cgal/orient2D.h>

bool is_face_flipped(const Eigen::Matrix<double,3,2>& T){
  double a[2] = {T(0,0),T(0,1)};
  double b[2] = {T(1,0),T(1,1)};
  double c[2] = {T(2,0),T(2,1)};
  return igl::copyleft::cgal::orient2D(a,b,c) <= 0;
}

void grad_to_eqtri(
  double target_area,
  Eigen::Matrix<double,2,3>& G_t
){
  double h = std::sqrt( target_area/sin(M_PI / 3.0));
  Eigen::Matrix<double,3,3> Tx;
  Tx<<0,0,0,h,0,0,h/2.,(std::sqrt(3)/2.)*h,0;
  Eigen::Matrix<double,3,3> gx;
  grad_operator(Tx,gx);
  G_t = gx.topRows(2);
}

void grad_operator(
    const Eigen::Matrix3d& T, // reference shape
    Eigen::Matrix3d& g
){
    Eigen::Matrix<double,1,3> v01 = T.row(1) - T.row(0);
    Eigen::Matrix<double,1,3> v12 = T.row(2) - T.row(1);
    Eigen::Matrix<double,1,3> v20 = T.row(0) - T.row(2);
    Eigen::Matrix<double, 1, 3> n = v01.cross(v12);
    // area of parallelogram is twice area of triangle
    // area of parallelogram is || v1 x v2 ||
    // This does correct l2 norm of rows, so that it contains #F list of twice
    // triangle areas
    double dblA = n.norm();
    Eigen::Matrix<double, 1, 3> u(0,0,1);
    // now normalize normals to get unit normals
    // u = n / dblA;

    Eigen::Matrix<double,1,3> eperp20,eperp01;
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

void flipped_elements(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::VectorXi& I
){
  I.setZero(F.rows());
  for(int i=0;i<F.rows();i++){
    double a[2] = {V(F(i,0),0),V(F(i,0),1)};
    double b[2] = {V(F(i,1),0),V(F(i,1),1)};
    double c[2] = {V(F(i,2),0),V(F(i,2),1)};
    if(igl::copyleft::cgal::orient2D(a,b,c)<=0) {
      I(i) = 1; // if cw or collinear, it's flipped
    }else
      I(i) = 0; // if ccw, it's not flipped
  }
}
