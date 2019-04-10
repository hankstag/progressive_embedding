#include "count_flipped_elements.h"
#include <igl/copyleft/cgal/orient2D.h>
#include <iostream>

void test_flip(
  const Eigen::MatrixXd& V, 
  const Eigen::MatrixXi& F, 
  std::vector<int>& fp
){
    int count = 0;
    fp.clear();
    for(int i=0;i<F.rows();i++){
        if(F.row(i).sum()==0) continue;
        double a[2] = {V(F(i,0),0),V(F(i,0),1)};
        double b[2] = {V(F(i,1),0),V(F(i,1),1)};
        double c[2] = {V(F(i,2),0),V(F(i,2),1)};
        if(igl::copyleft::cgal::orient2D(a,b,c)<=0) {
            count++;
            fp.push_back(i);
        }
    }
    std::cout<<"flipped :"<<count<<std::endl;
}

int test_flip(
  const Eigen::MatrixXd& V, 
  const Eigen::MatrixXi& F
){
    int count = 0;
    for(int i=0;i<F.rows();i++){
        if(F.row(i).sum()==0) continue;
        double a[2] = {V(F(i,0),0),V(F(i,0),1)};
        double b[2] = {V(F(i,1),0),V(F(i,1),1)};
        double c[2] = {V(F(i,2),0),V(F(i,2),1)};
        if(igl::copyleft::cgal::orient2D(a,b,c)<=0) {
            count++;
        }
    }
    return count;
}

void count_flipped_element(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::VectorXi& I
){
  I.setZero(F.rows());
  for(int i=0;i<F.rows();i++){
    if(F.row(i).minCoeff()<=0) continue;
    double a[2] = {V(F(i,0),0),V(F(i,0),1)};
    double b[2] = {V(F(i,1),0),V(F(i,1),1)};
    double c[2] = {V(F(i,2),0),V(F(i,2),1)};
    if(igl::copyleft::cgal::orient2D(a,b,c)<=0) {
      I(i) = 1;
    }else
      I(i) = 0;
  }
}