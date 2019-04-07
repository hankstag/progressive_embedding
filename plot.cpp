#include "plot.h"

void plot_polygon(
  igl::opengl::glfw::Viewer& viewer,
  const Eigen::MatrixXd& poly
){
  viewer.data().clear();
  for(int i=0;i<poly.rows();i++){
    int i_1 = (i+1) % poly.rows();
    viewer.data().add_edges(poly.row(i),poly.row(i_1),Eigen::RowVector3d(1,0,0));
  }
  viewer.launch();
}

void plot_mesh(
  igl::opengl::glfw::Viewer& viewer,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F
){
  viewer.data().clear();
  viewer.data().set_mesh(V,F);
  viewer.launch();
}