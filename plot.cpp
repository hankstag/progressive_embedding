#include "plot.h"
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

void plot_polygon(
  igl::opengl::glfw::Viewer& viewer,
  const Eigen::VectorXi& H,
  const Eigen::MatrixXd& poly
){
  viewer.data().clear();
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);
  viewer.core().align_camera_center(poly);
  for(int i=0;i<poly.rows();i++){
    int i_1 = (i+1) % poly.rows();
    viewer.data().add_label(poly.row(i),std::to_string(i));
    viewer.data().add_edges(poly.row(i),poly.row(i_1),Eigen::RowVector3d(1,0,0));
  }
  for(int i=0;i<H.rows();i++){
    if(H(i))
      viewer.data().add_points(poly.row(i),Eigen::RowVector3d(0,0,0));
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