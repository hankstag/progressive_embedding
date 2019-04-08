#ifndef PLOT
#define PLOT

#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>

void plot_polygon(
  igl::opengl::glfw::Viewer& viewer,
  const Eigen::VectorXi& H,
  const Eigen::MatrixXd& poly
);

void plot_mesh(
  igl::opengl::glfw::Viewer& viewer,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F
);

#endif