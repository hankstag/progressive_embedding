#include "plot.h"
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/boundary_loop.h>

void plot_polygon(
  igl::opengl::glfw::Viewer& viewer,
  const Eigen::VectorXi& H,
  const Eigen::MatrixXd& poly
){
  viewer.data().clear();
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
}

void plot_mesh(
  igl::opengl::glfw::Viewer& viewer,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const std::vector<int>& id,
  const Eigen::VectorXi& L, //hight faces
  bool show_boundary
){
  viewer.data().clear();
  viewer.data().set_mesh(V,F);
  viewer.core().align_camera_center(V,F);
  Eigen::MatrixXd C(F.rows(),3);
  C.setConstant(1);
  if(show_boundary){
    Eigen::VectorXi bd;
    igl::boundary_loop(F,bd);
    for(int i=0;i<bd.rows();i++){
      int i_1 = (i+1) % bd.rows();
      viewer.data().add_edges(V.row(bd(i)),V.row(bd(i_1)),Eigen::RowVector3d(0,1,0));
    }
  }
  for(int i: id){
    viewer.data().add_points(V.row(i),Eigen::RowVector3d(1,0,0));
  }
  for(int i=0;i<L.rows();i++){
    if(L(i)!=0)
      C.row(i) << 1,0,0;
  }
  if(L.sum()>0)
    viewer.data().set_colors(C);
}

// For several mesh/polygon at the same time
// using number key to switch between them
// usage example:  
// std::vector<Object> Os = {Object(V,F,OTYPE::MESH),
//                           Object(polygon,F,OTYPE::POLYGON),
//                           Object(uv,F,OTYPE::MESH)};
// plots(Os);

void plots(
  std::vector<Object>& data 
){
  igl::opengl::glfw::Viewer vr;
  int item_num = data.size();
  auto key_down = [&](
    igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier
  ){
    int id = (key - '0');
    std::cout<<"choose "<<id<<std::endl;
    if(id >= data.size() || id < 0) return false;
    Object obj = data[id];
    
    Eigen::MatrixXd V=data[id].V;
    Eigen::MatrixXi F=data[id].F;
    Eigen::VectorXi H;
    H.setZero(V.rows());
    switch(obj.get_type()){
      case OTYPE::POLYGON: plot_polygon(viewer,H,V);break;
      case OTYPE::MESH: plot_mesh(viewer,V,F,{},Eigen::VectorXi()); break;
      default: std::cout<<"other case"<<std::endl;
    }
    return false;
  };
  vr.callback_key_down = key_down;
  vr.launch();
}