#include <igl/opengl/glfw/Viewer.h>
#include <igl/matrix_to_list.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include "progressive_embedding.h"
#include "count_flipped_elements.h"
#include "plot.h"
#include "loader.h"
#include "argh.h"

int main(int argc, char *argv[])
{
  auto cmdl = argh::parser(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);
  if(cmdl[{"-h","-help"}]){
      std::cout<<"Usage: ./matchmaker_bin -options"<<std::endl;
      std::cout<<"-in: input model name"<<std::endl;
      std::cout<<"-l: # local iteration in progressive fix"<<std::endl;
      std::cout<<"-t: exponent of energy threshold for edge collapsing"<<std::endl;
      std::cout<<"-p: whether use progressive embedding"<<std::endl;
      std::cout<<"-b: whether use the boundary info in obj file"<<std::endl;
      exit(0);
  }

  int loop, threshold;
  std::string model,uvfile;
  cmdl("-in") >> model;
  
  Eigen::MatrixXd V,polygon,uv;
  Eigen::MatrixXi F;
  Eigen::VectorXi T,R,bd;
  load_model(model,V,uv,F,polygon,R,T);
  
  #define FIXING
  #ifndef FIXING
  F.conservativeResize(F.rows()-1,3);
  igl::boundary_loop(F,bd);
  Eigen::MatrixXd circle;
  igl::map_vertices_to_circle(V,bd,circle);
  Eigen::MatrixXd H;
  igl::harmonic(F,bd,circle,1,H);
  progressive_embedding(V,F,H,bd,circle,1e20);
  #else 
  Eigen::VectorXi b;
  Eigen::MatrixXd bc;
  igl::boundary_loop(F,b);
  igl::slice(uv,b,1,bc);
  Eigen::VectorXi I;
  count_flipped_element(uv,F,I);
  std::cout<<"count flipped: "<<I.sum()<<std::endl;
  bool x = progressive_embedding(V,F,uv,b,bc,1e20);
  #endif

  igl::opengl::glfw::Viewer vr;
  plot_mesh(vr,uv,F,{});

}
