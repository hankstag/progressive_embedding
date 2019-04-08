#include <igl/opengl/glfw/Viewer.h>
#include <igl/matrix_to_list.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include "progressive_embedding.h"
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
  
  F.conservativeResize(F.rows()-1,3);
  igl::boundary_loop(F,bd);
  Eigen::MatrixXd circle;
  igl::map_vertices_to_circle(V,bd,circle);
  Eigen::MatrixXd H;
  igl::harmonic(F,bd,circle,1,H);
  Eigen::VectorXi ci;
  progressive_fix(ci,bd,circle,V,F,H);

  igl::opengl::glfw::Viewer vr;
  plot_mesh(vr,H,F,{});

}
