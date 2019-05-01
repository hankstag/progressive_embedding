#include <igl/boundary_loop.h>
#include <igl/writeOBJ.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include "progressive_embedding.h"
#include "validity_check.h"
#include "plot.h"
#include "loader.h"
#include "argh.h"
#include <igl/Timer.h>

int main(int argc, char *argv[])
{
  auto cmdl = argh::parser(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);
  if(cmdl[{"-h","-help"}]){
      std::cout<<"Usage: ./ -options"<<std::endl;
      std::cout<<"-in: input model name"<<std::endl;
      exit(0);
  }

  int loop, threshold;
  std::string model,outfile;
  bool verbose;
  cmdl("-in") >> model;
  cmdl("-t", 20) >> threshold;
  cmdl("-vb",false) >> verbose;
  cmdl("-o", model+"_out.obj") >> outfile;
  double eps = std::pow(10,threshold);
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
  igl::Timer tm;
  tm.start();
  double t0 = tm.getElapsedTime();
  progressive_embedding(V,F,H,bd,circle,eps,true);
  double t1 = tm.getElapsedTime();
  uv = H;
  
  Eigen::VectorXi I;
  flipped_elements(uv,F,I);
  std::cout<<"# flips: "<<I.sum()<<std::endl;
  std::cout<<"Running time: ";
  printf("%.3f (s)\n",t1-t0);
  
  igl::opengl::glfw::Viewer vr;
  plot_mesh(vr,uv,F,{},Eigen::VectorXi());
  vr.launch();
  
  Eigen::MatrixXd CN;
  Eigen::MatrixXi FN;
  igl::writeOBJ(outfile,V,F,CN,FN,uv,F);

}
