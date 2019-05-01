#include <igl/boundary_loop.h>
#include <igl/Timer.h>
#include "progressive_embedding.h"
#include "validity_check.h"
#include "plot.h"
#include "loader.h"
#include "target_polygon.h"
#include "argh.h"
  
int main(int argc, char *argv[])
{
  auto cmdl = argh::parser(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);
  if(cmdl[{"-h","-help"}]){
      std::cout<<"Usage: ./matchmaker_bin -options"<<std::endl;
      std::cout<<"-in: input model name"<<std::endl;
      exit(0);
  }

  int threshold;
  std::string model,outfile;
  
  cmdl("-in") >> model;
  cmdl("-t",20) >> threshold;
  cmdl("-o", model+"_out.obj") >> outfile;
  double eps = std::pow(10,threshold);
  
  Eigen::MatrixXd V,polygon,uv;
  Eigen::MatrixXi F;
  Eigen::VectorXi T,R;
  load_model(model,V,uv,F,polygon,R,T);
  
  Eigen::VectorXi b;
  Eigen::MatrixXd bc;
  igl::boundary_loop(F,b);

  igl::Timer tm;
  tm.start();
  Eigen::VectorXi I;
  flipped_elements(uv,F,I);
  std::cout<<"#flips in initial map: "<<I.sum()<<std::endl;  
  
  double t0 = tm.getElapsedTime();
  bool x = progressive_embedding(V,F,uv,b,bc,eps,true);
  double t1 = tm.getElapsedTime();
  flipped_elements(uv,F,I);
  std::cout<<"\n#flips after PE: "<<I.sum()<<std::endl;
  std::cout<<"Running time: ";
  printf("%.3f (s)\n",t1-t0);
  
  igl::opengl::glfw::Viewer vr;
  plot_mesh(vr,uv,F,{},Eigen::VectorXi());
  vr.launch();
  
  Eigen::MatrixXd CN;
  Eigen::MatrixXi FN;
  igl::writeOBJ(outfile,V,F,CN,FN,uv,F);

}
