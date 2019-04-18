#include <igl/opengl/glfw/Viewer.h>
#include <igl/matrix_to_list.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include "progressive_embedding.h"
#include "validity_check.h"
#include "plot.h"
#include "loader.h"
#include "target_polygon.h"
#include "argh.h"
#include "convexation.h"
#include <igl/segment_segment_intersect.h>
#include <igl/copyleft/cgal/point_inside_polygon.h>
#include <igl/copyleft/cgal/segment_segment_intersect.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
  // auto strech = [](
  //   Eigen::RowVector2d p0,
  //   Eigen::RowVector2d v0,
  //   Eigen::RowVector2d v1,
  //   Eigen::RowVector2d p1,
  //   int N,
  //   double scale,
  //   Eigen::RowVector2d center,
  //   Eigen::MatrixXd& polyline
  // ){
    
  //   polyline.resize(N,2);
  //   for(int i=0;i<polyline.rows();i++){
      
  //   }
    
  //   Eigen::RowVector3d p0v0 = v0-p0;
  //   Eigen::RowVector3d p1v1 = v1-p1;
    
  //   double u,v;
  //   bool cross = igl::segment_segment_intersect(p0,p0v0,p1,p1v1,u,v,1e10);
  //   Eigen::RowVector3d joint;
  //   joint = p0 + p0v0*u;
  //   Eigen::RowVector3d v0v1;
  //   v0v1 << v1-v0;
  //   Eigen::RowVector3d nm;
  //   nm << 0,0,1;
  //   Eigen::RowVector3d dir = v0v1.cross(nm);
  //   double max_h = std::min((v0-joint).dot(dir),1.0);
  //   if(!cross) max_h = 1e-1;
  //   Eigen::RowVector3d mid_p = (v0+v1)*.5 + dir*max_h*.5;
  //   for(int i=0;i<polyline.rows();i++){
  //     double t = 1.0*i/N;
  //     polyline.row(i) << (1-t)*(1-t)*v0 + t*t*v1 + 2*t*(1-t)*mid_p;
  //   }
  //   // igl::opengl::glfw::Viewer viewer;
  //   // for(int i=0;i<polyline.rows();i++)
  //   //   viewer.data().add_points(polyline.row(i),Eigen::RowVector3d(1,0,0));
  //   // viewer.launch();
  //   std::cout<<" to sample: "<<N<<std::endl;
  // };
  
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

  int loop, threshold, size;
  std::string model,uvfile;
  cmdl("-in") >> model;
  cmdl("-s",0) >>size;
  
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
  std::map<int,int> mm;
  if(size!=0)
    background_mesh(size,bc,b,mm,uv,F);
  // std::vector<Object> O = {Object(uv,F,OTYPE::MESH)};
  // plots(O);
  
  igl::boundary_loop(F,b);
  igl::slice(uv,b,1,bc);  
  Eigen::MatrixXd H;
  
  //convexation(uv,F);
  Eigen::MatrixX2d bc2 = bc;
  Eigen::MatrixX2d bc2n;
  convexify(bc2,bc2n);
  igl::harmonic(F,b,bc2n,1,uv);
  
  Eigen::VectorXi I;
  flipped_elements(uv,F,I);
  std::cout<<"count flipped: "<<I.sum()<<std::endl;
  std::vector<Object> Os = {Object(uv,F,OTYPE::MESH)};
  plots(Os);
  
  bool x = progressive_embedding(V,F,uv,b,bc,1e20);
  #endif

  igl::opengl::glfw::Viewer vr;
  plot_mesh(vr,uv,F,{},Eigen::VectorXi());

}
