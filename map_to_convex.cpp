#include <igl/opengl/glfw/Viewer.h>
#include <igl/matrix_to_list.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include "progressive_embedding.h"
#include "validity_check.h"
#include "plot.h"
#include "loader.h"
#include "argh.h"
#include <igl/segment_segment_intersect.h>
#include <igl/copyleft/cgal/segment_segment_intersect.h>
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
  
void convexation(
  Eigen::MatrixXd& uv,
  const Eigen::MatrixXi& F
){
  Eigen::VectorXi bd;
  igl::boundary_loop(F,bd);
  Eigen::MatrixXd polygon(bd.rows(),3);
  auto rotation_matrix = [](
    double angle,
    Eigen::Matrix3d& m
  ){
    m<< std::cos(angle),std::sin(angle),0,
       -std::sin(angle),std::cos(angle),0,
       0,0,1;
  };
  
  auto strech = [&](
    Eigen::RowVector2d p0_2d,
    Eigen::RowVector2d v0_2d,
    Eigen::RowVector2d v1_2d,
    Eigen::RowVector2d p1_2d,
    int N,
    double scale,
    Eigen::RowVector2d center,
    Eigen::MatrixXd& polyline
  ){
    polyline.resize(N,3);
    Eigen::RowVector3d p0;
    p0 << p0_2d,0;
    Eigen::RowVector3d v0;
    v0 << v0_2d,0;
    Eigen::RowVector3d v1;
    v1 << v1_2d,0;
    Eigen::RowVector3d p1;
    p1 << p1_2d,0;
    Eigen::RowVector3d p0v0 = v0-p0;
    Eigen::RowVector3d p1v1 = v1-p1;
    double u,v;
    bool cross = igl::segment_segment_intersect(p0,p0v0,p1,p1v1,u,v,0);
    Eigen::RowVector3d joint;
    joint = p0 + p0v0*u;
    Eigen::RowVector3d jt_v0 = v0-joint;
    Eigen::RowVector3d jt_v1 = v1-joint;
    double angle_top = std::max(igl::PI*2/3,std::acos(jt_v0.dot(jt_v1)/(jt_v0.norm()*jt_v1.norm())));
    double angle_bottom = (igl::PI - angle_top)/2;
    double delta_angle = (igl::PI-angle_top) / (N-2);
    std::cout<<"angle top: "<<angle_top/igl::PI * 180<<std::endl;
    std::cout<<"angle bottom: "<<angle_bottom/igl::PI * 180<<std::endl;
    std::cout<<"angle delta: "<<delta_angle/igl::PI * 180<<std::endl;
    Eigen::Matrix3d R0, Rd;
    rotation_matrix(angle_bottom,R0);
    Eigen::RowVector3d v0v1 = v1-v0;
    Eigen::RowVector3d h0 = (R0*v0v1.transpose()).transpose().eval(); // reference vector
    
    double l = v0v1.norm();
    rotation_matrix(-delta_angle,Rd);
    
    Eigen::RowVector3d hp = h0;
    double t0 = v0v1.norm()*(l/(N-1))/h0.dot(v0v1);
    Eigen::RowVector3d vp = v0+t0*h0;
    polyline.row(0) << v0;
    if(N > 2)
      polyline.row(1) << vp;
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(uv,F);
    viewer.data().add_points(p0,Eigen::RowVector3d(1,0,0));
    viewer.data().add_points(p1,Eigen::RowVector3d(1,0,0));
    viewer.data().add_points(v0,Eigen::RowVector3d(1,1,0));
    viewer.data().add_points(v1,Eigen::RowVector3d(1,1,0));
    viewer.data().add_points(v0+t0*h0,Eigen::RowVector3d(0,1,0));
    std::cout<<"N = "<<N<<std::endl;
    for(int i=2;i<N-1;i++){
      Eigen::RowVector3d hi = (Rd * hp.transpose()).transpose().eval();
      double len = v0v1.norm()*(l/(N-1))/hi.dot(v0v1);
      Eigen::RowVector3d vi = vp + hi * len;
      polyline.row(i) << vi;
      viewer.data().add_points(vi,Eigen::RowVector3d(0,0,1));
      hp = hi;
      vp = vi;
    }
    polyline.bottomRows(1) << v1;
    
    // for(int i=1;i<N;i++){
    //   Eigen::RowVector3d h_i;
    //   Eigen::Matrix3d rotation;
    //   rotation<<std::cos(delta_angle),-std::sin(delta_angle),0,
    //             std::sin(delta_angle), std::cos(delta_angle),0,
    //             0,0,1;
    //   h_i = (rotation * h.transpose()).transpose();
    //viewer.data().add_points(h0+v0,Eigen::RowVector3d(0,0,1));
    //   h = h_i; 
    // }
    // viewer.data().add_points(joint,Eigen::RowVector3d(0,0,1));
    // // for(int i=1;i<polyline.rows()-1;i++)
    // //   viewer.data().add_points(polyline.row(i),Eigen::RowVector3d(1,0,0));
    viewer.launch();
  };
  auto is_corner = [](
    const Eigen::RowVector2d& vi,
    const Eigen::RowVector2d& vx,
    const Eigen::RowVector2d& vy
  ){    
    auto bc = vx-vi;
    auto ba = vy-vi;
    double angle = std::acos(bc.dot(ba)/(bc.norm()*ba.norm()));
    if(std::isnan(angle) || std::abs(angle-igl::PI)<1e-7) // straight
      return false;
    else {
      std::cout<<"corner "<<angle-igl::PI<<std::endl;
      return true;
    }
  };
  
  bool met_corner = false;
  int c0 = 0,c1;
  int begin = 0;
  // find the first corner
  Eigen::RowVector2d center;
  center.setZero();
  int count_corner = 0;
  for(int i=0;i<bd.rows();i++){
    int n = bd.rows();
    if(is_corner(uv.row(bd(i)),uv.row(bd((i-1+n)%n)),uv.row(bd((i+1)%n)))){
      begin = i;
      count_corner++;
      center += uv.row(bd(i));
    }
  }
  center /= count_corner;
  for(int j=0;j<bd.rows()+1;j++){
    int n = bd.rows();
    int i = (j+begin) % n;
    bool cc = is_corner(uv.row(bd(i%n)),uv.row(bd((i-1+n)%n)),uv.row(bd((i+1)%n)));
    if(cc && met_corner){
      std::cout<<" meet corner at "<<i%n<<std::endl;
      c1 = i%n;
      Eigen::MatrixXd polyline;
      // sample between c0 and c1
      Eigen::RowVector2d a,b,c,d;
      a = uv.row(bd((c0-1+n)%n));
      b = uv.row(bd(c0));
      c = uv.row(bd(c1));
      d = uv.row(bd((c1+1)%n));
      strech(a,b,c,d,(c1-c0+1+n)%n,0.01,center,polyline);
      for(int j=0;j<polyline.rows();j++)
        polygon.row((j+c0)%n) << polyline.row(j);
      c0 = i;
    }else if(cc){
      c0 = i;
      met_corner = true;
    }
  }
  igl::harmonic(F,bd,polygon,1,uv);
  Eigen::VectorXi T;
  flipped_elements(uv,F,T);
  std::cout<<"flipped "<<T.sum()<<std::endl;
  std::vector<Object> Os = {
    Object(uv,F,OTYPE::MESH)
  };
  plots(Os);
  // reverse scaling
  for(int i=0;i<uv.rows();i++){
    // find which sector does vertex vi belong to
    for(int j=0;j<bd.rows();j++){
      Eigen::MatrixXd T(3,2);
      T<< uv.row(bd(j)), uv.row(bd((j+1)%bd.rows())), center;
      if(igl::copyleft::cgal::point_inside_polygon(uv.row(i),T)){
        double u1,u2,t1,t2;
        bool cross1 = igl::segment_segment_intersect(p,r,q,s,u,t);
        bool cross2 = igl::segment_segment_intersect(p,r,q,s,u,t);
        double scaling = u2 / u1;
        uv.row(i) << scaling * (uv.row(i) - center) + center;
        break;
      }
    }
  }
}

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
  Eigen::MatrixXd H;
  igl::harmonic(F,b,bc,1,H);
  Eigen::VectorXi I;
  flipped_elements(H,F,I);
  convexation(uv,F);
  std::cout<<"count flipped: "<<I.sum()<<std::endl;
  bool x = progressive_embedding(V,F,uv,b,bc,1e20);
  #endif

  igl::opengl::glfw::Viewer vr;
  plot_mesh(vr,uv,F,{},Eigen::VectorXi());

}
