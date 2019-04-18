#include "convexation.h"
#include "decompose_polygon.h"
#include "plot.h"

std::vector<Eigen::MatrixXd> polylines;
double EPS=1e-11;
Eigen::RowVector2d bc;

void corners_of_polygon(
  const Eigen::Matrix<double,Eigen::Dynamic,2>& P,
  std::vector<int>& I
){
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
    else
      return true;
  };
  for(int i=0;i<P.rows();i++){
    int curr = i;
    int prev = (i-1+P.rows()) % P.rows();
    int next = (i+1) % P.rows();
    if(is_corner(P.row(curr),P.row(prev),P.row(next)))
      I.push_back(i);
  }
}

 void extend_mid_vertex(
    const int l,
    const int r,
    Eigen::MatrixX2d& pl,
    double l0,
    Eigen::VectorXi& A
  ){
    if(l >= r-1)
      return;
    int m = (l+r)/2;
    double li = (pl.row(l)-pl.row(r)).norm();
    double eps = EPS * li/l0;
    auto sample_between = [](
      const Eigen::RowVector2d& p1,
      const Eigen::RowVector2d& p2,
      Eigen::MatrixX2d& L
    ){
      //int n = (t-s+1+pl.rows()) % pl.rows();
      int n = L.rows();
      for(int i=0;i<n;i++)
        L.row(i) << (i*p2 + (n-1-i)*p1)/(n-1);
    };
    Eigen::MatrixX2d L1,L2;
    L1.resize(m-l+1,2);
    L2.resize(r-m+1,2);
    A(m) = 1;
    Eigen::MatrixX2d P0;
    P0.resize(A.sum(),2);
    do{
      auto pt = pl.row(m) + eps * (pl.row(m)-bc).normalized();
      sample_between(pl.row(l),pt,L1);
      sample_between(pt,pl.row(r),L2);
      eps /= 2.0;
      int c = 0;
      for(int i=0;i<A.rows();i++){
        if(A(i)){
          if(i == m){
            P0.row(c++) << pt;
          }else{
            P0.row(c++) << pl.row(i);
          }
        }
      }

    }while(!is_convex(P0));
    pl.block(l,0,m-l+1,2) = L1;
    pl.block(m,0,r-m+1,2) = L2;
    std::cout<<"l = "<<l<<std::endl;
    std::cout<<"m = "<<m<<std::endl;
    std::cout<<"r = "<<r<<std::endl;
    std::cout<<"is convex? "<<is_convex(P0)<<std::endl;
    extend_mid_vertex(l,m,pl,l0,A);
    extend_mid_vertex(m,r,pl,l0,A);

  }

void convexify(
  const Eigen::Matrix<double,Eigen::Dynamic,2>& P,
  Eigen::Matrix<double,Eigen::Dynamic,2>& Pn
){
  
  Pn.resizeLike(P);
  
  // collect corners    
  std::vector<int> I;
  corners_of_polygon(P,I);
  
  bc.setZero();
  for(int id: I)
    bc += P.row(id);
  bc /= I.size();
  
  for(int i=0;i<I.size();i++){
    int i_1 = (i+1) % I.size();
    int n = (I[i_1]-I[i]+1+P.rows()) % P.rows();
    assert(n >= 2 && "edge has at least two endpoints");
    Eigen::MatrixX2d pl;
    pl.resize(n,2);
    std::cout<<"start at "<<I[i]<<std::endl;
    for(int j=0;j<pl.rows();j++){
      pl.row(j) << P.row((I[i]+j)%P.rows());
    }
    igl::opengl::glfw::Viewer viewer;
    Eigen::VectorXi E;
    E.setConstant(pl.rows(),1);
    plot_polygon(viewer,E,pl);
    viewer.launch();
    Eigen::VectorXi A; // active vertices
    A.setZero(pl.rows());
    A(0) = 1;
    A(pl.rows()-1) = 1;
    double len = (P.row(I[i_1]) - P.row(I[i])).norm();
    EPS = len * 1e-1;
    extend_mid_vertex(0,n-1,pl,len,A);
 
    for(int j=0;j<pl.rows();j++){
      Pn.row((I[i]+j)%Pn.rows()) << pl.row(j);
    }

  }
  std::cout<<"check: "<<is_convex(Pn)<<std::endl;
  igl::opengl::glfw::Viewer viewer;
  Eigen::VectorXi E;
  E.setConstant(Pn.rows(),1);
  plot_polygon(viewer,E,Pn);
  viewer.launch();
  
}

void reverse(
  Eigen::MatrixXd& uv
){
  
}