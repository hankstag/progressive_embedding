#include <igl/opengl/glfw/Viewer.h>
#include <igl/matrix_to_list.h>
#include <igl/serialize.h>
#include "matchmaker.h"
#include "target_polygon.h"
#include "progressive_embedding.h"
#include "plot.h"
#include "loader.h"
#include "argh.h"
#include "cut_mesh/edge_flaps.h"
#include "local_smooth/local_smooth.h"
#include <igl/remove_unreferenced.h>
#include <igl/copyleft/cgal/orient2D.h>

// void test_local_smoothing(){
//   Eigen::MatrixXi F(13,3);
//   F<<0,1,2,
//      3,2,1,
//      2,4,5,
//      6,7,2,
//      8,6,9,
//      10,9,2,
//      2,11,0,
//      10,2,12,
//      12,2,5,
//      2,3,4,
//      2,7,11,
//      9,6,2,
//      8,7,6;
//   Eigen::MatrixXd V(13,2);
//   V<<                -39,            -9.8125,
//                 -39,              -9.75,
// -39.000000000000007,-9.9365594558871901,
//                 -39,            -9.6875,
//  -39.00001456174526,-9.9732687460454432,
//                 -39,-10.421052631578947,
// -39.000000000000007,-9.9365594558871901,
//                 -39,            -9.9375,
//                 -39,                -10,
//                 -39,-10.105263157894736,
//                 -39,-10.210526315789474,
//                 -39,             -9.875,
//                 -39,-10.315789473684211;
//   double target_area = 0.432473;
//   Eigen::VectorXi B(V.rows());
//   Eigen::MatrixXd V0 = V;
//   local_smoothing(V,F,B,V,10,1e10);
//   std::cout<<"diff: "<<(V-V0).norm()<<std::endl;
//   std::cout<<V<<std::endl;

// }

void quadratic_func(
    const Eigen::Vector2d& s,
    const Eigen::Vector2d& t,
    const Eigen::Vector2d& pos,
    Eigen::Vector2d& result,
    double eps
){
    Eigen::Vector2d direction = s-t;
    Eigen::MatrixXd rotation(2,2);
    rotation<<0,-1,
              1, 0;
    Eigen::Vector2d normal;
    normal = (rotation * direction).normalized();
    double length = (s-t).norm();
    double ratio = (s-pos).norm()/length;
    double a = -4*eps;
    double b = 4*eps;
    double offset = a*ratio*ratio+b*ratio;
    normal = normal.array()*offset;
    result = normal+pos;
}

void perturb_curve(
    const Eigen::MatrixXd& P,
    Eigen::MatrixXd& polygon_perturbed
){
    //for all consecutive colinear edges, map them to a bounding polynomial
    int s_id = 0;
    int t_id = 0;
    bool colinear = false;
    int count = 0;
    polygon_perturbed.resize(P.rows(),2);
    Eigen::VectorXi corners(P.rows());
    corners.setZero();
    // [should make sure vertex 0 is a corner point]
    for(int i=0;i<P.rows()+1;i++){
        int prev = (i-1+P.rows())%P.rows();
        int next = (i+1)%P.rows();
        int curr = i % P.rows();
        double a[2] = {P(prev,0),P(prev,1)};
        double b[2] = {P(curr,0),P(curr,1)};
        double c[2] = {P(next,0),P(next,1)};
        short r = igl::copyleft::cgal::orient2D(a,b,c);
        if(r == 0){ // colinear
            t_id = curr;
            colinear = true;
        }else{
            // if its end point of colinear line
            if(colinear){
                t_id = curr;
                Eigen::Vector2d s,t,result;
                s<< P(s_id,0),P(s_id,1);
                t<< P(t_id,0),P(t_id,1);
                int k=(s_id+1)%P.rows();
                while(k!=t_id){
                    Eigen::Vector2d pos;
                    pos<<P(k,0),P(k,1);
                    double eps;
                    quadratic_func(s,t,pos,result,eps);
                    //std::cout<<eps<<std::endl;
                    polygon_perturbed.row(count++)<<result(0),result(1);
                    k = (k+1)%P.rows();
                }
                colinear = false;
            }
            s_id = curr;
            if(i!=P.rows())
                polygon_perturbed.row(count++)<<P(i,0),P(i,1);
        }
    }
    polygon_perturbed.conservativeResize(count,2);
    std::vector<std::pair<int,int>> edges;
    for(int i=0;i<polygon_perturbed.rows();i++){
        edges.push_back(std::make_pair(i,(i+1)%polygon_perturbed.rows()));
    }
    std::cout<<"polygon perturbed size "<<polygon_perturbed.rows()<<std::endl;
    // std::cout<<std::setprecision(17)<<polygon_perturbed<<std::endl;
    //Eigen::MatrixXi F_;
    //display(polygon_perturbed,F_,edges,{},{},{});
    
}


void random_internal_vertices(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::VectorXi& pi
){
  // [ generate 3 random points ]
  std::vector<int> r;
  Eigen::VectorXi bd;
  std::vector<int> bd_list;
  igl::boundary_loop(F,bd);
  igl::matrix_to_list(bd,bd_list);
  int seed = 10;
  for(int i=0;i<3;i++){
    int t = -1;
    do{
        srand(seed);
        int z = rand();
        t = z%V.rows();
        seed=seed*4+i;
    }while(std::find(bd_list.begin(),bd_list.end(),t)!=bd_list.end() || 
           std::find(r.begin(),r.end(),t)!=r.end());
    r.push_back(t);
  }
  igl::list_to_matrix(r,pi);
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
  cmdl("-uv") >> uvfile;
  // cmdl("-l",20) >> loop;
  // cmdl("-t",20) >> threshold;
  #define QUAD
  #ifndef QUAD
  Eigen::MatrixXd V,polygon,uv;
  Eigen::MatrixXi F;
  Eigen::VectorXi T,R;
  load_model(model,V,uv,F,polygon,R,T);
  
  Eigen::MatrixXd c(3,2);
  c<<0,0,1,0,0,1;
  Eigen::VectorXi ci;
  random_internal_vertices(V,F,ci);
  
  std::vector<Eigen::MatrixXd> polys;
  target_polygon(V,F,c,ci,polys);
  
  R.setZero(polys[0].rows());
  match_maker(V,F,uv,c,ci,R,T,polys[0]);
  #else
  // load model and uv
  Eigen::MatrixXd V,uv,polygon;
  Eigen::MatrixXi F,Fuv;
  Eigen::VectorXi bd0,bd1;
  load_model_with_seam(model,V,F,polygon,bd0);
  
  std::pair<int,int> match;
  load_matching_info(uvfile,match);
  Eigen::MatrixXd _polygon;
  Eigen::VectorXi R;
  load_model_with_seam(uvfile,uv,Fuv,_polygon,bd1);
  uv.conservativeResize(uv.rows(),2);
  int id0 = -1,id1 = -1;
  for(int i=0;i<bd0.rows();i++){
      if(bd0(i) == match.first)
          id0 = i;
      if(bd1(i) == match.second)
          id1 = i;
  }
  int offset = (id1-id0+bd1.rows())%bd1.rows();
  std::cout<<"setting rotation index..."<<std::endl;
  set_rotation_index(uv,Fuv,R,offset);
  assert(bd0.rows()==bd1.rows());

  igl::opengl::glfw::Viewer viewer1;
  //plot_mesh(viewer1,V,F,{},true);

  Eigen::VectorXi ci;
  Eigen::MatrixXd c;
  
  #define SHORTCUT
  #ifdef SHORTCUT
  std::string serial_name = "carter_save_pb";
  igl::deserialize(V,"V",serial_name);
  igl::deserialize(uv,"uv",serial_name);
  igl::deserialize(F,"F",serial_name);
  igl::deserialize(ci,"bi",serial_name);
  igl::deserialize(c,"b",serial_name);
  // Eigen::VectorXi EMAP,EE;
  // Eigen::MatrixXi E,EF,EI;
  // Eigen::MatrixXi dEF,dEI,allE;
  // igl::edge_flaps(F,E,allE,EMAP,EF,EI,dEF,dEI,EE);
  // Eigen::VectorXi B;
  // B.setZero(V.rows());
  // for(int i=0;i<ci.rows();i++)
  //   B(ci(i)) = 1;
  // std::set<int> lengths;
  // for(int i=0;i<F.rows();i++){
  //   Eigen::MatrixXi local_F,NF;
  //   Eigen::MatrixXd local_V,NV,Nuv;
  //   expand_to_boundary(V,F,B,dEF,i,local_F);
  //   Eigen::VectorXi II;
  //   igl::remove_unreferenced(V,local_F,NV,NF,II);
  //   igl::remove_unreferenced(uv,local_F,Nuv,NF,II);
  //   if(lengths.find(NF.rows())==lengths.end()){
  //     Eigen::MatrixXd CN;
  //     Eigen::MatrixXi FN;
  //     igl::writeOBJ("patch"+std::to_string(i)+".obj",NV,NF,CN,FN,Nuv,NF);
  //     Eigen::VectorXi FLIP;
  //     count_flipped_element(Nuv,NF,FLIP);
  //     if(FLIP.sum()!=0){
  //       std::cout<<"patch "<<i<<" has #flipped: "<<FLIP.sum()<<std::endl;
  //       igl::opengl::glfw::Viewer viewer;
  //       viewer.data().set_mesh(NV,NF);
  //       Eigen::VectorXi bb;
  //       igl::boundary_loop(NF,bb);
  //       Eigen::MatrixXd poly,npoly;
  //       igl::slice(Nuv,bb,1,poly);
  //       perturb_curve(poly,npoly);
  //       igl::harmonic(NF ,bb,npoly,1,Nuv);
  //       count_flipped_element(Nuv,NF,FLIP);
  //       std::vector<Object> Os = {
  //         Object(V,F,OTYPE::MESH),
  //         Object(npoly,F,OTYPE::POLYGON)
  //       };
  //       plots(Os);
  //       // for(int i=0;i<bb.rows();i++)
  //       //   viewer.data().add_points(Nuv.row(bb(i)),Eigen::RowVector3d(1,0,0));
  //       for(int i=0;i<FLIP.rows();i++){
  //         if(FLIP(i)){
  //           for(int k=0;k<3;k++)
  //             viewer.data().add_points(NV.row(NF(i,k)),Eigen::RowVector3d(0,1,0));
  //           }
  //         }
  //         viewer.launch();
  //       }
  //       lengths.insert(NF.rows());
  //   }
  // }
  progressive_embedding(V,F,uv,ci,c,1e100);
  #endif

  match_maker(V,F,uv,c,ci,R,bd0,polygon);
  #endif
  igl::opengl::glfw::Viewer vr;
  plot_mesh(vr,uv,F,{});

}
