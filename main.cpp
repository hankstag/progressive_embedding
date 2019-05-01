#include <igl/opengl/glfw/Viewer.h>
#include <igl/matrix_to_list.h>
#include <igl/serialize.h>
#include "matchmaker.h"
#include "target_polygon.h"
#include "progressive_embedding.h"
#include "plot.h"
#include "loader.h"
#include "argh.h"
#include "slim/slim.h"

#include "validity_check.h"

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
  //test_local_smoothing();
  int loop, threshold, size;
  std::string model,uvfile;
  cmdl("-in") >> model;
  cmdl("-uv") >> uvfile;
  cmdl("-s", 0) >> size;
  cmdl("-t",20) >> threshold;
  
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
  
  igl::SLIMData sData;
  sData.slim_energy = igl::SLIMData::SYMMETRIC_DIRICHLET;
  igl::SLIMData::SLIM_ENERGY energy_type=igl::SLIMData::SYMMETRIC_DIRICHLET;
  //Eigen::SparseMatrix<double> Aeq;
  Eigen::VectorXd E;
  slim_precompute(V,F,uv,sData,igl::SLIMData::SYMMETRIC_DIRICHLET,ci,c,0,true,E,1.0);
  igl::opengl::glfw::Viewer vr;
  auto key_down = [&](
    igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier
  ){
    if (key == ' ') {
      slim_solve(sData,20,E);
      viewer.data().clear();
      viewer.data().set_mesh(V,F);
      viewer.data().set_uv(sData.V_o,F);
      viewer.data().show_texture = true;
    }
    if(key == '1'){
      slim_solve(sData,20,E);
      viewer.data().clear();
      viewer.data().set_mesh(sData.V_o,F);
      std::cout<<sData.V_o.row(ci(1))<<std::endl;
      viewer.data().show_texture = false;
    }
    return false;
  };
  vr.callback_key_down = key_down;
  //plot_mesh(vr,uv,F,{},Eigen::VectorXi());
  vr.launch();
}
