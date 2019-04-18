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
    Eigen::MatrixXd V,uv,polygon,c;
    Eigen::MatrixXi F;
    Eigen::VectorXi bd0,R,ci;
    
//#define LOADIN
#ifdef LOADIN
    load_in(model,V,F,polygon,bd0,R);
#else    
    Eigen::MatrixXd _polygon;
    Eigen::MatrixXi Fuv;
    Eigen::VectorXi bd1;
    std::pair<int,int> match;
    
    load_model_with_seam(model,V,F,polygon,bd0);
    load_matching_info(uvfile,match);
    load_model_with_seam(uvfile,uv,Fuv,_polygon,bd1);
    
    uv.conservativeResize(uv.rows(),2);
    int id0 = -1,id1 = -1;
    for(int i=0;i<bd0.rows();i++)
    {
        if(bd0(i) == match.first)
            id0 = i;
        if(bd1(i) == match.second)
            id1 = i;
    }
    int offset = (id1-id0+bd1.rows())%bd1.rows();
    std::cout<<"setting rotation index..."<<std::endl;
    set_rotation_index(uv,Fuv,R,offset);
    assert(bd0.rows()==bd1.rows());
#endif
    
    std::cout<<"rotation index:"<<std::endl;
    for (int i=0; i<R.rows(); i++)
    {
        if(R(i)!=0)
            std::cout<<i<<"-"<<bd0(i)<<":"<<R(i)<<std::endl;
    }
    
    //#define SHORTCUT
#ifdef SHORTCUT
    std::string serial_name = "local_save_pb";
    igl::deserialize(V,"V",serial_name);
    igl::deserialize(uv,"uv",serial_name);
    igl::deserialize(F,"F",serial_name);
    igl::deserialize(ci,"bi",serial_name);
    igl::deserialize(c,"b",serial_name);
    
    progressive_embedding(V,F,uv,ci,c,1e20);
#endif
    
    std::cout<<"match_maker para:"<<std::endl;
    std::cout<<"V:"<<V.rows()<<"*"<<V.cols()<<std::endl;
    std::cout<<"F:"<<F.rows()<<"*"<<F.cols()<<std::endl;
    std::cout<<"uv:"<<uv.rows()<<"*"<<uv.cols()<<std::endl;
    std::cout<<"c:"<<c.rows()<<"*"<<c.cols()<<std::endl;
    std::cout<<"ci:"<<ci.rows()<<"*"<<ci.cols()<<std::endl;
    std::cout<<"R:"<<R.rows()<<"*"<<R.cols()<<std::endl;
    std::cout<<"bd0:"<<bd0.rows()<<"*"<<bd0.cols()<<std::endl;
    std::cout<<"polygon:"<<polygon.rows()<<"*"<<polygon.cols()<<std::endl;
    for(int i=0;i<polygon.rows();i++)
        std::cout<<polygon.row(i)<<std::endl;
    
    match_maker(V,F,uv,c,ci,R,bd0,polygon);
#endif
    
    /*
     igl::SLIMData sData;
     sData.slim_energy = igl::SLIMData::SYMMETRIC_DIRICHLET;
     igl::SLIMData::SLIM_ENERGY energy_type=igl::SLIMData::SYMMETRIC_DIRICHLET;
     Eigen::SparseMatrix<double> Aeq;
     Aeq.resize(ci.rows()*2,uv.rows()*2);
     Eigen::VectorXd rb(ci.rows()*2);
     for(int i=0;i<ci.rows();i++){
     Aeq.coeffRef(2*i  ,ci(i))           = 1;
     Aeq.coeffRef(2*i+1,ci(i)+uv.rows()) = 1;
     rb(2*i) = uv(ci(i),0);
     rb(2*i+1) = uv(ci(i),1);
     }
     slim_precompute(V,F,uv,sData,igl::SLIMData::SYMMETRIC_DIRICHLET,ci,c,Aeq,rb,0);
     //slim_precompute(V,F,uv,sData,igl::SLIMData::SYMMETRIC_DIRICHLET,ci,c,0,true,E,Aeq,rb,1.0);
     */
    
    igl::opengl::glfw::Viewer vr;
    plot_mesh(vr,uv,F,{},Eigen::VectorXi());
    /*
     auto key_down = [&](
     igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier
     ){
     if (key == ' ') {
     slim_solve(sData,20);
     viewer.data().clear();
     viewer.data().set_mesh(V,F);
     viewer.data().set_uv(sData.V_o,F);
     viewer.core().align_camera_center(V);
     viewer.data().show_texture = true;
     }
     if(key == '1'){
     slim_solve(sData,20);
     viewer.data().clear();
     viewer.data().set_mesh(sData.V_o,F);
     viewer.core().align_camera_center(sData.V_o);
     viewer.data().show_texture = false;
     }
     return false;
     };
     
     vr.callback_key_down = key_down;
     */
    vr.launch();
    
}
