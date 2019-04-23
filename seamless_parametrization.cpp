#include <igl/opengl/glfw/Viewer.h>
#include "matchmaker.h"
#include "plot.h"
#include "loader.h"
#include "argh.h"
#include "slim/slim.h"
#include "process_ee.h"
#include "validity_check.h"

#include <time.h>

int main(int argc, char *argv[])
{
    clock_t now = clock();
    
    auto cmdl = argh::parser(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);
    
    if(cmdl[{"-h","-help"}])
    {
        std::cout<<"Usage: ./seamless_bin -options"<<std::endl;
        std::cout<<"-in: input model name"<<std::endl;
        exit(0);
    }
    
    std::string model;
    cmdl("-in") >> model;
    
    // init map
    Eigen::MatrixXd V,polygon,uv,c;
    Eigen::MatrixXi F,EE,cut;
    Eigen::VectorXi R,bd,ci;
    
    load_in(model,V,F,R,bd,polygon,EE);
    process_ee(EE,F,cut);
    match_maker(V,F,uv,c,ci,R,bd,polygon);

    // tmp
    igl::opengl::glfw::Viewer vr;
    plot_mesh(vr,uv,F,{},Eigen::VectorXi());
    vr.launch();
 /*
    // seamless optimization
    Eigen::SparseMatrix<double> Aeq;
    Eigen::VectorXd rb,E;
    igl::SLIMData sData;
    sData.slim_energy = igl::SLIMData::SYMMETRIC_DIRICHLET;
    igl::SLIMData::SLIM_ENERGY energy_type=igl::SLIMData::SYMMETRIC_DIRICHLET;
    
    buildAeq(cut,V,F,F,uv,Aeq);
    rb.setZero(Aeq.rows());
    slim_precompute(V,F,uv,sData,igl::SLIMData::SYMMETRIC_DIRICHLET,ci,c,Aeq,rb,0);
    //slim_solve(sData,500);
 
    igl::opengl::glfw::Viewer vr;
    double scale = 1.0;

    auto key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
    {
        
        if (key == '1')
        {
            slim_solve(sData,20);
            viewer.data().clear();
            viewer.data().set_mesh(sData.V_o,F);
            viewer.core().align_camera_center(sData.V_o);
            viewer.data().show_texture = false;
        }
        if (key == '2')
        {
            slim_solve(sData,20);
            viewer.data().clear();
            viewer.data().set_mesh(V,F);
            viewer.data().set_uv(sData.V_o,F);
            viewer.core().align_camera_center(V);
            viewer.data().show_texture = true;
        }
        if(key == 'a' || key == 'A')
        {
            viewer.data().clear();
            viewer.data().set_mesh(sData.V_o,F);
            viewer.core().align_camera_center(sData.V_o);
            viewer.data().show_texture = false;
        }
        if(key == 's' || key == 'S')
        {
            scale = 1.0;
            viewer.data().clear();
            viewer.data().set_mesh(V,F);
            viewer.data().set_uv(sData.V_o,F);
            viewer.core().align_camera_center(V);
            viewer.data().show_texture = true;
        }
        if(key == ',')
        {
            scale *= 2.0;
            std::cout<<"scale:"<<scale<<std::endl;
            viewer.data().clear();
            viewer.data().set_mesh(V,F);
            viewer.data().set_uv(sData.V_o*scale,F);
            viewer.core().align_camera_center(V);
            viewer.data().show_texture = true;
        }
        if(key == '.')
        {
            scale /= 2.0;
            std::cout<<"scale:"<<scale<<std::endl;
            viewer.data().clear();
            viewer.data().set_mesh(V,F);
            viewer.data().set_uv(sData.V_o*scale,F);
            viewer.core().align_camera_center(V);
            viewer.data().show_texture = true;
        }
        return false;
    };
    
    plot_mesh(vr,uv,F,{},Eigen::VectorXi());
    vr.callback_key_down = key_down;
    vr.launch();
*/
    
    std::cout<<"Total time:"<< (double)(clock()-now) / CLOCKS_PER_SEC<<"s"<<std::endl;
}
