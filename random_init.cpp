#include "local_smooth/auto_grad.hpp"
#include "validity_check.h"
#include "local_smooth/local_smooth.h"
#include "progressive_embedding.h"
#include <igl/harmonic.h>
#include <igl/boundary_loop.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/predicates/point_inside_convex_polygon.h>
#include <cstdlib>
#include "loader.h"
#include <ctime>
#include <igl/jet.h>
#include <igl/colormap.h>
#include "argh.h"

using namespace std;

void compuate_energy(
    const Eigen::MatrixXd& uv,
    const Eigen::MatrixXi& F,
    Eigen::VectorXd& Energy
){
    Eigen::VectorXd area;
    igl::doublearea(uv,F,area);
    int ctf=0;
    double target_area = 0;
    for(int k=0;k<area.rows();k++){
        if(!isnan(area(k)) && !isinf(area(k))){
            target_area += area(k);
            ctf++;
        }
    }
    auto dirichlet_energy = [](
        double target_area,
        const Eigen::MatrixXd& uv,
        int v0,
        int v1,
        int v2
    ){
        // calc gradient operator
        double h = std::sqrt(target_area/sin(M_PI / 3.0));
        Eigen::Matrix<double,3,3> Tx;
        Tx<<0,0,0,h,0,0,h/2.,(std::sqrt(3)/2.)*h,0;
        Eigen::Matrix<double,3,3> gx;
        grad_operator(Tx,gx);
        Eigen::Matrix<double,2,3> G_t = gx.topRows(2);
        Eigen::Matrix<double,3,2> Tuv;
        Tuv<<uv.row(v0),uv.row(v1),uv.row(v2);
        // calc energy
        return autogen::sd_energy(Tuv,G_t);
    };
    target_area /= ctf;
    Energy.resize(F.rows());
    for(int i=0;i<Energy.rows();i++)
        Energy(i) = dirichlet_energy(target_area,uv,F(i,0),F(i,1),F(i,2));
}
vector<double> randPoint(double x_c, double y_c, double r)
{
    auto fRand = [](double fMin, double fMax) {
        double f = (double)rand() / RAND_MAX;
        return fMin + f * (fMax - fMin);
    };
    auto inside_circle = [&](std::vector<double> &pt) {
        if (((pt[0] - x_c) * (pt[0] - x_c)) + (pt[1] - y_c) * (pt[1] - y_c) > r * r)
            return false;
        else
            return true;
    };
    std::vector<double> pt;
    bool accept = false;
    while (!accept)
    {
        double x = fRand(x_c - r, x_c + r);
        double y = fRand(y_c - r, y_c + r);
        pt = {x, y};
        if (inside_circle(pt))
            accept = true;
    }
    return pt;
}

// random_initialization
int main(int argc, char *argv[])
{
  auto cmdl = argh::parser(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);
  if(cmdl[{"-h","-help"}]){
    std::cout<<"Usage: ./rand_init_bin -options"<<std::endl;
    std::cout<<"-in: input model name"<<std::endl;
    std::cout<<"-o: output model name"<<std::endl;
    exit(0);
  }
  int threshold;
  std::string model,outfile;
  cmdl("-in") >> model;
  cmdl("-o", model+"_out.obj") >> outfile;
  cmdl("-t", 20) >> threshold;
  double eps = std::pow(10,threshold);

  Eigen::MatrixXd V;
  Eigen::MatrixXd uv;
  Eigen::MatrixXi F;
  Eigen::MatrixXd P;

  Eigen::VectorXi R, T, bd;
  
  load_model(model,V,uv,F,P,R,T);
  
  igl::boundary_loop(F, bd);
  std::set<int> bd_set;
  for (int i = 0; i < bd.rows(); i++){
    bd_set.insert(bd(i));
  }
  Eigen::MatrixXd circle;
  igl::map_vertices_to_circle(V, bd, circle);
  igl::harmonic(F, bd, circle, 1, uv);
  for (int i = 0; i < V.rows(); i++){
    srand(i * i * i + std::sin(i) + 5);
    if (bd_set.find(i) == bd_set.end()){
      auto pt = randPoint(0, 0, 1);
      uv.row(i) << pt[0], pt[1];
    }
  }
  Eigen::MatrixXd CN, TC;
  Eigen::MatrixXi FN, FTC;
  igl::writeOBJ(model + "_rand.obj", V, F, CN, FN, uv, F);

}
