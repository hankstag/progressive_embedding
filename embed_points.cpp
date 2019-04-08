
#include "embed_points.h"
#include <igl/AABB.h>
#include <igl/in_element.h>
#include <igl/triangle/triangulate.h>

// given a planar mesh and a bunch of points
// embed these points into the mesh
// if the points are located in the overlapping part
// then will embed in the triangle face with smaller index

void embed_points(
    const Eigen::MatrixXd& C,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F
){
    igl::AABB<Eigen::MatrixXd,2> aabb;
	aabb.init(V,F);
    Eigen::VectorXi I;
	igl::in_element(V,F,C,aabb,I);
    // gather info of where each constraint points land
    // since multiple points may land in the same triangle
    std::map<int,std::vector<int>> I_vec;
    for(int i=0;i<I.rows();i++){
        assert(I(i) != -1 && "constraints lying outside mesh");
        if(I_vec.find(I(i)) == I_vec.end())
            I_vec[I(i)] = {i};
        else
            I_vec[I(i)].push_back(i);
    }
    int nv0 = V.rows();
    V.conservativeResize(V.rows()+C.rows(),Eigen::NoChange);
    // triangulate individually
    for(auto p: I_vec){
        int fid = p.first;
        std::vector<int> v_list = p.second;
        // init
        Eigen::MatrixXd in_poly(3+v_list.size(),2);
        for(int i=0;i<3;i++)
            in_poly.row(i)<<V.row(F(fid,i));
        for(int k=0;k<v_list.size();k++){
            in_poly.row(3+k)<<C.row(v_list[k]);
        }
        Eigen::MatrixXi LE(3,2);
        LE<<0,1,1,2,2,0;
        Eigen::MatrixXd LV;
        Eigen::MatrixXi LF;
		igl::triangle::triangulate(in_poly,LE,Eigen::MatrixXd(),"YQ",LV,LF);
        for(int i=0;i<v_list.size();i++){
            V.row(nv0+v_list[i])<<C.row(v_list[i]);
        }
        Eigen::VectorXi M(LV.rows());
        M(0) = F(fid,0);
        M(1) = F(fid,1);
        M(2) = F(fid,2);
        for(int i=0;i<LV.rows()-3;i++){
            M(i+3) = nv0+v_list[i];
        }

        int nf = F.rows();
        F.conservativeResize(F.rows()+LF.rows(),3);
        for(int i=0;i<LF.rows();i++){
            for(int k=0;k<3;k++){
                F(nf+i,k) = M(LF(i,k));
            }
        }
        F.row(fid)<<-1,-1,-1;

    }
    auto Ft = F;
    int nFt = 0;
    for(int i=0;i<F.rows();i++){
        if(F.row(i).sum()!=-3){
            Ft.row(nFt++) = F.row(i);
        }
    }
    Ft.conservativeResize(nFt,3);
    F = Ft;
}