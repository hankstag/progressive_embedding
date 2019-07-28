#include "path_tracing.h"
#include <igl/boundary_loop.h>
#include <igl/matrix_to_list.h>
#include <igl/triangle_triangle_adjacency.h>
#include <iostream>
#include "dijkstra.h"
// #include <igl/opengl/glfw/imgui/ImGuiMenu.h>
// #include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
// #include "util.h"

// tracing through faces on dual graph of F
// generate a path
void direct_geodesic(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const std::vector<std::pair<int,int>>& M, // should also consider edges already cut()
    const std::set<int>& no_enter_f,
    const Eigen::MatrixXi& TT,
    int s,
    int t,
    Eigen::MatrixXd& Vn,
    Eigen::MatrixXi& Fn,
    std::vector<int>& L
){
    assert("should never be inside geodesic");
    std::cout<<"geodesic wrong"<<std::endl;
    exit(0);
    L.clear();
    std::vector<std::vector<int>> VF,VFi;
    igl::vertex_triangle_adjacency(V,F,VF,VFi);
    // Eigen::MatrixXi TT,TTi;
    // igl::triangle_triangle_adjacency(F,TT,TTi);

    // for(int f: no_enter_f){
    //     TT.row(f)<<-1,-1,-1;
    // }
    // for(int i=0;i<F.rows();i++){
    //     for(int k=0;k<3;k++)
    //         if(no_enter_f.find(TT(i,k))!=no_enter_f.end())
    //             TT(i,k) = -1;
    // }

    // for(int i=0;i<M.size();i++){
    //     int a = M[i].first;
    //     int b = M[i].second;
    //     std::vector<int> faces = VF[a];
    //     std::vector<int> adf;
    //     for(int j=0;j<faces.size();j++){ 
    //         for(int k=0;k<3;k++){
    //             int f = faces[j];
    //             if(F(f,k) == b)
    //                 adf.push_back(f);
    //         }
    //     }
    //     for(int j=0;j<3;j++){
    //         int f1 = adf[0];
    //         int f2 = adf[1];
    //         if(TT(f1,j) == f2)
    //             TT(f1,j) = -1;
    //         if(TT(f2,j) == f1)
    //             TT(f2,j) = -1;
    //     }
    // }
    //std::cout<<"difference "<<(TT-TT_x).norm()<<std::endl;
    std::vector<std::vector<int>> Du;
    igl::matrix_to_list(TT,Du);
    // drop the -1 terms
    for(auto& R: Du){
        for(int i=0;i<R.size();i++)
            if(R[i]==-1){
                R.erase(R.begin()+i);
                i--;
            }
    }
    std::cout<<"blocking faces"<<std::endl;
    for(int i: no_enter_f)
        std::cout<<i<<" ";
    std::cout<<std::endl;
    // drop faces in no_enter_f from VF[s] and VF[t]
    std::vector<int> vfs,vft;
    for(int i=0;i<VF[s].size();i++){
        if(no_enter_f.find(VF[s][i])==no_enter_f.end())
            vfs.push_back(VF[s][i]);
    }
    for(int i=0;i<VF[t].size();i++){
        if(no_enter_f.find(VF[t][i])==no_enter_f.end())
            vft.push_back(VF[t][i]);
    }
    // add s and t as virtual faces into Du
    Du.push_back(vfs);
    Du.push_back(vft);
    for(int f: vfs){
        Du[f].push_back(Du.size()-2);
    }
    for(int f: vft)
        Du[f].push_back(Du.size()-1);

    // [find shortest face sequence]
    std::set<int> T = {int(Du.size()-1)};
    Eigen::VectorXd minDistance;
    Eigen::VectorXi previous;
    igl::dijkstra(int(Du.size()-2),T,Du,minDistance,previous);
    igl::dijkstra(int(Du.size()-1),previous,L);
    // if(L.size()==1){
    //     std::cout<<"did not find path"<<std::endl;
    //     std::vector<int> P = {s,t};
    //     std::vector<std::pair<int,int>> E;
    //     std::vector<std::vector<int>> Q;
    //     std::vector<int> H(no_enter_f.begin(),no_enter_f.end());
    //     display(V,F,M,P,H,Q);
    //     exit(0);
    // }
    // drop the fake face
    L.erase(L.begin());
    L.resize(L.size()-1);

    // [connect the barycenters of faces in list L]
    // [collect all barycenter and midpoints]

    int nv = Vn.rows();
    int nf = Fn.rows();
    std::vector<int> P = {t,s};
    std::vector<int> Pt(2*L.size()-1);
    std::iota(Pt.begin(),Pt.end(),Vn.rows());
    P.insert(P.begin()+1,Pt.begin(),Pt.end());
    Vn.conservativeResize(Vn.rows()+2*L.size()-1,Eigen::NoChange);
    Fn.conservativeResize(Fn.rows()+5*L.size()-2,3);
    int v10,v20,v1,v2,v3;
    for(int i=0;i<L.size();i++){
        int f = L[i];
        int id = -1; // id of the vertex corresponds to the common edge
        if(i==L.size()-1)
            id = 0;
        else{
            for(int k=0;k<3;k++){
               if(TT(f,k)==L[i+1]){
                    id = k;
                    break;
                }
            }
        }
        v1 = F(f,id);
        v2 = F(f,(id+1)%3);
        v3 = F(f,(id+2)%3);
        Vn.row(nv)<<(Vn.row(v1)+Vn.row(v2)+Vn.row(v3)) / 3;
        if(i == 0){
            Fn.block(nf,0,4,3)<<v1,nv+1,nv,nv+1,v2,nv,v2,v3,nv,v3,v1,nv;
            Vn.row(nv+1)<<(Vn.row(v1)+Vn.row(v2)) / 2;
            nf += 4;
            nv += 2;
        }else if(i == L.size()-1){
            Fn.block(nf,0,4,3)<<s,v20,nv,v20,nv-1,nv,nv-1,v10,nv,v10,s,nv;
            nf += 4;
        }else{
            if(v10 != v1){
                Fn.block(nf,0,5,3)<<v1,nv+1,nv,nv+1,v2,nv,v2,nv-1,nv,nv-1,v10,nv,v10,v1,nv;
            }else
                Fn.block(nf,0,5,3)<<v1,nv+1,nv,nv+1,v2,nv,v2,v20,nv,v20,nv-1,nv,nv-1,v10,nv;
            Vn.row(nv+1)<<(Vn.row(v1)+Vn.row(v2)) / 2;
            nf += 5;
            nv += 2;
        }
        Fn.row(f)<<-1,-1,-1;
        v10 = v1;
        v20 = v2;
    }
    L = P;
    // drop the (-1,-1,-1) rows
    Eigen::MatrixXi tF=Fn;
	int k=0;
	for(int i=0;i<Fn.rows();i++)
		if(Fn.row(i).sum()!=-3){
			tF.row(k++)<<Fn.row(i);
        }
	tF.conservativeResize(k,3);
	Fn = tF;
}

bool path_tracing(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const std::pair<int,int>& T,
    const std::set<int>& A, // vertices needs to avoid
    const std::vector<std::pair<int,int>>& Mask,
    const std::set<int>& no_enter_f,
    const Eigen::MatrixXi& TT,
    Eigen::MatrixXd& Vn,
    Eigen::MatrixXi& Fn,
    std::vector<int>& E
){
    Vn = V;
    Fn = F;
    E.clear();
    std::vector<std::vector<int>> J;
    igl::adjacency_list(F,J);

    auto cut_off_vertices = [&Fn](
        const std::pair<int,int>& T,
        const std::set<int>& A, 
        std::vector<std::vector<int>>& M, 
        const std::set<int>& no_enter_f
    ){
        if(A.empty()) return;
        for(int v=0;v<M.size();v++){
            if(A.find(v)==A.end()) continue;
            for(int neighbor : M[v]){
                auto ptr = std::find(M[neighbor].begin(),M[neighbor].end(),v);
                if(ptr!=M[neighbor].end())
                    M[neighbor].erase(ptr);
            }
            M[v].clear();
        }
        for(int f: no_enter_f){
            for(int i=0;i<3;i++){
                if(Fn(f,i) == T.first || Fn(f,i) == T.second){ // suppse it's c
                    int a = Fn(f,(i+2)%3);
                    int b = Fn(f,(i+1)%3);
                    int c = Fn(f,i);
                    // break (c,a) and (c,b)
                    // std::cout<<"break connection of "<<a<<","<<c<<std::endl;
                    // std::cout<<"break connection of "<<b<<","<<c<<std::endl;
                    auto it1 = std::find(M[c].begin(),M[c].end(),a);
                    if(it1!=M[c].end())
                        M[c].erase(it1);
                    auto it2 = std::find(M[c].begin(),M[c].end(),b);
                    if(it2!=M[c].end())
                        M[c].erase(it2);
                    auto it3 = std::find(M[a].begin(),M[a].end(),c);
                    if(it3!=M[a].end())
                        M[a].erase(it3);
                    auto it4 = std::find(M[b].begin(),M[b].end(),c);
                    if(it4!=M[b].end())
                        M[b].erase(it4);
                }
            }
        }
    };

    // cut off the vertices in A from adj_map
    cut_off_vertices(T,A,J,no_enter_f);
    Eigen::VectorXd min_d_;
    Eigen::VectorXi prev_;
    std::set<int> x_;
    igl::dijkstra(V,J,T.first,{T.second},min_d_,prev_);
    igl::dijkstra(T.second, prev_, E);
    if(E.size() == 1){
        return false;
    }
    return true;
}
