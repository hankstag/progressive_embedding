#include "matchmaker.h"
#include "decompose_polygon.h"
#include "mst.h"
#include "cut_mesh/HalfEdgeIterator.h"
#include "edge_split.h"
#include "plot.h"
#include "validity_check.h"

#include <igl/vertex_triangle_adjacency.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/boundary_loop.h>
#include <igl/matrix_to_list.h>
// #include <igl/opengl/glfw/imgui/ImGuiMenu.h>
// #include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include <unordered_set>
#include <unordered_map>
#include <igl/writeOBJ.h>
#include <igl/remove_unreferenced.h>
#include "cut_mesh/edge_flaps.h"
#include "progressive_embedding.h"

#include <igl/facet_components.h>
#include <igl/copyleft/cgal/point_inside_polygon.h>

void boundary_straightening(
  Eigen::MatrixXd& P
){
  int n = P.rows();
  // find first non-horizontal, non-vertical corner point
  int i=0;
  int offset=-1;
  while(i<n){
    int prev = (i+n-1)%n;
    int curr = i;
    int next = (i+1) % n;
    bool vertical   = (std::abs(P(curr,0) - P(prev,0))<1e-15 && std::abs(P(curr,0) - P(next,0))<1e-15);
    bool horizontal = (std::abs(P(curr,1) - P(prev,1))<1e-15 && std::abs(P(curr,1) - P(next,1))<1e-15);
    if(!vertical && !horizontal)
      offset=i;
    i++;
  }
  bool colinear;
  double x_mem = P(offset,0);
  double y_mem = P(offset,1);
  for(int d=0;d<n;d++){
    int i = (d+offset) % n;
    int prev = (i+n-1)%n;
    int curr = i;
    int next = (i+1) % n;
    bool vertical   = (std::abs(P(curr,0) - P(prev,0))<1e-15 && std::abs(P(curr,0) - P(next,0))<1e-15);
    bool horizontal = (std::abs(P(curr,1) - P(prev,1))<1e-15 && std::abs(P(curr,1) - P(next,1))<1e-15);
    if(vertical)
      P(curr,0) = x_mem;
    if(horizontal)
      P(curr,1) = y_mem;
    if(!horizontal && !vertical){
      x_mem = P(curr,0);
      y_mem = P(curr,1);
    }
  }

}

// precondition of Tutte's embedding
// make sure the mesh is 3-connected
// by spliting edges that connects any
// non-consecutive boundary vertices
void remove_ears(
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F
){
    // get rid of ears
    Eigen::VectorXi bd;
    igl::boundary_loop(F,bd);
    std::unordered_map<int,int> mm;
    std::unordered_set<int> is_constrained;
    for(int i=0;i<bd.rows();i++){
        is_constrained.insert(bd(i));
        mm[bd(i)] = i;
    }
    
    Eigen::MatrixXi TT,TTi;
    igl::triangle_triangle_adjacency(F,TT,TTi);
    for(int i=0;i<F.rows();i++){
        bool all_constrained = true;
        for(int k=0;k<3;k++){
            bool a1 = (is_constrained.find(F(i,k))==is_constrained.end());
            bool a2 = (is_constrained.find(F(i,(k+1)%3))==is_constrained.end());
            if(!a1 && !a2 && std::abs(mm[F(i,k)]-mm[F(i,(k+1)%3)])!=1){
                edge_split(V,F,TT,TTi,i,k);
            }
        }
    }
}

// based on the polygon list build a graph
void build_graph(
    int n, // size of boundary
    const Eigen::MatrixXd& V,
    const std::vector<std::vector<int>>& L,
    Eigen::SparseMatrix<int>& graph
){
    int nm = 0;
    for(auto z: L)
        nm += z.size();
    
    // initialize graph
    typedef Eigen::Triplet<int> TP;
    std::vector<TP> tripletList;
    tripletList.reserve(nm);
    bool drop = false; // drop one edge on the boundary
    for(auto l: L){
        for(int j=0;j<l.size();j++){
            int j_1 = (j+1) % l.size();
            if(std::abs(l[j]-l[j_1])==n-1){
                continue;
            }
            double d = (V.row(l[j])-V.row(l[j_1])).norm();
            int val = (std::abs(l[j]-l[j_1]) == 1) ? -1 : 1;
            tripletList.push_back(TP(l[j],l[j_1],val));
            tripletList.push_back(TP(l[j_1],l[j],val));
        }
    }
    graph.resize(V.rows(),V.rows());
    graph.setFromTriplets(tripletList.begin(), tripletList.end());
}

void neighbor_sector_bar(
    const std::vector<std::vector<int>>& L, // polygon list
    const Eigen::MatrixXd& V2,
    const Eigen::MatrixXi& F2,
    const std::set<std::pair<int,int>>& sector_bar,
    int s,
    int t,
    std::pair<int,int>& b00,
    std::pair<int,int>& b01
){
    // build a fake mesh around s
    Eigen::MatrixXi F_fk(L.size(),3);
    int n = 0;
    for(int i=0;i<L.size();i++){
        if(L[i].empty()) continue;
        for(int j=0;j<L[i].size();j++){
            if(L[i][j] == s){
                int j_prev = (j-1+L[i].size())%L[i].size();
                int j_next = (j+1)%L[i].size();
                F_fk.row(n++)<<L[i][j_prev],s,L[i][j_next];
            }
        }
    }
    F_fk.conservativeResize(n,3);
    Eigen::VectorXi b_fk;
    igl::boundary_loop(F_fk,b_fk);
    int location_of_t = -1;
    for(int i=0;i<b_fk.rows();i++){
        if(b_fk(i) == t){
            location_of_t = i;
            break;
        }
    }
    assert(location_of_t != -1);
    bool prev_found = false;
    for(int i=0; i<b_fk.size(); i++) {
        // find immediate prev in sector
        int prev = b_fk[(location_of_t-i-1+b_fk.size())%b_fk.size()];
        auto pair1 = std::make_pair(prev, s);
        auto pair2 = std::make_pair(s, prev);
        if (sector_bar.find(pair1)!=sector_bar.end()){
            b00 = pair1;
            prev_found = true;
            break;
        }
        if(sector_bar.find(pair2)!=sector_bar.end()){
            b00 = pair2;
            prev_found = true;
            break;
        }
    }
    bool next_found = false;
    for(int i=0; i<b_fk.size(); i++) {
        // find immediate next in sector
        int next = b_fk[(location_of_t+i)%b_fk.size()];
        auto pair1 = std::make_pair(next, s);
        auto pair2 = std::make_pair(s, next);
        if (sector_bar.find(pair1)!=sector_bar.end()){
            b01 = pair1;
            next_found = true;
            break;
        }
        if(sector_bar.find(pair2)!=sector_bar.end()){
            b01 = pair2;
            next_found = true;
            break;
        }
    }
    // if not prev and not next
    if(!prev_found && !next_found){
        b00.first = -1;
        b00.second = -1;
        b01.first = -1;
        b01.second = -1;
    }
    
}

bool match_sector_bar(
    std::map<std::pair<int,int>,std::vector<int>>& splits,
    const Eigen::VectorXi& T,
    std::pair<int,int>& b,
    std::pair<int,int>& c,
    int center
){
    if(b.first == -1) return false; // invalid b
    auto b_r = std::make_pair(b.second,b.first);
    if(splits.find(b)!=splits.end()){ // this sector bar is a trace
        if(b.first == center)
            c = std::make_pair(splits[b][0],
                                splits[b][1]);
        else if(b.second == center){
            int n = splits[b].size();
            c = std::make_pair(splits[b][n-2],splits[b][n-1]);
        }
    }else if(splits.find(b_r)!=splits.end()){
        if(b_r.first == center)
            c = std::make_pair(splits[b_r][0],
                                splits[b_r][1]);
        else if(b_r.second == center){
            int n = splits[b_r].size();
            c = std::make_pair(splits[b_r][n-2],
                                splits[b_r][n-1]);
        }
    }else{ // this sector bar is a boundary edge
        c = std::make_pair(T(b.first),T(b.second));
    }
    return true;
}

// -- face block --
void collect_blocked_face(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    std::set<int>& no_enter_f,
    const std::pair<int,int>& c00,
    const std::pair<int,int>& c01
){
    // sweep from c01 to c00
    int c = c00.first; // center point
    if(c != c01.first && c != c01.second)
        c = c00.second;
    // find neighbor of c
    std::vector<int> l;
    std::vector<int> vl; // valid list
    int f_start,e_start,f_end,e_end;
    int v1 = c00.first+c00.second-c; // the other end point of c00
    int v2 = c01.first+c01.second-c; // ...                    c01
    for(int i=0;i<F.rows();i++){
        for(int j=0;j<3;j++){
            if(F(i,j) == c){
                l.push_back(i);
                if(F(i,(j+1)%3) == v1)
                    f_start = l.size()-1, e_start = j;
                if(F(i,(j+2)%3) == v2)
                    f_end = l.size()-1, e_end = (j+1)%3; 
            }
        }
    }
    Eigen::MatrixXi F_local(l.size(),3);
    for(int i=0;i<l.size();i++){
        F_local.row(i)<<F.row(l[i]);
    }
    Eigen::MatrixXi TT,TTi;
    igl::triangle_triangle_adjacency(F_local,TT,TTi);
    igl::HalfEdgeIterator<Eigen::MatrixXi> heIter(F_local.derived(),TT,TTi,f_start,e_start);
    int current_f = f_start;
    if(heIter.Vi1() != c){
        heIter.flipE();
        vl.push_back(l[current_f]);
    }
    int count = 0;
    while(current_f != f_end){
        count++;
        if(heIter.Fif() == -1)
            break;
        heIter.nextFE();
        
        current_f = heIter.getState().fi;
        vl.push_back(l[current_f]);
    }
    std::sort(l.begin(),l.end());
    std::sort(vl.begin(),vl.end());

    std::vector<int> diff;
    std::set_difference(l.begin(), l.end(), vl.begin(), vl.end(), 
                    std::inserter(diff, diff.begin()));

    std::vector<int> P;
    for(int i: diff)
        no_enter_f.insert(i);
}

// given a pair of trace [a,b]
// mark faces and vertices as impassible
// so that new tracing path won't go through them
void mark_impassible(
    const std::vector<std::vector<int>>& L,
    const Eigen::MatrixXd& V2,
    const Eigen::MatrixXi& F2,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::VectorXi& T, // match 2d to 3d
    const std::set<std::pair<int,int>>& sector_bar,
    const std::pair<int,int>& t, // trace [a,b]
    std::map<std::pair<int,int>,std::vector<int>>& splits, // path traced in 3d
    std::set<int>& no_enter_f
){
    int end_a = t.first;
    int end_b = t.second;
    // [ pick the trace intersect with a and b ]
    std::set<std::pair<int,int>> Ta,Tb;
    for(auto e: sector_bar){
        if(e.first == end_a || e.second == end_a)
            Ta.insert(e);
        if(e.first == end_b || e.second == end_b)
            Tb.insert(e);
    }
    std::pair<int,int> b00,b01,b10,b11;
    neighbor_sector_bar(L,V2,F2,Ta,t.first,t.second,b00,b01); // around t.first
    neighbor_sector_bar(L,V2,F2,Tb,t.second,t.first,b10,b11); // around t.second
    
    std::vector<std::pair<int,int>> E; // highlight edges
    std::vector<int> P; // highlight points
    std::vector<int> H;    
    std::vector<std::vector<int>> Q;
    //display(V2,F2,{b00,b01,b10,b11},P,H,Q);
    // find the match for these four sector bars   
    std::pair<int,int> c00,c01,c10,c11;
    int x1 = b00.first; // common point of b00,b01
    if(x1 != b01.first && x1 != b01.second)
        x1 = b00.second;
    int x2 = b10.first;// common point of b10,b11
    if(x2 != b11.first && x2 != b11.second)
        x2 = b10.second;

    if(match_sector_bar(splits,T,b00,c00,end_a) && match_sector_bar(splits,T,b01,c01,end_a))
        collect_blocked_face(V,F,no_enter_f,c00,c01);
    if(match_sector_bar(splits,T,b10,c10,end_b) && match_sector_bar(splits,T,b11,c11,end_b))
        collect_blocked_face(V,F,no_enter_f,c10,c11);
    // std::cout<<"==> b00 "<<b00.first<<","<<b00.second<<std::endl;
    // std::cout<<"==> b01 "<<b01.first<<","<<b01.second<<std::endl;
    // std::cout<<"==> b10 "<<b10.first<<","<<b10.second<<std::endl;
    // std::cout<<"==> b11 "<<b11.first<<","<<b11.second<<std::endl;
    // std::cout<<"==> c00 "<<c00.first<<","<<c00.second<<std::endl;
    // std::cout<<"==> c01 "<<c01.first<<","<<c01.second<<std::endl;
    // std::cout<<"==> c10 "<<c10.first<<","<<c10.second<<std::endl;
    // std::cout<<"==> c11 "<<c11.first<<","<<c11.second<<std::endl;
    for(auto p: splits)
        Q.push_back(p.second);
    //display(V,F,{c00,c01,c10,c11},P,H,Q);

}

void prapare_TT(
    const Eigen::MatrixXi& F,
    const std::set<int>& no_enter_f,
    const std::vector<std::pair<int,int>>& M,
    const std::vector<std::vector<int>>& VF,
    Eigen::MatrixXi& TT,
    Eigen::MatrixXi& TTi
){

    for(int f: no_enter_f){
        assert(f<TT.rows());
        TT.row(f)<<-1,-1,-1;
        TTi.row(f)<<-1,-1,-1;
    }
    for(int i=0;i<F.rows();i++){
        for(int k=0;k<3;k++)
            if(no_enter_f.find(TT(i,k))!=no_enter_f.end()){
                assert(i<TT.rows());
                TT(i,k) = -1;
                TTi(i,k) = -1;
            }
    }
    for(int i=0;i<M.size();i++){
        int a = M[i].first;
        int b = M[i].second;
        std::vector<int> faces = VF[a];
        std::vector<int> adf;
        for(int j=0;j<faces.size();j++){ 
            for(int k=0;k<3;k++){
                int f = faces[j];
                if(F(f,k) == b)
                    adf.push_back(f);
            }
        }
        assert(adf.size()!=0);
        for(int j=0;j<3;j++){
            int f1 = adf[0];
            int f2 = adf[1];
            if(TT(f1,j) == f2)
                TT(f1,j) = -1;
            if(TT(f2,j) == f1)
                TT(f2,j) = -1;
        }
    }
}
void test_decompose(){
  Eigen::MatrixXd P(9,2);
  P << 0,0,2,0,4,0,4,2,1,2,1,1,2,1,2,4,0,4;
  Eigen::VectorXi R(8);
  R.setZero();
  R(4) = 1;
  Eigen::MatrixXd c;
  Eigen::MatrixXd V2;
  Eigen::MatrixXi F2;
  std::vector<std::vector<int>> L;

  decompose_polygon(P,R,c,V2,F2,L);
}
void match_maker(
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  Eigen::MatrixXd& uv,
  const Eigen::MatrixXd& c,
  const Eigen::VectorXi& ci,
  const Eigen::VectorXi& R,
  const Eigen::VectorXi& T_s,
  const Eigen::MatrixXd& P_s
){
  //test_decompose();
  Eigen::VectorXi T = T_s;
  Eigen::MatrixXd P = P_s;
  boundary_straightening(P);
  Eigen::VectorXi H(P.rows());
  H.setZero();
  //igl::opengl::glfw::Viewer vr;
  //plot_polygon(vr,H,P);
  remove_ears(V,F);
  // [ use Shor to get list of polygons ]
  Eigen::MatrixXd V2;
  Eigen::MatrixXi F2;
  std::vector<std::vector<int>> L;

  decompose_polygon(P,R,c,V2,F2,L);

  Eigen::VectorXi bound;
  igl::boundary_loop(F,bound);
    int bod_num = T.rows();
    // update corresponding map
    T.conservativeResize(T.rows()+c.rows());
    T.bottomRows(c.rows())=ci;

    // [ build a minimum spanning tree ]
    std::vector<int> parent(T.rows(),0);
    std::set<int> mst_set;
    Eigen::SparseMatrix<int> graph;
    build_graph(bod_num,V2,L,graph);

    for(int i=0;i<bod_num;i++){
        mst_set.insert(i);
        parent[i] = i-1;
    }
    mst(V2,F2,T.rows()-1,parent,mst_set,graph);

    // [ by the merit of shor algorithm, we know the ]
    // [ range of boundary vertices are 0: bod_num-1 ]
    std::set<std::pair<int,int>> to_trace;
    for(auto l: L){
        for(int j=0;j<l.size();j++){
            int j_1 = (j+1)%l.size();
            bool is_internal_edge = (l[j] > bod_num - 1 || l[j_1] > bod_num -1);
            bool is_adjacent = (std::abs(l[j]-l[j_1])==1) || (std::abs(l[j]-l[j_1])==bod_num-1);
            bool is_boundary_edge = !is_internal_edge && is_adjacent;
            if(!is_boundary_edge){
                // edge (j,j_1) at least one internal vertex
                auto p = std::make_pair(std::min(l[j],l[j_1]),std::max(l[j],l[j_1]));
                to_trace.insert(p);
            }
        }
    }
    // std::cout<<"to trace number "<<to_trace.size()<<std::endl;
    std::set<int> no_enter; // set of vertices that should avoid
    std::vector<std::pair<int,int>> mask_e;
    std::map<std::pair<int,int>,std::vector<int>> splits;
    // std::cout<<"# to trace "<<to_trace.size()<<std::endl;
    Eigen::VectorXi bi;
    igl::boundary_loop(F,bi);
    for(int i=0;i<bi.rows();i++)
        no_enter.insert(bi(i));

    // [ traces already in the mesh together with 
    //   boundary edges forms a set of sectors ]
    std::set<std::pair<int,int>> sector_bar;
    for(auto e: to_trace){
        for(int i : {e.first,e.second}){
            if(i<bod_num){
                int i_prev = (i-1+bod_num)%bod_num;
                int i_next = (i+1)%bod_num;
                sector_bar.insert(std::make_pair(i,i_prev));
                sector_bar.insert(std::make_pair(i_next,i));
            }
        }
    }

    // randomize to_trace
    std::vector<std::pair<int,int>> to_trace_rand(to_trace.begin(),to_trace.end());
    std::random_shuffle(to_trace_rand.begin(), to_trace_rand.end());
    std::vector<std::pair<int,int>> to_trace_ordered;
    std::vector<std::pair<int,int>> trace_mst;
    int trace_in_mst = 0;
    for(auto e: to_trace_rand){
        int l = e.first;
        int r = e.second;
        // if this edge is within the mst tree
        if((parent[l] == r && graph.coeff(r,l)) ||
           (parent[r] == l && graph.coeff(l,r))
        ){ 
            to_trace_ordered.insert(to_trace_ordered.begin(),e);
            trace_mst.insert(trace_mst.begin(),e);
        }else
            to_trace_ordered.push_back(e);
    }

    std::set<std::pair<int,int>> on_path;
    std::set<int> on_path_v;
    for(int i=0;i<bound.rows();i++){
        int a = std::min(bound(i),bound((i+1)%bound.rows()));
        int b = std::max(bound(i),bound((i+1)%bound.rows()));
        on_path.insert(std::make_pair(a,b));
        on_path_v.insert(bound(i));
    }
    Eigen::MatrixXi TT,TTi;
    igl::triangle_triangle_adjacency(F,TT,TTi);
    std::vector<std::vector<int>> VF,VFi;
    igl::vertex_triangle_adjacency(V,F,VF,VFi);
    int xc = 0;
    for(auto e: to_trace_ordered){
        int l = e.first;
        int r = e.second;
        // std::cout<<"traced 2d "<<l<<" to "<<r<<std::endl;        
        // std::cout<<"traced 3d "<<T(l)<<" to "<<T(r)<<std::endl;
    }
    for(auto e: to_trace_ordered){
        int l = e.first;
        int r = e.second;
        std::cout<<"traced 2d "<<l<<" to "<<r<<std::endl;        
        std::cout<<"traced 3d "<<T(l)<<" to "<<T(r)<<std::endl;
        for(int i=0;i<ci.rows();i++){
            if(ci(i) != T(l) && ci(i) != T(r))
                no_enter.insert(ci(i));
        }
        // start and end point should not be avoided
        for(int x: {T(l),T(r)}){
            auto it = no_enter.find(x);
            if(it!=no_enter.end()){
                no_enter.erase(it);
            }
        }
        std::set<int> no_enter_f;
        mark_impassible(L,V2,F2,V,F,T,sector_bar,e,splits,no_enter_f);
        auto Vn = V;
        auto Fn = F;
        std::vector<int> E; // traced path
        path_tracing(V,F,std::make_pair(T(l),T(r)),no_enter,mask_e,no_enter_f,TT,Vn,Fn,E);
        std::reverse(E.begin(),E.end());
        splits[e] = E;
        assert(E.size()>=2);
        assert(T(e.first) == E[0] && T(e.second) == E.back());
        // update impassible edges info
        for(int i=0;i<E.size();i++){
            no_enter.insert(E[i]);
            if(i == E.size()-1) break;
            int a = E[i];
            int b = E[i+1];
            mask_e.push_back(std::make_pair(a,b));
            on_path.insert(std::make_pair(std::min(a,b),std::max(a,b)));
            on_path_v.insert(a);
            on_path_v.insert(b);
        }
        sector_bar.insert(e);
        //igl::triangle_triangle_adjacency(F,TT,TTi);
        igl::vertex_triangle_adjacency(V,F,VF,VFi);
        no_enter_f.clear();
        prapare_TT(F,no_enter_f,mask_e,VF,TT,TTi);
        // split edge whose both end points are on boundary
        std::vector<std::vector<double>> V_vec;
        std::vector<std::vector<int>> F_vec;
        for(int i=0;i<V.rows();i++)
            V_vec.push_back({V(i,0),V(i,1),V(i,2)});
        for(int i=0;i<F.rows();i++)
            F_vec.push_back({F(i,0),F(i,1),F(i,2)});
        std::vector<std::vector<int>> TT_vec;
        std::vector<std::vector<int>> TTi_vec;
        igl::matrix_to_list(TT,TT_vec);
        igl::matrix_to_list(TTi,TTi_vec);
        for(int i=0;i<F_vec.size();i++){
            for(int k=0;k<3;k++){
                int a1 = std::min(F_vec[i][k],F_vec[i][(k+1)%3]);
                int a2 = std::max(F_vec[i][k],F_vec[i][(k+1)%3]);
                auto t = std::make_pair(a1,a2);
                if(on_path.find(t) != on_path.end()) continue;
                if(on_path_v.find(F_vec[i][k])!=on_path_v.end() &&
                   on_path_v.find(F_vec[i][(k+1)%3])!=on_path_v.end()){
                    bool status = edge_split(V_vec,F_vec,TT_vec,TTi_vec,i,k);
                    assert(status!=false && "split boundary");
                }
            }
        }
        igl::list_to_matrix(TT_vec,TT);
        igl::list_to_matrix(TTi_vec,TTi);
        Eigen::MatrixXd V0;
        Eigen::MatrixXi F0;
        igl::list_to_matrix(F_vec,F);
        igl::list_to_matrix(V_vec,V);
        std::cout<<"iteration "<<xc++<<"/"<<to_trace_ordered.size()<<std::endl;
        std::cout<<"current #F "<<F.rows()<<std::endl;
    }
    Eigen::MatrixXd CN,TC;
    Eigen::MatrixXi FN,FTC;

    // sample on to_trace
    Eigen::MatrixXd known;
    std::vector<std::vector<double>> known_vec;
    Eigen::VectorXi ki;
    std::vector<int> ki_vec;
    for(auto e: to_trace_ordered){
        int a = e.first;
        int b = e.second;
        // on 2d [a     -     b]
        // on 3d [T(a) --- T(b)]
        std::vector<int> path = splits[e];
        if(path.size()>2){
            Eigen::VectorXd ratio;
            ratio.setZero(path.size()-2);
            double total_len = 0.0;
            for(int j=0;j<path.size()-1;j++){
                total_len += (V.row(path[j])-V.row(path[(j+1)%path.size()])).norm();
            }
            double temp_len = 0.0;
            for(int j=1;j<path.size()-1;j++){
                temp_len+=(V.row(path[j])-V.row(path[j-1])).norm();
                ratio(j-1) = 1.0*j / (path.size()-1);
            }
            for(int k=0;k<ratio.size();k++){
                Eigen::RowVector2d t = ratio(k)*(V2.row(b)-V2.row(a))+V2.row(a);
                known_vec.push_back({t(0),t(1)});
                ki_vec.push_back(path[k+1]);
            }
        }
    }
    Eigen::VectorXi bdr;
    igl::boundary_loop(F,bdr);
    igl::list_to_matrix(known_vec,known);
    V2.conservativeResize(V2.rows()+known.rows(),V2.cols());
    V2.bottomRows(known.rows())<<known;
    igl::list_to_matrix(ki_vec,ki);
    T.conservativeResize(T.rows()+ki.rows());
    T.bottomRows(ki.rows())<<ki;
    Eigen::MatrixXd H1;
    Eigen::VectorXi nb=T;
    Eigen::MatrixXd nbc=V2;
    igl::harmonic(F ,nb,nbc,1,H1);
    H1.conservativeResize(H1.rows(),2);
    uv = H1;
    igl::opengl::glfw::Viewer vr;
    plot_mesh(vr,H1,F,{},Eigen::VectorXi());
    vr.launch();
    // decompose mesh into patches
    Eigen::VectorXi group;
    Eigen::VectorXi counts;
    std::vector<std::vector<std::vector<int>>> TTv;
    for(int i=0;i<TT.rows();i++){
      std::vector<int> x,y,z;
      if(TT(i,0)>=0) x.push_back(TT(i,0));
      if(TT(i,1)>=0) x.push_back(TT(i,1));
      if(TT(i,2)>=0) x.push_back(TT(i,2));
      TTv.push_back({x,y,z});
    }
    igl::facet_components(TTv,group,counts);
    std::cout<<counts<<std::endl;
    
    int patch_num = group.maxCoeff()+1;
    std::vector<int> fn(patch_num,0);
    for(int i=0;i<group.rows();i++){
      fn[group(i)]++;
    }
    // build the patches
    std::vector<Eigen::MatrixXd> PV; // PV[i] is vertex positions of ith patch
    std::vector<Eigen::MatrixXd> Puv; // PV[i] is vertex positions of ith patch
    std::vector<Eigen::MatrixXi> PF; // PF[i] is faces of ith patch
    std::vector<Eigen::VectorXi> PM; // PM[i] is the correspondence for each patch to original mesh
    std::vector<int> ptr(patch_num,0); // ptr[i] points to the last valid line in PF[i]
    for(int k=0;k<patch_num;k++){ // initialize size of PF
      Eigen::MatrixXi kF(fn[k],3);
      Eigen::MatrixXd kuv, kV;
      PV.push_back(kV);
      Puv.push_back(kuv);
      std::cout<<"size of patch "<<k<<" is "<<fn[k]<<std::endl;
      PF.push_back(kF);
    }
    for(int f=0;f<F.rows();f++){
      int g = group(f);
      PF[g].row(ptr[g]++) << F.row(f);
    }
    for(int k=0;k<patch_num;k++){
      Eigen::MatrixXi kF;
      Eigen::VectorXi II;
      igl::remove_unreferenced(V, PF[k],PV[k] ,kF,II); // PF[k] to V
      igl::remove_unreferenced(H1,PF[k],Puv[k],kF,II);
      Eigen::VectorXi fix;
      Eigen::MatrixXd fix_pos;
      igl::boundary_loop(kF,fix);
      igl::slice(Puv[k],fix,1,fix_pos);
      //igl::slice(nb,II,1,fix);
      //igl::slice(nbc,II,1,fix_pos);
      Eigen::VectorXi fl;
      flipped_elements(Puv[k],kF,fl);
      if(fl.sum() > 0){ // Tutte generates flips
        std::cout<<"flips "<<fl.sum()<<std::endl;
        Eigen::RowVector2d bc;
        bc.setZero();
        for(int r=0;r<fix_pos.rows();r++)
          bc += fix_pos.row(r);
        bc /= fix_pos.rows();
        
        Eigen::VectorXi BI(Puv[k].rows());
        BI.setZero();
        for(int i=0;i<fix.rows();i++)
          BI(fix(i)) = 1;
        for(int id = 0;id<Puv[k].rows();id++){
          Eigen::RowVector2d pt = Puv[k].row(id);
          if(!BI(id) && !igl::copyleft::cgal::point_inside_polygon(fix_pos,pt))
            Puv[k].row(id) << bc;
        }
        std::vector<Object> Os = {Object(PV[k],kF,OTYPE::MESH),
                                  Object(Puv[k],kF,OTYPE::MESH)};
        plots(Os);
        // Eigen::MatrixXd CN;
        // Eigen::MatrixXi FN;
        //igl::writeOBJ("filigree_patch.obj",PV[k],kF,CN,FN,Puv[k],kF);
        progressive_embedding(PV[k],kF,Puv[k],fix,fix_pos,1e15);
        // copy new positions back
        for(int i=0;i<uv.rows();i++){
          if(II(i) != -1)
            uv.row(i) << Puv[k].row(II(i));
        }
      }

      
    }
    Eigen::VectorXi I;
    flipped_elements(uv,F,I);
    std::cout<<"check flips: "<<I.sum()<<std::endl;
    igl::opengl::glfw::Viewer viewer;
    viewer.data().clear();
    viewer.data().set_mesh(uv,F);
    //viewer.data().add_points(bc,Eigen::RowVector3d(1,0,0));
    viewer.launch();
}
