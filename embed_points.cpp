
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


// #include <igl/boundary_loop.h>
// #include "shor.h"
// #include "is_inside_polygon.h"
// #include "util.h"
// #include "path_tracing.h"
// #include "cut_mesh/cut_mesh_simple.h"
// #include <igl/on_boundary.h>
// #include <igl/opengl/glfw/Viewer.h>
// #include <igl/opengl/glfw/imgui/ImGuiMenu.h>
// #include <igl/opengl/glfw/imgui/ImGuiHelpers.h>


// // global variable
// Eigen::MatrixXi EM; // indicating edge is original/new

// void sequence(
//     const Eigen::MatrixXd& V,
//     const Eigen::MatrixXi& F,
//     const Eigen::VectorXi& I,
//     const Eigen::MatrixXi& E,
//     const std::vector<Eigen::MatrixXd>& P_set,
//     Eigen::VectorXi& Seq,
//     Eigen::VectorXi& mask
// ){
//     // generate the ordered sequence to connect
//     // find closest point in T and S
//     auto match = [](const Eigen::MatrixXd& S, const Eigen::MatrixXd& T, 
//                 int last_i){

//         int step_s = S.rows() / 20 + 1;
//         int step_t = T.rows() / 20 + 1;
//         double md = std::numeric_limits<double>::max();
//         int a, b;
//         for(int i=0;i<S.rows();i+=step_s){
//             for(int j=0;j<T.rows();j+=step_t){
//                 if(i == last_i) continue;
//                 if((S.row(i)-T.row(j)).norm() < md){
//                     a = i;
//                     b = j;
//                     md = (S.row(i)-T.row(j)).norm();
//                 }
//             }
//         }
//         return std::make_pair(a,b);
//     };
//     Eigen::VectorXi B;
//     igl::boundary_loop(F,B);

//     if(I.rows()==0 && P_set.size() == 1) return;
//     // [ideally should use a minimium spanning tree]
//     Eigen::VectorXi K;
//     int np = P_set.size()-1;
//     K.setZero(I.rows() + np);
//     std::vector<int> Seqv;
//     std::vector<std::vector<int>> bd;
//     igl::boundary_loop(F,bd);
//     for(int i=0;i<bd.size();i++){
//         if(bd[i].size() == B.rows() && bd[i][0] == B(0)){
//             bd.erase(bd.begin()+i);
//             i--;
//         }
//     }

//     // calculate all the center point
//     Eigen::MatrixXd center(K.rows(),V.cols());
//     for(int i=0;i<I.rows();i++)
//         center.row(i) << V.row(I(i));
//     int ct = 0;
//     for(int i=0;i<bd.size();i++){
//         Eigen::RowVector3d barycenter;
//         barycenter.setZero();
//         for(int j=0;j<bd[i].size();j++){
//             barycenter += V.row(bd[i][j]);
//         }
//         center.row(ct+I.rows())<<barycenter.array() / bd[i].size();
//         ct++;
//     }

//     // connect center(0) to boundary
//     Eigen::VectorXi opi;
//     Eigen::MatrixXd Vb,op;
    
//     igl::slice(V,B,1,Vb);
//     if(I.rows()==0){
//         igl::list_to_matrix(bd[0],opi);
//         igl::slice(V,opi,1,op);
//     }else{
//         op = V.row(I(0));
//     }
//     auto m = match(Vb,op,-1);
//     int last_i = -1; // last connected boundary's connector index
//     std::vector<int> maskv={1}; // first two must be traced
//     if(I.rows()!=0){
//         Seqv = {B(m.first),I(0)};
//         maskv.push_back(1);
//     }else{
//         Seqv = {B(m.first),bd[0][m.second]};
//         last_i = m.second;
//         maskv.push_back(-1);
//     }
//     // std::cout<<Seqv[0]<<","<<Seqv[1]<<std::endl;

//     int c = 0;
//     K(0) = 1;
//     int t = c;
    
//     while(K.sum()<K.rows()){
//         double x = std::numeric_limits<double>::max();
//         for(int i=0;i<K.rows();i++){
//             if(K(i) == 1) continue;
//             if((center.row(c)-center.row(i)).norm() < x){
//                 x = (center.row(c)-center.row(i)).norm();
//                 t = i;
//             }
//         }
//         Eigen::MatrixXd op1, op2;
//         Eigen::VectorXi bt1, bt2;
//         // after I.rows meaning it's a boundary loop
//         if(c < I.rows()){
//             op1 = V.row(I(c));
//         }else{
//             igl::list_to_matrix(bd[c-I.rows()],bt1);
//             igl::slice(V,bt1,1,op1);
//         }
//         if(t < I.rows()){
//             op2 = V.row(I(t));
//         }else{
//             igl::list_to_matrix(bd[t-I.rows()],bt2);
//             igl::slice(V,bt2,1,op2);
//         }
//         auto p = match(op1,op2,last_i);
//         if(op1.rows()>1){ // meaning op is a boundary loop
//             Seqv.push_back(bt1(p.first));
//             maskv.push_back(1);
//         }
//         if(op2.rows()>1){
//             last_i = p.second;
//             Seqv.push_back(bt2(p.second));
//             maskv.push_back(-1);
//         }else{
//             last_i = -1;
//             Seqv.push_back(I(t));
//             maskv.push_back(1);
//         }
//         K(t) = 1;
//         c = t;
//     }
//     igl::list_to_matrix(Seqv,Seq);
//     igl::list_to_matrix(maskv,mask);
// }

// void split(
//     Eigen::MatrixXd& V,
//     Eigen::MatrixXi& F,
//     const int N,
//     std::vector<int>& path,
//     std::unordered_map<int,int>& m
// ){
//     int num = V.rows();
//     // [make sure path is correct direction]
//     assert(path.size() > 1 && "path invalid");
//     Eigen::VectorXi B; 
//     igl::boundary_loop(F,B);
//     int n_f = F.rows();
//     int n_v = V.rows();
//     int dim = V.cols();

//     int sz,x,n2,n1;
//     // number of edges available
//     sz = path.size()-1;
//     // number of maximum numbers for the first (front) edges
//     x = N / sz + 1;
//     n2 = x*sz - N;
//     n1 = sz - n2; // front * x + back * (x-1) = ptAdd
    
//     V.conservativeResize(V.rows()+N*2,dim);
//     F.conservativeResize(F.rows()+2*(N+sz),3);
//     if(dim == 3)
//         EM.conservativeResize(EM.rows()+2*(N+sz),3);
//     std::vector<int> n_path;
//     // remember the points coord for pts already inserted
//     std::map<std::pair<int,int>,Eigen::MatrixXd> points_added;
//     std::map<std::pair<int,int>,std::vector<int>> points_added_i;

//     for(int i=0;i<path.size()-1;i++){
//         int a = path[i];
//         int b = path[i+1];
//         int t1,t2,t3;
//         n_path.push_back(a);
//         for(int j=0;j<n_f;j++){
//             int index = -1;
//             for(int k=0;k<3;k++){
//                 if((F(j,k)==a    && F(j,(k+1)%3)==b) ||
//                    (F(j,k)==b    && F(j,(k+1)%3)==a) ||
//                    (F(j,k)==m[b] && F(j,(k+1)%3)==m[a]) ||
//                    (F(j,k)==m[a] && F(j,(k+1)%3)==m[b]) ||
//                    (F(j,k)==a    && F(j,(k+1)%3)==m[b]) ||
//                    (F(j,k)==m[b] && F(j,(k+1)%3)==a) ||
//                    (F(j,k)==b    && F(j,(k+1)%3)==m[a]) ||
//                    (F(j,k)==m[a] && F(j,(k+1)%3)==b)){
//                        index = k;
//                        break;
//                 }
//             }
//             if(index == -1) continue;
//             t1 = F(j,index);
//             t2 = F(j,(index+1)%3);
//             t3 = F.row(j).sum()-t1-t2;
//             int n = i>=n1? x-1:x;
//             if(n==0)
//                 break;
//             int start = n_v;
//             std::pair<int,int> c(m[t2],m[t1]);
//             if(points_added.find(c)!=points_added.end()){
//                 for(int k=n-1;k>=0;k--) {
//                     m[n_v] = points_added_i[c][k];
//                     m[points_added_i[c][k]] = n_v;
//                     V.row(n_v++) << points_added[c].row(k);
//                 }
//             }else{
//                 c = std::make_pair(t1,t2);
//                 points_added[c]=Eigen::MatrixXd(n,dim);
//                 for(int k=0;k<n;k++){
//                     V.row(n_v)=(V.row(t2)-V.row(t1)).array()*(k+1)/(n+1)+V.row(t1).array();
//                     points_added[c].row(k)=V.row(n_v);
//                     points_added_i[c].push_back(n_v);
//                     n_path.push_back(n_v);
//                     n_v++;
//                 }
//             }
//             // add face
//             std::vector<int> vl{t1}; // list of vertex to make the faces
//             for(int p=0;p<n;p++){
//                 vl.push_back(start+p);
//             } // t1 x1 x2 x3 x4 t2
//             vl.push_back(t2);
//             for(int k=0;k<n+1;k++){
//                 F.row(n_f)<<t3,vl[k],vl[k+1];
//                 if(dim == 3){
//                     if(k==0)
//                         EM.row(n_f++)<<0,0,1;
//                     else if(k==n)
//                         EM.row(n_f++)<<1,0,0;
//                     else
//                         EM.row(n_f++)<<1,0,1;
//                 }else
//                     n_f++;
//             }
//             F.block(j,0,F.rows()-j-1,3) = F.bottomRows(F.rows()-j-1);
//             F.conservativeResize(F.rows()-1,3);
//             if(dim == 3){
//                 EM.block(j,0,EM.rows()-j-1,3) = EM.bottomRows(EM.rows()-j-1);
//                 EM.conservativeResize(EM.rows()-1,3);
//             }
//             n_f-=1;
//             j-=1;
//         }
//     }
//     n_path.push_back(path.back());
//     path = n_path;
//     V.conservativeResize(n_v,dim);
//     F.conservativeResize(n_f,3); 
//     if(dim == 3)
//         EM.conservativeResize(n_f,3); 
//     assert(n_v-num == 2*N && "split wrong");
// }

// void cut_mesh_to_disk(
//     const Eigen::MatrixXd& V,
//     const Eigen::MatrixXi& F,
//     const Eigen::VectorXi& Seq,
//     const Eigen::VectorXi& mask,
//     Eigen::MatrixXd& Vn,
//     Eigen::MatrixXi& Fn,
//     std::unordered_map<int,int>& m,
//     std::vector<std::vector<int>>& cuts
// ){
//     Fn = F;
//     Vn = V;
//     // connect in order of the index sequence [seq]
//     std::vector<int> E;
//     std::vector<std::vector<int>> Av;

//     // initialize A to be its boundary
//     for(int i=0;i<Seq.rows()-1;i++){
//         if(mask(i)==-1) continue;
//         igl::boundary_loop(Fn,Av);
//         std::vector<int> A;
//         for(int i=0;i<Av.size();i++){
//             A.insert(A.end(),Av[i].begin(),Av[i].end());
//         }
//         for(int j=0;j<A.size();j++){
//             if(A[j] == Seq(i) || A[j] == Seq(i+1)){
//                 A.erase(A.begin()+j);
//                 j--;
//             }
//         }
//         // remove boundary constraints from the avoid list first (add back later)
//         std::vector<int> SQ;
//         igl::matrix_to_list(Seq,SQ);
//         std::vector<int> diff;
//         std::set_difference(A.begin(), A.end(), SQ.begin(), SQ.end(), 
//                         std::inserter(diff, diff.begin()));
//         A = diff;
//         for(int j=i+2;j<Seq.rows();j++) // avoid unconnected constraint points
//             A.push_back(Seq(j));
//         std::cout<<"tracing "<<Seq(i)<<" to "<<Seq(i+1)<<std::endl;
//         std::vector<std::pair<int,int>> Mask;
//         std::set<int> no_enter_f;
//         Eigen::MatrixXi TT;
//         std::set<int> s(A.begin(), A.end());
//         path_tracing(Vn,Fn,std::make_pair(Seq(i),Seq(i+1)),s,Mask,no_enter_f,TT,Vn,Fn,E);
//         std::vector<std::vector<int>> cv;
//         std::vector<std::vector<int>> ck;   
//         igl::cut_mesh(Vn,Fn,{E},cv,ck);
//         for(int g = 0;g<cv.size();g++){
//             if(cv[g].size()==1){
//                 m.insert(std::make_pair(cv[g][0], cv[g][0]));
//             }else {
//                 if(m.find(cv[g][0])!=m.end() || m.find(cv[g][1])!=m.end()){
//                     m[cv[g][0]]=cv[g][1];
//                     m[cv[g][1]]=cv[g][0];
//                 }else{
//                     m.insert(std::make_pair(cv[g][0], cv[g][1]));
//                     m.insert(std::make_pair(cv[g][1], cv[g][0]));
//                 }
//             }
//         }
//         cuts.push_back(E);
//     }
// }

// void embed_constraints(
//     Eigen::MatrixXd& V,
//     Eigen::MatrixXi& F,
//     std::vector<Eigen::MatrixXd>& P_set,
//     const Eigen::MatrixXd& c,
//     const Eigen::VectorXi& ci,
//     const Eigen::MatrixXi& E,
//     std::vector<Eigen::VectorXi>& T,
//     Eigen::VectorXi& T0,
//     Eigen::VectorXi& R,
//     Eigen::MatrixXd& Pn,
//     Eigen::MatrixXi& En
// ){  
//     Eigen::MatrixXi F2;
//     Eigen::MatrixXd V2;
//     std::unordered_map<int,int> match;
//     std::cout<<"polygon is simple? "<<is_simple_polygon(P_set[0])<<std::endl;
//     //std::cout<<P_set[0]<<std::endl;
//     if(P_set.size()>1 || (P_set.size()==1 && is_simple_polygon(P_set[0]))){
//         int max_len = -1;
//         int idx = -1;
//         Eigen::MatrixXd PV;
//         Eigen::MatrixXi PE;
//         Eigen::VectorXi Tn;
//         for(int i=0;i<P_set.size();i++){
//             if(P_set[i].rows()>max_len){
//                 max_len = P_set[i].rows();
//                 idx = i;
//             }
//         }
//         for(int i=0;i<P_set.size();i++){
//             // check orientation of polygon set
//             Eigen::MatrixXd tmpP = P_set[i];
//             if(i != idx && orientation(P_set[i])>0){
//                 for(int j=P_set[i].rows()-1;j>=0;j--){
//                     tmpP.row(tmpP.rows()-j-1)<<P_set[i].row(j);
//                 }
//                 P_set[i] = tmpP;
//                 T[i].reverse();
//             }
//             PV.conservativeResize(PV.rows()+P_set[i].rows(),2);
//             Tn.conservativeResize(Tn.rows()+P_set[i].rows());
//             Tn.bottomRows(P_set[i].rows())<<T[i];
//             PV.bottomRows(P_set[i].rows())<<P_set[i];
//             int o_pe = PE.rows();
//             PE.conservativeResize(PE.rows()+P_set[i].rows(),2);
//             PE.bottomRows(P_set[i].rows())<<Eigen::VectorXi::LinSpaced(P_set[i].rows(),o_pe,o_pe+P_set[i].rows()-1),
//                                             Eigen::VectorXi::LinSpaced(P_set[i].rows(),o_pe+1,o_pe+P_set[i].rows());
            
//             PE.bottomRows(1)<<o_pe+P_set[i].rows()-1,o_pe;
//         }
//         if(c.rows()>0){
//             PV.conservativeResize(PV.rows()+c.rows(),Eigen::NoChange);
//             PV.bottomRows(c.rows())<<c;
//         }
//         for(int i=0;i<Tn.rows();i++){
//             match[Tn(i)] = i;
//         }
//         for(int i=0;i<ci.rows();i++){
//             match[ci(i)] = PV.rows()-c.rows()+i;
//         }
//         // [triangulate P_set]
//         Eigen::MatrixXd H(P_set.size()-1,2);
//         Eigen::MatrixXd jc(P_set.size()-1,2);
//         int ct=0;
//         for(int i=0;i<P_set.size();i++){
//             if(i == idx) continue;
//             // get a position inside P_set[i]
//             Eigen::RowVector3d n;
//             n<<0,0,1;
//             Eigen::RowVector3d eg;
//             eg << (P_set[i].row(1) - P_set[i].row(0)),0;
//             Eigen::RowVector3d perp = eg.cross(n);
//             Eigen::RowVector2d perp2d;
//             perp2d << perp(0),perp(1);
//             double alpha = 1.0;
//             Eigen::RowVector2d x = perp2d * alpha +(P_set[i].row(0)+P_set[i].row(1))/2;
//             while(!point_in_poly(P_set[i],x)){
//                 alpha /= 2;
//                 x = perp2d * alpha +(P_set[i].row(0)+P_set[i].row(1))/2;
//             }
//             jc.row(ct)<<(P_set[i].row(0)+P_set[i].row(1))/2;
//             H.row(ct++)<<x;
//         }
//         igl::triangle::triangulate(PV,PE,H,"YQq33",V2,F2);
//         // igl::opengl::glfw::Viewer viewer;
//         // igl::opengl::glfw::imgui::ImGuiMenu menu;
//         // viewer.plugins.push_back(&menu);
//         // viewer.data().set_mesh(V2,F2);
//         // viewer.data().add_points(PV,Eigen::RowVector3d(0,0,0));
//         // viewer.launch();
//     }else{
//         if(!T.empty()){
//             for(int i=0;i<T[0].rows();i++){
//                 match[T[0](i)] = i;
//             }    
//         }
//         for(int i=0;i<ci.rows();i++){
//             match[ci(i)] = P_set[0].rows()+i;
//         }

//         Eigen::VectorXi cp(c.rows());
//         cp.setConstant(-1);
//         std::vector<std::vector<int>> _L;
//         Shor_van_wyck(P_set[0],R,c,cp,E,"YQq33",V2,F2,_L);
//         // igl::opengl::glfw::Viewer viewer;
//         // igl::opengl::glfw::imgui::ImGuiMenu menu;
//         // viewer.plugins.push_back(&menu);
//         // viewer.data().set_mesh(V2,F2);
//         // viewer.launch();
//     }

//     // [connect objects(point,end points of segment,hole)]
//     Eigen::VectorXi Seq, Seq2;
//     Eigen::VectorXi mk; // whether tracing a path between Seq(i) and Seq(i+1)
//     sequence(V,F,ci,E,P_set,Seq,mk);
//     Seq2 = Seq;
//     std::for_each(Seq2.data(),Seq2.data()+Seq2.size(),[&match](int & s){s=match[s];});

//     Eigen::MatrixXd Vn,V2n;
//     Eigen::MatrixXi Fn,F2n;
//     std::vector<std::vector<int>> cuts1,cuts2;
//     std::unordered_map<int,int> m1,m2;
//     cut_mesh_to_disk(V,F,Seq,mk,Vn,Fn,m1,cuts1);
//     cut_mesh_to_disk(V2,F2,Seq2,mk,V2n,F2n,m2,cuts2);
//     EM.setZero(Fn.rows(),Fn.cols());
//     std::cout<<"finish cut"<<std::endl;
//     //show_bd(V2n,F2n);
//     // [split edges]
//     std::vector<int> drop_extra; // [vertices in between end points of an constrained edge]
//     for(int i=0;i<cuts1.size();i++){
//         if(cuts1[i].size()>cuts2[i].size()){
//             int a = cuts1[i].size()-cuts2[i].size();
//             split(V2n,F2n,a,cuts2[i],m2);
//         }else if(cuts1[i].size()<cuts2[i].size()){
//             //if(cuts1[i].size()!=2){
//                 int a = cuts2[i].size()-cuts1[i].size();
//                 split(Vn,Fn,a,cuts1[i],m1);
//             // }else{
//             //     drop_extra.insert(drop_extra.end(),cuts2[i].begin()+1,cuts2[i].end()-1);
//             // }
//         }
//     }
//     // [extract boundary as new polygon]
//     Eigen::VectorXi B,D;
//     // bool ok1 = is_mesh_valid(Vn,Fn);
//     // bool ok2 = is_mesh_valid(V2n,F2n);
//     //if(!ok1 || !ok2) exit(0);
//     get_ri(V2n,F2n,R);
//     std::cout<<"rotation index sum "<<R.sum()<<std::endl;
//     igl::boundary_loop(F2n,B);
//     Eigen::VectorXi ext;
//     std::vector<int> extv;
//     for(int i=0;i<B.rows();i++){
//         if(std::find(drop_extra.begin(),drop_extra.end(),B(i))==drop_extra.end() &&
//            std::find(drop_extra.begin(),drop_extra.end(),m2[B(i)])==drop_extra.end() )
//             extv.push_back(i);
//     }
//     igl::list_to_matrix(extv,ext);
//     Eigen::VectorXi tB;
//     igl::slice(B,ext,1,tB);
//     igl::boundary_loop(Fn,D);
//     igl::slice(V2n,B,1,Pn);
//     if(Seq2.rows()>0 && P_set.size() == 1){
//         std::cout<<"the last one "<<Seq2.bottomRows(1)<<std::endl;
//         for(int i=0;i<tB.rows();i++)
//             if(tB(i) == Seq2.bottomRows(1)(0))
//                 R(i) = 1;
//     }
//     // [matching two boundaries]
//     //show_bd(Vn,Fn);
//     T0.setZero(Pn.rows());
//     int a = -1, b = -1;
//     if(Seq.rows()>0){
//         for(int i=0;i<tB.rows();i++){
//             if(tB(i) == Seq2(1)){
//                 a = i;
//             }
//             if(D(i) == Seq(1)){
//                 b = i;
//             }
//         }
//         for(int i=0;i<T0.rows();i++){
//             T0((a+i)%T0.rows()) = D((b+i)%D.rows());
//         }
//     }else{
//         if(!T.empty())
//             T0 = T[0];
//         else
//             T0 = D;
//     }
//     //show_bd(V2n,F2n);
//     F = Fn;
//     V = Vn;
//     En = EM;
// }