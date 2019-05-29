#include "shor.h"
#include <igl/PI.h>
#include <igl/doublearea.h>
#include <igl/boundary_loop.h>
#include <igl/adjacency_list.h>

#include <igl/copyleft/cgal/orient2D.h>
#include <igl/copyleft/cgal/ear_clipping.h>
#include "is_simple_polygon.h"
#include <iostream>

using namespace std;

short orientation(const Eigen::Matrix<double,3,2>& P){
    double a[2] = {P(0, 0), P(0, 1)};
    double b[2] = {P(1, 0), P(1, 1)};
    double c[2] = {P(2, 0), P(2, 1)};
    return igl::copyleft::cgal::orient2D(a, b, c);
}

void set_rotation_index(
   const Eigen::MatrixXd& uv,
   const Eigen::MatrixXi& F,
   Eigen::VectorXi& R,
   int offset
){
   Eigen::VectorXi bd;
   igl::boundary_loop(F,bd);
   R.setZero(bd.rows());
   std::vector<std::vector<int>> A;
   igl::adjacency_list(F,A,true);

   // for every boundary vertex, update its rotation index
   for(int i=0;i<bd.rows();i++){
       int v = bd((i+offset)%bd.rows());
       Angle sum = Angle();
       std::reverse(A[v].begin(),A[v].end());
       for(int j=0;j<A[v].size()-1;j++){
           int id = A[v][j];
           int id_1 = A[v][j+1];
           Angle I(uv.row(id),uv.row(v),uv.row(id_1)); // w1, v, w2
           sum = (j==0)? I: sum+I;
       }
       R(i) = sum.r;
   }
}

void add_triangle(
  Eigen::MatrixXi& F, 
  int i, int j, 
  std::vector<std::vector<int>>& K,
  std::vector<std::vector<int>>& Q
){
  if(K[i][j]==-1) return;
  int k = K[i][j];
  F.conservativeResize(F.rows()+1,3);
  F.bottomRows(1) << i,k,j;
  add_triangle(F,i,k,K,Q);
  add_triangle(F,k,j,K,Q);
}

bool weakly_self_overlapping(
  const Eigen::MatrixXd& P,
  const Eigen::VectorXi& R,
  Eigen::MatrixXi& F
){
  std::cout<<"weakly self overlapping test"<<std::endl;
  auto add_to_table = [](
    const Eigen::MatrixXd& P,
    const Eigen::VectorXi& R,
    std::vector<std::vector<Angle>>& FA,
    std::vector<std::vector<Angle>>& LA,
    int i, int k, int j
  ){
    const int N = P.rows();
    Eigen::Matrix<double,3,2> tri;
    tri<<P.row(i),P.row(k),P.row(j);
    if(orientation(tri)<=0) return false;
    
    Angle I(tri.row(2),tri.row(0),tri.row(1)); // J I K
    Angle J(tri.row(1),tri.row(2),tri.row(0)); // K J I
    
    Angle k1 = LA[i][k];        
    Angle k2(P.row(i),P.row(k),P.row(j));
    Angle k3 = FA[k][j];
    Angle Fr = I + FA[i][k];
    Angle La = LA[k][j] + J;
    if((Fr.r <= R(i)) && (La.r <= R(j)) && (k1+k2+k3).r == R(k)){
      FA[i][j] = Fr;
      LA[i][j] = La;
      return true;
    }else
      return false;
  };

  const int N = P.rows();
  std::vector<std::vector<int>> K(N, std::vector<int>(N,-1));
  std::vector<std::vector<int>> Q(N, std::vector<int>(N,0));
  std::vector<std::vector<Angle>> FA(N, std::vector<Angle>(N));
  std::vector<std::vector<Angle>> LA(N, std::vector<Angle>(N));
  bool succ = false;
  int i;
  for(int i=0;i<N;i++){
    Q[i][(i+1)%N]=1;
    FA[i][(i+1)%N] = Angle(P.row((i+1)%N),P.row(i),P.row((i+1)%N));
    LA[i][(i+1)%N] = Angle(P.row(i),P.row((i+1)%N),P.row(i));
  }
    int h = -1;
    std::vector<std::vector<double>> A(N, std::vector<double>(N,-1));
  for(int d=2;d<N && !succ;d++){
    for(i=0;i<N && !succ;i++){
      int j = (i + d)%N;
      int k = (i + 1)%N;
      double min_q = -1.0f;
      while (k!=j){
        if (Q[i][k]==1 && Q[k][j]==1){
          if(add_to_table(P,R,FA,LA,i,k,j)){
            Q[i][j] = 1;
            K[i][j] = k;
            if(i == (j+1)%N){
              succ = true;
              h = i;
              break;
            }
          }
        }
        k = (k + 1)%N;
      }
    }
  }
  std::cout<<"test done"<<std::endl; 
  if(h == -1) return false;
  add_triangle(F,h,(h-1+N)%N,K,Q);
  return true;
}

void drop_colinear(
  const Eigen::MatrixXd& P,
  const Eigen::VectorXi& R,
  Eigen::VectorXi& B,
  Eigen::MatrixXd& mP,
  Eigen::VectorXi& mR
){
  int dropped = 0;
  int N = P.rows();
  std::vector<int> Bv;
  for(int i=0;i<N;i++){
    int a = (i-1+N)%N;
    int v = i;
    int b = (i+1)%N;
    Eigen::Matrix<double,3,2> T;
    T<<P.row(a),P.row(v),P.row(b);
    if(R(v)!=0 || orientation(T)!=0){ // non-colinear or rotate index nonzero
      Bv.push_back(v);
    }else
      dropped++;
  }
  igl::list_to_matrix(Bv,B);
  igl::slice(P,B,1,mP);
  igl::slice(R,B,mR);
}

void add_colinear(
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXi& nF,
  const Eigen::VectorXi& B,
  Eigen::MatrixXi& F
){
  int added = 0;
  int sN = nF.maxCoeff()+1; // size of simplified polygon
  F.resizeLike(nF);
  for(int i=0;i<F.rows();i++)
    F.row(i)<<B(nF(i,0)),B(nF(i,1)),B(nF(i,2));
  typedef std::pair<int,int> Edge;
  std::vector<Edge> stk;
  std::set<std::pair<int,int>> sid; // edge vertex pairs needs split
  for(int i=0;i<nF.rows();i++){
    bool split = false;
    for(int k=0;k<3;k++){
      int a = nF(i,k);
      int b = nF(i,(k+1)%3);
      if(  (b-a==1 || a-b == sN-1) &&
        (B(b)-B(a)!=1 && B(a)-B(b)!=P.rows()-1)){ // edge need split
        if(!split){ // do not split a face twice (do it later)  
          stk.push_back(Edge(i,k));
          split = true;
        }
        sid.insert(std::make_pair(B(a),B(b)));
      }
    }
  }
  while(!stk.empty()){
    Edge e = stk.back();
    stk.pop_back();
    int a = F(e.first,e.second);
    int b = F(e.first,(e.second+1)%3);
    int c = F(e.first,(e.second+2)%3);
    int rg = a < b ? b-a : b+P.rows()-a;
    int r = F.rows();
    F.conservativeResize(F.rows()+rg,3);
    added += (rg-1);
    // from r+0 to r+rg-1
    for(int i=0;i<rg;i++){
      F.row(r+i)<<(a+i)%P.rows(),(a+i+1)%P.rows(),c;
    }
    // if bc needs split
    if(sid.find(std::make_pair(b,c))!=sid.end()){
      stk.push_back(Edge(r+rg-1,1)); // its now 2nd edge of face r+rg-1
    }
    // if ac needs split
    if(sid.find(std::make_pair(c,a))!=sid.end()){
      stk.push_back(Edge(r,2)); // its now 3rd edge of face r
    }
    F.row(e.first) << -1,-1,-1; // delete face
  }
  Eigen::MatrixXi tF=F;
  // drop the (-1,-1,-1) rows
  int k=0;
  for(int i=0;i<F.rows();i++)
    if(F.row(i).sum()!=-3)
      tF.row(k++)<<F.row(i);
  tF.conservativeResize(k,3);
  F = tF;
}

void subdivide_polygon(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  std::vector<std::vector<int>>& L
){
  // collect halfedge and corresponding face id
  // initialize L to be collection of all faces(polygons)
  std::map<std::pair<int,int>, int> H;
  Eigen::VectorXi G(F.rows());
  Eigen::MatrixXi FF,FFI;
  igl::triangle_triangle_adjacency(F,FF,FFI);
  L.resize(F.rows());
  for(int i=0;i<F.rows();i++){
    G(i) = i; // every face belongs to itself
    for(int k=0;k<3;k++){
      L[i].push_back(F(i,k));
      H[std::make_pair(F(i,k),F(i,(k+1)%3))] = FF(i,k);
    }
  }

  // traverse all the halfedges
  for(auto h: H){ // [he, pid]
    auto he = h.first;
    auto rhe = std::make_pair(he.second,he.first);
    int p1 = H[he]; // a -> b
    int p2 = H[rhe]; // b -> a
    if(p1 == -1 || p2 == -1) continue;
    
    // up to group root
    while(p1!=G[p1])
      p1 = G[p1];
    while(p2!=G[p2])
      p2 = G[p2];

    // combine p1 and p2
    Eigen::MatrixXd poly(L[p1].size()+L[p2].size()-2,2);
    auto a = std::find(L[p1].begin(),L[p1].end(),he.first);
    auto b = std::find(L[p2].begin(),L[p2].end(),he.second);

    std::vector<int> L1(L[p1].begin(),a);
    std::vector<int> R1(a,L[p1].end());
    std::vector<int> L2(L[p2].begin(),b);
    std::vector<int> R2(b,L[p2].end());
    
    std::vector<int> S;
    S = L1;
    auto c = R2.empty() ? R2.begin() : R2.begin()+1;
    auto d = R1.empty() ? R1.begin() : R1.begin()+1;
    S.insert(S.end(),c,R2.end());
    S.insert(S.end(),L2.begin(),L2.end());
    S.insert(S.end(),d,R1.end());
      
    for(int i=0;i<poly.rows();i++){
      poly.row(i)<<V.row(S[i]);
    }
    
    // if merged polygon is simple, drop edge he/rhe
    // erase L[p2], add to L[p1]
    if(is_simple_polygon(poly)){
      H[he]=-1;
      H[rhe]=-1;
      G[p2] = p1; // p2 belongs to p1 now
      L[p1] = S;
      L[p2].clear();
    }
  }
}

void simplify_triangulation(
  const Eigen::MatrixXd& V_i,
  const std::vector<std::vector<int>>& L_i,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F
){
  std::vector<std::vector<int>> L;
  V = V_i;
  F.resize(0,0);
  std::map<std::pair<int,int>, std::vector<int>> E;
  // calculate the average edge length
  double avl = 0.0f;
  for(int i=0;i<V_i.rows();i++){
    avl += (V_i.row((i+1)%V_i.rows())-V_i.row(i)).norm();
  }
  avl /= V_i.rows();
  // traverse L, split the internal edges
  for(int i=0;i<L_i.size();i++){
    if(!L_i[i].empty()){
      L.resize(L.size()+1); // conservativeResize
    }
    for(int j=0;j<L_i[i].size();j++){
      int a = L_i[i][j];
      int b = L_i[i][(j+1)%L_i[i].size()];
      L.back().push_back(a);
      if(E.find(std::make_pair(b,a))!=E.end()){
        auto tv = E[std::make_pair(b,a)]; // tmp vec
        std::reverse(tv.begin(),tv.end());
        E[std::make_pair(a,b)] = tv;
        L.back().insert(L.back().end(),tv.begin(),tv.end());
        continue;
      }
      int n = std::max((V.row(b)-V.row(a)).norm() / avl,1.0);
      int vn = V.rows();
      if(b-a!=1 && a-b!=V_i.rows()-1){
        V.conservativeResize(V.rows()+n-1,2);
        // split to n edges
        for(int k=0;k<n-1;k++){
          V.row(vn+k)<<V.row(a) + (V.row(b)-V.row(a))*(k+1)/n;
          E[std::make_pair(a,b)].push_back(vn+k);
          L.back().push_back(vn+k);
        }
      }
    }
  }
  // triangulate subpolygons
  Eigen::MatrixXd LP; // local polygon 
  Eigen::MatrixXi LF; // local faces
  Eigen::MatrixXd LV;
  for(int i=0;i<L.size();i++){
    LP.resize(L[i].size(),2);
    for(int j=0;j<LP.rows();j++){
      LP.row(j)<<V.row(L[i][j]);
    }
    //display(LP,0,LP.rows()-1);
    Eigen::MatrixXi LE(LP.rows(),2);
    LE<<Eigen::VectorXi::LinSpaced(LP.rows(),0,LP.rows()-1),
        Eigen::VectorXi::LinSpaced(LP.rows(),1,LP.rows());
    LE(LE.rows()-1,1) = 0;
    igl::triangle::triangulate(LP,LE,Eigen::MatrixXd(),"YQq33",LV,LF);    
    Eigen::MatrixXd nV = LV.bottomRows(LV.rows()-LP.rows());
    int n_o = V.rows();
    V.conservativeResize(V.rows()+nV.rows(),2);
    V.bottomRows(nV.rows())<<nV;
    for(int f=0;f<LF.rows();f++){
      for(int k=0;k<3;k++){
        if(LF(f,k) >= L[i].size())
          LF(f,k) += (n_o-LP.rows());
        else
          LF(f,k) = L[i][LF(f,k)];
      }
    }
    F.conservativeResize(F.rows()+LF.rows(),3);
    F.bottomRows(LF.rows())<<LF;
  }
}

// the re-implementation of Shor algorithm
bool Shor_van_wyck(
  const Eigen::MatrixXd& P,
  const Eigen::VectorXi& R,
  const std::string flags,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  bool do_refine
){

  // [drop colinear points]
  Eigen::VectorXi B; // remaining vertices: non-colinear/rotateindex!=0
  Eigen::MatrixXd mP; // P \ colinear vertices
  Eigen::VectorXi mR;
  drop_colinear(P,R,B,mP,mR);

  // [ear clipping]
  Eigen::VectorXi D;
  Eigen::MatrixXi eF;
  Eigen::MatrixXd nP;
  Eigen::VectorXi nR;
  igl::copyleft::cgal::ear_clipping(mP,mR,D,eF,nP);
  igl::slice(mR,D,1,nR);

  // [weakly-self-overlapping test]
  Eigen::MatrixXi nF;
  bool succ = (nP.rows()==0) || weakly_self_overlapping(nP,nR,nF);
  if(!succ){
    std::cout<<"shor failed"<<std::endl;
    exit(0);
    return false;
  }
  // [map simplified index to initial polygon]
  for(int i=0;i<nF.rows();i++){
    nF.row(i) << D(nF(i,0)),D(nF(i,1)),D(nF(i,2));
  }
  if(eF.rows()>0){
    nF.conservativeResize(nF.rows()+eF.rows(),3);
    nF.block(nF.rows()-eF.rows(),0,eF.rows(),3) = eF;
  }

  // [add back colinear vertices by spliting boundary edges]
  add_colinear(P,nF,B,F);
  if(!do_refine){
    V = P;
    return true;
  }
  V = P;

  // [simplify mesh (subdivide into small polygons)]
  std::vector<std::vector<int>> L;
  subdivide_polygon(V,F,L);

  // [refine each small polygon]
  Eigen::MatrixXd V0 = V;
  simplify_triangulation(V0,L,V,F);
  return true;
}
