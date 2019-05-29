#include "decompose_polygon.h"
#include "is_simple_polygon.h"
#include "embed_points.h"
#include <igl/triangle_triangle_adjacency.h>
#include <igl/copyleft/cgal/orient2D.h>

bool is_convex(
  const Eigen::MatrixXd& P
){
  for(int i=0;i<P.rows();i++){
    int prev = (i-1+P.rows())%P.rows();
    int next = (i+1)%P.rows();
    double a[2] = {P(prev,0),P(prev,1)};
    double b[2] = {P(i,0),P(i,1)};
    double c[2] = {P(next,0),P(next,1)};
    short r = igl::copyleft::cgal::orient2D(a,b,c);
    if(r < 0)
        return false;
  }
  return true;
}

void merge_triangles(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  std::vector<std::vector<int>>& L
){

  // use the greedy method for now
  // could be improved

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
    if(is_simple_polygon(poly) && is_convex(poly)){
      H[he]=-1;
      H[rhe]=-1;
      G[p2] = p1; // p2 belongs to p1 now
      L[p1] = S;
      L[p2].clear();
    }
  }
}

void decompose_polygon(
  const Eigen::MatrixXd& P,
  const Eigen::VectorXi& R,
  const Eigen::MatrixXd& C,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F, 
  std::vector<std::vector<int>>& L
){
  bool succ = Shor_van_wyck(P,R,"",V,F,false);
  assert(succ && "Shor failed");
  embed_points(C,V,F);
  merge_triangles(V,F,L);
}
