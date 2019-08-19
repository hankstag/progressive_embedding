#include "progressive_embedding.h"
#include "local_operation.h"
#include "validity_check.h"
#include "local_smooth/local_smooth.h"
#include "local_smooth/auto_grad.hpp"
#include <igl/boundary_loop.h>
#include <igl/slice.h>

#include "plot.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/Timer.h>

#include <limits>
#include <Eigen/Dense>

double domain_area = 0.0;
int n_face = 0;
bool vb = false;

bool is_face_valid(
  const Eigen::Matrix<double,2,3>& G_t,
  const Eigen::Matrix<double,3,2>& T,
  const double threshold
){
  bool flipped = is_face_flipped(T);
  double e = autogen::sd_energy(T,G_t);
  return std::isfinite(e) && (e < threshold) && !flipped;
}

void move_to_center(
  Eigen::MatrixXd& uv,
  const Eigen::MatrixXi& F,
  int i
){
  Eigen::VectorXi bd;
  igl::boundary_loop(F,bd);
  Eigen::RowVector2d center;
  center.setZero();
  for(int i=0;i<bd.rows();i++){
    center += uv.row(bd(i));
  }
  center /= bd.rows();
  uv.row(i) << center;
}

int expand_to_boundary(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::VectorXi& B,
  const Eigen::MatrixXi& dEF_s,
  const int fid,
  Eigen::VectorXi& X
){
  Eigen::MatrixXi dEF = dEF_s;
  std::set<int> N;
  std::deque<int> Q;
  for(int i=0;i<F.rows();i++){
    for(int k=0;k<3;k++){
      int e = F.rows()*((k+2)%3)+i;
      int k_1 = (k+1)%3;
      if(B(F(i,k)) && B(F(i,k_1)))
        dEF(e,1) = -1;
    }
  }
  if(fid != -1){
    Q.push_back(fid);
    N.insert(fid);
    while(!Q.empty()){
      int f = Q.front();
      Q.pop_front();
      for(int i=0;i<3;i++){
        int e = F.rows()*((i+2)%3)+f;
        if(dEF(e,1)!=-1 && N.find(dEF(e,1))==N.end()){
          Q.push_back(dEF(e,1));
          N.insert(dEF(e,1));
        }
      }
    }
  }
  X.resize(N.size());
  int i=0;
  for(int n: N)
    X(i++) = n;
  Eigen::MatrixXi local_F;
  igl::slice(F,X,1,local_F);
  std::cout<<"count interior points for "<<std::endl;
  int in_id= 0;
  int count=0;
  std::set<int> interior;
  for(int i=0;i<local_F.rows();i++){
    for(int k=0;k<3;k++){
      int id = local_F(i,k);
      if(!B(id)){
        in_id = id;
        interior.insert(in_id);
      }
    }
  }
  std::cout<<"#interior point "<<interior.size()<<std::endl;
  return interior.size()==1?in_id:-1;
}

void neighbor_k_ring(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXi& dEF,
    int fid,
    int layer,
    std::set<int>& N // fid collection of neighbors
){
    std::deque<std::pair<int,int>> Q;
    std::pair<int,int> p(fid,0);
    Q.push_back(p);
    N.insert(fid);
    int l = 0;
    while(l<layer && !Q.empty()){
        int f = Q.front().first;
        l = Q.front().second;
        if(l == layer) break;
        Q.pop_front();
        for(int i=0;i<3;i++){
          int e = F.rows()*((i+2)%3)+f;
          if(dEF(e,1)!=-1 && N.find(dEF(e,1))==N.end()){
                std::pair<int,int> p(dEF(e,1),l+1);
                Q.push_back(p);
                N.insert(dEF(e,1));
          }
        }
    }
}

void find_candidate_positions(
  const Eigen::MatrixXd& uv,
  const Eigen::MatrixXi& ring,
  int v,
  int avoid,
  Eigen::MatrixXd& pos
){
  int pi = 0;
  pos.resize(ring.rows()*10,2);
  for(int i=0;i<ring.rows();i++){
    bool common = false;
    for(int j=0;j<3;j++)
      if(ring(i,j) == avoid){
        common = true;
      }
    if(common) continue;
    for(int k=0;k<3;k++){
      if(ring(i,k) == v){
        int x = ring(i,(k+1)%3);
        int y = ring(i,(k+2)%3);
        pos.row(pi++) << uv.row(x);
        pos.row(pi++) << uv.row(y);
        pos.row(pi++) << (uv.row(x) + uv.row(y))/2;
      }
    }
  }
  pos.conservativeResize(pi,2);

}

double calculate_angle_sum(
  Eigen::MatrixXi& F,
  Eigen::MatrixXd& V,
  int v1,
  int avoid
){
  double angle_sum = 0;
  for(int i=0;i<F.rows();i++){
    bool common = false;
    for(int j=0;j<3;j++)
      if(F(i,j) == avoid)
          common = true;
    if(common) continue;
    for(int k=0;k<3;k++){
      if(F(i,k) == v1){
        Eigen::Vector2d a = V.row(F(i,(k+2)%3));
        Eigen::Vector2d b = V.row(v1);
        Eigen::Vector2d c = V.row(F(i,(k+1)%3));
        auto bc = c - b;
        auto ba = a - b;
        angle_sum += std::acos(bc.dot(ba)/(bc.norm()*ba.norm()));
        break;
      }
    }
  }
  return angle_sum;
}

// find a position for vertex v1 along the line (v0 --- t) s.t. the 
// maximum energy of 1-ring neighbor of v1 is minimized and no flips
std::pair<bool,double> flip_avoid_line_search(
  Eigen::MatrixXd& uv,
  const Eigen::MatrixXi& F, // the ring
  int a0,
  int b,
  const Eigen::RowVector2d& x, // initial guess for b is (x+a)/2
  Eigen::RowVector2d& y,       // output position
  double target_area,
  double eps
){
  double t = 1.0;
  int MAX_IT = 75;
  bool valid = true;
  double max_energy = std::numeric_limits<double>::max();
  std::vector<int> neighbor_b;
  for(int i=0;i<F.rows();i++){
    for(int k=0;k<3;k++){
      if(F(i,k) == b){
        neighbor_b.push_back(i);
      }
    }
  }
  Eigen::Matrix<double,2,3> G_t;
  grad_to_eqtri(target_area,G_t);

  for(int it = 0;it<MAX_IT;it++){
    Eigen::RowVector2d pt;
    pt(0) = (1-t)*uv(a0,0)+t*x(0);
    pt(1) = (1-t)*uv(a0,1)+t*x(1);
    t *= 0.8;
    uv.row(b) << pt;
    // calculate the energy of neighbors faces of b
    std::vector<double> E; // one ring energy
    valid = true;
    for(int f: neighbor_b){
        Eigen::Matrix<double,3,2> Tuv;
        Tuv<<uv.row(F(f,0)),uv.row(F(f,1)),uv.row(F(f,2));
        double e = autogen::sd_energy(Tuv,G_t);
        valid = (std::isfinite(e) && !is_face_flipped(Tuv));
        if(!valid) break;
        E.push_back(e);
    }
    if(!valid) continue;
    double em = *(std::max_element(E.begin(),E.end()));
    if(em < max_energy){
        max_energy = em;
        y = pt;
    }
  }
  
  // set invalid elements upper-bound energy 1e20
  if(max_energy > 1e20)  valid = false;
  return std::make_pair(valid,max_energy);
}

using Action = std::tuple<int,int,Eigen::MatrixXi,std::vector<int>>;

void collapse_invalid_elements(
  const Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  Eigen::MatrixXd& uv,
  Eigen::VectorXi& I, //invalid
  const Eigen::VectorXi& B, //boundary
  const double eps,
  const double avg,
  std::vector<Action>& L
){
  std::cout<<"collapsing ... "<<std::endl;
  n_face = F.rows();
  // Build adjacency info
  Eigen::VectorXi EMAP,EE;
  Eigen::MatrixXi E,EF,EI;
  Eigen::MatrixXi dEF,dEI,allE;
  igl::edge_flaps(F,E,allE,EMAP,EF,EI,dEF,dEI,EE);
  
  // Use the same size reference shape through out collapsing
  Eigen::Matrix<double,2,3> G;
  grad_to_eqtri(avg,G);

  // Try collapsing an edge 
  // mark the 1-ring faces as `updated` if succ
  auto collapse_if_valid = [&](
    int f, int k, Eigen::VectorXi& updated
  ){
    int e = F.rows()*((k+2)%3)+f;
    int a = F(f,k), b = F(f,(k+1)%3);
    if(B(a) == 1 || B(b) == 1) return false;
    std::vector<int> N; // 1-ring neighbor exclude the collapsed
    if(edge_collapse_is_valid(e,uv.row(b),uv,F,dEF,dEI,EE,allE,N)){
      // unmark the collapsed faces
      I(f)=0;
      I(dEF(e,1))=0;
      Eigen::MatrixXi ring(N.size(),3);
      for(int i=0;i<N.size();i++){
        ring.row(i) << F.row(N[i]);
        updated(N[i]) = 1;
      }
      updated(f)=0;
      updated(dEF(e,1))=0;
      if(a > b) std::swap(a,b);
      L.push_back(Action(b,a,ring,N));
      collapse_edge(e,uv.row(a),uv,F,dEF,dEI,EE,allE);
      n_face -= 2;
      return true;
    }
    return false;
  };

  int num_invalid = I.sum();
  int layer = 0;
  while(num_invalid!=0){
    bool do_collapse = false;
    Eigen::VectorXi updated = Eigen::VectorXi::Zero(F.rows());
    for(int f=0;f<I.rows();f++){
      if(I(f)==0) continue;
      for(int k=0;k<3;k++){
        do_collapse = (do_collapse || collapse_if_valid(f,k,updated));
      }
    }
    // check and mark invalid elements
    // for updated elements
    for(int i=0;i<I.rows();i++){
      if(updated(i)){
        Eigen::Matrix<double,3,2> T;
        T<<uv.row(F(i,0)),uv.row(F(i,1)),uv.row(F(i,2));
        I(i) = !is_face_valid(G,T,eps);
      }
    }
    num_invalid = I.sum();
    // check_result(uv,F,G);
    // std::cout<<"invalid size "<<num_invalid<<std::endl;

    layer = do_collapse ? 0 : layer+1;
    if(do_collapse) continue;
    
    // if did not collapse anything
    // also collapse its layer 1...k neighbors
    std::set<int> W;
    for(int i=0;i<I.rows();i++){
      if(I(i)==0) continue;
      std::set<int> N;
      neighbor_k_ring(uv,F,dEF,i,layer,N);
      for(int n: N)
        W.insert(n);
    }
    for(int s: W){
      if(B(F(s,0))+B(F(s,1))+B(F(s,2)) < 2)
        I(s) = 1;
    }
  }
  std::cout<<"collapsing done"<<std::endl;
}
bool insert_vertex_back(
  const std::vector<Action>& L,
  const Eigen::VectorXi& B,
  const Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  Eigen::MatrixXd& uv,
  double eps
){
  // for every action in the list
  igl::Timer timer;
  timer.start();
  double total_time = 0.0;
  int local_itr = 0;
  for(int ii=0;ii<L.size();ii++){
    
    double avg = domain_area / n_face;
    
    double time1 = timer.getElapsedTime();
    std::cerr<<"insert back "<<ii<<"/"<<L.size()<<" ";
    //std::cout<<std::get<0>(L[ii])<<"(n) <-> "<<std::get<1>(L[ii])<<"(o)"<<std::endl;
    
    Action ac = L[ii];
    Eigen::MatrixXi ring = std::get<2>(ac);
    std::vector<int> nbs = std::get<3>(ac);
    
    bool found = false;
    // position candidates
    Eigen::RowVector2d pos;
    double max_energy = std::numeric_limits<double>::max();

    int v0 = std::get<1>(ac); // v0 -> existing vertex
    int v1 = std::get<0>(ac); // v1 -> new vertex
    
    // save faces for insertion failure
    auto uv_store = uv;
    auto F_store = F;

    // restore collapsed
    uv.row(v1) = uv.row(v0);
    for(int j=0;j<nbs.size();j++)
      F.row(nbs[j])<<ring.row(j);
    // pick the valid sector with angle less than PI
    double angle_sum_of_v0 = calculate_angle_sum(ring,uv,v0,v1);
    if(angle_sum_of_v0 < igl::PI)
      std::swap(v0,v1);
    
    Eigen::MatrixXd cd;
    find_candidate_positions(uv,ring,v1,v0,cd);
    //std::cout<<"try pos: "<<cd.rows()<<std::endl;
    for(int j=0;j<cd.rows();j++){
        Eigen::RowVector2d q;
        Eigen::RowVector2d x;
        x << cd(j,0),cd(j,1);
        auto r = flip_avoid_line_search(uv,ring,v0,v1,x,q,avg,eps);
        bool succ  = std::get<0>(r);
        double e_m = std::get<1>(r);
        found = (found || succ);
        if(succ && max_energy > e_m){
          max_energy = e_m;
          pos = q;
          break;
        }
    }

    auto drop_empty_faces = [](const Eigen::MatrixXi& F, Eigen::MatrixXi& Fn){
      Fn = F;
      int k=0;
      for(int i=0;i<F.rows();i++){
        if(F.row(i).sum()!=0)
          Fn.row(k++) << F.row(i);
      }
      Fn.conservativeResize(k,3);
    };

    if(found){
      n_face += 2;
      if(vb)
        std::cout<<"current face size "<<n_face<<"/"<<F.rows()<<std::endl;
      uv.row(v1) << pos;
      Eigen::MatrixXi Ft;
      drop_empty_faces(F,Ft);
      local_smoothing(V,Ft,B,uv,10,1e10,avg);
      local_itr = 0; // reset record
    }else{
      F = F_store;
      uv = uv_store;
      Eigen::MatrixXi Ft;
      drop_empty_faces(F,Ft);
      local_smoothing(V,Ft,B,uv,100,1e10,avg);
      ii--;
      local_itr++;
    }

    double time2 = timer.getElapsedTime();
    total_time += (time2 - time1);
    double expect_total_time = (total_time / (std::max(1,ii+1))) * L.size();
    int minute = int((expect_total_time-total_time)/60.0);
    int second = (expect_total_time-total_time) - 60*minute;
    std::cerr<<"expect time: "<<minute<<" mins "<<second<<" seconds";
    std::cerr<<'\r';
    if(local_itr > 50)
      return false;
  }
  std::cout<<std::endl;
  return true;
}

void check_result(
  const Eigen::MatrixXd& uv,
  const Eigen::MatrixXi& F,
  const Eigen::Matrix<double,2,3>& G
){
  // [ check result ]
  Eigen::VectorXd A;
  Eigen::MatrixXi Fn=F;
  int k=0;
  for(int i=0;i<Fn.rows();i++){
    if(F.row(i).sum()!=0)
      Fn.row(k++) << F.row(i);
  }
  Fn.conservativeResize(k,3);
  igl::doublearea(uv,Fn,A);
  std::cout<<"max area: "<<A.maxCoeff()<<std::endl;
  std::cout<<"min area: "<<A.minCoeff()<<std::endl;
  Eigen::VectorXd E;
  E.setZero(Fn.rows());
  for(int i=0;i<Fn.rows();i++){
    Eigen::Matrix<double,3,2> T;
    T<<uv.row(Fn(i,0)),uv.row(Fn(i,1)),uv.row(Fn(i,2));
    E(i) = autogen::sd_energy(T,G);
  }
  std::cout<<"max energy "<<E.maxCoeff()<<std::endl;
  Eigen::VectorXi T;
  flipped_elements(uv,Fn,T);
  std::cout<<"flipped: "<<T.sum()<<std::endl;
}

bool progressive_embedding(
  const Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  Eigen::MatrixXd& uv,
  const Eigen::VectorXi& bi,
  const Eigen::MatrixXd& b,
  double eps,
  bool verbose
){

  // [ Reference shape area for Dirichlet Energy ]
  double area_total = .0;
  Eigen::VectorXd area;
  igl::doublearea(uv,F,area);
  double avg = 0;
  for(int k=0;k<area.rows();k++){
    if(std::isfinite(area(k))){
      area_total += area(k);
    }
  }
  domain_area = area_total;
  avg = area_total / F.rows();

  // [ Mark boundary vertices ]
  Eigen::VectorXi B;
  B.setZero(V.rows());
  for(int i=0;i<bi.rows();i++)
    B(bi(i)) = 1;

  // [ collect invalid elements ]
  Eigen::VectorXd E(F.rows());
  Eigen::Matrix<double,2,3> G;
  grad_to_eqtri(avg,G);
  for(int i=0;i<F.rows();i++){
    Eigen::Matrix<double,3,2> T;
    T<<uv.row(F(i,0)),uv.row(F(i,1)),uv.row(F(i,2));
    E(i) = autogen::sd_energy(T,G);
  }
  Eigen::VectorXi I;
  flipped_elements(uv,F,I);
  for(int i=0;i<F.rows();i++){
    if(!std::isfinite(E(i))||E(i)>eps)
      I(i)=1;
  }
  std::cerr<<"#total invalid: "<<I.sum()<<std::endl;
  std::vector<Action> L; // list of collapse operation stored
  collapse_invalid_elements(V,F,uv,I,B,eps,avg,L);

  check_result(uv,F,G);

  // [ invert vertex in reverse order of L ]
  std::reverse(L.begin(),L.end());
  bool succ = insert_vertex_back(L,B,V,F,uv,eps);
  return succ;
}
