#include "progressive_embedding.h"
#include <Eigen/Dense>

#include "local_operation.h"
#include <queue>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/cgal/orient2D.h>
#include "plot.h"

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include "local_smooth/local_smooth.h"
#include <math.h>
#include <igl/slim.h>
#include <igl/writeOFF.h>
#include <igl/Timer.h>
#include <igl/linprog.h>
#include <limits>
#include "local_smooth/auto_grad.hpp"
#define LOG

int expand_to_boundary(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::VectorXi& I,
  const Eigen::MatrixXi& FF_s,
  const int fid,
  Eigen::MatrixXi& local_F
){
  Eigen::MatrixXi FF = FF_s;
  std::set<int> N;
  std::deque<int> Q;
  N.clear();
  int level = 0;
  for(int i=0;i<FF.rows();i++){
    for(int k=0;k<3;k++){
      int k_1 = (k+1) % 3;
      if(I(F(i,k)) && I(F(i,k_1)))
        FF(i,k) = -1;
    }
  }
  if(fid != -1){
      Q.push_back(fid);
      N.insert(fid);
      while(!Q.empty()){
          int f = Q.front();
          Q.pop_front();
          for(int i=0;i<3;i++){
            if(FF(f,i)!=-1 && N.find(FF(f,i))==N.end()){
              Q.push_back(FF(f,i));
              N.insert(FF(f,i));
            }
          }
      }
  }
  Q.clear();
  local_F.resize(N.size(),3);
  int i=0;
  for(int x: N){    
      local_F.row(i++)<<F.row(x);
  }
  std::cout<<"count interior points"<<std::endl;
  int my_id = 0;
  for(int i=0;i<local_F.rows();i++){
    for(int k=0;k<3;k++){
      int id = local_F(i,k);
      if(!I(id)){
        my_id = id;
        //std::cout<<my_id<<std::endl;
      }
    }
  }
  return my_id;
}

void get_neighbor_at_level(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXi& FF,
    int fid,
    int neighbor,
    Eigen::VectorXi& nbs,
    Eigen::MatrixXi& region
){
    //std::cout << "Hello\n";
    std::set<int> N;
    std::deque<std::pair<int,int>> Q;
    N.clear();
    int level = 0;
    if(fid != -1){
        std::pair<int,int> p(fid,0);
        Q.push_back(p);
        N.insert(fid);
        int l = 0;
        while(l<neighbor && !Q.empty()){
            int f = Q.front().first;
            l = Q.front().second;
            if(l == neighbor) break;
            Q.pop_front();
            for(int i=0;i<3;i++){
                if(FF(f,i)!=-1 && N.find(FF(f,i))==N.end()){
                    std::pair<int,int> p(FF(f,i),l+1);
                    Q.push_back(p);
                    N.insert(FF(f,i));
                }
            }
        }
    }
    Q.clear();
    region.resize(N.size(),3);
    int i=0;
    nbs.resize(N.size());
    //std::cout<<"neighbor num "<<region.rows()<<std::endl;
    for(int x: N){
        nbs.row(i)<<x; 
        region.row(i++)<<F.row(x);
    }
}

template <typename DerivedT>
double doublearea(const Eigen::MatrixBase<DerivedT>& T){
    if(T.cols()==2){
        Eigen::Matrix<double,1,2> v01 = T.row(1) - T.row(0);
        Eigen::Matrix<double,1,2> v12 = T.row(2) - T.row(1);
        Eigen::Matrix<double,1,2> v20 = T.row(0) - T.row(2);
        double a = v01.norm();
        double b = v12.norm();
        double c = v20.norm();
        double s = (a+b+c)/2.;
        return std::sqrt(double(s*(s-a)*(s-b)*(s-c)));
    }else{
        Eigen::Matrix<double,1,3> v01,v12,v20;
        v01 << T(1,0)-T(0,0),T(1,1)-T(0,1),0;
        v12 << T(2,0)-T(1,0),T(2,1)-T(1,1),0;
        v20 << T(0,0)-T(2,0),T(0,1)-T(2,1),0;
        Eigen::Matrix<double,1,3> n = v01.cross(v12);
        return n.norm();
    }
}
bool is_flipped(const Eigen::MatrixXd& V3d, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const std::vector<int>& I){
    // check flip
    int flipped = 0;
    for(int i: I){
        int x = F(i,0);
        int y = F(i,1);
        int z = F(i,2);
        double a[2] = {double(V(x, 0)), double(V(x, 1))};
        double b[2] = {double(V(y, 0)), double(V(y, 1))};
        double c[2] = {double(V(z, 0)), double(V(z, 1))};
        short k = igl::copyleft::cgal::orient2D(a, b, c);
        if(k <= 0) {
            flipped++;
        }
    }
    if(flipped>0)
        return true;
    else
        return false;
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
            if(F(i,j) == avoid){
                common = true;
            }
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

// find a position for b that maximize the minarea in F
std::pair<bool,double> flip_avoid_line_search(
    const Eigen::MatrixXd& V,
    Eigen::MatrixXd& uv,
    const Eigen::MatrixXi& F,
    int a0,
    int b,
    const Eigen::RowVectorXd& x, // initial guess for b is (x+a)/2
    Eigen::RowVector2d& y, // output position
    double target_area
){
    double t = 1.0;
    int MAX_IT = 75;
    bool valid = false;
    double maxenergy = std::numeric_limits<double>::max();
    //Eigen::VectorXi dirty;
    std::vector<int> dirty;
    for(int i=0;i<F.rows();i++){
        for(int k=0;k<3;k++){
            if(F(i,k) == b){
                dirty.push_back(i);
            }
        }
    }
    std::vector<double> one_ring_chosen;
    double h = std::sqrt( target_area/sin(M_PI / 3.0));
    Eigen::Matrix<double,3,3> Tx;
    Tx<<0,0,0,h,0,0,h/2.,(std::sqrt(3)/2.)*h,0;
    Eigen::Matrix<double,3,3> gx;
    grad_operator(Tx,gx);
    Eigen::Matrix<double,2,3> G_t = gx.topRows(2);
    
    for(int it = 0;it<MAX_IT;it++){
        Eigen::RowVector2d pt;
        pt(0) = (1-t)*uv(a0,0)+t*x(0);
        pt(1) = (1-t)*uv(a0,1)+t*x(1);
        t *= 0.8;
        uv.row(b) << pt;
        // calculate the energy of dirty faces
        std::vector<double> one_ring_energy;
        for(int f: dirty){
            Eigen::Matrix<double,3,3> T;
            Eigen::Matrix<double,3,2> Tuv;
            for(int t=0;t<3;t++){
                // make sure vertex v is always the first in T
                Tuv.row(t) << uv.row(F(f,t));
                T.row(t) << V.row(F(f,t));
            }
            one_ring_energy.push_back(autogen::sd_energy(Tuv,G_t));
        }

        // }
        
        double max_one_ring_energy = 0.0;
        bool bad_element = false;
        for(double z: one_ring_energy){
            if(std::isnan(z) || std::isinf(z))
                bad_element = true;
            if(z > max_one_ring_energy)
                max_one_ring_energy = z;
        }
        if(is_flipped(V,uv,F,dirty) || bad_element) continue;
        if(std::isnan(max_one_ring_energy) || std::isinf(max_one_ring_energy))
            valid = false;
        else 
            valid = true;
        // // get max energy
        // if(std::isnan(max_one_ring_energy)){
        //     std::cout<<"its nan\n";
        //     exit(0);
        // }
        // if(std::isinf(max_one_ring_energy)){
        //     std::cout<<"its inf\n";
        //     exit(0);
        // }
        if(max_one_ring_energy < maxenergy){
            maxenergy = max_one_ring_energy;
            one_ring_chosen = one_ring_energy;
            y = pt;
        }
        // Eigen::VectorXd A;
        // igl::doublearea(uv,F,A);
        // if(minarea<A.minCoeff()){
        //     minarea = A.minCoeff();
        //     
        // }
    }
    return std::make_pair(valid,maxenergy);
}

using Action = std::tuple<int,int,Eigen::MatrixXi,Eigen::VectorXi,bool>;

void collect_invalid_elements(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& uv,
  const double eps,
  const double avg,
  Eigen::VectorXi& I
){
  // [init Dirichlet energy]
  std::cout<<"init dirichlet energy ... "<<std::endl;
  Eigen::VectorXd Energy;
  std::vector<int> FL(F.rows());
  std::iota(FL.begin(),FL.end(),0);
  compute_dirichlet_energy(avg,uv,F,FL,eps,Energy);
  std::cout<<"done"<<std::endl;
  
  // [collect invalid elements]
  std::cout<<"collect invalid elements"<<std::endl;
  I.setZero(F.rows());
  count_flipped_element(uv,F,I);
  for(int i=0;i<F.rows();i++){
    if(I(i)) continue;
    if(!std::isfinite(Energy(i))||Energy(i)>eps)
      I(i) = 1;
  }
}

void collapse_invalid_elements(
  const Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  Eigen::MatrixXd& uv,
  const Eigen::VectorXi& I,
  const Eigen::VectorXi& B,
  const double eps,
  const double avg,
  std::vector<Action>& L
){
  std::vector<int> FL(F.rows());
  std::iota(FL.begin(),FL.end(),0);
  Eigen::VectorXi EMAP,EE;
  Eigen::MatrixXi E,EF,EI;
  Eigen::MatrixXi dEF,dEI,allE;
  igl::edge_flaps(F,E,allE,EMAP,EF,EI,dEF,dEI,EE);
  std::vector<int> S;
  for(int i=0;i<I.rows();i++){
    if(I(i))
      S.push_back(i);
  }
  while(S.size()!=0){
    bool do_collapse = false;
    for(int f: S){
      typedef std::tuple<double,int,bool> edge_info;
      std::vector<edge_info> f_info;
      for(int k=0;k<3;k++){
        double l = (uv.row(F(f,k))-uv.row(F(f,(k+1)%3))).norm();
        bool on_bd = (B(F(f,k) == 1) || B(F(f,(k+1)%3)) == 1);
        f_info.push_back(edge_info(l,k,on_bd));
      }
      std::sort(f_info.begin(),f_info.end());
      // collapse it if valid
      for(auto info: f_info){
        int k = std::get<1>(info);
        int e = F.rows()*((k+2)%3)+f;
        int a = F(f,k), b = F(f,(k+1)%3);
        std::vector<int> N;
        //std::cout<<"try collpasing "<<a<<" and "<<b<<std::endl;
        if(std::get<2>(info)) {
            //std::cout<<"its on bd\n";
            continue;
        }
        if(edge_collapse_is_valid(e,uv.row(b),uv,F,dEF,dEI,EE,allE,N)){
            std::sort(N.begin(), N.end());
            N.erase( std::unique( N.begin(), N.end() ), N.end() );
            Eigen::VectorXi nbi;
            igl::list_to_matrix(N,nbi);
            Eigen::MatrixXd nb(N.size()*3,2);
            Eigen::MatrixXi F_store(nbi.rows(),3);
            // collect all the positions beside a and b
            int fn = 0;
            for(int fi: N){
                if(F.row(fi).sum()==0) continue;
                F_store.row(fn++)<<F.row(fi);
            }
            F_store.conservativeResize(fn,3);
            Eigen::RowVector2d pt;
            if(B(a) == 1){
                // if a is boundary collapse to a
                pt = uv.row(a);
                L.push_back(Action(a,b,F_store,nbi,std::get<2>(info)));
            }else if(B(b) == 1){
                // if b is boundary collapse to b
                pt = uv.row(b);
                L.push_back(Action(b,a,F_store,nbi,std::get<2>(info)));
            }else{
                // otherwise collapse to a small index
                pt = uv.row(std::min(a,b));
                L.push_back(Action(std::max(a,b),std::min(a,b),F_store,nbi,std::get<2>(info)));
            }
            //std::cout<<uv.row(std::get<0>(L.back()))<<std::endl;
            //std::cout<<uv.row(std::get<1>(L.back()))<<std::endl;
            collapse_edge(e,pt,uv,F,dEF,dEI,EE,allE,{});
            do_collapse = true;
        }
      }
    }
    S.clear();
    test_flip(uv, F, S);
    Eigen::VectorXd Energy;
    compute_dirichlet_energy(avg,uv,F,FL,eps,Energy);
    int valid_face = 0;
    for(int i=0;i<F.rows();i++){
      if(F.row(i).sum() == 0) continue;
      valid_face++;
      if(!std::isfinite(Energy(i)) || Energy(i)>eps){
          if(std::find(S.begin(),S.end(),i) == S.end()){
              S.push_back(i);
          }
      }
    }
    int neighbor = 0;
    if(!do_collapse){
      Eigen::MatrixXi FF,FFI; // adjacency map
      igl::triangle_triangle_adjacency(F,FF,FFI);
      std::cout<<"move to barycenter"<<std::endl;
      int n = 0;
      for(int fid: S){
        std::cout<<n++<<"/"<<S.size()<<std::endl;
        Eigen::MatrixXi local_F;
        int interior_id=expand_to_boundary(V,F,B,FF,fid,local_F);
        //igl::opengl::glfw::Viewer vr;
        Eigen::RowVector2d the_center;
        for(int i=0;i<local_F.rows();i++){
          if(local_F.row(i).sum()!=0){
            Eigen::RowVector2d barycenter = (uv.row(local_F(i,0)) + uv.row(local_F(i,1)) + uv.row(local_F(i,2))) / 3;
            the_center += barycenter;
          }
        }
        the_center /= local_F.rows();
        uv.row(interior_id) << the_center;
        break;
        //plot_mesh(vr,uv,local_F,{});
      }
      neighbor++;
    }else
      neighbor = 0;
    int invalid_size = S.size();
    std::cout<<"invalid size "<<invalid_size<<std::endl;
    std::cout<<"valid face "<<valid_face<<std::endl;
    // also collapse its level-k neighbors
    std::set<int> S_aux;
    Eigen::MatrixXi FF,FFI; // adjacency map
    igl::triangle_triangle_adjacency(F,FF,FFI);
    Eigen::MatrixXi region;
    for(int i=0;i<S.size() && neighbor!=0;i++){
        Eigen::VectorXi nbs;
        get_neighbor_at_level(uv,F,FF,S[i],neighbor,nbs,region);
        //S_aux.insert(S[i]);
        for(int j=0;j<nbs.rows();j++)
            S_aux.insert(nbs(j));   
    }
    // for(int s: S_aux){
    //     int count = 0;
    //     for(int k=0;k<3;k++){
    //         if(B(F(f,k)))
    //             count++;
    //         }
    //     }
    //     if(count == 2)
    //         continue;
    //     if(std::find(S.begin(),S.end(),s) == S.end())
    //         S.push_back(s);
    // }
    // if((valid_face == invalid_size || !do_collapse) && flipped == 0)
    //     invalid_size = 0;
    // }
  }
}
bool insert_vertex_back(
  const std::vector<Action>& L,
  const Eigen::VectorXi& B,
  const Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  Eigen::MatrixXd& uv
){
      // for every action in the list
    igl::Timer timer;
    timer.start();
    bool rollback = false;
    int degree = 2;
    double total_time = 0.0;
    int ii = 0;
    for(;ii<L.size();ii++){
        double time1 = timer.getElapsedTime();
        std::cout<<"rollback "<<ii<<"/"<<L.size()<<std::endl;
        std::cout<<std::get<0>(L[ii])<<" <-> "<<std::get<1>(L[ii])<<std::endl;
        auto a = L[ii];
        Eigen::MatrixXi ring = std::get<2>(a);
        Eigen::VectorXi nbs = std::get<3>(a);
        bool both_on_bd = std::get<4>(a);
        // position candidates
        //#define BRUTAL
        bool found = false;
        Eigen::RowVectorXd pos(2);
        double minarea = std::numeric_limits<double>::min();
        double maxenergy = std::numeric_limits<double>::max();
        Eigen::MatrixXd cd;
        int v0 = std::get<1>(a);
        int v1 = std::get<0>(a);
        auto uv_store = uv;
        uv.row(v1) = uv.row(v0);
        auto F_store = F;
        for(int j=0;j<nbs.rows();j++)
            F.row(nbs(j))<<ring.row(j);
        double angle_sum_of_v0 = calculate_angle_sum(ring,uv,v0,v1);
        if(angle_sum_of_v0 < igl::PI)
            std::swap(v0,v1);
        find_candidate_positions(uv,ring,v1,v0,cd); 
        // update average area
        Eigen::VectorXd area;
        double avg = 0.0;
        igl::doublearea(uv,F,area);
        int ctf=0;
        for(int k=0;k<area.rows();k++){
            if(!isnan(area(k)) && !isinf(area(k))){
                avg += area(k);
                if(F.row(k).sum()!=0)
                    ctf++;
            }
        }
        //std::cout<<avg<<" divide into "<<ctf<<std::endl;
        avg /= ctf;

        std::cout<<"try pos: "<<cd.rows()<<std::endl;
        //std::cout<<"use target area insertion "<<avg<<std::endl;
        for(int j=0;j<cd.rows();j++){
            Eigen::RowVector2d y;
            std::pair<bool,double> z;
            z = flip_avoid_line_search(V,uv,ring,v0,v1,cd.row(j),y,avg);
            //std::cout<<"position "<<j<<": "<<std::get<0>(z)<<","<<std::get<1>(z)<<std::endl;
            int count = 0;
            if(std::get<0>(z) == true && std::get<1>(z)<maxenergy){
                found = true;
                pos << y(0),y(1);
                maxenergy = std::get<1>(z);
            }
        }
        if(found){
            auto F_t = F;
            int d = 0;
            for(int i=0;i<F_t.rows();i++){
                if(F.row(i).sum() != 0)
                    F_t.row(d++)<<F.row(i);
            }
            F_t.conservativeResize(d,3);
            std::cout<<"current face size "<<d<<"/"<<F.rows()<<std::endl;
            uv.row(v1) << pos;
            local_smoothing(V,F_t,B,uv,10,1e10);
        }else{
            F = F_store;
            uv = uv_store;
            // show the ring
            Eigen::MatrixXi F_show(ring.rows(),3);
            for(int i=0;i<ring.rows();i++){
                F_show.row(i)<<ring.row(i);
            }
            // std::cout<<AAk<<std::endl;
            uv.row(v1) = uv.row(v0);
            for(int r=0;r<ring.rows();r++){
                for(int k=0;k<3;k++){
                    std::cout<<ring(r,k)<<" ";
                }
                std::cout<<std::endl;
            }
            std::map<int,int> ring_to_dump;
            int c = 0;
            Eigen::MatrixXi local_F = ring;
            for(int t=0;t<ring.rows();t++){
                for(int k=0;k<3;k++)
                    if(ring_to_dump.find(ring(t,k))==ring_to_dump.end())
                        ring_to_dump[ring(t,k)] = c++;
            }
            
            for(int i=0;i<ring.rows();i++){
                for(int k=0;k<3;k++){
                    local_F(i,k)=ring_to_dump[ring(i,k)];
                    std::cout<<ring_to_dump[ring(i,k)]<<" ";
                }
                std::cout<<std::endl;
            }
            Eigen::MatrixXd local_v(ring_to_dump.size(),2);
            Eigen::MatrixXd local_v3(ring_to_dump.size(),3);
            for(auto p: ring_to_dump){
                local_v3.row(p.second)<<V.row(p.first);
                local_v.row(p.second)<<uv.row(p.first);
            }
            std::vector<int> xx;
            xx.resize(ring_to_dump.size());
            Eigen::VectorXi xxt;
            for(auto t: ring_to_dump)
                xx[t.second] = t.first;
            igl::list_to_matrix(xx,xxt);
            Eigen::MatrixXd CN;
            Eigen::MatrixXi FN;
            igl::writeOBJ("small_region.obj",local_v3,local_F,CN,FN,local_v,local_F);
            //igl::writeOFF("small_region.off",local_v,local_F);
            std::streamsize ss = std::cout.precision();
            std::cout<<std::setprecision(17)<<local_v<<std::endl;
            std::cout << std::setprecision(6);
            // std::cout<<"the collapse before this"<<std::endl;
            Eigen::MatrixXi ring_prev = std::get<2>(L[ii+1]);
            std::cout<<ring_prev<<std::endl;
            uv = uv_store;
            auto F_t = F;
            int d = 0;
            for(int i=0;i<F_t.rows();i++){
                if(F.row(i).sum() != 0)
                    F_t.row(d++)<<F.row(i);
            }
            F_t.conservativeResize(d,3);
            auto uv_o = uv;
            local_smoothing(V,F_t,B,uv,100,1e10);
            ii--;
        }
        double time2 = timer.getElapsedTime();
        total_time += (time2 - time1);
        double expect_total_time = (total_time / (ii+1)) * L.size();
        std::cout<<"expect time: "<<(expect_total_time-total_time)/60.0<<" mins "<<std::endl;
    } 
    return true;
}

bool progressive_embedding(
  const Eigen::VectorXi& bi,
  const Eigen::MatrixXd& b,
  const Eigen::MatrixXd& V,
  double eps,
  Eigen::MatrixXi& F,
  Eigen::MatrixXd& uv
){

  Eigen::VectorXd area;
  igl::doublearea(uv,F,area);
  double avg = 0;
  for(int k=0;k<area.rows();k++){
      if(std::isfinite(area(k))){
          avg += area(k);
      }
  }
  avg /= area.rows();
  // avg = area.array().isFinite().select(0,area)

  // Mark boundary vertices
  Eigen::VectorXi B;
  B.setZero(V.rows());
  for(int i=0;i<bi.rows();i++)
    B(bi(i)) = 1;

  Eigen::VectorXi I;
  collect_invalid_elements(V,F,uv,eps,avg,I);
  
  std::vector<Action> L;
  collapse_invalid_elements(V,F,uv,I,B,eps,avg,L);

  // invert vertex in reverse order of L
  std::reverse(L.begin(),L.end());
  insert_vertex_back(L,B,V,F,uv);
  Eigen::VectorXi T;
  count_flipped_element(uv,F,T);
  std::cout<<"flipped: "<<T.sum()<<std::endl;
  igl::opengl::glfw::Viewer vr;
  vr.data().set_mesh(uv,F);
  vr.launch();
  return true;
}

bool progressive_fix(
    const Eigen::VectorXi& cs,
    Eigen::VectorXi& bi,
    Eigen::MatrixXd& b,
    const Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::MatrixXd& uv
){
    double eps = 1e20;
    Eigen::MatrixXi F_test=F;
    Eigen::MatrixXd uv_test=uv;
    progressive_embedding(bi,b,V,eps,F_test,uv_test);
    
    igl::Timer timer;
    timer.start();
    double t0 = timer.getElapsedTime();

    // [calculate reference triangle area for local_smoothing]
    Eigen::VectorXd area;
    igl::doublearea(uv,F,area);
    double avg = 0;
    for(int k=0;k<area.rows();k++){
        if(std::isfinite(area(k))){
            avg += area(k);
        }
    }
    avg /= area.rows();
    std::cout<<"averge area is "<<avg<<std::endl;
    Eigen::VectorXi B;
    B.setZero(V.rows());
    for(int i=0;i<bi.rows();i++)
      B(bi(i)) = 1;
    // [init Dirichlet energy]
    std::cout<<"init dirichlet energy ... "<<std::endl;
    Eigen::VectorXd Energy;
    std::vector<int> FL(F.rows());
    std::iota(FL.begin(),FL.end(),0);
    compute_dirichlet_energy(avg,uv,F,FL,eps,Energy);
    std::cout<<"done"<<std::endl;

    // [collect invalid elements]
    std::cout<<"collect invalid elements"<<std::endl;
    Eigen::VectorXi I;
    I.setZero(F.rows());
    std::vector<int> S;
    test_flip(uv, F, S);
    count_flipped_element(uv,F,I);
    std::cout<<"flipped: "<<I.sum()<<std::endl;
    for(int i=0;i<F.rows();i++){
      if(I(i)) continue;
      if(!std::isfinite(Energy(i))||Energy(i)>eps){
        I(i) = 1;
        std::cout<<area(i)<<std::endl;
      }
    }
    std::cout<<"invalid: "<<I.sum()<<std::endl;
    std::cout<<"done"<<std::endl;

    //std::vector<int> S;
    S.clear();
    for(int i=0;i<I.rows();i++){
      if(I(i)) S.push_back(i);
    }

    int flipped = I.sum();
    Eigen::MatrixXi F_original = F;
    Eigen::MatrixXd uv_original = uv;
    Eigen::VectorXi bdr;
    igl::boundary_loop(F,bdr);
    std::vector<int> bdr_v;
    igl::matrix_to_list(bi,bdr_v);
    std::set<int> bdr_s(bdr_v.begin(),bdr_v.end());

    Eigen::VectorXi M; // should not be moved by local smoothing
    M.setZero(V.rows());
    for(int i=0;i<bdr.rows();i++)
        M(bdr(i)) = 1;
    for(int i=0;i<cs.rows();i++)
        M(cs(i)) = 1;
    
    std::vector<Action> L; //action_list
    Eigen::VectorXi EMAP,EE;
    Eigen::MatrixXi E,EF,EI;
    Eigen::MatrixXi dEF,dEI,allE;
    igl::edge_flaps(F,E,allE,EMAP,EF,EI,dEF,dEI,EE);
    std::set<int> candidates(S.begin(),S.end());
    int invalid_size = S.size();
    int neighbor = 0;
    int valid_faces = 0;
    while(invalid_size!=0){
        bool do_collapse = false;
        for(int f: S){
            typedef std::tuple<double,int,bool> edge_info;
            std::vector<edge_info> f_info;
            for(int k=0;k<3;k++){
                double l = (uv.row(F(f,k))-uv.row(F(f,(k+1)%3))).norm();
                bool on_bd = (bdr_s.find(F(f,k))!=bdr_s.end()) || (bdr_s.find(F(f,(k+1)%3))!=bdr_s.end());
                f_info.push_back(edge_info(l,k,on_bd));
            }
            std::sort(f_info.begin(),f_info.end());
            // collapse it if valid
            for(auto info: f_info){
                int k = std::get<1>(info);
                int e = F.rows()*((k+2)%3)+f;
                int a = F(f,k), b = F(f,(k+1)%3);
                std::vector<int> N;
                //std::cout<<"try collpasing "<<a<<" and "<<b<<std::endl;
                if(std::get<2>(info)) {
                    //std::cout<<"its on bd\n";
                    continue;
                }
                if(edge_collapse_is_valid(e,uv.row(b),uv,F,dEF,dEI,EE,allE,N)){
                    std::sort(N.begin(), N.end());
                    N.erase( std::unique( N.begin(), N.end() ), N.end() );
                    Eigen::VectorXi nbi;
                    igl::list_to_matrix(N,nbi);
                    Eigen::MatrixXd nb(N.size()*3,2);
                    Eigen::MatrixXi F_store(nbi.rows(),3);
                    // collect all the positions beside a and b
                    int fn = 0;
                    for(int fi: N){
                        if(F.row(fi).sum()==0) continue;
                        F_store.row(fn++)<<F.row(fi);
                    }
                    F_store.conservativeResize(fn,3);
                    Eigen::RowVector2d pt;
                    if(bdr_s.find(a) != bdr_s.end()){
                        // if a is boundary collapse to a
                        pt = uv.row(a);
                        L.push_back(Action(a,b,F_store,nbi,std::get<2>(info)));
                    }else if(bdr_s.find(b) != bdr_s.end()){
                        // if b is boundary collapse to b
                        pt = uv.row(b);
                        L.push_back(Action(b,a,F_store,nbi,std::get<2>(info)));
                    }else{
                        // otherwise collapse to a small index
                        pt = uv.row(std::min(a,b));
                        L.push_back(Action(std::max(a,b),std::min(a,b),F_store,nbi,std::get<2>(info)));
                    }
                    bool is_bd1 = (bdr_s.find(std::get<0>(L.back()))!=bdr_s.end());
                    bool is_bd2 = (bdr_s.find(std::get<1>(L.back()))!=bdr_s.end());
                    //std::cout<<uv.row(std::get<0>(L.back()))<<std::endl;
                    //std::cout<<uv.row(std::get<1>(L.back()))<<std::endl;
                    collapse_edge(e,pt,uv,F,dEF,dEI,EE,allE,bdr_s);
                    do_collapse = true;
                    //std::cout<<"collapse "<<std::get<0>(L.back())<<"("<<is_bd1<<") to "<<std::get<1>(L.back())<<"("<<is_bd2<<")"<<std::endl;
                    for(int i=0;i<nbi.rows();i++)
                        candidates.insert(nbi(i));
                }else{
                    std::sort(N.begin(), N.end());
                    N.erase( std::unique( N.begin(), N.end() ), N.end() );
                    Eigen::VectorXi nbi;
                    igl::list_to_matrix(N,nbi);
                    for(int i=0;i<nbi.rows();i++)
                        candidates.insert(nbi(i));
                }
            }
        }
        S.clear();
        // Eigen::MatrixXi cF(candidates.size(),3);
        // int i=0;
        // for(auto s: candidates)
        //     cF.row(i++)<<F.row(s);
        test_flip(uv, F, S);
        flipped = S.size();
        compute_dirichlet_energy(avg,uv,F,FL,eps,Energy);
        int valid_face = 0;
        for(int i=0;i<F.rows();i++){
            if(F.row(i).sum() == 0) continue;
            valid_face++;
            if(!std::isfinite(Energy(i)) || Energy(i)>eps){
                if(std::find(S.begin(),S.end(),i) == S.end()){
                    // std::cout<<Energy(i)<<std::endl;
                    S.push_back(i);
                }
            }
        }
        if(!do_collapse){
          Eigen::MatrixXi FF,FFI; // adjacency map
          igl::triangle_triangle_adjacency(F,FF,FFI);
          std::cout<<"move to barycenter"<<std::endl;
          int n = 0;
          for(int fid: S){
            std::cout<<n++<<"/"<<S.size()<<std::endl;
            Eigen::MatrixXi local_F;
            int interior_id=expand_to_boundary(V,F,B,FF,fid,local_F);
            //igl::opengl::glfw::Viewer vr;
            Eigen::RowVector2d the_center;
            for(int i=0;i<local_F.rows();i++){
              if(local_F.row(i).sum()!=0){
                Eigen::RowVector2d barycenter = (uv.row(local_F(i,0)) + uv.row(local_F(i,1)) + uv.row(local_F(i,2))) / 3;
                the_center += barycenter;
              }
            }
            the_center /= local_F.rows();
            uv.row(interior_id) << the_center;
            break;
            //plot_mesh(vr,uv,local_F,{});
          }
          neighbor++;
        }
        else
            neighbor = 0;
        invalid_size = S.size();
        std::cout<<"invalid size "<<invalid_size<<std::endl;
        std::cout<<"valid face "<<valid_face<<std::endl;
        // also collapse its level-k neighbors
        std::set<int> S_aux;
        Eigen::MatrixXi FF,FFI; // adjacency map
        igl::triangle_triangle_adjacency(F,FF,FFI);
        Eigen::MatrixXi region;
        for(int i=0;i<S.size() && neighbor!=0;i++){
            Eigen::VectorXi nbs;
            get_neighbor_at_level(uv,F,FF,S[i],neighbor,nbs,region);
            //S_aux.insert(S[i]);
            for(int j=0;j<nbs.rows();j++)
                S_aux.insert(nbs(j));   
        }
        for(int s: S_aux){
            int count = 0;
            for(int k=0;k<3;k++){
                if(bdr_s.find(F(s,k)) != bdr_s.end()){
                    count++;
                }
            }
            if(count == 2)
                continue;
            if(std::find(S.begin(),S.end(),s) == S.end())
                S.push_back(s);
        }
        if((valid_face == invalid_size || !do_collapse) && flipped == 0)
            invalid_size = 0;
    }
    std::cout<<"L.size() is "<<L.size()<<std::endl;
    
    std::reverse(L.begin(),L.end());
    // for every action in the list
    bool rollback = false;
    int degree = 2;
    double total_time = 0.0;
    int ii = 0;
    for(;ii<L.size();ii++){
        double time1 = timer.getElapsedTime();
        std::cout<<"rollback "<<ii<<"/"<<L.size()<<std::endl;
        bool is_bd1 = (bdr_s.find(std::get<0>(L[ii]))!=bdr_s.end());
        bool is_bd2 = (bdr_s.find(std::get<1>(L[ii]))!=bdr_s.end());
        std::cout<<std::get<0>(L[ii])<<"("<<is_bd1<<") <-> "<<std::get<1>(L[ii])<<"("<<is_bd2<<")"<<std::endl;
        auto a = L[ii];
        Eigen::MatrixXi ring = std::get<2>(a);
        Eigen::VectorXi nbs = std::get<3>(a);
        bool both_on_bd = std::get<4>(a);
        // position candidates
        //#define BRUTAL
        bool found = false;
        Eigen::RowVectorXd pos(2);
        double minarea = std::numeric_limits<double>::min();
        double maxenergy = std::numeric_limits<double>::max();
        Eigen::MatrixXd cd;
        int v0 = std::get<1>(a);
        int v1 = std::get<0>(a);
        auto uv_store = uv;
        uv.row(v1) = uv.row(v0);
        auto F_store = F;
        for(int j=0;j<nbs.rows();j++)
            F.row(nbs(j))<<ring.row(j);
        double angle_sum_of_v0 = calculate_angle_sum(ring,uv,v0,v1);
        if(angle_sum_of_v0 < igl::PI)
            std::swap(v0,v1);
        find_candidate_positions(uv,ring,v1,v0,cd); 
        // calc average area
        Eigen::VectorXd area;
        double avg = 0.0;
        igl::doublearea(uv,F,area);
        int ctf=0;
        for(int k=0;k<area.rows();k++){
            if(!isnan(area(k)) && !isinf(area(k))){
                avg += area(k);
                if(F.row(k).sum()!=0)
                    ctf++;
            }
        }
        //std::cout<<avg<<" divide into "<<ctf<<std::endl;
        avg /= ctf;

        std::cout<<"try pos: "<<cd.rows()<<std::endl;
        //std::cout<<"use target area insertion "<<avg<<std::endl;
        for(int j=0;j<cd.rows();j++){
            Eigen::RowVector2d y;
            std::pair<bool,double> z;
            z = flip_avoid_line_search(V,uv,ring,v0,v1,cd.row(j),y,avg);
            //std::cout<<"position "<<j<<": "<<std::get<0>(z)<<","<<std::get<1>(z)<<std::endl;
            int count = 0;
            if(std::get<0>(z) == true && std::get<1>(z)<maxenergy){
                found = true;
                pos << y(0),y(1);
                maxenergy = std::get<1>(z);
                //std::cout<<"getting "<<std::get<1>(z)<<std::endl;
                //std::cout<<"position is "<<pos<<std::endl;
            }
        }
        // if(found){
        //     min_area = std::min(min_area,minarea);
        // }
        // std::cout<<"min area is now: "<<min_area<<std::endl;
        #ifdef SERIAL
        if(ii%4 == 0){
            // serialization
            // - F
            // - uv
            // - ii
            // model_name + eps + progressive_bin
            std::string serial_name = model_name + "_eps" + std::to_string(threshold) + "_pb";
            igl::serialize(ii,"ii",serial_name,true);
            igl::serialize(uv,"uv",serial_name);
            igl::serialize(F,"F",serial_name);
            std::cout<<"searialize "<<ii<<std::endl;
        }
        #endif
        if(found){
            auto F_t = F;
            int d = 0;
            for(int i=0;i<F_t.rows();i++){
                if(F.row(i).sum() != 0)
                    F_t.row(d++)<<F.row(i);
            }
            F_t.conservativeResize(d,3);
            std::cout<<"current face size "<<d<<"/"<<F.rows()<<std::endl;
            uv.row(v1) << pos;
            local_smoothing(V,F_t,M,uv,10,1e10);
        }else{
            F = F_store;
            uv = uv_store;
            // show the ring
            Eigen::MatrixXi F_show(ring.rows(),3);
            for(int i=0;i<ring.rows();i++){
                F_show.row(i)<<ring.row(i);
            }
            // std::cout<<AAk<<std::endl;
            uv.row(v1) = uv.row(v0);
            for(int r=0;r<ring.rows();r++){
                for(int k=0;k<3;k++){
                    bool is_bdi = (bdr_s.find(ring(r,k)) != bdr_s.end());
                    std::cout<<ring(r,k)<<"("<<is_bdi<<")"<<" ";
                }
                std::cout<<std::endl;
            }
            std::map<int,int> ring_to_dump;
            int c = 0;
            Eigen::MatrixXi local_F = ring;
            for(int t=0;t<ring.rows();t++){
                for(int k=0;k<3;k++)
                    if(ring_to_dump.find(ring(t,k))==ring_to_dump.end())
                        ring_to_dump[ring(t,k)] = c++;
            }
            
            for(int i=0;i<ring.rows();i++){
                for(int k=0;k<3;k++){
                    local_F(i,k)=ring_to_dump[ring(i,k)];
                    std::cout<<ring_to_dump[ring(i,k)]<<" ";
                }
                std::cout<<std::endl;
            }
            Eigen::MatrixXd local_v(ring_to_dump.size(),2);
            Eigen::MatrixXd local_v3(ring_to_dump.size(),3);
            for(auto p: ring_to_dump){
                local_v3.row(p.second)<<V.row(p.first);
                local_v.row(p.second)<<uv.row(p.first);
            }
            std::vector<int> xx;
            xx.resize(ring_to_dump.size());
            Eigen::VectorXi xxt;
            for(auto t: ring_to_dump)
                xx[t.second] = t.first;
            igl::list_to_matrix(xx,xxt);
            Eigen::MatrixXd CN;
            Eigen::MatrixXi FN;
            igl::writeOBJ("small_region.obj",local_v3,local_F,CN,FN,local_v,local_F);
            //igl::writeOFF("small_region.off",local_v,local_F);
            std::streamsize ss = std::cout.precision();
            std::cout<<std::setprecision(17)<<local_v<<std::endl;
            std::cout << std::setprecision(6);
            // std::cout<<"the collapse before this"<<std::endl;
            Eigen::MatrixXi ring_prev = std::get<2>(L[ii+1]);
            std::cout<<ring_prev<<std::endl;
            is_bd1 = (bdr_s.find(std::get<0>(L[ii+1]))!=bdr_s.end());
            is_bd2 = (bdr_s.find(std::get<1>(L[ii+1]))!=bdr_s.end());
            uv = uv_store;
            auto F_t = F;
            int d = 0;
            for(int i=0;i<F_t.rows();i++){
                if(F.row(i).sum() != 0)
                    F_t.row(d++)<<F.row(i);
            }
            F_t.conservativeResize(d,3);   
            // igl::opengl::glfw::Viewer vv;
            // vv.data().set_mesh(uv,F);
            // vv.launch();
            // optimze and try again
            auto uv_o = uv;
            local_smoothing(V,F_t,M,uv,100,1e10);
            // if((uv_o-uv).norm()<1e-20) // negligible change by local smooth
            //     break;
            ii--;
        }
        double time2 = timer.getElapsedTime();
        total_time += (time2 - time1);
        double expect_total_time = (total_time / (ii+1)) * L.size();
        std::cout<<"expect time: "<<(expect_total_time-total_time)/60.0<<" mins "<<std::endl;
    }  
    // Eigen::MatrixXd bdr_pos1;
    // igl::slice(uv_original,bdr,1,bdr_pos1);
    // igl::slice(uv,bdr,1,bdr_pos2);
    // auto dd = bdr_pos2-bdr_pos1;
    // std::cout<<"check constraints "<<std::endl;
    // for(int i=0;i<dd.rows();i++){
    //     if(dd(i)!=0){
    //         std::cout<<">>>>"<<std::setprecision(17)<<bdr_pos1.row(i)<<std::endl;
    //         std::cout<<"<<<<"<<std::setprecision(17)<<bdr_pos2.row(i)<<std::endl;
    //     }
    // }
    // std::cout<<"check constraints done"<<std::endl;
    // double tn = timer.getElapsedTime();
    // std::ofstream log;
    // log.open("progressive_embedding_log.txt",std::ios_base::app);
    // log << "modelname, #V, #F, #invalid, #flip, progressive embedding time\n";
    // log << model_name <<","<< V.rows() << ","
    //        << F.rows() << "," 
    //        << invalid_all << ","
    //        << flip_all << ","
    //        << stat_progressive_time <<"\n";
    return true;
}