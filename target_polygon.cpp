#include "target_polygon.h"
#include "path_tracing.h"
#include "cut_mesh/cut_mesh_simple.h"
#include <igl/setdiff.h>
#include <igl/list_to_matrix.h>
#include <igl/triangle/triangulate.h>
#include <igl/sortrows.h>
#include <iostream>
#include <deque>
#include <igl/opengl/glfw/Viewer.h>

using namespace std;

// generate background mesh, with all boudnary constraints embeded
// p: boundary constraint points coords
void background_mesh(
    int N,                     
    const Eigen::MatrixXd& c,
    const Eigen::VectorXi& ci,
    map<int,int>& m,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F
){

    // calculate the barycenter b of boudnary constraints
    // calculate range_x
    // calculate range_y
    // bounding box centered b, with edge width max(range_x,range_y)
    
    double range_x,range_y;
    Eigen::RowVector2d barycenter;
    if(c.rows()==0){
        range_x = 10.0f;
        range_y = 10.0f;
        barycenter << 0,0;
    }else{
        range_x = c.col(0).maxCoeff() - c.col(0).minCoeff();
        range_y = c.col(1).maxCoeff() - c.col(1).minCoeff();
        barycenter = c.colwise().sum() / c.rows();
    }
    double width = std::max(range_x, range_y)*1.5;
    double e = width / N;
    Eigen::MatrixXd P(4*N,2);
    Eigen::MatrixXi E(4*N,2);
    Eigen::RowVector2d Z;
    Z<<barycenter(0)-width/2,barycenter(1)+width/2;
    //#define counterclock
    #ifdef counterclock
    P.block(0  ,0,N,1) = Eigen::VectorXd::LinSpaced(N,Z(0),      Z(0));
    P.block(0  ,1,N,1) = Eigen::VectorXd::LinSpaced(N,Z(1),      Z(1)-width+e);
    P.block(N  ,0,N,1) = Eigen::VectorXd::LinSpaced(N,Z(0),      Z(0)+width-e);
    P.block(N  ,1,N,1) = Eigen::VectorXd::LinSpaced(N,Z(1)-width,Z(1)-width);
    P.block(2*N,0,N,1) = Eigen::VectorXd::LinSpaced(N,Z(0)+width,Z(0)+width);
    P.block(2*N,1,N,1) = Eigen::VectorXd::LinSpaced(N,Z(1)-width,Z(1)-e);
    P.block(3*N,0,N,1) = Eigen::VectorXd::LinSpaced(N,Z(0)+width,Z(0)+e);
    P.block(3*N,1,N,1) = Eigen::VectorXd::LinSpaced(N,Z(1),      Z(1));
    #else
    Z << barycenter(0)-width/2,barycenter(1)-width/2;
    P.block(0  ,0,N,1) = Eigen::VectorXd::LinSpaced(N,Z(0),Z(0)+width-e);
    P.block(0  ,1,N,1) = Eigen::VectorXd::LinSpaced(N,Z(1),      Z(1));
    P.block(N  ,0,N,1) = Eigen::VectorXd::LinSpaced(N,Z(0)+width,Z(0)+width);
    P.block(N  ,1,N,1) = Eigen::VectorXd::LinSpaced(N,Z(1), Z(1)+width-e);
    P.block(2*N,0,N,1) = Eigen::VectorXd::LinSpaced(N,Z(0)+width,Z(0)+e);
    P.block(2*N,1,N,1) = Eigen::VectorXd::LinSpaced(N,Z(1)+width,Z(1)+width);
    P.block(3*N,0,N,1) = Eigen::VectorXd::LinSpaced(N,Z(0),Z(0));
    P.block(3*N,1,N,1) = Eigen::VectorXd::LinSpaced(N,Z(1)+width,Z(1)+e);
    #endif
    E<<Eigen::VectorXi::LinSpaced(4*N,0,4*N),
        Eigen::VectorXi::LinSpaced(4*N,1,4*N+1);
    E(4*N-1,1)=0;
    
    // embed hard constraints
    P.conservativeResize(P.rows()+c.rows(),2);
    if(c.rows()>0)
        P.bottomRows(c.rows())<<c;
    for(int i=0;i<ci.rows();i++){
        m[ci(i)] = 4*N+i;
    }
    // std::cout<<P.rows()<<std::endl;
    // igl::opengl::glfw::Viewer vr;
    // igl::opengl::glfw::imgui::ImGuiMenu menu;
    // vr.plugins.push_back(&menu);
    // //for(int i=0;i<P.rows();i++){
    // vr.data().add_points(P.row(0),Eigen::RowVector3d(0,0,0));
    // vr.data().add_points(P.row(N),Eigen::RowVector3d(0,0,0));
    // vr.data().add_points(P.row(N*2),Eigen::RowVector3d(0,0,0));
    // vr.data().add_points(P.row(N*3),Eigen::RowVector3d(0,0,0));
    // for(int i=0;i<P.rows();i++){
    //     //vr.data().add_points(P.row(i),Eigen::RowVector3d(1,0,0));
    //     //vr.data().add_label(P.row(i),std::to_string(i));
    // }
    igl::opengl::glfw::Viewer viewer;
    // for(int i=0;i<P.rows();i++){
    //     viewer.data().add_points(P.row(i),Eigen::RowVector3d(1,0,0));
    // }
    // viewer.core().align_camera_center(P);
    // viewer.launch();
    std::vector<std::pair<int,int>> E_pair;
    for(int i=0;i<E.rows();i++){
        E_pair.push_back(std::make_pair(E(i,0),E(i,1)));
    }
    //display(P,F,E_pair,{},{},{});
    igl::triangle::triangulate(P,E,Eigen::MatrixXd(),"YQq33",V,F);
    //display(V,F,{},{},{},{});
    // viewer.data().set_mesh(V,F);
    // viewer.launch();
}

void generate_hilbert(
    int order,
    const Eigen::MatrixXd& c,
    const Eigen::VectorXi& ci,
    map<int,int>& m,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F
){
    using namespace std; 
    auto cum_sum = [&](Eigen::MatrixXd& X)->Eigen::MatrixXd{
        Eigen::MatrixXd Y(X.rows(),X.cols());
        for(int i=0;i<X.rows();i++){
            Eigen::RowVectorXd sum_i(X.cols());
            sum_i.setZero();
            for(int j=0;j<=i;j++){
                sum_i += X.row(j);
            }
            Y.row(i) << sum_i;
        }
        return Y;
    };
    
    auto Hilbert_curve = [&]()->Eigen::MatrixXd{
        Eigen::MatrixXd A(0,2);
        Eigen::MatrixXd B(0,2);
        Eigen::MatrixXd C(0,2);
        Eigen::MatrixXd D(0,2);
        Eigen::MatrixXd north(1,2);
        Eigen::MatrixXd east(1,2);
        Eigen::MatrixXd south(1,2);
        Eigen::MatrixXd west(1,2);
        north<<0,1;
        east<<1,0;
        south<<0,-1;
        west<<-1,0;
        A.setZero();
        B.setZero();
        C.setZero();
        D.setZero();
        int order = 3;
        for(int i=0;i<order;i++){
            int row_AA = static_cast<int>(B.rows() + A.rows()*2 + C.rows() + 3);
            int row_BB = static_cast<int>(A.rows() + B.rows()*2 + D.rows() + 3);
            int row_CC = static_cast<int>(D.rows() + C.rows()*2 + A.rows() + 3);
            int row_DD = static_cast<int>(C.rows() + D.rows()*2 + B.rows() + 3);
            Eigen::MatrixXd AA(row_AA,2);
            Eigen::MatrixXd BB(row_BB,2);
            Eigen::MatrixXd CC(row_CC,2);
            Eigen::MatrixXd DD(row_DD,2);
            AA.block(0,0,B.rows(),2) = B;
            AA.block(B.rows(),0,1,2) = north;
            AA.block(B.rows()+1,0,A.rows(),2) = A;
            AA.block(B.rows()+1+A.rows(),0,1,2) = east;
            AA.block(B.rows()+1+A.rows()+1,0,A.rows(),2) = A;
            AA.block(B.rows()+1+A.rows()+1+A.rows(),0,1,2) = south;
            AA.block(B.rows()+1+A.rows()+1+A.rows()+1,0,C.rows(),2) = C;
            // BB
            BB.block(0,0,A.rows(),2) = A;
            BB.block(A.rows(),0,1,2) = east;
            BB.block(A.rows()+1,0,B.rows(),2) = B;
            BB.block(A.rows()+1+B.rows(),0,1,2) = north;
            BB.block(A.rows()+1+B.rows()+1,0,B.rows(),2) = B;
            BB.block(A.rows()+1+B.rows()+1+B.rows(),0,1,2) = west;
            BB.block(A.rows()+1+B.rows()+1+B.rows()+1,0,D.rows(),2) = D;
            // CC
            CC.block(0,0,B.rows(),2) = D;
            CC.block(D.rows(),0,1,2) = west;
            CC.block(D.rows()+1,0,C.rows(),2) = C;
            CC.block(D.rows()+1+C.rows(),0,1,2) = south;
            CC.block(D.rows()+1+C.rows()+1,0,C.rows(),2) = C;
            CC.block(D.rows()+1+C.rows()+1+C.rows(),0,1,2) = east;
            CC.block(D.rows()+1+C.rows()+1+C.rows()+1,0,C.rows(),2) = A;
            // DD
            DD.block(0,0,C.rows(),2) = C;
            DD.block(C.rows(),0,1,2) = south;
            DD.block(C.rows()+1,0,D.rows(),2) = D;
            DD.block(C.rows()+1+D.rows(),0,1,2) = west;
            DD.block(C.rows()+1+D.rows()+1,0,D.rows(),2) = D;
            DD.block(C.rows()+1+D.rows()+1+D.rows(),0,1,2) = north;
            DD.block(C.rows()+1+D.rows()+1+D.rows()+1,0,C.rows(),2) = B;
            A = AA;
            B = BB;
            C = CC;
            D = DD;
        }
        Eigen::MatrixXd Y(A.rows()+1,2);
        Y.setZero();
        Y<<0,0,cum_sum(A);
        return Y;
    };
    
    auto extend_hilbert_to_poly = [&](Eigen::MatrixXd& Y)->deque<complex<double>>{
        deque<complex<double>> poly;
        vector<int> angles;
        angles.resize(Y.rows());
        // extend the curve to a polygon
        for(int i=1;i<Y.rows()-1;i++){
            // do not care about the first and last one (for now)
            Eigen::Vector3d a;
            a<<(Y.row(i) - Y.row(i-1)).transpose(),0;
            Eigen::Vector3d b;
            b<<(Y.row(i) - Y.row(i+1)).transpose(),0;
            if((a.cross(b))(2) > 0){
                angles[i] = 1;
            }else if((a.cross(b))(2) < 0){
                angles[i] = -1;
            }else{ // == 0
                angles[i] = 0;
            }
        }
        if(Y(0,1) == Y(1,1)){
            double l = (Y(0,0)-Y(1,0))/4;
            poly.push_back(complex<double>(Y(0,0)+l,Y(0,1)+l));
            poly.push_front(complex<double>(Y(0,0)+l,Y(0,1)-l));
        }else if(Y(0,0) == Y(1,0)){
            double l = (Y(1,1)-Y(0,1))/4;
            poly.push_back(complex<double>(Y(0,0)+l,Y(0,1)-l));
            poly.push_front(complex<double>(Y(0,0)-l,Y(0,1)-l));
        }
        for(int i=1;i<angles.size()-1;i++){
            Eigen::Vector2d p_a;
            p_a << (((Y.row(i-1)+Y.row(i+1))/2+Y.row(i))/2).transpose();
            Eigen::Vector2d p_b;
            Eigen::Vector2d t;
            t<<(2*Y.row(i)).transpose();
            p_b << t-p_a;
            if(angles[i] == 1){
                poly.push_back(complex<double>(p_a(0),p_a(1)));
                poly.push_front(complex<double>(p_b(0),p_b(1)));
            }else if(angles[i] == -1){
                poly.push_back(complex<double>(p_b(0),p_b(1)));
                poly.push_front(complex<double>(p_a(0),p_a(1)));
            }else{
                if(poly[poly.size()-1].imag() + poly[0].imag() == Y(i,1)*2){
                    poly.push_back(complex<double>(Y(i,0),poly[poly.size()-1].imag()));
                    poly.push_front(complex<double>(Y(i,0),poly[0].imag()));
                }else{
                    poly.push_back(complex<double>(poly[poly.size()-1].real(),Y(i,1)));
                    poly.push_front(complex<double>(poly[0].real(),Y(i,1)));
                }
            }
        }
        int N = Y.rows();
        if(Y(N-1,1) == Y(N-2,1)){
            double l = (Y(N-1,0)-Y(N-2,0))/4;
            poly.push_front(complex<double>(Y(N-1,0)+l,Y(N-1,1)+l));
            poly.push_back(complex<double>(Y(N-1,0)+l,Y(N-1,1)-l));
        }else if(Y(N-1,0) == Y(N-2,0)){
            double l = (Y(N-2,1)-Y(N-1,1))/4;
            poly.push_back(complex<double>(Y(N-1,0)-l,Y(N-1,1)-l));
            poly.push_front(complex<double>(Y(N-1,0)+l,Y(N-1,1)-l));
        }
        return poly;
    };
    auto Y = Hilbert_curve();
    auto poly = extend_hilbert_to_poly(Y);

    Eigen::MatrixXd p(poly.size(),2);
    Eigen::MatrixXi e(p.rows(),2);
    for(int i=0;i<poly.size();i++){
        p.row(i) << poly[i].real(),poly[i].imag();
    }
    e.resize(p.rows(),2);
    e << Eigen::VectorXi::LinSpaced(p.rows(),0,p.rows()-1),
         Eigen::VectorXi::LinSpaced(p.rows(),1,p.rows());
    e(p.rows()-1,1) = 0;
    // embed hard constraints
    int init=p.rows();
    p.conservativeResize(p.rows()+c.rows(),2);
    if(c.rows()>0)
        p.bottomRows(c.rows())<<c;
    for(int i=0;i<ci.rows();i++){
        m[ci(i)] = init+i;
        std::cout<<ci(i)<<" <=== "<<init+i<<std::endl;
    }
    igl::triangle::triangulate(p,e,Eigen::MatrixXd(),"YQq33a0.1",V,F);
    

}

void line_ratio(
    const Eigen::MatrixXd& py_line,
    Eigen::VectorXd& rt
){
    // get total length of polyline
    double l = 0.0f;
    for(int i=0;i<py_line.rows()-1;i++){
        l += (py_line.row(i) - py_line.row(i+1)).norm();
    }

    // calculate ratio
    rt.setZero(py_line.rows());
    for(int i=1;i<py_line.rows()-1;i++){
        rt[i] = rt[i-1] + (py_line.row(i)-py_line.row(i+1)).norm() / l;
    }
    rt(rt.rows()-1) = 1.0f;
}

void sample_polyline(
    const Eigen::MatrixXd& L,
    const Eigen::VectorXd& rt,
    Eigen::MatrixXd& SL
){
    Eigen::VectorXd rt0;
    line_ratio(L,rt0);
    int n = L.rows();
    SL.resize(rt.rows(),L.cols());
    SL.setZero();
    SL.row(0) << L.row(0);
    SL.bottomRows(1) << L.bottomRows(1);
    for(int k=1;k<rt.size()-1;k++){
        for(int l=0;l<rt0.size()-1;l++){
            if(rt[k]>=rt0[l] && rt[k]<=rt0[l+1]){
                double r = (rt0[l+1]-rt[k])/(rt0[l+1]-rt0[l]);
                double a = L(l,0)*r+L((l+1)%n,0)*(1-r);
                if(L(l,0) == L((l+1)%n,0))
                    a = L(l,0);
                double b = L(l,1)*r+L((l+1)%n,1)*(1-r);
                if(L(l,1) == L((l+1)%n,1))
                    b = L(l,1);
                SL(k,0) = a;
                SL(k,1) = b;
                break;
            }
        }
    }
}


// connect all outer-bd constraints to 
// bd of boundary and extract polygon
// input:
//       c, index of constraint points on 3d mesh
//       m, [index 3d -> index 2d]
//       bnd, boundary vertex positions of 3d mesh
//       bndi, boundary vertices index
//       (V2, F2), 2d background mesh
// ouput:
//       P0, extracted polygon
//       T,  correspondence of 3d boundary and extracted polygon
void outer_boundary(
    const Eigen::VectorXi& c,
    std::map<int,int>& m,
    const Eigen::MatrixXd& bnd,
    const Eigen::VectorXi& bndi,
    Eigen::MatrixXd& V2,
    Eigen::MatrixXi& F2,
    Eigen::MatrixXd& P0,
    Eigen::VectorXi& T
){
    P0 = Eigen::MatrixXd();
    T = Eigen::VectorXi();
    // [pick and sort ci that's in bndi]
    Eigen::VectorXi _,cdf,ci;
    igl::setdiff(c,bndi,cdf,_);
    if(cdf.rows()!=0){
        igl::setdiff(c,cdf,ci,_);
    }else
        ci = c;

    auto locate_closest = [](const Eigen::RowVectorXd& p, const Eigen::MatrixXd& rg){
        double dt = std::numeric_limits<double>::max();
        int x = 0;
        for(int i=0;i<rg.rows();i++){
            double d = (p-rg.row(i)).norm();
            if(d < dt){
                dt = d;
                x = i;
            }
        }
        return x;
    };
    Eigen::MatrixXd bnd2;
    Eigen::VectorXi bnd2i;
    igl::boundary_loop(F2,bnd2i);
    igl::slice(V2,bnd2i,1,bnd2);
    int pos = 0;
    Eigen::VectorXi rgi = bnd2i;
    Eigen::MatrixXd rg = bnd2;
    std::vector<int> K;
    std::map<int,int> D; //duplicated vertices
    for(int i=0;i<ci.rows();i++){
        // [connect m[ci(i)] to outer boundary -> path_i]
        pos = locate_closest(V2.row(m[ci(i)]),rg);
        K.push_back(rgi(pos));
        Eigen::VectorXi ti;
        Eigen::MatrixXd t;
        int N = rg.rows();
        if(i == 0){
            t.resize(N-1,2);
            ti.resize(N-1,1);
            t<<rg.bottomRows(N-pos-1),rg.topRows(pos);
            ti<<rgi.bottomRows(N-pos-1),rgi.topRows(pos);
        }else{
            t.resize(N-pos-1,2);
            ti.resize(N-pos-1,1);
            t=rg.bottomRows(N-pos-1);
            ti=rgi.bottomRows(N-pos-1);
        }
        rg = t;
        rgi = ti;
        std::vector<int> E;
        std::set<int> A;
        Eigen::MatrixXd V2n;
        Eigen::MatrixXi F2n;
        std::vector<std::pair<int,int>> Mask;
        std::set<int> no_enter_f;
        Eigen::MatrixXi TT;
        path_tracing(V2,F2,std::make_pair(m[ci(i)],K.back()),A,Mask,no_enter_f,TT,V2n,F2n,E);
        std::vector<std::vector<int>> cv;
        std::vector<std::vector<int>> ck;   
        igl::cut_mesh(V2n,F2n,{E},cv,ck);
        for(int g = 0;g<cv.size();g++){
            if(cv[g].size()==1){
                D.insert(std::make_pair(cv[g][0], cv[g][0]));
            }else {
                if(D.find(cv[g][0])!=D.end() || D.find(cv[g][1])!=D.end()){
                    D[cv[g][0]]=cv[g][1];
                    D[cv[g][1]]=cv[g][0];
                }else{
                    D.insert(std::make_pair(cv[g][0], cv[g][1]));
                    D.insert(std::make_pair(cv[g][1], cv[g][0]));
                }
            }
        }
        V2 = V2n;
        F2 = F2n;
    }
    igl::boundary_loop(F2,bnd2i);
    igl::slice(V2,bnd2i,1,bnd2);
    // sample polygon
    Eigen::VectorXi bi(ci.rows()),b2i(ci.rows()),ki(ci.rows());
    for(int i=0;i<ci.rows();i++){
        for(int j=0;j<bndi.rows();j++)
            if(bndi(j) == ci(i))
                bi(i) = j;
        for(int j=0;j<bnd2i.rows();j++){
            if(bnd2i(j) == K[i])
                b2i(i) = j;
            if(bnd2i(j) == D[K[i]])
                ki(i) = j; 
        }
    }
    for(int k=0;k<ci.rows();k++){
        int i1 = bi(k);
        int i2 = bi((k+1)%bi.rows());
        Eigen::VectorXi r(bnd.rows());
        int c = 0;
        for(int z = i1;z%bnd.rows() != i2;z++){
            r(z-i1) = z%bnd.rows();
            c++;
        }
        Eigen::VectorXi ti;
        r.conservativeResize(c,1);
        igl::slice(bndi,r,1,ti);
        T.conservativeResize(T.rows()+ti.rows(),1);
        T.bottomRows(ti.rows())<<ti;

        r.conservativeResize(c-1,1);
        Eigen::MatrixXd L1,L2;
        igl::slice(bnd,r,1,L1);

        Eigen::VectorXd rt;
        line_ratio(L1,rt);
        int j1 = ki(k);
        int j2 = b2i((k+1)%b2i.rows());
        Eigen::VectorXi r2(bnd2.rows());
        c = 0;
        for(int z=j1;z%bnd2.rows() != j2;z++){
            r2(z-j1) = z%bnd2.rows();
            c++;
        }
        r2.conservativeResize(c,1);
        igl::slice(bnd2,r2,1,L2);
        Eigen::MatrixXd SL;
        sample_polyline(L2,rt,SL);            
        P0.conservativeResize(P0.rows()+SL.rows()+1,2);
        P0.bottomRows(SL.rows()+1)<<V2.row(m[ci(k)]),SL;
    }
    // match polygon to bnd -> T
    
    // igl::opengl::glfw::Viewer viewer;
    // viewer.data().set_mesh(V2,F2);
    // for(int i=0;i<P0.rows();i++){
    //     viewer.data().add_points(P0.row(i),Eigen::RowVector3d(0,0,0));
    //     viewer.data().add_edges(P0.row(i),P0.row((i+1)%P0.rows()),Eigen::RowVector3d(0,0,0));
    // }
    // igl::opengl::glfw::imgui::ImGuiMenu menu;
    // viewer.plugins.push_back(&menu);
    // viewer.launch();

}

// generate boundary polygon based on constraint point position
// p only constains boundary constraints
// c(i) stores index of boundary that p_i is on
void target_polygon(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& p, 
    const Eigen::VectorXi& pi,
    std::vector<Eigen::MatrixXd>& poly
){
    cout<<"target polygon"<<endl;
    // generate background mesh [1. rec; 2.hilbert curve]
    Eigen::MatrixXd VV;
    Eigen::MatrixXi FF;
    map<int,int> m; // map from 3d to 2d
    //#define HILBERT
    #ifndef HILBERT
    background_mesh(REC_EDGE_NUM/4,p,pi,m,VV,FF);
    #else
    generate_hilbert(3,p,pi,m,VV,FF);
    #endif
    //Eigen::VectorXi b(pi.rows());
    //std::vector<std::vector<int>> bdt,bds;
    Eigen::VectorXi bdt,bds;
    igl::boundary_loop(FF,bdt);
    igl::boundary_loop(F,bds);
    // sample the boundary of FF
    
    // Eigen::MatrixXd bds_pos;
    // igl::slice(V,bds,1,bds_pos);
    // outer_boundary(pi,m,bds_pos,bds,VV,FF,P0,T);
    // display_poly(P0);
    
    int s1 = bds(0);
    int s2 = bds(0);
    int t1 = bdt(0);
    int t2 = bdt(0);
    // need also insert the first and last point (closed loop)
    Eigen::MatrixXd L1,L2;
    igl::slice(V ,bds,1,L1);
    igl::slice(VV,bdt,1,L2);
    Eigen::VectorXd rt,rs;
    //#define PRESERVE_CORNER
    #ifndef HILBERT
    L1.conservativeResize(L1.rows()+1,3);
    L1.bottomRows(1)<<V.row(bds(0));
    L2.conservativeResize(L2.rows()+1,Eigen::NoChange);
    L2.bottomRows(1)<<VV.row(bdt(0));
    Eigen::MatrixXd PTS;
    line_ratio(L1,rt);
    sample_polyline(L2,rt,PTS);
    PTS.conservativeResize(PTS.rows()-1,Eigen::NoChange);
    poly.push_back(PTS);
    #else
    // evenly sample between every point of initial space filling curve
    #ifdef EVEN
    int nv_to_insert = L1.rows()-L2.rows();
    Eigen::VectorXi nv_each_edge(L2.rows());
    nv_each_edge.setConstant(nv_to_insert/L2.rows());
    int nv_left = nv_to_insert % L2.rows();
    for(int i=0;i<nv_left;i++)
        nv_each_edge(i) += 1;
    if(nv_each_edge.sum()!=L1.rows()-L2.rows()){
        std::cout<<"sample wrong"<<std::endl;
        exit(0);
    }
    int s = 0;
    Eigen::MatrixXd sampled(L1.rows()+1,2);
    for(int k=0;k<L2.rows();k++){
        // subdivide edge (k,k+1)
        Eigen::VectorXd rs(nv_each_edge(k)+2);
        for(int j=0;j<rs.rows();j++){
            rs(j) = 1.0*j/(rs.rows()-1);
        }
        Eigen::MatrixXd local_L2;
        Eigen::MatrixXd L2_t(2,2);
        L2_t<<L2.row(k),L2.row((k+1)%L2.rows());
        sample_polyline(L2_t,rs,local_L2);
        sampled.block(s,0,local_L2.rows(),2) = local_L2;
        s += (nv_each_edge(k)+1);
    }
    #else
    std::cout<<"sample goal:"<<L1.rows()<<std::endl;
    L1.conservativeResize(L1.rows()+1,3);
    L1.bottomRows(1)<<V.row(bds(0));
    L2.conservativeResize(L2.rows()+1,Eigen::NoChange);
    L2.bottomRows(1)<<VV.row(bdt(0));
    line_ratio(L1,rs);
    line_ratio(L2,rt);
    Eigen::MatrixXd sampled(L1.rows(),2);
    int last_s = 0;
    for(int k=0;k<L2.rows()-1;k++){
        double rk = rt(k);
        double rk1 = rt(k+1);
        int s1=last_s;
        int s2=-1;
        for(int j=0;j<L1.rows()-1;j++){
            // if(rk >= rs(j) && rk < rs(j+1))
            //     s1 = j;
            if(rk1 >= rs(j) && rk1 < rs(j+1))
                s2 = j;
        }
        if(s1 >= s2 || s2 == -1)
            s2 = (s1+1)%L1.rows();
        if(k == L2.rows()-2)
            s2 = L1.rows()-1;
        last_s = s2;
        //std::cout<<s1<<","<<s2<<std::endl;
        // sample between s1, s2
        Eigen::VectorXd local_rs;
        Eigen::MatrixXd local_L2;
        int l = std::max(s2-s1+1,2);
        line_ratio(L1.block(s1,0,l,3),local_rs);
        sample_polyline(L2.block(k,0,2,2),local_rs,local_L2);
        sampled.block(s1,0,local_L2.rows(),2) = local_L2;
    }
    #endif
    sampled.conservativeResize(sampled.rows()-1,Eigen::NoChange);
    Eigen::MatrixXi F_;
    std::vector<std::pair<int,int>> edges;
    std::vector<int> points;
    for(int i=0;i<sampled.rows();i++){
        points.push_back(i);
        edges.push_back(std::make_pair(i,(i+1)%sampled.rows()));
    }
    //display(sampled,F_,edges,points,{},{});
    poly.push_back(sampled);
    #endif
    std::cout<<"sample end"<<std::endl;
}
