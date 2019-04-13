#include "local_operation.h"
#include <igl/copyleft/cgal/orient2D.h>
#include <igl/sortrows.h>
#include <unordered_set>

void remove_empty_faces(Eigen::MatrixXi& F){
    Eigen::MatrixXi F2(F.rows(),3);
    int m = 0;
    for(size_t f = 0;f<F.rows();f++)
    {
        if(F(f,0) != IGL_COLLAPSE_EDGE_NULL ||
           F(f,1) != IGL_COLLAPSE_EDGE_NULL ||
           F(f,2) != IGL_COLLAPSE_EDGE_NULL){
            F2.row(m) = F.row(f);
            m++;
        }
    }
    F2.conservativeResize(m,F2.cols());
    F = F2;
}

std::vector<int> circulation(
    const int e,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi& dEF,
    const Eigen::MatrixXi& dEI,
    const Eigen::VectorXi& EE,
    const Eigen::MatrixXi& allE
){
    using namespace std;
    // prepare output
    std::vector<int> N;
    N.reserve(6);
    const int m = F.rows();

    const auto & step = [&](
            const int e,
            const int ff,
            int & ne,
            int & nf,
            bool flip)
    {
        const int nv = dEI(e,1);
        // get next face
        nf = dEF(e,1);
        if(nf == -1){
            ne = -1;
            return;
        }
        // get next edge
        int dir = flip? -1 : 1;
        ne = nf+m*((nv+3+dir)%3);
    };

    const auto & cycle = [&](const int e, bool flip)->bool{
        const int f0 = dEF(e,0);
        int fi = f0;
        int ei = e;
        while(true) {
            step(ei, fi, ei, fi, flip);
            if (fi == -1){
                return false;
            } // reach boundary
            N.push_back(fi);
            if (fi == f0) {
                assert(ei == e);
                return true;
            }
        }

    };
    // cycle option: true -> cycle around d
    //               false-> cycle around s
    bool succ = cycle(e,false);
    if(!succ) {
        cycle(EE(e), true);
    }
    return N;
}

void neighbor(const int e, bool op, 
    const Eigen::MatrixXi& F, 
    const Eigen::MatrixXi& dEF,
    const Eigen::MatrixXi& dEI,
    const Eigen::VectorXi& EE,
    const Eigen::MatrixXi& allE, 
    std::vector<int>& NF
){
    // e = (s->d)
    // op indicate the neighbor of s(true) or d(false)
    auto switch_edge = [&](const int e,int& ne, bool flip){
        const int cf = dEF(e,0);
        const int cv = dEI(e,0);
        const int dir = flip? -1: 1;
        int pe = cf+F.rows()*((cv+dir+3)%3);
        ne = EE(pe);
    };
    int ee=op?e:EE(e);
    if(dEF(ee,1)==-1){
        switch_edge(e,ee,op);
        if(!op) ee = EE(ee);
    }
    NF = circulation(ee,F,dEF,dEI,EE,allE);
    if(NF.size()==0)
        NF.emplace_back(dEF(ee,0));  
};

bool edge_collapse_is_valid(
    const int ae,
    const Eigen::RowVectorXd & p,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXi& dEF,
    const Eigen::MatrixXi& dEI,
    const Eigen::VectorXi& EE,
    const Eigen::MatrixXi& allE,
    std::vector<int>& NF)
{
    using namespace Eigen;
    using namespace std;
    // For consistency with collapse_edge.cpp, let's determine edge flipness
    // (though not needed to check validity)
    // source and destination
    const int s = allE(ae,0);
    const int d = allE(ae,1);

    if(s == IGL_COLLAPSE_EDGE_NULL && d == IGL_COLLAPSE_EDGE_NULL)
    {
        return false;
    }

    auto collect_vertices = [&](vector<int>& NF, VectorXi& v){
        vector<int> NV,uNV;
        for(auto f : NF) {
            NV.push_back(F(f,0));
            NV.push_back(F(f,1));
            NV.push_back(F(f,2));
        }
        vector<size_t> _1,_2;
        igl::unique(NV,uNV,_1,_2);
        igl::list_to_matrix(uNV,v);
    };

    // check flip
    bool flipped = false;
    // check connectivity
    bool valid = false;
    // get neighbor faces for s and d
    // base case: if ae has a neighbor edge is also bd edge
    vector<int> Ns,Nd;
    VectorXi Ns_v,Nd_v;
    
    neighbor(ae,true ,F,dEF,dEI,EE,allE,Ns);
    neighbor(ae,false,F,dEF,dEI,EE,allE,Nd);
    collect_vertices(Ns,Ns_v);
    collect_vertices(Nd,Nd_v);
    VectorXi Nint = igl::intersect(Ns_v,Nd_v);
    if((dEF(ae,1) == -1 && Nint.size() == 3) ||
       (dEF(ae,1) != -1 && Nint.size() == 4))
        valid = true;
    Ns.insert(Ns.begin(),Nd.begin(),Nd.end());
    NF = Ns;
    std::sort(NF.begin(),NF.end());
    NF.erase( std::unique(NF.begin(),NF.end()),NF.end());
    int f0 = dEF(ae,0);
    int f1 = dEF(ae,1);
    auto loc = std::find(NF.begin(),NF.end(),f0);
    NF.erase(loc);
    loc = std::find(NF.begin(),NF.end(),f1);
    NF.erase(loc);
    // Eigen::VectorXi nf,nc(2),lf,IA_;
    // nc<<dEF(ae,0),dEF(ae,1);
    // igl::list_to_matrix(NF,nf);
    // igl::setdiff(nf,nc,lf,IA_);
    // Eigen::MatrixXd LV = V;
    // LV.row(s) << p;
    // LV.row(d) << p;
    // Eigen::MatrixXi LF(lf.rows(),3);
    // for(int i=0;i<LF.rows();i++)
    //     LF.row(i)<<F.row(lf(i));
    // auto test_flip = [](const Eigen::MatrixXd& V, const Eigen::MatrixXi& F){
    //     bool flipped = false;
    //     int count = 0;
    //     for(int i=0;i<F.rows();i++){
    //         double a[2] = {V(F(i,0),0),V(F(i,0),1)};
    //         double b[2] = {V(F(i,1),0),V(F(i,1),1)};
    //         double c[2] = {V(F(i,2),0),V(F(i,2),1)};
    //         if(igl::copyleft::cgal::orient2D(a, b, c) < 0){
    //             count++;
    //             flipped = true;
    //         }
    //     }
    //     return flipped;
    // };
    //flipped = test_flip(LV,LF);
    //flipped = false;
    if(!valid) {
        return false;
    }
    return true;
}

bool collapse_edge(
    const int ae,
    const Eigen::RowVectorXd & p,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::MatrixXi& dEF,
    Eigen::MatrixXi& dEI,
    Eigen::VectorXi& EE,
    Eigen::MatrixXi& allE
){
    // Assign this to 0 rather than, say, -1 so that deleted elements will get
    // draw as degenerate elements at vertex 0 (which should always exist and
    // never get collapsed to anything else since it is the smallest index)
    using namespace Eigen;
    using namespace std;
    int s = allE(ae,0);
    int d = allE(ae,1);
    // Important to grab neighbors of d before monkeying with edges

    // if non is control points, then pick the smaller
    bool choose_d = d < s ? true : false;
    int to_drop = choose_d ? s : d;
    int to_keep = choose_d ? d : s;
    // std::cout<<to_drop<<" is dropped "<<std::endl;
    // std::cout<<to_keep<<" is kept "<<std::endl;
    vector<int> neighbor_of_d;
    neighbor(ae,choose_d,F,dEF,dEI,EE,allE,neighbor_of_d);

    // move source and destination to midpoint
    V.row(s) = p;
    V.row(d) = p;

    const auto & kill_face = [&](const int f){
        for(int i=0;i<3;i++){
            int e = f+F.rows()*i;
            for(int j=0;j<2;j++){
                dEF(e,j) = IGL_COLLAPSE_EDGE_NULL;
                dEI(e,j) = IGL_COLLAPSE_EDGE_NULL;
                allE(e,j) = IGL_COLLAPSE_EDGE_NULL;
            }
        }
        F(f,0) = IGL_COLLAPSE_EDGE_NULL;
        F(f,1) = IGL_COLLAPSE_EDGE_NULL;
        F(f,2) = IGL_COLLAPSE_EDGE_NULL;
    };

    // update edge info
    // for each flap
    const int m = F.rows();
    vector<int> edges={ae};
    if(dEF(ae,1)!=-1){
        edges.emplace_back(EE(ae));
    }
    for(int e: edges)
    {
        const int f = dEF(e,0);
        const int v = dEI(e,0);
        // next edge emanating from d
        const int ae1 = f+m*((v+1)%3);
        // prev edge pointing to s
        const int ae2 = f+m*((v+2)%3);
        // face adjacent to f on e1, also incident on d
        const int f1 = dEF(ae1,1); // -1
        // across from which vertex of f1 does e1 appear?
        const int v1 = dEI(ae1,1); // -1
        const int f2 = dEF(ae2,1); 
        const int v2 = dEI(ae2,1);

        int oe1,oe2;
        oe1 = EE(ae1);
        oe2 = EE(ae2);
        dEF(oe1,1) = f2;
        dEF(oe2,1) = f1;
        dEI(oe1,1) = v2;
        dEI(oe2,1) = v1;
        
        // if ae1 is bd, then EE(ae1) will be rubbish
        if(oe2 != ae2)
            EE(oe1) = oe2; 
        else
            EE(oe1) = oe1;
        
        // if ae2 is bd, then EE(ae2) will be rubbish
        if(oe1 != ae1)
            EE(oe2) = oe1;
        else
            EE(oe2) = oe2;
        // kill everything inside f
        kill_face(f);
    }

    for(auto f : neighbor_of_d)
    {
        for(int v = 0;v<3;v++)
        {
            // this edge should be
            if(F(f,v) == to_drop)
            {
                int e1=f+m*((v+1)%3);
                int e2=f+m*((v+2)%3);
                allE(e1,1) = to_keep;
                allE(e2,0) = to_keep;
                F(f,v) = to_keep;
                break;
            }
        }
    }

    return true;
}