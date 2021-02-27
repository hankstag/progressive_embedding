#include "decompose_polygon.h"
#include "is_simple_polygon.h"
#include "embed_points.h"
#include "plot.h"
#include <igl/triangle_triangle_adjacency.h>
#include <igl/copyleft/cgal/orient2D.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Partition_traits_2<K> Traits;
typedef Traits::Point_2 Point;
typedef Traits::Polygon_2 Polygon_2;

int find_next_mark(const Eigen::VectorXi &mark, int start)
{
  int size = mark.rows();
  int next = (start + 1) % size;

  while (mark(next) != 1)
    next = (next + 1) % size;

  return next;
}

bool is_convex(
    const Eigen::MatrixXd &P)
{
  for (int i = 0; i < P.rows(); i++)
  {
    int prev = (i - 1 + P.rows()) % P.rows();
    int next = (i + 1) % P.rows();
    double a[2] = {P(prev, 0), P(prev, 1)};
    double b[2] = {P(i, 0), P(i, 1)};
    double c[2] = {P(next, 0), P(next, 1)};
    short r = igl::copyleft::cgal::orient2D(a, b, c);
    if (r < 0)
      return false;
  }
  return true;
}

void merge_triangles(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const Eigen::VectorXi &mark,
    std::vector<std::vector<int>> &L)
{

  // use the greedy method for now
  // could be improved

  // collect halfedge and corresponding face id
  // initialize L to be collection of all faces(polygons)
  std::map<std::pair<int, int>, int> H;
  Eigen::VectorXi G(F.rows());
  Eigen::MatrixXi FF, FFI;
  igl::triangle_triangle_adjacency(F, FF, FFI);
  L.resize(F.rows());
  for (int i = 0; i < F.rows(); i++)
  {
    G(i) = i; // every face belongs to itself
    for (int k = 0; k < 3; k++)
    {
      L[i].push_back(F(i, k));
      H[std::make_pair(F(i, k), F(i, (k + 1) % 3))] = FF(i, k);
    }
  }

  // traverse all the halfedges
  for (auto h : H)
  { // [he, pid]
    auto he = h.first;
    auto rhe = std::make_pair(he.second, he.first);
    int p1 = H[he];  // a -> b
    int p2 = H[rhe]; // b -> a
    if (p1 == -1 || p2 == -1)
      continue;

    // up to group root
    while (p1 != G[p1])
      p1 = G[p1];
    while (p2 != G[p2])
      p2 = G[p2];

    // combine p1 and p2
    Eigen::MatrixXd poly(L[p1].size() + L[p2].size() - 2, 2);
    auto a = std::find(L[p1].begin(), L[p1].end(), he.first);
    auto b = std::find(L[p2].begin(), L[p2].end(), he.second);

    std::vector<int> L1(L[p1].begin(), a);
    std::vector<int> R1(a, L[p1].end());
    std::vector<int> L2(L[p2].begin(), b);
    std::vector<int> R2(b, L[p2].end());

    std::vector<int> S;
    S = L1;
    auto c = R2.empty() ? R2.begin() : R2.begin() + 1;
    auto d = R1.empty() ? R1.begin() : R1.begin() + 1;
    S.insert(S.end(), c, R2.end());
    S.insert(S.end(), L2.begin(), L2.end());
    S.insert(S.end(), d, R1.end());

    for (int i = 0; i < poly.rows(); i++)
    {
      poly.row(i) << V.row(S[i]);
    }

    // if merged polygon is simple, drop edge he/rhe
    // erase L[p2], add to L[p1]
    Polygon_2 pgn;
    for (int i = 0; i < poly.rows(); i++)
      pgn.push_back(Point(poly(i, 0), poly(i, 1)));

    // if(is_simple_polygon(poly) && is_convex(poly)){
    if (pgn.is_simple())
    {
      H[he] = -1;
      H[rhe] = -1;
      G[p2] = p1; // p2 belongs to p1 now
      L[p1] = S;
      L[p2].clear();
    }
  }

  std::cout << "TODO here" << std::endl;

  std::cout << "L.size() = " << L.size() << std::endl;
  // TODO: get simplified L
  std::vector<std::vector<int>> L_sim;
  std::vector<int> L_size_full; // for partition check
  std::vector<std::vector<int>> L_full;
  for (int ii = 0; ii < L.size(); ii++)
  {
    if (L[ii].size() == 0)
      continue;
    std::vector<int> poly_sim; // simplified polygon
    std::cout << "L[" << ii << "].size() = " << L[ii].size() << std::endl;
    Eigen::MatrixXd PL, PL_sim;
    for (int j = 0; j < L[ii].size(); j++)
    {
      int v_id = L[ii][j], v_prev = L[ii][(j + L[ii].size() - 1) % L[ii].size()], v_next = L[ii][(j + 1) % L[ii].size()];
      if (mark(v_id) == 1) // non_subdivide node
      {
        poly_sim.push_back(v_id);
      }
      else 
      {
        int index0 = (v_next - v_id + V.rows()) % V.rows(), index1 = (find_next_mark(mark, v_id) - v_id + V.rows()) % V.rows();
        if (index0 > index1) 
          poly_sim.push_back(v_id);
        else
        {
          index0 = (v_id - v_prev + V.rows()) % V.rows(), index1 = (find_next_mark(mark, v_prev) - v_prev + V.rows()) % V.rows();
          if (index0 > index1)
            poly_sim.push_back(v_id);
        }
        // double l1 = (V.row(v_id) - V.row(v_prev)).norm();
        // double l2 = (V.row(v_id) - V.row(v_next)).norm();
        // double l3 = (V.row(v_next) - V.row(v_prev)).norm();
        // double cos_a = (l1 * l1 + l2 * l2 - l3 * l3) / (2 * l1 * l2);
        // if (fabs(cos_a + 1) > 0.01)
        // {
        //   poly_sim.push_back(v_id);
        //   std::cout << "keep subdivide node " << v_id << ", angle " << std::acos(cos_a) / igl::PI * 180 << std::endl;
        // }
      }
      PL.conservativeResize(PL.rows() + 1, 2);
      PL.row(PL.rows() - 1) = V.row(v_id);
    }
    // igl::opengl::glfw::Viewer viewer;
    // Eigen::VectorXi TL(PL.rows());
    // TL.setConstant(1);
    // viewer.data().set_mesh(V, F);
    // plot_polygon(viewer, TL, PL);
    // viewer.launch();
    std::cout << "size of simplified polygon: " << poly_sim.size() << std::endl;
    for (int v_id : poly_sim)
    {
      PL_sim.conservativeResize(PL_sim.rows() + 1, 2);
      PL_sim.row(PL_sim.rows() - 1) = V.row(v_id);
    }
    // igl::opengl::glfw::Viewer viewer_sim;
    // Eigen::VectorXi TL_sim(PL_sim.rows());
    // TL_sim.setConstant(1);
    // viewer_sim.data().set_mesh(V, F);
    // plot_polygon(viewer_sim, TL_sim, PL_sim);
    // viewer_sim.launch();
    L_sim.push_back(poly_sim);
    L_size_full.push_back(L[ii].size());
    L_full.push_back(L[ii]);
  }

  L = L_sim; // use the simplified polygon

  std::vector<std::vector<int>> L_new;
  for (int ii = 0; ii < L.size(); ii++)
  {
    std::cout << "\nL[" << ii << "].size() = " << L[ii].size() << std::endl;
    if (L[ii].size() == 0) // not likely to happen for L_sim
      continue;
    // construct pgn for cgal partition
    Polygon_2 pgn;
    for (int v_id : L[ii])
    {
      pgn.push_back(Point(V(v_id, 0), V(v_id, 1)));
    }
    std::vector<Polygon_2> partition_polys;
    CGAL::approx_convex_partition_2(pgn.vertices_begin(),
                                    pgn.vertices_end(),
                                    std::back_inserter(partition_polys));
    int cnt_sim = 0, cnt_full = 0;
    for (int i = 0; i < partition_polys.size(); i++)
    {
      auto poly = partition_polys[i];
      std::vector<int> poly_id;
      Eigen::MatrixXd P;
      int index = 0;
      for (auto it = poly.vertices_begin(); it != poly.vertices_end(); it++)
      {
        P.conservativeResize(P.rows() + 1, 2);
        P.row(index) << it->x(), it->y();
        index++;
      }

      for (int row_iter = 0; row_iter < P.rows(); row_iter++)
      {
        for (int v_id : L[ii])
        {
          if (V.row(v_id) == P.row(row_iter))
          {
            poly_id.push_back(v_id);
          }
        }
      }

      std::cout << "poly_id.size() = " << poly_id.size() << "\tP.rows() = " << P.rows() << std::endl;
      if (poly_id.size() != P.rows())
        std::cout << "find v_id error" << std::endl;
      if (poly_id.size() < 3)
        std::cout << "degenerate polygon error" << std::endl;
      // std::cout << "P:" << std::endl
      //           << std::setprecision(17) << P << std::endl;

      // subdivide poly_id
      std::vector<int> poly_sub;
      for (int k = 0; k < poly_id.size(); k++)
      {
        int v0 = poly_id[k], v1 = poly_id[(k + 1) % poly_id.size()];
        poly_sub.push_back(v0);
        int index0 = (v1 - v0 + V.rows()) % V.rows(), index1 = (find_next_mark(mark, v0) - v0 + V.rows()) % V.rows();
        std::cout << "v0 = " << v0 << "," << mark(v0) << "\tv1 = " << v1 << "," << mark(v1) << std::endl; 
        std::cout << "v0_to_v1 = " << index0 << " v0_to_next_mark = " << index1 << std::endl;
        
        if (index0 <= index1)
        {
          int kk = (v0 + 1) % V.rows();
          while (kk != v1)
          {
            poly_sub.push_back(kk);
            kk = (kk + 1) % V.rows();
          }
        }
      }
      L_new.push_back(poly_sub);
      std::cout << "subdivide poly size: " << poly_sub.size() << std::endl;
      cnt_sim += poly_id.size();
      cnt_full += poly_sub.size();
      // igl::opengl::glfw::Viewer viewer;
      // Eigen::VectorXi T(P.rows());
      // T.setConstant(1);
      // viewer.data().set_mesh(V, F);
      // plot_polygon(viewer, T, P);
      // viewer.launch();

      Eigen::MatrixXd Pfull(poly_sub.size(), 2);
      for (int k = 0; k < poly_sub.size(); k++)
      {
        Pfull.row(k) = V.row(poly_sub[k]);
      }
      // igl::opengl::glfw::Viewer viewerfull;
      // Eigen::VectorXi Tfull(Pfull.rows());
      // Tfull.setConstant(1);
      // viewerfull.data().set_mesh(V, F);
      // plot_polygon(viewerfull, Tfull, Pfull);
      // viewerfull.launch();
      for (auto v : poly_id) std::cout << v << " ";
      std::cout << std::endl;
      for (auto v : poly_sub) std::cout << v << " ";
      std::cout << std::endl;
    }
    std::cout << "check subdivide:" << L_size_full[ii] - cnt_full << " " <<(int)L[ii].size() - cnt_sim << std::endl;
    for (auto v : L_full[ii]) std::cout << v << " ";
    std::cout << std::endl;
    for (auto v : L[ii]) std::cout << v << " ";
    std::cout << std::endl;
  } 
  L = L_new;
}

void decompose_polygon(
    const Eigen::MatrixXd &P,
    const Eigen::VectorXi &R,
    const Eigen::MatrixXd &C,
    const Eigen::VectorXi &mark,
    Eigen::MatrixXd &V,
    Eigen::MatrixXi &F,
    std::vector<std::vector<int>> &L)
{

  bool succ = Shor_van_wyck_v2(P, R, mark, "", V, F, false);
  assert(succ && "Shor failed");
  embed_points(C, V, F);
  igl::opengl::glfw::Viewer vr;
  std::cout << V.rows() << " " << mark.rows() << std::endl;
  merge_triangles(V, F, mark, L);
  std::cout << "end of decompose polygon" << std::endl;
}
