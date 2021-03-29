#include "shor.h"
#include "plot.h"
#include <igl/boundary_loop.h>
#include <igl/adjacency_list.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/predicates/predicates.h>
#include "is_simple_polygon.h"
#include <iostream>
#include <unsupported/Eigen/MPRealSupport>

using namespace std;

short orientation_flt(const Eigen::Matrix<double, 3, 2> &P)
{
  Eigen::Matrix<double, 1, 2> a_db(P(0, 0), P(0, 1));
  Eigen::Matrix<double, 1, 2> b_db(P(1, 0), P(1, 1));
  Eigen::Matrix<double, 1, 2> c_db(P(2, 0), P(2, 1));
  return short(igl::predicates::orient2d(a_db, b_db, c_db));
}

short orientation_mpf(const Eigen::Matrix<mpfr::mpreal, 3, 2> &P)
{
  using mpfr::mpreal;
  mpfr_prec_t n_bits = mpreal::get_default_prec();
  mpreal::set_default_prec(n_bits * 8);
  mpreal a0, a1, b0, b1, c0, c1;
  a0 = P(0, 0);
  a1 = P(0, 1);
  b0 = P(1, 0);
  b1 = P(1, 1);
  c0 = P(2, 0);
  c1 = P(2, 1);
  a0.set_prec(n_bits * 8);
  a1.set_prec(n_bits * 8);
  b0.set_prec(n_bits * 8);
  b1.set_prec(n_bits * 8);
  c0.set_prec(n_bits * 8);
  c1.set_prec(n_bits * 8);
  mpreal signed_area = a0 * b1 - b0 * a1 +
                       b0 * c1 - c0 * b1 +
                       c0 * a1 - a0 * c1;
  short code = -9999;
  if (signed_area < 0.0)
  {
    code = -1;
  }
  else if (signed_area == 0.0)
  {
    code = 0;
  }
  else
  {
    code = 1;
  }
  mpreal::set_default_prec(n_bits);
  assert(code != -9999);
  return code;
}

template <typename Scalar>
short orientation(const Eigen::Matrix<Scalar, 3, 2> &P)
{
  if (std::is_same<Scalar, mpfr::mpreal>::value)
  {
    Eigen::Matrix<mpfr::mpreal, 3, 2> P_mpf = P.template cast<mpfr::mpreal>();
    return orientation_mpf(P_mpf);
  }
  else
  {
    Eigen::Matrix<double, 3, 2> P_flt = P.template cast<double>();
    return orientation_flt(P_flt);
  }
}

template <typename Scalar>
Scalar signed_area(const Eigen::Matrix<Scalar, 3, 2> &P)
{
  if (std::is_same<Scalar, mpfr::mpreal>::value)
  {
    using mpfr::mpreal;
    mpfr_prec_t n_bits = mpreal::get_default_prec();
    mpreal::set_default_prec(n_bits * 8);
    mpreal a0, a1, b0, b1, c0, c1;
    a0 = mpfr::mpreal(P(0, 0));
    a1 = mpfr::mpreal(P(0, 1));
    b0 = mpfr::mpreal(P(1, 0));
    b1 = mpfr::mpreal(P(1, 1));
    c0 = mpfr::mpreal(P(2, 0));
    c1 = mpfr::mpreal(P(2, 1));
    a0.set_prec(n_bits * 8);
    a1.set_prec(n_bits * 8);
    b0.set_prec(n_bits * 8);
    b1.set_prec(n_bits * 8);
    c0.set_prec(n_bits * 8);
    c1.set_prec(n_bits * 8);
    mpfr::mpreal area = (a0 * b1 - b0 * a1 +
                         b0 * c1 - c0 * b1 +
                         c0 * a1 - a0 * c1);
    mpreal::set_default_prec(n_bits);
    return Scalar(area);
  }
  else
  {
    double a0, a1, b0, b1, c0, c1;
    a0 = double(P(0, 0));
    a1 = double(P(0, 1));
    b0 = double(P(1, 0));
    b1 = double(P(1, 1));
    c0 = double(P(2, 0));
    c1 = double(P(2, 1));
    return a0 * b1 - b0 * a1 +
           b0 * c1 - c0 * b1 +
           c0 * a1 - a0 * c1;
  }
}

template <typename Scalar>
bool segment_segment_intersect(
    const Eigen::Matrix<Scalar, 1, 2> &a,
    const Eigen::Matrix<Scalar, 1, 2> &b,
    const Eigen::Matrix<Scalar, 1, 2> &c,
    const Eigen::Matrix<Scalar, 1, 2> &d)
{
  Eigen::Matrix<Scalar, 3, 2> T1, T2, T3, T4;
  T1 << a, b, c;
  T2 << b, c, d;
  T3 << a, b, d;
  T4 << a, c, d;

  auto t1 = orientation(T1);
  auto t2 = orientation(T2);
  auto t3 = orientation(T3);
  auto t4 = orientation(T4);

  // assume m,n,p are colinear, check whether p is in range [m,n]
  auto on_segment = [](
                        const Eigen::Matrix<Scalar, 1, 2> &m,
                        const Eigen::Matrix<Scalar, 1, 2> &n,
                        const Eigen::Matrix<Scalar, 1, 2> &p) {
    return ((p(0) >= std::min(m(0), n(0))) &&
            (p(0) <= std::max(m(0), n(0))) &&
            (p(1) >= std::min(m(1), n(1))) &&
            (p(1) <= std::max(m(1), n(1))));
  };

  // colinear case
  if ((t1 == 0 && on_segment(a, b, c)) ||
      (t2 == 0 && on_segment(c, d, b)) ||
      (t3 == 0 && on_segment(a, b, d)) ||
      (t4 == 0 && on_segment(c, d, a)))
    return true;

  // ordinary case
  return (t1 != t3 && t2 != t4);
}

// a wrapper for ear clipping functionalities that supports both mpf and floats
template <typename Scalar>
void ear_clipping(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &P,
    const Eigen::VectorXi &RT,
    Eigen::VectorXi &I,
    Eigen::MatrixXi &eF,
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &nP)
{
  // check whether vertex i is an ear
  auto is_ear = [](
                    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &P,
                    const Eigen::VectorXi &RT,
                    const Eigen::VectorXi &L,
                    const Eigen::VectorXi &R,
                    const int i) {
    int a = L(i), b = R(i);
    if (RT(i) != 0 || RT(a) != 0 || RT(b) != 0)
      return false;
    Eigen::Matrix<Scalar, 1, 2> pa = P.row(a);
    Eigen::Matrix<Scalar, 1, 2> pb = P.row(b);
    Eigen::Matrix<Scalar, 1, 2> pi = P.row(i);
    Eigen::Matrix<Scalar, 3, 2> T;
    T << pa, pi, pb;
    short r = orientation(T);
    if (r != 1)
      return false;

    int k = R(b);
    int l = R(k);
    while (k != a)
    {
      Eigen::Matrix<Scalar, 1, 2> q = P.row(k);
      Eigen::Matrix<Scalar, 1, 2> p = P.row(l);
      if (segment_segment_intersect(pa, pb, q, p) || segment_segment_intersect(pa, pi, q, p) || segment_segment_intersect(pi, pb, q, p))
      {
        return false;
      }
      k = l;
      l = R(k);
    }

    return true;
  };

  Eigen::VectorXi L(P.rows());
  Eigen::VectorXi R(P.rows());
  for (int i = 0; i < P.rows(); i++)
  {
    L(i) = int((i - 1 + P.rows()) % P.rows());
    R(i) = int((i + 1) % P.rows());
  }

  Eigen::VectorXi ears; // mark ears
  Eigen::VectorXi X;    // clipped vertices
  ears.setZero(P.rows());
  X.setZero(P.rows());

  // initialize ears
  for (int i = 0; i < P.rows(); i++)
  {
    ears(i) = is_ear(P, RT, L, R, i);
  }

  // clip ears until none left
  while (ears.maxCoeff() == 1)
  {

    // find the first ear
    int e = 0;
    while (e < ears.rows() && ears(e) != 1)
      e++;

    // find valid neighbors
    int a = L(e), b = R(e);
    if (a == b)
      break;

    // clip ear and store face
    eF.conservativeResize(eF.rows() + 1, 3);
    eF.bottomRows(1) << a, e, b;
    L(b) = a;
    L(e) = -1;
    R(a) = b;
    R(e) = -1;
    ears(e) = 0; // mark vertex e as non-ear

    // update neighbor's ear status
    ears(a) = is_ear(P, RT, L, R, a);
    ears(b) = is_ear(P, RT, L, R, b);
    X(e) = 1;

    // when only one edge left
    // mark the endpoints as clipped
    if (L(a) == b && R(b) == a)
    {
      X(a) = 1;
      X(b) = 1;
    }
  }

  // collect remaining vertices if theres any
  for (int i = 0; i < X.rows(); i++)
    X(i) = 1 - X(i);
  I.resize(X.sum());
  int j = 0;
  for (int i = 0; i < X.rows(); i++)
    if (X(i) == 1)
    {
      I(j++) = i;
    }
  igl::slice(P, I, 1, nP);
}

template <typename Scalar>
void set_rotation_index(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &uv,
    const Eigen::MatrixXi &F,
    Eigen::VectorXi &R,
    int offset)
{
  Eigen::VectorXi bd;
  igl::boundary_loop(F, bd);
  R.setZero(bd.rows());
  std::vector<std::vector<int>> A;
  igl::adjacency_list(F, A, true);

  // for every boundary vertex, update its rotation index
  for (int i = 0; i < bd.rows(); i++)
  {
    int v = bd((i + offset) % bd.rows());
    Angle<Scalar> sum = Angle<Scalar>();
    std::reverse(A[v].begin(), A[v].end());
    for (int j = 0; j < A[v].size() - 1; j++)
    {
      int id = A[v][j];
      int id_1 = A[v][j + 1];
      Angle<Scalar> I(uv.row(id), uv.row(v), uv.row(id_1)); // w1, v, w2
      sum = (j == 0) ? I : sum + I;
    }
    R(i) = sum.r;
  }
}

Eigen::VectorXi set_ri(
    const Eigen::MatrixXd &uv,
    const Eigen::MatrixXi &F)
{
  // std::cout<<F<<std::endl;
  int offset = 0;
  Eigen::VectorXi R;
  Eigen::VectorXi bd;
  igl::boundary_loop(F, bd);
  R.setZero(bd.rows());
  std::cout << "boundary lengths: " << bd.rows() << std::endl;
  std::vector<std::vector<int>> A;
  igl::adjacency_list(F, A, true);

  // for every boundary vertex, update its rotation index
  for (int i = 0; i < bd.rows(); i++)
  {
    int v = bd((i + offset) % bd.rows());
    Angle<double> sum = Angle<double>();
    std::reverse(A[v].begin(), A[v].end());
    for (int j = 0; j < A[v].size() - 1; j++)
    {
      int id = A[v][j];
      int id_1 = A[v][j + 1];
      Angle<double> I(uv.row(id), uv.row(v), uv.row(id_1)); // w1, v, w2
      sum = (j == 0) ? I : sum + I;
    }
    R(i) = sum.r;
  }
  return R;
}

template <typename Scalar>
Eigen::VectorXi set_all_ri(
    const Eigen::Matrix<Scalar, -1, -1> &uv,
    const Eigen::MatrixXi &F)
{
  std::vector<std::vector<int>> A;
  igl::adjacency_list(F, A, true);

  int offset = 0;
  std::vector<std::vector<int>> bds;
  igl::boundary_loop(F, bds);

  Eigen::VectorXi R;
  R.setZero(uv.rows());

  for (auto bd : bds)
  {
    std::cout << "boundary lengths: " << bd.size() << std::endl;
    // for every boundary vertex, update its rotation index
    for (int i = 0; i < bd.size(); i++)
    {
      int v = bd[(i + offset) % bd.size()];
      Angle<Scalar> sum = Angle<Scalar>();
      std::reverse(A[v].begin(), A[v].end());
      for (int j = 0; j < A[v].size() - 1; j++)
      {
        int id = A[v][j];
        int id_1 = A[v][j + 1];
        Angle<Scalar> I(uv.row(id), uv.row(v), uv.row(id_1)); // w1, v, w2
        sum = (j == 0) ? I : sum + I;
      }
      R(v) = sum.r;
    }
  }
  return R;
}

Eigen::VectorXi set_all_ri_str(
    const std::vector<std::vector<std::string>> &uv_str,
    const Eigen::MatrixXi &F)
{

  auto str_to_num = [](const std::vector<std::vector<std::string>> &uv_str, Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> &uv, int n_digits){

    using mpfr::mpreal;
    uv.resize(uv_str.size(), 2);

    mpreal::set_default_prec(mpfr::digits2bits(n_digits));
    std::cout << " set digits: " << n_digits << std::endl;

  #ifdef CHECK_MPF_CONFIG
    using namespace std;
    const mpreal one = 1.0;
    const mpreal zero = 0.0;
    const mpreal eps = std::numeric_limits<mpreal>::epsilon();
    const int base = std::numeric_limits<mpreal>::radix;
    const mpreal prec = eps * base;
    const int bindigits = std::numeric_limits<mpreal>::digits(); // eqv. to mpfr::mpreal::get_default_prec();
    const mpreal rnd = std::numeric_limits<mpreal>::round_error();
    const mpreal maxval = std::numeric_limits<mpreal>::max();
    const mpreal minval = std::numeric_limits<mpreal>::min();
    const mpreal small = one / maxval;
    const mpreal sfmin = (small > minval) ? small * (one + eps) : minval;
    const mpreal round = std::numeric_limits<mpreal>::round_style();
    const int min_exp = std::numeric_limits<mpreal>::min_exponent;
    const mpreal underflow = std::numeric_limits<mpreal>::min();
    const int max_exp = std::numeric_limits<mpreal>::max_exponent;
    const mpreal overflow = std::numeric_limits<mpreal>::max();

    // Additionally compute pi with required accuracy - just for fun :)
    const mpreal pi = mpfr::const_pi();

    cout.precision(100); // Show all the digits
    cout << "pi         =    " << pi << endl;
    cout << "eps        =    " << eps << endl;
    cout << "base       =    " << base << endl;
    cout << "prec       =    " << prec << endl;
    cout << "b.digits   =    " << bindigits << endl;
    cout << "rnd        =    " << rnd << endl;
    cout << "maxval     =    " << maxval << endl;
    cout << "minval     =    " << minval << endl;
    cout << "small      =    " << small << endl;
    cout << "sfmin      =    " << sfmin << endl;
    cout << "1/sfmin    =    " << 1 / sfmin << endl;
    cout << "round      =    " << round << endl;
    cout << "max_exp    =    " << max_exp << endl;
    cout << "min_exp    =    " << min_exp << endl;
    cout << "underflow  =    " << underflow << endl;
    cout << "overflow   =    " << overflow << endl;
  #endif
    for (int i = 0; i < uv_str.size(); i++)
    {
      for (int k = 0; k < 2; k++)
      {
        uv(i, k) = mpfr::mpreal(uv_str[i][k]);
      }
    }
    std::cout << "prec passed in from c++: " << uv(0, 0).get_prec() << std::endl;
  };
  int n_digits = uv_str[0][0].size();
  Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> uv;
  str_to_num(uv_str, uv, n_digits);
  std::vector<std::vector<int>> A;
  igl::adjacency_list(F, A, true);

  int offset = 0;
  std::vector<std::vector<int>> bds;
  igl::boundary_loop(F, bds);

  Eigen::VectorXi R;
  R.setZero(uv.rows());

  for (auto bd : bds)
  {
    std::cout << "boundary lengths: " << bd.size() << std::endl;
    // for every boundary vertex, update its rotation index
    for (int i = 0; i < bd.size(); i++)
    {
      int v = bd[(i + offset) % bd.size()];
      Angle<mpfr::mpreal> sum = Angle<mpfr::mpreal>();
      std::reverse(A[v].begin(), A[v].end());
      for (int j = 0; j < A[v].size() - 1; j++)
      {
        int id = A[v][j];
        int id_1 = A[v][j + 1];
        Angle<mpfr::mpreal> I(uv.row(id), uv.row(v), uv.row(id_1)); // w1, v, w2
        sum = (j == 0) ? I : sum + I;
      }
      R(v) = sum.r;
    }
  }
  return R;
}

void add_triangle(
    Eigen::MatrixXi &F,
    int i, int j,
    std::vector<std::vector<int>> &K,
    std::vector<std::vector<int>> &Q)
{
  if (K[i][j] == -1)
    return;
  int k = K[i][j];
  F.conservativeResize(F.rows() + 1, 3);
  F.bottomRows(1) << i, k, j;
  add_triangle(F, i, k, K, Q);
  add_triangle(F, k, j, K, Q);
}

void inspect_table(
    std::vector<std::vector<int>> &K,
    std::vector<std::vector<int>> &Q)
{

  // find the max distance between (j < i) for every i
  int max_size = 0;
  std::pair<int, int> max_ij;
  int N = K.size();
  for (int i = 0; i < N; i++)
  {
    for (int j = (i - 1 + N) % N; j != i; j = (j - 1 + N) % N)
    {
      if (Q[i][j] == 1)
      {
        int p_size = i > j ? (N + 1 - i + j) : (j - i + 1);
        if (max_size < p_size)
        {
          max_ij = std::make_pair(i, j);
          max_size = p_size;
        }
      }
    }
  }
  std::cout << "max polygon size: " << max_size << "/" << N << std::endl;
  std::cout << "range = [" << max_ij.first << "," << max_ij.second << "]" << std::endl;
}

template <typename Scalar>
bool add_to_table(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &P,
    const Eigen::VectorXi &R,
    std::vector<std::vector<Angle<Scalar>>> &FA,
    std::vector<std::vector<Angle<Scalar>>> &LA,
    int i, int k, int j)
{

  const int N = P.rows();
  Eigen::Matrix<Scalar, 1, 2> pi = P.row(i);
  Eigen::Matrix<Scalar, 1, 2> pk = P.row(k);
  Eigen::Matrix<Scalar, 1, 2> pj = P.row(j);

  Eigen::Matrix<Scalar, 3, 2> T;
  T << pi, pk, pj;
  if (orientation(T) != 1)
  {
    return false;

  }

  Angle<Scalar> I(pj, pi, pk); // J I K
  Angle<Scalar> J(pk, pj, pi); // K J I

  Angle<Scalar> k1 = LA[i][k];
  Angle<Scalar> k2(P.row(i), P.row(k), P.row(j));
  Angle<Scalar> k3 = FA[k][j];
  Angle<Scalar> Fr = I + FA[i][k];
  Angle<Scalar> La = LA[k][j] + J;
  bool ij_2gon = (j + 1) % P.rows() == i; // whether Pji is 2gon
  bool cond_i = ij_2gon ? (Fr.r == R(i)) : (Fr.r <= R(i));
  bool cond_j = ij_2gon ? (La.r == R(j)) : (La.r <= R(j));
  if (cond_i && cond_j && (k1 + k2 + k3).r == R(k))
  {
    // // update this outside
    // FA[i][j] = Fr;
    // LA[i][j] = La;

    return true;
  }
  else
  {

    return false;
  }
}

bool weakly_self_overlapping_str(
    const std::vector<std::vector<std::string>> &P_str,
    const Eigen::VectorXi &R,
    Eigen::MatrixXi &F)
{
  using mpfr::mpreal;
  int dps = P_str[0][0].size();
  std::cout << "dps inside wso check: " << dps << std::endl;
  mpreal::set_default_prec(mpfr::digits2bits(dps));
  Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic> P(P_str.size(), 2);
  std::cout << "actual prec: " << P(0, 0).get_prec() << std::endl;
  for (int i = 0; i < P_str.size(); i++)
  {
    for (int k = 0; k < 2; k++)
    {
      P(i, k) = mpreal(P_str[i][k]);
    }
  }
  return weakly_self_overlapping(P, R, F);
}

template <typename Scalar>
bool weakly_self_overlapping(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &P,
    const Eigen::VectorXi &R,
    Eigen::MatrixXi &F)
{
  std::cout << "weakly self overlapping test" << std::endl;

  const int N = P.rows();
  std::vector<std::vector<int>> K(N, std::vector<int>(N, -1));
  std::vector<std::vector<int>> Q(N, std::vector<int>(N, 0));
  std::vector<std::vector<Angle<Scalar>>> FA(N, std::vector<Angle<Scalar>>(N));
  std::vector<std::vector<Angle<Scalar>>> LA(N, std::vector<Angle<Scalar>>(N));

  // TODO: select the best triangulation
  Eigen::Matrix<Scalar, -1, -1> min_A(N, N);
  min_A.setConstant(1e10); // INF=1e10

  bool succ = false;
  // int i;
  for (int i = 0; i < N; i++)
  {
    Q[i][(i + 1) % N] = 1;
    FA[i][(i + 1) % N] = Angle<Scalar>(P.row((i + 1) % N), P.row(i), P.row((i + 1) % N));
    LA[i][(i + 1) % N] = Angle<Scalar>(P.row(i), P.row((i + 1) % N), P.row(i));
  }
  int h = -1;

  Scalar best_min_A = -1.0;
  // std::vector<std::vector<Scalar>> A(N, std::vector<Scalar>(N, -1));
  // for (int d = 2; d < N && !succ; d++)
  for (int d = 2; d < N; d++)
  {
    // std::cout << "d = " << d << "/" << N - 1 << std::endl;
    // for (i = 0; i < N && !succ; i++)
    for (int i = 0; i < N; i++)
    {
      // std::cout << "i = " << i << "/" << N - 1 << std::endl;
      int j = (i + d) % N;
      int k = (i + 1) % N;
      // Scalar min_q = -1.0f;
      while (k != j)
      {
        // std::cout << "k=" << k << "\tj=" << j << std::endl;
        if (Q[i][k] == 1 && Q[k][j] == 1)
        {
          if (add_to_table(P, R, FA, LA, i, k, j))
          {

            // Q[i][j] = 1;
            // K[i][j] = k;
            // seg min area
            Eigen::Matrix<Scalar, 1, 2> pi = P.row(i);
            Eigen::Matrix<Scalar, 1, 2> pk = P.row(k);
            Eigen::Matrix<Scalar, 1, 2> pj = P.row(j);
            Eigen::Matrix<Scalar, 3, 2> T;
            T << pi, pk, pj;
            Scalar tmp = signed_area(T);
            Scalar tmp1;
            tmp1 = min_A(i, k) < min_A(k, j) ? min_A(i, k) : min_A(k, j);
            tmp1 = tmp < tmp1 ? tmp : tmp1;

            if (Q[i][j] == 0 || tmp1 > min_A(i, j))
            {
              Q[i][j] = 1;
              K[i][j] = k;
              min_A(i, j) = tmp1;

              // update FA and LA (instead of updating this in add_to_table)
              Angle<Scalar> I(pj, pi, pk); // J I K
              Angle<Scalar> J(pk, pj, pi); // K J I
              FA[i][j] = I + FA[i][k];
              LA[i][j] = LA[k][j] + J;
            }

            if (i == (j + 1) % N)
            {
              // succ = true;
              // std::cout << "h=" << i << "successful with area_min = " << min_A(i, j) << std::endl;
              if (min_A(i, j) > best_min_A)
              {
                std::cout << "new best" << std::endl;
                best_min_A = min_A(i, j);
                h = i;
              }
              // break;
            }

          }
        }

        k = (k + 1) % N;
      }
    }
  }
  std::cout << "test done, h = " << h << "\nmin_area = " << best_min_A << std::endl;
  inspect_table(K, Q);
  if (h == -1)
    return false;
  add_triangle(F, h, (h - 1 + N) % N, K, Q);
  return true;
}

void drop_colinear(
    const Eigen::MatrixXd &P,
    const Eigen::VectorXi &R,
    Eigen::VectorXi &B,
    Eigen::MatrixXd &mP,
    Eigen::VectorXi &mR)
{
  int dropped = 0;
  int N = P.rows();
  std::vector<int> Bv;
  for (int i = 0; i < N; i++)
  {
    int a = (i - 1 + N) % N;
    int v = i;
    int b = (i + 1) % N;
    Eigen::Matrix<double, 3, 2> T;
    T << P.row(a), P.row(v), P.row(b);
    if (R(v) != 0 || orientation(T) != 0)
    { // non-colinear or rotate index nonzero
      Bv.push_back(v);
    }
    else
      dropped++;
  }
  igl::list_to_matrix(Bv, B);
  igl::slice(P, B, 1, mP);
  igl::slice(R, B, mR);
}

template <typename Scalar>
void drop_colinear_v2(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &P,
    const Eigen::VectorXi &R,
    const Eigen::VectorXi &mark,
    Eigen::VectorXi &B,
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &mP,
    Eigen::VectorXi &mR)
{
  int dropped = 0;
  int N = P.rows();
  std::vector<int> Bv;
  for (int i = 0; i < N; i++)
  {
    int a = (i - 1 + N) % N;
    int v = i;
    int b = (i + 1) % N;
    if (R(v) != 0 || mark(v) != 0)
    { // non-colinear or rotate index nonzero
      Bv.push_back(v);
    }
    else
      dropped++;
  }
  igl::list_to_matrix(Bv, B);
  igl::slice(P, B, 1, mP);
  igl::slice(R, B, mR);
}

template <typename Scalar>
void add_colinear(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &P,
    const Eigen::MatrixXi &nF,
    const Eigen::VectorXi &B,
    Eigen::MatrixXi &F)
{
  int added = 0;
  int sN = nF.maxCoeff() + 1; // size of simplified polygon
  F.resizeLike(nF);
  for (int i = 0; i < F.rows(); i++)
    F.row(i) << B(nF(i, 0)), B(nF(i, 1)), B(nF(i, 2));
  typedef std::pair<int, int> Edge;
  std::vector<Edge> stk;
  std::set<std::pair<int, int>> sid; // edge vertex pairs needs split
  for (int i = 0; i < nF.rows(); i++)
  {
    bool split = false;
    for (int k = 0; k < 3; k++)
    {
      int a = nF(i, k);
      int b = nF(i, (k + 1) % 3);
      if ((b - a == 1 || a - b == sN - 1) &&
          (B(b) - B(a) != 1 && B(a) - B(b) != P.rows() - 1))
      { // edge need split
        if (!split)
        { // do not split a face twice (do it later)
          stk.push_back(Edge(i, k));
          split = true;
        }
        sid.insert(std::make_pair(B(a), B(b)));
      }
    }
  }
  while (!stk.empty())
  {
    Edge e = stk.back();
    stk.pop_back();
    int a = F(e.first, e.second);
    int b = F(e.first, (e.second + 1) % 3);
    int c = F(e.first, (e.second + 2) % 3);
    int rg = a < b ? b - a : b + P.rows() - a;
    int r = F.rows();
    F.conservativeResize(F.rows() + rg, 3);
    added += (rg - 1);
    // from r+0 to r+rg-1
    for (int i = 0; i < rg; i++)
    {
      F.row(r + i) << (a + i) % P.rows(), (a + i + 1) % P.rows(), c;
    }
    // if bc needs split
    if (sid.find(std::make_pair(b, c)) != sid.end())
    {
      stk.push_back(Edge(r + rg - 1, 1)); // its now 2nd edge of face r+rg-1
    }
    // if ac needs split
    if (sid.find(std::make_pair(c, a)) != sid.end())
    {
      stk.push_back(Edge(r, 2)); // its now 3rd edge of face r
    }
    F.row(e.first) << -1, -1, -1; // delete face
  }
  Eigen::MatrixXi tF = F;
  // drop the (-1,-1,-1) rows
  int k = 0;
  for (int i = 0; i < F.rows(); i++)
    if (F.row(i).sum() != -3)
      tF.row(k++) << F.row(i);
  tF.conservativeResize(k, 3);
  F = tF;
}

void subdivide_polygon(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    std::vector<std::vector<int>> &L)
{
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
    if (is_simple_polygon(poly))
    {
      H[he] = -1;
      H[rhe] = -1;
      G[p2] = p1; // p2 belongs to p1 now
      L[p1] = S;
      L[p2].clear();
    }
  }
}

void simplify_triangulation(
    const Eigen::MatrixXd &V_i,
    const std::vector<std::vector<int>> &L_i,
    Eigen::MatrixXd &V,
    Eigen::MatrixXi &F)
{
  std::vector<std::vector<int>> L;
  V = V_i;
  F.resize(0, 0);
  std::map<std::pair<int, int>, std::vector<int>> E;
  // calculate the average edge length
  double avl = 0.0f;
  for (int i = 0; i < V_i.rows(); i++)
  {
    avl += (V_i.row((i + 1) % V_i.rows()) - V_i.row(i)).norm();
  }
  avl /= V_i.rows();
  // traverse L, split the internal edges
  for (int i = 0; i < L_i.size(); i++)
  {
    if (!L_i[i].empty())
    {
      L.resize(L.size() + 1); // conservativeResize
    }
    for (int j = 0; j < L_i[i].size(); j++)
    {
      int a = L_i[i][j];
      int b = L_i[i][(j + 1) % L_i[i].size()];
      L.back().push_back(a);
      if (E.find(std::make_pair(b, a)) != E.end())
      {
        auto tv = E[std::make_pair(b, a)]; // tmp vec
        std::reverse(tv.begin(), tv.end());
        E[std::make_pair(a, b)] = tv;
        L.back().insert(L.back().end(), tv.begin(), tv.end());
        continue;
      }
      int n = std::max((V.row(b) - V.row(a)).norm() / avl, 1.0);
      int vn = V.rows();
      if (b - a != 1 && a - b != V_i.rows() - 1)
      {
        V.conservativeResize(V.rows() + n - 1, 2);
        // split to n edges
        for (int k = 0; k < n - 1; k++)
        {
          V.row(vn + k) << V.row(a) + (V.row(b) - V.row(a)) * (k + 1) / n;
          E[std::make_pair(a, b)].push_back(vn + k);
          L.back().push_back(vn + k);
        }
      }
    }
  }
  // triangulate subpolygons
  Eigen::MatrixXd LP; // local polygon
  Eigen::MatrixXi LF; // local faces
  Eigen::MatrixXd LV;
  for (int i = 0; i < L.size(); i++)
  {
    LP.resize(L[i].size(), 2);
    for (int j = 0; j < LP.rows(); j++)
    {
      LP.row(j) << V.row(L[i][j]);
    }
    //display(LP,0,LP.rows()-1);
    Eigen::MatrixXi LE(LP.rows(), 2);
    LE << Eigen::VectorXi::LinSpaced(LP.rows(), 0, LP.rows() - 1),
        Eigen::VectorXi::LinSpaced(LP.rows(), 1, LP.rows());
    LE(LE.rows() - 1, 1) = 0;
    igl::triangle::triangulate(LP, LE, Eigen::MatrixXd(), "YQq33", LV, LF);
    Eigen::MatrixXd nV = LV.bottomRows(LV.rows() - LP.rows());
    int n_o = V.rows();
    V.conservativeResize(V.rows() + nV.rows(), 2);
    V.bottomRows(nV.rows()) << nV;
    for (int f = 0; f < LF.rows(); f++)
    {
      for (int k = 0; k < 3; k++)
      {
        if (LF(f, k) >= L[i].size())
          LF(f, k) += (n_o - LP.rows());
        else
          LF(f, k) = L[i][LF(f, k)];
      }
    }
    F.conservativeResize(F.rows() + LF.rows(), 3);
    F.bottomRows(LF.rows()) << LF;
  }
}

// the re-implementation of Shor algorithm
bool Shor_van_wyck(
    const Eigen::MatrixXd &P,
    const Eigen::VectorXi &R,
    const std::string flags,
    Eigen::MatrixXd &V,
    Eigen::MatrixXi &F,
    bool do_refine)
{

  // [drop colinear points]
  Eigen::VectorXi B;  // remaining vertices: non-colinear/rotateindex!=0
  Eigen::MatrixXd mP; // P \ colinear vertices
  Eigen::VectorXi mR;
  drop_colinear(P, R, B, mP, mR);

  // [ear clipping]
  Eigen::VectorXi D;
  Eigen::MatrixXi eF;
  Eigen::MatrixXd nP;
  Eigen::VectorXi nR;
  ear_clipping(mP, mR, D, eF, nP);
  igl::slice(mR, D, 1, nR);

  // [weakly-self-overlapping test]
  Eigen::MatrixXi nF;
  bool succ = (nP.rows() == 0) || weakly_self_overlapping(nP, nR, nF);
  if (!succ)
  {
    std::cout << "shor failed" << std::endl;
    exit(0);
    return false;
  }
  // [map simplified index to initial polygon]
  for (int i = 0; i < nF.rows(); i++)
  {
    nF.row(i) << D(nF(i, 0)), D(nF(i, 1)), D(nF(i, 2));
  }
  if (eF.rows() > 0)
  {
    nF.conservativeResize(nF.rows() + eF.rows(), 3);
    nF.block(nF.rows() - eF.rows(), 0, eF.rows(), 3) = eF;
  }

  // [add back colinear vertices by spliting boundary edges]
  add_colinear(P, nF, B, F);
  if (!do_refine)
  {
    V = P;
    return true;
  }
  V = P;

  // [simplify mesh (subdivide into small polygons)]
  std::vector<std::vector<int>> L;
  subdivide_polygon(V, F, L);

  // [refine each small polygon]
  Eigen::MatrixXd V0 = V;
  simplify_triangulation(V0, L, V, F);
  return true;
}

// the re-implementation of Shor algorithm
template <typename Scalar>
bool Shor_van_wyck_v2(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &P,
    const Eigen::VectorXi &R,
    const Eigen::VectorXi &mark,
    const std::string flags,
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &V,
    Eigen::MatrixXi &F,
    Eigen::MatrixXi &Fn,
    bool do_refine)
{

  // [drop colinear points]
  Eigen::VectorXi B;                                        // remaining vertices: non-colinear/rotateindex!=0
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> mP; // P \ colinear vertices
  Eigen::VectorXi mR;
  drop_colinear_v2(P, R, mark, B, mP, mR);
  std::cout << "after removing colinear #P: " << mP.rows() << std::endl;

  // [ear clipping]
  Eigen::VectorXi D;
  Eigen::MatrixXi eF;
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> nP;
  Eigen::VectorXi nR;
  ear_clipping(mP, mR, D, eF, nP);
  igl::slice(mR, D, 1, nR);
  std::cout << "after ear clipping #P: " << nP.rows() << std::endl;

  // [weakly-self-overlapping test]
  Eigen::MatrixXi nF;
  bool succ = (nP.rows() == 0) || weakly_self_overlapping(nP, nR, nF);
  if (!succ)
  {
    std::cout << "shor failed" << std::endl;
    return false;
  }
  else
  {
    // verification
    Eigen::VectorXi R_verify;
    set_rotation_index(nP, nF, R_verify);
    bool verify_succ = true;
    for (int i = 0; i < R_verify.rows(); i++)
    {
      if (R_verify(i) != nR(i))
      {
        // std::cout<<i<<":"<<R_verify(i)<<" --- "<< nR(i)<<std::endl;
        verify_succ = false;
      }
    }
    if (!verify_succ)
    {
      std::cout << "shor failed in verification" << std::endl;
      return false;
    }
  }
  // [map simplified index to initial polygon]
  // for(int i=0;i<nF.rows();i++){
  //   nF.row(i) << D(nF(i,0)),D(nF(i,1)),D(nF(i,2));
  // }

  if (eF.rows() > 0)
  {
    nF.conservativeResize(nF.rows() + eF.rows(), 3);
    nF.block(nF.rows() - eF.rows(), 0, eF.rows(), 3) = eF;
  }

  // [add back colinear vertices by spliting boundary edges]
  std::cout << "add back vertices\n";
  add_colinear(P, nF, B, F);

  V = P;
  Fn = F;
  Eigen::VectorXi bd;
  igl::boundary_loop(Fn, bd);
  std::cout << "#P size: " << V.rows() << std::endl;
  std::cout << "#bd: " << bd.rows() << std::endl;

  return true;
}

template bool add_to_table<mpfr::mpreal>(
    const Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> &,
    const Eigen::VectorXi &,
    std::vector<std::vector<Angle<mpfr::mpreal>>> &,
    std::vector<std::vector<Angle<mpfr::mpreal>>> &,
    int i, int k, int j);
template bool add_to_table<double>(
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &,
    const Eigen::VectorXi &,
    std::vector<std::vector<Angle<double>>> &,
    std::vector<std::vector<Angle<double>>> &,
    int i, int k, int j);

template Eigen::VectorXi set_all_ri<mpfr::mpreal>(
    const Eigen::Matrix<mpfr::mpreal, -1, -1> &,
    const Eigen::MatrixXi &F);
template Eigen::VectorXi set_all_ri<double>(
    const Eigen::Matrix<double, -1, -1> &,
    const Eigen::MatrixXi &F);
template short orientation<mpfr::mpreal>(const Eigen::Matrix<mpfr::mpreal, 3, 2> &P);
template short orientation<double>(const Eigen::Matrix<double, 3, 2> &P);
template void drop_colinear_v2<mpfr::mpreal>(const Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> &, const Eigen::Matrix<int, -1, 1> &, const Eigen::Matrix<int, -1, 1> &, Eigen::Matrix<int, -1, 1> &, Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> &, Eigen::Matrix<int, -1, 1> &);
template void drop_colinear_v2<double>(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &, const Eigen::Matrix<int, -1, 1> &, const Eigen::Matrix<int, -1, 1> &, Eigen::Matrix<int, -1, 1> &, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &, Eigen::Matrix<int, -1, 1> &);
template void add_colinear<double>(const Eigen::Matrix<double, -1, -1> &, const Eigen::Matrix<int, -1, -1> &nF, const Eigen::Matrix<int, -1, 1> &B, Eigen::Matrix<int, -1, -1> &F);
template void set_rotation_index<double>(const Eigen::Matrix<double, -1, -1> &, const Eigen::Matrix<int, -1, -1> &, Eigen::Matrix<int, -1, 1> &, int);
template void set_rotation_index<mpfr::mpreal>(const Eigen::Matrix<mpfr::mpreal, -1, -1> &, const Eigen::Matrix<int, -1, -1> &, Eigen::Matrix<int, -1, 1> &, int);
template bool Shor_van_wyck_v2<double>(const Eigen::Matrix<double, -1, -1> &, const Eigen::Matrix<int, -1, 1> &, const Eigen::Matrix<int, -1, 1> &, const std::string, Eigen::Matrix<double, -1, -1> &, Eigen::Matrix<int, -1, -1> &, Eigen::Matrix<int, -1, -1> &, bool);
template bool Shor_van_wyck_v2<mpfr::mpreal>(const Eigen::Matrix<mpfr::mpreal, -1, -1> &, const Eigen::Matrix<int, -1, 1> &, const Eigen::Matrix<int, -1, 1> &, const std::string, Eigen::Matrix<mpfr::mpreal, -1, -1> &, Eigen::Matrix<int, -1, -1> &, Eigen::Matrix<int, -1, -1> &, bool);