// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Michael Rabinovich
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "line_search.h"
#include <igl/copyleft/cgal/orient2D.h>

IGL_INLINE double igl::line_search(
  Eigen::MatrixXd& x,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& d,
  double step_size,
  std::function<double(Eigen::MatrixXd&)> energy,
  double cur_energy)
{
  double old_energy;
  if (cur_energy > 0)
  {
    old_energy = cur_energy;
  }
  else
  {
    old_energy = energy(x); // no energy was given -> need to compute the current energy
  }
  double new_energy = old_energy;
  int cur_iter = 0; int MAX_STEP_SIZE_ITER = 12;

  while (new_energy >= old_energy && cur_iter < MAX_STEP_SIZE_ITER)
  {
    Eigen::MatrixXd new_x = x + step_size * d;
    // for 2D case
    bool flipped = false;
    if(x.cols()==2) {
      for (int i = 0; i < F.rows(); i++) {
        double a[2] = {x(F(i,0), 0), x(F(i,0), 1)};
        double b[2] = {x(F(i,1), 0), x(F(i,1), 1)};
        double c[2] = {x(F(i,2), 0), x(F(i,2), 1)};
        if(igl::copyleft::cgal::orient2D(a, b, c)<=0)
            flipped = true;
      }
    }
    double cur_e = energy(new_x);
    if (std::isnan(cur_e) || std::isinf(cur_e) || cur_e >= old_energy || flipped)
    {
      step_size /= 2;
    }
    else
    {
      x = new_x;
      new_energy = cur_e;
    }
    cur_iter++;
  }
  return new_energy;
}


#ifdef IGL_STATIC_LIBRARY
#endif
