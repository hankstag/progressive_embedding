// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "edge_flaps.h"
#include <igl/unique_edge_map.h>
#include <vector>
#include <cassert>
#include <iostream>
IGL_INLINE void igl::edge_flaps(
        const Eigen::MatrixXi & F,
        const Eigen::MatrixXi & E,
        const Eigen::MatrixXi & AllE,
        const Eigen::VectorXi & EMAP,
        Eigen::MatrixXi & EF,
        Eigen::MatrixXi & EI,
        Eigen::MatrixXi & dEF,
        Eigen::MatrixXi & dEI
){
    // Initialize to boundary value
    EF.setConstant(E.rows(),2,-1);
    EI.setConstant(E.rows(),2,-1);
    // loop over all faces
    for(int f = 0;f<F.rows();f++)
    {
        // loop over edges across from corners
        for(int v = 0;v<3;v++)
        {
            // get edge id
            const int e = EMAP(v*F.rows()+f);
            // See if this is left or right flap w.r.t. edge orientation
            if( F(f,(v+1)%3) == E(e,0) && F(f,(v+2)%3) == E(e,1))
            {
                EF(e,0) = f;
                EI(e,0) = v;
            }else
            {
                assert(F(f,(v+1)%3) == E(e,1) && F(f,(v+2)%3) == E(e,0));
                EF(e,1) = f;
                EI(e,1) = v;
            }
        }
    }
    // ===  keep the same === //
    dEF.setConstant(AllE.rows(),2,-1);
    dEI.setConstant(AllE.rows(),2,-1);
    for(int f = 0;f<F.rows();f++){
        for(int v = 0;v < 3;v++){
            // directed edge
            const int e = v*F.rows()+f;
            const int ue = EMAP(v*F.rows()+f);
            if(F(f,(v+1)%3) == E(ue,0) && F(f,(v+2)%3) == E(ue,1)){
                dEF(e,0) = EF(ue,0);
                dEF(e,1) = EF(ue,1);
                dEI(e,0) = EI(ue,0);
                dEI(e,1) = EI(ue,1);
            }else{
                dEF(e,0) = EF(ue,1);
                dEF(e,1) = EF(ue,0);
                dEI(e,0) = EI(ue,1);
                dEI(e,1) = EI(ue,0);
            }
        }
    }
}

IGL_INLINE void igl::edge_flaps(
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXi & E,
  const Eigen::VectorXi & EMAP,
  Eigen::MatrixXi & EF,
  Eigen::MatrixXi & EI)
{
  // Initialize to boundary value
  EF.setConstant(E.rows(),2,-1);
  EI.setConstant(E.rows(),2,-1);
  // loop over all faces
  for(int f = 0;f<F.rows();f++)
  {
    // loop over edges across from corners
    for(int v = 0;v<3;v++)
    {
      // get edge id
      const int e = EMAP(v*F.rows()+f);
      // See if this is left or right flap w.r.t. edge orientation
      if( F(f,(v+1)%3) == E(e,0) && F(f,(v+2)%3) == E(e,1))
      {
        EF(e,0) = f;
        EI(e,0) = v;
      }else
      {
        assert(F(f,(v+1)%3) == E(e,1) && F(f,(v+2)%3) == E(e,0));
        EF(e,1) = f;
        EI(e,1) = v;
      }
    }
  }
}

IGL_INLINE void igl::edge_flaps(
  const Eigen::MatrixXi & F,
  Eigen::MatrixXi & E,
  Eigen::VectorXi & EMAP,
  Eigen::MatrixXi & EF,
  Eigen::MatrixXi & EI)
{
  Eigen::MatrixXi allE;
  std::vector<std::vector<int> > uE2E;
  igl::unique_edge_map(F,allE,E,EMAP,uE2E);
  // Const-ify to call overload
  const auto & cE = E;
  const auto & cEMAP = EMAP;
  return edge_flaps(F,cE, cEMAP,EF,EI);
}

IGL_INLINE void igl::edge_flaps(
        const Eigen::MatrixXi & F,
        Eigen::MatrixXi & E,
        Eigen::MatrixXi & allE,
        Eigen::VectorXi & EMAP,
        Eigen::MatrixXi & EF,
        Eigen::MatrixXi & EI,
        Eigen::MatrixXi & dEF,
        Eigen::MatrixXi & dEI,
        Eigen::VectorXi & EE)
{
    std::vector<std::vector<int> > uE2E;
    igl::unique_edge_map(F,allE,E,EMAP,uE2E);
    // Const-ify to call overload
    const auto & cE = E;
    const auto & cEMAP = EMAP;
    EE.setConstant(allE.rows(),-1);
    for(int i=0;i<uE2E.size();i++){
        for(int j=0;j<uE2E[i].size();j++){
            // if directed edge is boundary, it maps to itself
            EE(uE2E[i][j]) = uE2E[i][(j+1)%uE2E[i].size()];
        }
    }
    return edge_flaps(F,cE, allE, cEMAP,EF,EI, dEF, dEI);
}