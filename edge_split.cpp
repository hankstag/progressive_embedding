#include "edge_split.h"
#include <iostream>

bool edge_split_m2(
  std::vector<std::vector<double>>& V0,
  std::vector<std::vector<int>>& F0,
  std::vector<std::vector<int>>& FF_vec,
  std::vector<std::vector<int>>& FFi_vec,
  std::vector<Property>& props,
  int f0,
  int e0
){
  int f0_curr_layer = f0;
  int f0_next_layer = f0%2==0? (f0/2)*2+1 : (f0/2)*2;
  
  int f1_curr_layer = FF_vec[f0_curr_layer][e0];
  int f1_next_layer = FF_vec[f0_next_layer][e0];
  
  // corner id should be independent with the layer
  int e1 = int(FFi_vec[f0_curr_layer][e0]);
  int e01 = (e0 + 1) % 3;
  int e02 = (e0 + 2) % 3;
  int e11 = (e1 + 1) % 3;
  int e12 = (e1 + 2) % 3;
  
  int f01_curr_layer = int(FF_vec[f0_curr_layer][e01]);
  int f02_curr_layer = int(FF_vec[f0_curr_layer][e02]);
  int f11_curr_layer = int(FF_vec[f1_curr_layer][e11]);
  int f12_curr_layer = int(FF_vec[f1_curr_layer][e12]);
  
  int f01_next_layer = int(FF_vec[f0_next_layer][e01]);
  int f02_next_layer = int(FF_vec[f0_next_layer][e02]);
  int f11_next_layer = int(FF_vec[f1_next_layer][e11]);
  int f12_next_layer = int(FF_vec[f1_next_layer][e12]);
  
  int u1_curr_layer = F0[f0_curr_layer][e01];
  int u0_curr_layer = F0[f1_curr_layer][e11];
  
  // pushback middle currrent layer
  V0.push_back({(0.6*V0[u0_curr_layer][0]+0.4*V0[u1_curr_layer][0]),
                (0.6*V0[u0_curr_layer][1]+0.4*V0[u1_curr_layer][1]),
                (0.6*V0[u0_curr_layer][2]+0.4*V0[u1_curr_layer][2])
              });
  
  int u1_next_layer = F0[f0_next_layer][e01];
  int u0_next_layer = F0[f1_next_layer][e11];
  
  // pushback middle next layer
  V0.push_back({(0.6*V0[u0_next_layer][0]+0.4*V0[u1_next_layer][0]),
                (0.6*V0[u0_next_layer][1]+0.4*V0[u1_next_layer][1]),
                (0.6*V0[u0_next_layer][2]+0.4*V0[u1_next_layer][2])
              });
  
  F0.push_back(F0[f0_curr_layer]);
  F0.push_back(F0[f0_next_layer]);
  F0.push_back(F0[f1_curr_layer]);
  F0.push_back(F0[f1_next_layer]);
  
  int ux_curr_layer = V0.size()-2;
  int ux_next_layer = V0.size()-1;
  
  F0[f0_curr_layer][e0] = ux_curr_layer;
  F0[f0_next_layer][e0] = ux_next_layer;
  F0[f1_curr_layer][e1] = ux_curr_layer;
  F0[f1_next_layer][e1] = ux_next_layer;
  
  int fx1_next_layer = F0.size()-1;
  int fx1_curr_layer = F0.size()-2;
  int fx0_next_layer = F0.size()-3;
  int fx0_curr_layer = F0.size()-4;
  
  F0[fx1_next_layer][e11] = ux_next_layer;
  F0[fx0_next_layer][e01] = ux_next_layer;
  F0[fx1_curr_layer][e11] = ux_curr_layer;
  F0[fx0_curr_layer][e01] = ux_curr_layer;
  
  if(f12_curr_layer != -1){
    FF_vec[f12_curr_layer][FFi_vec[f1_curr_layer][e12]] = fx1_curr_layer;
  }
  if(f12_next_layer != -1){
    FF_vec[f12_next_layer][FFi_vec[f1_next_layer][e12]] = fx1_next_layer;
  }
  if(f02_curr_layer != -1){
    FF_vec[f02_curr_layer][FFi_vec[f0_curr_layer][e02]] = fx0_curr_layer;
  }
  if(f02_next_layer != -1){
    FF_vec[f02_next_layer][FFi_vec[f0_next_layer][e02]] = fx0_next_layer;
  }
  
  FF_vec.push_back(FF_vec[f0_curr_layer]);
  FF_vec.push_back(FF_vec[f0_next_layer]);
  FF_vec.push_back(FF_vec[f1_curr_layer]);
  FF_vec.push_back(FF_vec[f1_next_layer]);
  
  FF_vec[f0_curr_layer][e02] = fx0_curr_layer;
  FF_vec[f0_next_layer][e02] = fx0_next_layer;
  
  FF_vec[f0_curr_layer][e0] = fx1_curr_layer;
  FF_vec[f0_next_layer][e0] = fx1_next_layer;
  
  FF_vec[f1_curr_layer][e12] = fx1_curr_layer;
  FF_vec[f1_next_layer][e12] = fx1_next_layer;
  
  FF_vec[f1_curr_layer][e1] = fx0_curr_layer;
  FF_vec[f1_next_layer][e1] = fx0_next_layer;
  
  FF_vec[fx1_curr_layer][e11] = f1_curr_layer;
  FF_vec[fx1_next_layer][e11] = f1_next_layer;
  
  FF_vec[fx0_curr_layer][e01] = f0_curr_layer;
  FF_vec[fx0_next_layer][e01] = f0_next_layer;
  
  // update FFi
  FFi_vec.push_back(FFi_vec[f0_curr_layer]);
  FFi_vec.push_back(FFi_vec[f0_next_layer]);
  FFi_vec.push_back(FFi_vec[f1_curr_layer]);
  FFi_vec.push_back(FFi_vec[f1_next_layer]);
  FFi_vec[f0_curr_layer][e02] = e01;
  FFi_vec[f0_next_layer][e02] = e01;
  FFi_vec[f1_curr_layer][e12] = e11;
  FFi_vec[f1_next_layer][e12] = e11;
  FFi_vec[fx0_curr_layer][e01] = e02;
  FFi_vec[fx0_next_layer][e01] = e02;
  FFi_vec[fx1_curr_layer][e11] = e12;
  FFi_vec[fx1_next_layer][e11] = e12;
    
  // update properties
  props.push_back(props[f0_curr_layer]);
  props.push_back(props[f0_next_layer]);
  props.push_back(props[f1_curr_layer]);
  props.push_back(props[f1_next_layer]);
  
  props[f0_curr_layer].pj((e0+2)%3) = 0;
  props[f0_next_layer].pj((e0+2)%3) = 0;
  props[f0_curr_layer].kappa((e0+2)%3) = 0;
  props[f0_next_layer].kappa((e0+2)%3) = 0;
  
  props[f1_curr_layer].pj((e1+2)%3) = 0;
  props[f1_next_layer].pj((e1+2)%3) = 0;
  props[f1_curr_layer].kappa((e1+2)%3) = 0;
  props[f1_next_layer].kappa((e1+2)%3) = 0;
  
  props[fx0_curr_layer].pj((e0+1)%3) = 0;
  props[fx0_next_layer].pj((e0+1)%3) = 0;
  props[fx0_curr_layer].kappa((e0+1)%3) = 0;
  props[fx0_next_layer].kappa((e0+1)%3) = 0;
  
  props[fx1_curr_layer].pj((e1+1)%3) = 0;
  props[fx1_next_layer].pj((e1+1)%3) = 0;
  props[fx1_curr_layer].kappa((e1+1)%3) = 0;
  props[fx1_next_layer].kappa((e1+1)%3) = 0;
  
  return true;
  
}

bool edge_split(
  std::vector<std::vector<double>>& V0,
  std::vector<std::vector<int>>& F0,
  std::vector<std::vector<int>>& FF_vec,
  std::vector<std::vector<int>>& FFi_vec,
  std::vector<Property>& props,
  int f0,
  int e0
){
  int f1 = FF_vec[f0][e0];
  if(f1 == -1) return false;
  int e1 = int(FFi_vec[f0][e0]);
  int e01 = (e0 + 1) % 3;
  int e02 = (e0 + 2) % 3;
  int e11 = (e1 + 1) % 3;
  int e12 = (e1 + 2) % 3;
  int f01 = int(FF_vec[f0][e01]);
  int f02 = int(FF_vec[f0][e02]);
  int f11 = int(FF_vec[f1][e11]);
  int f12 = int(FF_vec[f1][e12]);

  int u1 = F0[f0][e01];
  int u0 = F0[f1][e11];
  // pushback middle
  V0.push_back({(0.6*V0[u0][0]+0.4*V0[u1][0]),
                (0.6*V0[u0][1]+0.4*V0[u1][1]),
                (0.6*V0[u0][2]+0.4*V0[u1][2])
               });
  F0.push_back(F0[f0]);
  F0.push_back(F0[f1]);
  int ux = V0.size()-1;

  F0[f1][e1] = ux;
  F0[f0][e0] = ux;
  int fx1 = F0.size()-1;
  int fx0 = F0.size()-2;
  F0[fx1][e11] = ux;
  F0[fx0][e01] = ux;

  if(f12 != -1){
    FF_vec[f12][FFi_vec[f1][e12]] = fx1;
  }
  if(f02 != -1){
    FF_vec[f02][FFi_vec[f0][e02]] = fx0;
  }

  FF_vec.push_back(FF_vec[f0]);
  FF_vec.push_back(FF_vec[f1]);
  FF_vec[f0][e02] = fx0;
  FF_vec[f0][e0] = fx1;
  FF_vec[f1][e12] = fx1;
  FF_vec[f1][e1] = fx0;
  FF_vec[fx1][e11] = f1;
  FF_vec[fx0][e01] = f0;

  FFi_vec.push_back(FFi_vec[f0]);
  FFi_vec.push_back(FFi_vec[f1]);
  FFi_vec[f0][e02] = e01;
  FFi_vec[f1][e12] = e11;
  FFi_vec[fx0][e01] = e02;
  FFi_vec[fx1][e11] = e12;
  
  // update properties
  props.push_back(props[f0]);
  props.push_back(props[f1]);
  
  props[f0].pj((e0+2)%3) = 0;
  props[f0].kappa((e0+2)%3) = 0;
  
  props[f1].pj((e1+2)%3) = 0;
  props[f1].kappa((e1+2)%3) = 0;
  
  props[fx0].pj((e0+1)%3) = 0;
  props[fx0].kappa((e0+1)%3) = 0;
  
  props[fx1].pj((e1+1)%3) = 0;
  props[fx1].kappa((e1+1)%3) = 0;

  return true;

}
