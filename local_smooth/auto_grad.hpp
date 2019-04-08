#pragma once

#include <Eigen/Dense>

namespace autogen{
template<typename Scalar>
Scalar sd_energy(const Eigen::Matrix<Scalar,3,2>& V, const Eigen::Matrix<Scalar,2,3>& G)  {
Scalar result_0;
const auto helper_0 = G(0,0)*V(0,0) + G(0,1)*V(1,0) + G(0,2)*V(2,0);
const auto helper_1 = pow(helper_0, 2);
const auto helper_2 = G(0,0)*V(0,1) + G(0,1)*V(1,1) + G(0,2)*V(2,1);
const auto helper_3 = pow(helper_2, 2);
const auto helper_4 = G(1,0)*V(0,0) + G(1,1)*V(1,0) + G(1,2)*V(2,0);
const auto helper_5 = pow(helper_4, 2);
const auto helper_6 = G(1,0)*V(0,1) + G(1,1)*V(1,1) + G(1,2)*V(2,1);
const auto helper_7 = pow(helper_6, 2);
const auto helper_8 = pow(helper_0*helper_6 - helper_2*helper_4, -2);
result_0 = helper_1*helper_8 + helper_1 + helper_3*helper_8 + helper_3 + helper_5*helper_8 + helper_5 + helper_7*helper_8 + helper_7;
return result_0;
}
template<typename Scalar>
void sd_grad(const Eigen::Matrix<Scalar,3,2>& V, const Eigen::Matrix<Scalar,2,3>& G, Eigen::Matrix<Scalar,1,2>& res)  {
const auto helper_0 = G(0,0)*V(0,0) + G(0,1)*V(1,0) + G(0,2)*V(2,0);
const auto helper_1 = G(0,0)*helper_0;
const auto helper_2 = G(1,0)*V(0,0) + G(1,1)*V(1,0) + G(1,2)*V(2,0);
const auto helper_3 = G(1,0)*helper_2;
const auto helper_4 = G(1,0)*V(0,1) + G(1,1)*V(1,1) + G(1,2)*V(2,1);
const auto helper_5 = G(0,0)*V(0,1) + G(0,1)*V(1,1) + G(0,2)*V(2,1);
const auto helper_6 = helper_0*helper_4 - helper_2*helper_5;
const auto helper_7 = pow(helper_6, -2);
const auto helper_8 = pow(helper_0, 2);
const auto helper_9 = pow(helper_6, -3);
const auto helper_10 = helper_9*(G(0,0)*helper_4 - G(1,0)*helper_5);
const auto helper_11 = pow(helper_5, 2);
const auto helper_12 = pow(helper_2, 2);
const auto helper_13 = pow(helper_4, 2);
const auto helper_14 = G(0,0)*helper_5;
const auto helper_15 = G(1,0)*helper_4;
const auto helper_16 = helper_9*(G(0,0)*helper_2 - G(1,0)*helper_0);
res(0) = 2*helper_1*helper_7 + 2*helper_1 - 2*helper_10*helper_11 - 2*helper_10*helper_12 - 2*helper_10*helper_13 - 2*helper_10*helper_8 + 2*helper_3*helper_7 + 2*helper_3;
res(1) = 2*helper_11*helper_16 + 2*helper_12*helper_16 + 2*helper_13*helper_16 + 2*helper_14*helper_7 + 2*helper_14 + 2*helper_15*helper_7 + 2*helper_15 + 2*helper_16*helper_8;
}
template<typename Scalar>
void sd_hess(const Eigen::Matrix<Scalar,3,2>& V, const Eigen::Matrix<Scalar,2,3>& G, Eigen::Matrix<Scalar,2,2>& res)  {
const auto helper_0 = pow(G(0,0), 2);
const auto helper_1 = pow(G(1,0), 2);
const auto helper_2 = G(0,0)*V(0,0) + G(0,1)*V(1,0) + G(0,2)*V(2,0);
const auto helper_3 = G(1,0)*V(0,1) + G(1,1)*V(1,1) + G(1,2)*V(2,1);
const auto helper_4 = G(0,0)*V(0,1) + G(0,1)*V(1,1) + G(0,2)*V(2,1);
const auto helper_5 = G(1,0)*V(0,0) + G(1,1)*V(1,0) + G(1,2)*V(2,0);
const auto helper_6 = helper_2*helper_3 - helper_4*helper_5;
const auto helper_7 = pow(helper_6, -2);
const auto helper_8 = helper_0*helper_7 + helper_0 + helper_1*helper_7 + helper_1;
const auto helper_9 = G(0,0)*helper_3 - G(1,0)*helper_4;
const auto helper_10 = pow(helper_6, -3);
const auto helper_11 = 4*G(0,0)*helper_10;
const auto helper_12 = 4*G(1,0)*helper_10;
const auto helper_13 = pow(helper_2, 2);
const auto helper_14 = pow(helper_6, -4);
const auto helper_15 = 3*helper_14*pow(helper_9, 2);
const auto helper_16 = pow(helper_4, 2);
const auto helper_17 = pow(helper_5, 2);
const auto helper_18 = pow(helper_3, 2);
const auto helper_19 = 2*G(0,0);
const auto helper_20 = G(0,0)*helper_5 - G(1,0)*helper_2;
const auto helper_21 = 2*G(1,0);
const auto helper_22 = 3*helper_20*helper_9/helper_6;
const auto helper_23 = 2*helper_10*(-helper_13*helper_22 - helper_16*helper_22 - helper_17*helper_22 - helper_18*helper_22 + helper_19*helper_2*helper_20 - helper_19*helper_4*helper_9 + helper_20*helper_21*helper_5 - helper_21*helper_3*helper_9);
const auto helper_24 = 3*helper_14*pow(helper_20, 2);
res(0) = -2*helper_11*helper_2*helper_9 - 2*helper_12*helper_5*helper_9 + 2*helper_13*helper_15 + 2*helper_15*helper_16 + 2*helper_15*helper_17 + 2*helper_15*helper_18 + 2*helper_8;
res(1) = helper_23;
res(2) = helper_23;
res(3) = 2*helper_11*helper_20*helper_4 + 2*helper_12*helper_20*helper_3 + 2*helper_13*helper_24 + 2*helper_16*helper_24 + 2*helper_17*helper_24 + 2*helper_18*helper_24 + 2*helper_8;
}

}
