#include "Property.h"

Property::Property(){
  pj.setZero();
  kappa.setZero();
  TP.setZero();
  angle = 0;
}

Property::Property(Eigen::Vector3i _pj, Eigen::Vector3d _kappa, Eigen::MatrixXd _TP, double r){
  pj = _pj;
  kappa = _kappa;
  TP = _TP;
  angle = r;
}
