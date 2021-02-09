#ifndef PROP_H
#define PROP_H

class fData{

public:
  fData(Eigen::Vector3i _pj, Eigen::Vector3d _k, Eigen::MatrixXd _TP, double r): pj(_pj), kappa(_k), TP(_TP), angle(r){};
  ~fData(){};

  // period jumps
  Eigen::Vector3i pj;
  
  // diff between neighbor frames
  Eigen::Vector3d kappa;
  
  // reference frame (2 principle directions)
  Eigen::MatrixXd TP;
  
  double angle;

};

#endif