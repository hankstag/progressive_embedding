#ifndef PROPERTY_H
#define PROPERTY_H

#include <Eigen/Core>
#include <vector>

// collection of properties defined per face
class Property{

public:
  
  Property();
  Property(Eigen::Vector3i, Eigen::Vector3d, Eigen::MatrixXd, double);
  ~Property(){};

  // period jumps
  Eigen::Vector3i pj;
  
  // diff between neighbor frames
  Eigen::Vector3d kappa;
  
  // reference frame (2 principle directions)
  Eigen::MatrixXd TP;

  // angle of the field
  double angle;
  
};

#endif
