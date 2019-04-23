
#include "process_ee.h"
#include <igl/PI.h>
#include <iostream>

bool read_cone_ee(
  std::string &fname, 
  Eigen::MatrixXi& RI,
  Eigen::MatrixXi& EE
){
  // Open file, and check for error
  FILE * obj_file = fopen(fname.c_str(),"r");
  if(NULL==obj_file){
    fprintf(stderr,"IOError: %s could not be opened...\n",
            fname.c_str());
    return false;
  }
  #define LINE_MAX 2048
  char line[LINE_MAX];
  int line_no = 1;
  EE.setZero();
  RI.setZero();
  while (fgets(line, LINE_MAX, obj_file) != NULL){
    char type[LINE_MAX];
    // Read first word containing type
    if(sscanf(line, "%s",type) == 1){
      char * l = &line[strlen(type)];
      if(strlen(type) >= 1){
        if (type[0] == 'c'){
          int x[2];
          int count = sscanf(l,"%d %d\n",&x[0],&x[1]);
          RI.conservativeResize(RI.rows()+1,2);
          RI.bottomRows(1)<<x[0]-1,x[1];
        }else if(type[0] == 'e'){
          int x[6];
          int count = sscanf(l,"%d %d %d %d %d %d\n",&x[0],&x[1],&x[2],&x[3],&x[4],&x[5]);
          EE.conservativeResize(EE.rows()+1,6);
          EE.bottomRows(1)<<x[0]-1,x[1]-1,x[2]-1,x[3]-1,x[4],x[5];
        }
      }
    }else{
      // ignore empty line
    }   
    line_no++;
  }
  fclose(obj_file);
  return true;
}

void process_ee(
  const Eigen::MatrixXi& ee,
  const Eigen::MatrixXi& Fuv,
  Eigen::MatrixXi& cut
){
  int c = 0;
  cut.resize(ee.rows(),4);
  for(int i=0;i<ee.rows();i++){
    if(ee(i,0) > ee(i,2) || ee(i,5)==1) continue;
    int A2,B2,C2,D2;
    for(int j=0;j<3;j++){
        if(Fuv(ee(i,0),j) == ee(i,1)){
            A2 = ee(i,1);
            B2 = Fuv(ee(i,0),(j+1)%3);
        }
        if(Fuv(ee(i,2),j) == ee(i,3)){
            C2 = ee(i,3);
            D2 = Fuv(ee(i,2),(j+1)%3);
        }
    }
    cut.row(c++)<<A2,B2,C2,D2;
  }
  cut.conservativeResize(c,4);
}

void buildAeq(
  const Eigen::MatrixXi& cut,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& Fuv,
  const Eigen::MatrixXd& uv,
  Eigen::SparseMatrix<double>& Aeq
){
  Eigen::VectorXd tail;
  int N = V.rows();
  int c = 0;
  int m = cut.rows();
  Aeq.resize(2*m,V.rows()*2);
  int A,B,C,D,A2,B2,C2,D2;
   for(int i=0;i<cut.rows();i++){
     int A2=cut(i,0);
    int B2=cut(i,1);
    int C2=cut(i,2);
    int D2=cut(i,3);
      

    std::complex<double> l0,l1,r0,r1;
    l0 = std::complex<double>(uv(A2,0),uv(A2,1));
    l1 = std::complex<double>(uv(B2,0),uv(B2,1));
    r0 = std::complex<double>(uv(C2,0),uv(C2,1));
    r1 = std::complex<double>(uv(D2,0),uv(D2,1));
  

    int r = std::round(2.0 * std::log((l0 - l1) / (r0 - r1)).imag() / igl::PI);
    r = ((r % 4) + 4) % 4; // ensure that r is between 0 and 3
     switch(r){
      case 0:
        Aeq.coeffRef(c,A2) += 1;
        Aeq.coeffRef(c,B2) += -1;
        Aeq.coeffRef(c,C2) += -1;
        Aeq.coeffRef(c,D2) += 1;
        Aeq.coeffRef(c+1,A2+N) += 1;
        Aeq.coeffRef(c+1,B2+N) += -1;
        Aeq.coeffRef(c+1,C2+N) += -1;
        Aeq.coeffRef(c+1,D2+N) += 1;
        c = c+2;
        break;
      case 1:
        Aeq.coeffRef(c,A2) += 1;
        Aeq.coeffRef(c,B2) += -1;
        Aeq.coeffRef(c,C2+N) += 1;
        Aeq.coeffRef(c,D2+N) += -1;
        Aeq.coeffRef(c+1,C2) += 1;
        Aeq.coeffRef(c+1,D2) += -1;
        Aeq.coeffRef(c+1,A2+N) += -1;
        Aeq.coeffRef(c+1,B2+N) += 1;   
        c = c+2;               
        break;
      case 2:
        Aeq.coeffRef(c,A2) += 1;
        Aeq.coeffRef(c,B2) += -1;
        Aeq.coeffRef(c,C2) += 1;
        Aeq.coeffRef(c,D2) += -1;
        Aeq.coeffRef(c+1,A2+N) += 1;
        Aeq.coeffRef(c+1,B2+N) += -1;
        Aeq.coeffRef(c+1,C2+N) += 1;
        Aeq.coeffRef(c+1,D2+N) += -1;
        c = c+2;
        break;
      case 3:
        Aeq.coeffRef(c,A2) += 1;
        Aeq.coeffRef(c,B2) += -1;
        Aeq.coeffRef(c,C2+N) += -1;
        Aeq.coeffRef(c,D2+N) += 1;
        Aeq.coeffRef(c+1,C2) += 1;
        Aeq.coeffRef(c+1,D2) += -1;
        Aeq.coeffRef(c+1,A2+N) += 1;
        Aeq.coeffRef(c+1,B2+N) += -1;
        c = c+2;
        break;
    }
  }
    
  Aeq.conservativeResize(c,V.rows()*2);
  Aeq.makeCompressed();
  std::cout<<"Aeq size "<<Aeq.rows()<<","<<Aeq.cols()<<std::endl;
  // test initial violation 
  Eigen::VectorXd UV(uv.rows()*2);
  UV<<uv.col(0),uv.col(1);
  Eigen::SparseMatrix<double> t = UV.sparseView();
  t.makeCompressed();
  Eigen::SparseMatrix<double> mm = Aeq*t;
  Eigen::VectorXd z = Eigen::VectorXd(mm);
  if(z.rows()>0)
  std::cout<<"max violation "<<z.cwiseAbs().maxCoeff()<<std::endl;
}
