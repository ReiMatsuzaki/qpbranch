#include <iostream>

#include <pwgto.hpp>
#include <operator.hpp>
#include <eigenplus.hpp>

using namespace std;
using namespace qpbranch;

int main () {

  // system parameter
  double m = 2000.0;        // mass
  complex<double> b(1.3);   // potential width
  complex<double> v0(1.2);  // potential height
  complex<double> q0(0.0);  // potential position

  // operators
  auto op_id = new OperatorId();
  auto op_p1 = new OperatorPn(1);
  auto op_p2 = new OperatorPn(2);
  auto op_v  = new OperatorGausspot(v0, b, q0);
  
  // basis function
  int num = 2;
  VectorXi ns(num); ns << 0, 2;
  vector<Operator*> ops{op_id, op_p1, op_p2, op_v};
  PlaneWaveGto *basis = new PlaneWaveGto(ns, ops);
  basis->gs_ << VectorXcd::Constant(num, 1.0);
  basis->Rs_ << VectorXd::Constant(num, 0.0);
  basis->Ps_ << VectorXd::Constant(num, 0.0);
  basis->setup();

  // run on grid
  /*int nR = 200;*/
int nR = 2;
  VectorXd Rs = VectorXd::LinSpaced(nR, -10.0, 10.0);
  for(int iR = 0; iR < nR; iR++) {
    cout << iR << "/" << nR << endl;
    
    basis->Rs_ = VectorXd::Constant(num, Rs[iR]);
    basis->setup();

    MatrixXcd S(num, num), H(num, num), P2(num,num), P1(num, num);
    basis->matrix(op_id, op_id, &S);
    basis->matrix(op_id, op_p1, &P1);
    basis->matrix(op_id, op_p2, &P2);
    basis->matrix(op_id, op_v,  &H);
    H += P2/(2*m);

    MatrixXcd U(num,num);
    VectorXd w(num);
    zhegv(H, S, &U, &w);
    cout << "w:" << w << endl;
  }
  
  delete basis;
}


