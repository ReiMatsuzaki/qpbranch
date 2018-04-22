#include <iostream>
#include <gtest/gtest.h>
#include "pwgto.hpp"

using namespace std;
using namespace qpbranch;

TEST(utest_pwgto, overlap) {
  int num = 4;
  VectorXi ns(num); ns << 2, 3, 4, 2;
  vector<Operator> ops; ops.push_back(kOp0);
  PlaneWaveGTO *basis = new PlaneWaveGTO(ns, ops);
  basis->gs_ << 1.0,  1.0, 1.0, 1.0;
  basis->Rs_ << 10.0, 0.0, 0.1, 0.1;
  basis->Ps_ << 1.0,  0.0, 0.1, 0.1;
  basis->ops_[0] = kOp0;
  basis->setup();
  MatrixXcd S(num, num);
  basis->overlap(kOp0, kOp0, &S);

  for(int A = 1; A<num; A++) {
    ASSERT_DOUBLE_EQ(1.0, S(A,A).real());
    ASSERT_DOUBLE_EQ(0.0, S(A,A).imag());
  }
  
  delete basis;
}
TEST(utest_pwgto, multipole0) {
  int num = 4;
  vector<Operator> ops; ops.push_back(kOp0);
  VectorXi ns(num); ns << 0, 0, 2, 2;
  PlaneWaveGTO *pw_basis = new PlaneWaveGTO(ns, ops);
  pw_basis->gs_ << 1.1, 1.1, 1.1, 1.1;
  pw_basis->Rs_ << 0.0, 0.1, 0.0, 0.0;
  pw_basis->Ps_ = VectorXd::Zero(num);
  
  GaussBasis *basis = pw_basis;
  basis->setup();

  MatrixXcd S(num, num);
  basis->overlap(kOp0, kOp0, &S);

}

