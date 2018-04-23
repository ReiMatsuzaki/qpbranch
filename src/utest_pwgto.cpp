#include <iostream>
#include <gtest/gtest.h>
#include "mathplus.hpp"
#include "pwgto.hpp"

using namespace std;
using namespace qpbranch;

TEST(utest_pwgto, overlap) {
  int num = 4;
  VectorXi ns(num); ns << 0, 0, 2, 2;
  vector<Operator> ops; ops.push_back(kOp0);
  PlaneWaveGTO *basis = new PlaneWaveGTO(ns, ops);
  basis->gs_ << 1.0,  complex<double>(0.9,-0.8), 1.0, complex<double>(0.2, 0.1);
  basis->Rs_ << 0.0, 0.0, 0.2, 0.2;
  basis->Ps_ << 1.0, 0.0, 0.5, 0.5;
  
  basis->setup();
  MatrixXcd S(num, num);
  basis->overlap(kOp0, kOp0, &S);

  for(int A = 0; A<num; A++) {
    ASSERT_DOUBLE_EQ(1.0, S(A,A).real()) << "A:" << A << endl;
    ASSERT_DOUBLE_EQ(0.0, S(A,A).imag()) << "A:" << A << endl;
  }

  for(int A = 0; A<num; A++) {
    for(int B = 0; B<A; B++) {
      ASSERT_DOUBLE_EQ(real(S(A,B)), +real(S(B,A))) << "A:" << A << endl << "B:" << B << endl;
      ASSERT_DOUBLE_EQ(imag(S(A,B)), -imag(S(B,A))) << "A:" << A << endl << "B:" << B << endl;
    }
  }

  int A = 2;
  int B = 3;
  auto calc = S(A, B)/(basis->Ns_(A)*basis->Ns_(B));
  VectorXcd gg(9+1);
  gtoint2n(9, conj(basis->gs_(A))+basis->gs_(B), &gg);
  auto ref = gg((basis->ns_(A) + basis->ns_(B))/2);
  ASSERT_DOUBLE_EQ(real(ref), real(calc)) << "A:" << A << endl << "B:" << B << endl;
  ASSERT_DOUBLE_EQ(imag(ref), imag(calc)) << "A:" << A << endl << "B:" << B << endl;
  
}
TEST(utest_pwgto, multipole0) {
  int num = 5;
  vector<Operator> ops = {kOp0, kOp1, kOp2};
  VectorXi ns(num); ns << 0, 0, 2, 1, 0;
  PlaneWaveGTO *pw_basis = new PlaneWaveGTO(ns, ops);
  pw_basis->gs_ << 1.1, 1.1, 1.1, 0.4, 0.8;
  pw_basis->Rs_ << 0.0, 0.1, 0.1, 0.1, 0.3;
  pw_basis->Ps_ << 0.0, 0.0, 0.3, 0.0, 0.1;
  
  GaussBasis *basis = pw_basis;
  basis->setup();

  MatrixXcd S(num, num);
  basis->overlap(kOp0, kOp0, &S);

  MatrixXcd M01, M10;
  basis->overlap(kOp0, kOp1, &M01);
  basis->overlap(kOp1, kOp0, &M10);

  for(int A = 0; A < num; A++) {
    for(int B = 0; B < num; B++) {
      ASSERT_DOUBLE_EQ(real(M01(A,B)), real(M10(A,B)));
      ASSERT_DOUBLE_EQ(imag(M01(A,B)), imag(M10(A,B)));
    }
  }

  /*
  int A = 1;
  int B = 1;
  VectorXcd gg(3);
  gtoint2n(2, conj(basis->gs_(A))+basis_->gs(B), &gg);
  auto ref  = gg((basis->ns_[A]+basis->ns_[B])/2) * basis->Ns_[A] * basis_->Ns[B];
  auto calc = 
  */
  delete basis;
}

