#include <iostream>
#include <Eigen/Eigenvalues>
#include <gtest/gtest.h>
#include "mathplus.hpp"
#include "pwgto.hpp"

using namespace std;
using namespace qpbranch;
using namespace boost;

TEST(utest_pwgto, test_coef_d) {
  complex<double> gAB(1, 0.2);
  complex<double> wAB(1.2, 0);
  double RA(0.3), RB(0.3);
  const int maxnA = 3;
  const int maxnB = 2;
  multi_array<complex<double>, 3> d(extents[maxnA+1][maxnB+1][maxnA+maxnB+1]);;
  hermite_coef_d(gAB, wAB, RA, RB, maxnA, maxnB, &d);
  
  for(int nA=0; nA<=maxnA; nA++) {
    for(int nB=0; nB<=maxnB; nB++) {
      for(int Nk=0; Nk<=maxnA+maxnB; Nk++) {
	auto ref = hermite_coef_d_0(gAB, wAB, RA, RB, nA, nB, Nk);
	ASSERT_DOUBLE_EQ(real(ref), real(d[nA][nB][Nk])) << "nA:" << nA << endl;
	ASSERT_DOUBLE_EQ(imag(ref), imag(d[nA][nB][Nk])) << "nA:" << nA << endl;	
      }
    }
  }

}
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
  pw_basis->gs_ << 1.1, 1.1, 1.2, 0.4, 0.8;
  pw_basis->Rs_ << 0.0, 0.1, 0.1, 0.1, 0.3;  
  pw_basis->Ps_ << 0.0, 0.0, 0.3, 0.0, 0.1;
  
  GaussBasis *basis = pw_basis;
  basis->setup();

  MatrixXcd S(num, num);
  basis->overlap(kOp0, kOp0, &S);
  
  MatrixXcd M01(num, num);
  basis->overlap(kOp0, kOp1, &M01);

  MatrixXcd M10(num, num);
  basis->overlap(kOp1, kOp0, &M10);

  MatrixXcd M02(num, num);
  basis->overlap(kOp0, kOp2, &M02);

  MatrixXcd M11(num, num);
  basis->overlap(kOp1, kOp1, &M11);

  MatrixXcd M20(num, num);
  basis->overlap(kOp2, kOp0, &M20);
  
  for(int A = 0; A < num; A++) {
    auto g = conj(basis->gs_(A)) + basis->gs_(A);
    int n = basis->ns_(A);
    VectorXcd gg(n+1+1);
    gtoint2n(n+1, g, &gg);
    auto ref = gg(n+1) * pow(basis->Ns_(A), 2);
    ASSERT_DOUBLE_EQ(real(ref), real(M20(A,A)));
    ASSERT_DOUBLE_EQ(real(ref), real(M02(A,A)));
  }

  int A = 1;
  int B = 2;
  ASSERT_DOUBLE_EQ(real(M01(A,B)), real(M10(A,B))) << "A:"<<A<<endl<<"B:"<<B<<endl;
  ASSERT_DOUBLE_EQ(imag(M01(A,B)), imag(M10(A,B))) << "A:"<<A<<endl<<"B:"<<B<<endl;
  
  ASSERT_DOUBLE_EQ(real(M02(A,B)), real(M11(A,B))) << "A:"<<A<<endl<<"B:"<<B<<endl;
  ASSERT_DOUBLE_EQ(imag(M02(A,B)), imag(M11(A,B))) << "A:"<<A<<endl<<"B:"<<B<<endl;
  
  ASSERT_DOUBLE_EQ(real(M20(A,B)), real(M11(A,B))) << "A:"<<A<<endl<<"B:"<<B<<endl;
  ASSERT_DOUBLE_EQ(imag(M20(A,B)), imag(M11(A,B))) << "A:"<<A<<endl<<"B:"<<B<<endl;

  delete basis;
}
TEST(utest_pwgto, test_harmonic) {
  int num = 4;
  double x0 = 0.3d;
  double k = 0.39;
  double m = 2.0;
  double w = sqrt(k/m);

  VectorXi ns(num); ns << 0, 1, 2, 3;
  vector<Operator> ops = {kOp0, kOp2, kOpP2};
  PlaneWaveGTO *basis = new PlaneWaveGTO(ns, ops);

  complex<double> g = m*w/2;
  basis->gs_ = VectorXcd::Ones(num)*g;
  basis->Rs_ = VectorXd::Ones(num )*x0;
  basis->Ps_ = VectorXd::Zero(num);
  basis->setup();

  MatrixXcd P2(num,num), R2(num,num), H(num,num), S(num,num);

  basis->overlap(kOp0, kOp0,  &S);
  basis->overlap(kOp0, kOp2,  &P2);
  basis->overlap(kOp0, kOpP2, &R2);
  H = P2/(2*m) + k/2*R2;

  cout << "H:" << endl;
  cout << H << endl;
  cout << "w:" << endl;
  cout << w /2 << endl;

  GeneralizedSelfAdjointEigenSolver<MatrixXcd> es;
  es.compute(H, S);

  cout << "Es:" << endl;
  cout << es.eigenvalues() << endl;

  ASSERT_DOUBLE_EQ(w/2, es.eigenvalues()(0));
  
  delete basis;  
}
