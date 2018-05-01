#include <iostream>

#include <gtest/gtest.h>

#include <qpbranch/mathplus.hpp>
#include <qpbranch/eigenplus.hpp>
#include <qpbranch/pwgto.hpp>
#include <qpbranch/pwgto_buf.hpp>

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
  //  OperatorId opopid;  
  //  auto *op_id = &opopid;
  auto *op_id = new OperatorId();
  
  vector<Operator*> ops = {op_id};
  auto *basis = new Pwgto(ns, ops);
  basis->gs_ << 1.0,  complex<double>(0.9,-0.8), 1.0, complex<double>(0.2, 0.1);
  basis->Rs_ << 0.0, 0.0, 0.2, 0.2;
  basis->Ps_ << 0.0, 0.0, 0.5, 0.5;
  basis->setup();

  MatrixXcd S(num, num);
  basis->matrix(op_id, op_id, &S);

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
TEST(utest_pwgto, multipole) {

  int num = 5;
  auto id = new OperatorId();
  auto R1 = new OperatorRn(1);
  auto R2 = new OperatorRn(2);
  vector<Operator*> ops = {id, R1, R2};
  VectorXi ns(num); ns << 0, 0, 2, 1, 0;
  auto *basis = new Pwgto(ns, ops);
  basis->gs_ << 1.1, 1.1, 1.2, 0.4, 0.8;
  basis->Rs_ << 0.0, 0.1, 0.1, 0.1, 0.3;  
  basis->Ps_ << 0.0, 0.0, 0.3, 0.0, 0.1;
  basis->setup();

  MatrixXcd S(num, num);
  basis->matrix(id, id, &S);
  
  MatrixXcd M01(num, num);
  basis->matrix(id, R1, &M01);

  MatrixXcd M10(num, num);
  basis->matrix(R1, id, &M10);

  MatrixXcd M02(num, num);
  basis->matrix(id, R2, &M02);

  MatrixXcd M11(num, num);
  basis->matrix(R1, R1, &M11);

  MatrixXcd M20(num, num);
  basis->matrix(R2, id, &M20);
  
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
TEST(utest_pwgto, pn) {
  
  int num = 5;
  auto id = new OperatorId();
  auto P1 = new OperatorPn(1);
  auto P2 = new OperatorPn(2);
  vector<Operator*> ops = {id, P1, P2};
  VectorXi ns(num); ns << 0, 0, 2, 1, 0;
  auto *basis = new Pwgto(ns, ops);
  basis->gs_ << 1.1, 1.1, 1.2, 0.4, 0.8;
  basis->Rs_ << 0.0, 0.1, 0.1, 0.1, 0.3;  
  basis->Ps_ << 0.0, 0.0, 0.3, 0.0, 0.1;
  basis->setup();

  MatrixXcd M01(num, num);
  basis->matrix(id, P1, &M01);

  MatrixXcd M10(num, num);
  basis->matrix(P1, id, &M10);

  MatrixXcd M02(num, num);
  basis->matrix(id, P2, &M02);

  MatrixXcd M11(num, num);
  basis->matrix(P1, P1, &M11);

  MatrixXcd M20(num, num);
  basis->matrix(P2, id, &M20);
  
  for(int A = 0; A < num; A++) {
    for(int B = 0; B < num; B++) {
      double tol = pow(10.0, -10);
      ASSERT_NEAR(real(M01(A,B)), real(M10(A,B)), tol) <<"A:"<<A<<endl<<"B:"<<B;
      ASSERT_NEAR(real(M11(A,B)), real(M02(A,B)), tol) <<"A:"<<A<<endl<<"B:"<<B;
      ASSERT_NEAR(real(M11(A,B)), real(M20(A,B)), tol) <<"A:"<<A<<endl<<"B:"<<B;
    }
  }

  delete basis;
}
TEST(utest_pwgto, da) {

  int num = 4;
  
  complex<double> b(1.3);
  complex<double> v0(1.2);
  complex<double> q0(0.0);    

  double dx = 0.01;
  double R0 = 0.1;
  double R1 = 0.2;
  double P0 = 0.3;
  double P1 = 0.2;
  complex<double> g0 = complex<double>(1.0, 0.1);
  complex<double> g1 = complex<double>(1.3, 0.15);
  
  for(int n = 0; n < 2; n++) {
    VectorXi ns = VectorXi::Constant(num, n);
    vector<int> ops_id = {kIdDR, kIdDP, kIdDgr, kIdDgi};  
    for(auto it = ops_id.begin(); it!=ops_id.end(); ++it) {
      
      int kId = *it;
      auto op_id = new OperatorId();
      auto op_da = new OperatorDa(kId);
      auto op_pot = new OperatorGausspot(v0, b, q0);
      vector<Operator*> ops = {op_id, op_da, op_pot};
      
      auto *basis = new Pwgto(ns, ops);
      basis->Rs_ << R0, R0, R0, R1;
      basis->Ps_ << P0, P0, P0, P1;
      basis->gs_ << g0, g0, g0, g1;
      if(kId==kIdDR) {
	basis->Rs_(1) += dx;
	basis->Rs_(2) -= dx;
      } else if(kId==kIdDP) {
	basis->Ps_(1) += dx;
	basis->Ps_(2) -= dx;
      } else if(kId==kIdDgr) {
	basis->gs_(1) += complex<double>(dx, 0.0);
	basis->gs_(2) -= complex<double>(dx, 0.0);
      } else if(kId==kIdDgi) {
	basis->gs_(1) += complex<double>(0.0, dx);
	basis->gs_(2) -= complex<double>(0.0, dx);
      }
      basis->setup();
      
      MatrixXcd M1(num,num), M2(num,num);
      basis->matrix(op_id,  op_da, &M1);
      basis->matrix(op_id,  op_id, &M2);
      EXPECT_NEAR(M1(3, 0).real(), ((M2(3,1)-M2(3,2))/(2*dx)).real(), dx*dx) <<
	"kId:" << kId << endl <<
	"n  :" << n   << endl;
      
      basis->matrix(op_da, op_pot, &M1);
      basis->matrix(op_id, op_pot, &M2);
      EXPECT_NEAR(M1(0, 3).real(), ((M2(1,3)-M2(2,3))/(2*dx)).real(), dx*dx) <<
	"kId:" << kId << endl <<
	"n  :" << n   << endl;
    }
  }
}
TEST(utest_pwgto, test_harmonic) {

  int num = 4;
  double x0 = 0.3d;
  double k = 0.39;
  double m = 2.0;
  double w = sqrt(k/m);

  VectorXi ns(num); ns << 0, 1, 2, 3;
  auto opid = new OperatorId();
  auto opR2 = new OperatorRn(2);
  auto opP2 = new OperatorPn(2);  
  vector<Operator*> ops = {opid, opR2, opP2};
  auto *basis = new Pwgto(ns, ops);

  complex<double> g = m*w/2;
  basis->gs_ = VectorXcd::Ones(num)*g;
  basis->Rs_ = VectorXd::Ones(num )*x0;
  basis->Ps_ = VectorXd::Zero(num);
  basis->setup();

  MatrixXcd P2(num,num), R2(num,num), H(num,num), S(num,num);
  basis->matrix(opid, opid, &S);
  basis->matrix(opid, opR2, &R2);
  basis->matrix(opid, opP2, &P2);
  H = P2/(2*m) + k/2*R2;

  VectorXd eigs;
  MatrixXcd U;
  zhegv(H, S, &U, &eigs);

  int nx(3);
  VectorXd xs(nx);
  xs << 0.2, 0.25, 0.3;

  for(int n = 0; n < num-1; n++) {
    ASSERT_DOUBLE_EQ(w*(n+0.5), eigs(n)) << "n: " << n;
    VectorXcd ys(nx);
    basis->at(opid, U.col(n), xs, &ys);
  }
  
  delete basis;  
  
}
TEST(utest_pwgto, test_gausspot) {

  int num = 1;
  complex<double> b(1.3);
  complex<double> v0(1.2);
  complex<double> q0(0.0);
  
  VectorXi ns(num); ns << 0;
  auto opid = new OperatorId();
  auto opdR = new OperatorDa(kIdDR);
  auto opV  = new OperatorGausspot(v0, b, q0);
  vector<Operator*> ops = {opid, opdR, opV};
  auto *basis = new Pwgto(ns, ops);
  basis->gs_ << 1000.0;
  basis->Rs_ << 0.3;
  basis->setup();
  
  MatrixXcd V(num,num);
  basis->matrix(opid, opV, &V);
  double R0 = basis->Rs_(0);

  double dR = 0.001;
  VectorXd Rs(3); Rs << R0, R0-dR, R0+dR;
  VectorXcd refs(3);
  opV->at(Rs, &refs);
  ASSERT_NEAR(real(refs(0)), real(V(0,0)), 3.0*pow(10.0, -3));

  basis->matrix(opdR, opV, &V);
  complex<double> ref = (refs[2]-refs[1])/(2.0*dR);
  ASSERT_NEAR(real(ref), 2.0*real(V(0,0)), pow(10.0, -3));

  delete basis;
  
}

