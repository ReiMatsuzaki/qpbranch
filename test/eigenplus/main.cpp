#include <iostream>
#include <gtest/gtest.h>
#include <qpbranch/gtestplus.hpp>
#include <qpbranch/eigenplus.hpp>
#include <qpbranch/con.hpp>

using namespace std;
using namespace qpbranch;

TEST(TestEigenplus, TestZhegv) {

  int n = 3;
  complex<double> i(0,1);
  MatrixXcd H(n,n), S(n,n), U(n,n);
  VectorXd w(n);
  H <<
    1.0,       0.5-0.3*i, 0.2+0.1*i,
    0.5+0.3*i, 1.5,       0.1+0.3*i,
    0.2-0.1*i, 0.1-0.3*i, 0.8;
  S <<
    1.0,       0.0,       0.4+0.5*i,
    0.0,       1.0,       0.1+0.3*i,
    0.4-0.5*i, 0.1-0.3*i, 1.0;
  Zhegv(H, S, &U, &w);

  MatrixXcd UHU = U.adjoint() * H * U;
  MatrixXcd DiagW = w.asDiagonal();

  MatrixXcd USU = U.adjoint() * S * U;
  MatrixXcd id = MatrixXcd::Identity(n, n);
  
  EXPECT_MATRIXXCD_EQ(UHU, DiagW);  
  EXPECT_MATRIXXCD_EQ(USU,  id);
}
TEST(TestEigenplus, TestZgesv) {
  int n = 3;
  complex<double> i(0,1);
  MatrixXcd H(n,n), S(n,n), U(n,n);
  VectorXcd b(n), x(n);
  H <<
    1.0,       0.5-0.3*i, 0.2+0.1*i,
    0.5+0.3*i, 1.5,       0.1+0.3*i,
    0.2-0.1*i, 0.1-0.3*i, 0.8;
  b <<
    1.0,       0.0,       0.4+0.5*i;

  Zgesv(H, b, &x);

  MatrixXcd Hx = H*x;
  EXPECT_VECTORXCD_EQ(Hx, b);

}
TEST(TestEigenplus, TestIntetDiag) {
  complex<double> i(0,1);
  int n = 3;
  MatrixXcd H(n,n);

  H <<
    1.0,       0.5-0.3*i, 0.2+0.1*i,
    0.5+0.3*i, 1.5,       0.1+0.3*i,
    0.2-0.1*i, 0.1-0.3*i, 0.8;
    
  VectorXcd c0(n); c0 << 1.0, 2.0, 0.3;
  c0 /= sqrt(c0.dot(c0));

  // check EOM
  double dt = 0.001;
  VectorXcd cp = c0;
  VectorXcd cm = c0;
  IntetDiag(H, +dt, &cp);
  IntetDiag(H, -dt, &cm);
  VectorXcd idc1 = i*(cp-cm)/(2*dt);
  VectorXcd Hc0 = H*c0;
  EXPECT_VECTORXCD_NEAR(idc1, Hc0, 3*pow(10.0, -6));

  // check norm conversion for long integration
  IntetDiag(H, 11.4, &cp);
  auto norm2 = real(cp.dot(cp));
  EXPECT_NEAR(1.0, norm2, pow(10.0, -10.0));
  
}
TEST(TestEigenplus, TestIntetGdiag) {
  complex<double> i(0,1);
  int n = 3;
  MatrixXcd H(n,n), S(n,n);
  /*
  H <<
    1.0,       0.5-0.3*i, 0.2+0.1*i,
    0.5+0.3*i, 1.5,       0.1+0.3*i,
    0.2-0.1*i, 0.1-0.3*i, 0.8;
  S <<
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0;
  */
  
  H <<
    1.0, 0.0, 0.0,
    0.0, 1.5, 0.0,
    0.0, 0.0, 2.0;
  S <<
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0;

  double dt = 5.0;
  VectorXcd c0(n); c0 << 1.0, 0.0, 0.0;
  VectorXcd c = c0;
  IntetGdiag(H, S, dt, &c);
  auto norm2 = real(c.dot(S*c));
  EXPECT_NEAR(1.0, norm2, pow(10.0, -10.0));
  EXPECT_DOUBLE_EQ(real(c(0)), cos(1.0*dt));
  
}
TEST(TestEigenplus, TestIntetGdiagRabi) {
  // see ${QPBRANCH}/doc/two_state/
  // WARNING:
  // The result in Tannor's book may be wrong.

  mangan4::Con& con = mangan4::Con::getInstance();
  con.set_root("_con_rabi");

  // target system
  complex<double> i(0,1);
  int n = 2;
  MatrixXcd H(n,n), S(n,n);
  auto Ea = 0.0;
  auto Eb = 1.0;
  auto x = 1.0;
  H <<
    Ea,  -x/2,
    -x/2, Eb;
  S <<
    1.0, 0.0, 
    0.0, 1.0;

  // time propagation
  auto dt = 0.1;
  auto nt = 100;
  VectorXcd c0(n); c0 << 1.0, 0.0;
  VectorXcd c(n), ce(n);
  c = c0;
  for(int it = 0; it < nt; it++) {

    // time
    auto t = it*dt;

    // analytic solution
    auto w0 = Eb-Ea;
    auto Omega = sqrt(w0*w0 + x*x);
    auto arg = Omega*t/2;
    ce(0) = exp(-i*Ea*t) * exp(-i*w0*t/2.0) * (cos(arg) + i*w0/Omega*sin(arg));
    ce(1) = exp(-i*Eb*t) * exp(+i*w0*t/2.0) * (x/(2*Omega)) * 2.0*i * sin(arg);

    // dump
    con.write_f("t", it, t);
    con.write_f1("exact_cr", it, ce.real());
    con.write_f1("exact_ci", it, ce.imag());
    con.write_f1("cr",       it, c.real());
    con.write_f1("ci",       it, c.imag());

    // update calculation results
    if(it != nt-1)
      IntetGdiag(H, S, dt, &c);
  }

  auto tol = pow(10.0, -7);

  EXPECT_VECTORXCD_NEAR(c, ce, tol);
}
TEST(TestEigenplus, TestIntetGdiag_nume) {
  complex<double> i(0,1);
  int n = 2;
  MatrixXcd H(n,n), S(n,n);

  auto Ea = 0.2;
  auto Eb = 1.3;
  auto x = 1.4;
  //  auto mu_eps = 0.0;
  H <<
    Ea,  -x/2,
    -x/2, Eb;
  S <<
    1.0, 0.3, 
    0.3, 1.0;

  double dt = 0.001;
  VectorXcd c0(n); c0 << 1.0/sqrt(3.0), sqrt(2.0/3.0);
  VectorXcd cp = c0;
  VectorXcd cm = c0;
  IntetGdiag(H, S, +dt, &cp);
  IntetGdiag(H, S, -dt, &cm);

  VectorXcd iSdc1 = i*S*(cp-cm)/(2*dt);
  VectorXcd Hc0 = H*c0;

  EXPECT_VECTORXCD_NEAR(iSdc1, Hc0, 2*pow(10.0, -6));

}
