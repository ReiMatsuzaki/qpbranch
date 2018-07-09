#include <iostream>
#include <gtest/gtest.h>
#include <qpbranch/gtestplus.hpp>
#include <qpbranch/dy_mono.hpp>
#include <qpbranch/dy_nac.hpp>
#include <qpbranch/con.hpp>

using namespace std;
using namespace qpbranch;
using namespace mangan4;

class TestMonoState : public ::testing::Test {
protected:
  int num;
  DyNac *dy_nac;
  DySetPoly *dy_old;
  void SetUp() {
    
    // basic info
    num = 1;
    VectorXi ns(num); ns << 0;
    
    // syste
    int numI = 1;
    complex<double> v0(0.225, 0.0);
    complex<double> b  = 1.0;
    complex<double> q0  = 0.0;
    auto v = new OperatorGausspot(v0, b, q0);
    double m = 2000.0;
    OpMat opHeIJ(boost::extents[numI][numI]);
    OpMat opXkIJP(boost::extents[numI][numI]);
    opHeIJ[0][0] = v;
    opXkIJP[0][0] = nullptr;
    
    // dynamics for Non Adiabatic Coupled 
    dy_nac = new DyNac(opHeIJ, opXkIJP, ns, "thawed", 2, 0.01);
    dy_nac->gamma0_ = complex<double>(1.0, 0.0);
    dy_nac->q0_ = -10.0;
    dy_nac->p0_ = 21.2312; // = sqrt(real(v0) *2.0*m)
    dy_nac->m_ = m;    
    dy_nac->SetUp();
    dy_nac->UpdateBasis();
    
    // Old
    dy_old = new DySetPoly(v, ns, "thawed", 2, 0.01);
    dy_old->m_   = dy_nac->m_;
    dy_old->gr0_ = dy_nac->gamma0_.real();
    dy_old->gi0_ = dy_nac->gamma0_.imag();
    dy_old->q0_  = dy_nac->q0_;
    dy_old->p0_  = dy_nac->p0_;
    dy_old->m_   = dy_nac->m_;
    dy_old->SetUp();
    dy_old->UpdateBasis();
  }
  void TearDown() {
    delete dy_nac;
    delete dy_old;
  }  
};
TEST_F(TestMonoState, At) {

  int nx = 2;
  VectorXd xs(nx); xs << -10.3, -9.6;
  VectorXcd y1s(nx), y2s(nx);
  dy_nac->At( 0, xs, &y1s);
  dy_old->At( xs, &y2s);
  EXPECT_VECTORXCD_EQ(y1s, y2s);
  
}
TEST_F(TestMonoState, Hamiltonian) {

  MatrixXcd H1(num,num), H2(num,num);
  dy_nac->Hamiltonian(dy_nac->id_, &H1);
  dy_old->Hamiltonian(dy_old->id_, &H2);
  EXPECT_MATRIXXCD_EQ(H1, H2);

  dy_nac->Hamiltonian(dy_nac->DR_, &H1);
  dy_old->Hamiltonian(dy_old->DR_, &H2);
  EXPECT_MATRIXXCD_EQ(H1, H2);

  dy_nac->Hamiltonian(dy_nac->DP_, &H1);
  dy_old->Hamiltonian(dy_old->DP_, &H2);
  EXPECT_MATRIXXCD_EQ(H1, H2);

  dy_nac->Hamiltonian(dy_nac->Dgr_, &H1);
  dy_old->Hamiltonian(dy_old->Dgr_, &H2);
  EXPECT_MATRIXXCD_EQ(H1, H2);

  dy_nac->Hamiltonian(dy_nac->Dgi_, &H1);
  dy_old->Hamiltonian(dy_old->Dgi_, &H2);
  EXPECT_MATRIXXCD_EQ(H1, H2);
  
}
TEST_F(TestMonoState, Dotx) {
  VectorXd dotx1(4), dotx2(4);
  dy_nac->Dotx(&dotx1);
  dy_old->Dotx(&dotx2);
  double tol = pow(10.0, -10.0);
  EXPECT_VECTORXD_NEAR(dotx1, dotx2, tol);
}

TEST(TestDelta, TestTully1) {
  
  Con& con = Con::getInstance();
  con.set_root("_con_delta_tully1");
  
  // potential
  int numI = 2;
  OpSpMat opHeIJ(boost::extents[numI][numI]);
  VectorXd xs = VectorXd::LinSpaced(500, -10.0, 10.0);
  double A = 0.01;
  double B = 1.6;
  double C = 0.005;
  double D = 1.0;
  VectorXd y00s(xs.size());
  VectorXd y01s(xs.size());
  VectorXd y10s(xs.size());
  VectorXd y11s(xs.size());
  for(int i = 0; i < xs.size(); i++) {
    auto x = xs(i);
    if(x>0)
      y00s(i) = +A*(1-exp(-B*x));
    else
      y00s(i) = -A*(1-exp(+B*x));
    y01s(i) = C*exp(-D*x*x);
    y10s(i) = y01s(i);
    y11s(i) = -y00s(i);
  }
  opHeIJ[0][0] = new OperatorSpline(xs, y00s);
  opHeIJ[0][1] = new OperatorSpline(xs, y01s);
  opHeIJ[1][0] = new OperatorSpline(xs, y10s);
  opHeIJ[1][1] = new OperatorSpline(xs, y11s);

  // dynamics
  auto dy = new DyNacDelta(opHeIJ);
  dy->type_intenuc_ = "RK4";
  dy->q0_ = -7.0;
  //  dy->q0_ = -1.0;
  dy->p0_ = 15.0;
  dy->m_  = 2000.0;
  dy->SetUp();

  // time
  auto dt = 20.0;
  //  auto nt = 2;
  auto nt = 100;

  // calculate
  for(int it = 0; it < nt; it++) {
    con.write_f("t", it, it*dt);
    dy->DumpCon(it);
    if(it!=nt-1)
      dy->Update(dt);
  }
  
  delete dy;
  
}
TEST(TestDelta, TestBasic) {

 Con& con = Con::getInstance();
  con.set_root("_con_delta_tully1");
  
  // potential
  int numI = 2;
  OpSpMat opHeIJ(boost::extents[numI][numI]);
  VectorXd xs = VectorXd::LinSpaced(500, -10.0, 10.0);
  double A = 0.01;
  double B = 1.6;
  double C = 0.005;
  double D = 1.0;
  VectorXd y00s(xs.size());
  VectorXd y01s(xs.size());
  VectorXd y10s(xs.size());
  VectorXd y11s(xs.size());
  for(int i = 0; i < xs.size(); i++) {
    auto x = xs(i);
    if(x>0)
      y00s(i) = +A*(1-exp(-B*x));
    else
      y00s(i) = -A*(1-exp(+B*x));
    y01s(i) = C*exp(-D*x*x);
    y10s(i) = y01s(i);
    y11s(i) = -y00s(i);
  }
  opHeIJ[0][0] = new OperatorSpline(xs, y00s);
  opHeIJ[0][1] = new OperatorSpline(xs, y01s);
  opHeIJ[1][0] = new OperatorSpline(xs, y10s);
  opHeIJ[1][1] = new OperatorSpline(xs, y11s);

  // dynamics
  auto dy = new DyNacDelta(opHeIJ);
  dy->type_intenuc_ = "RK4";
  dy->q0_ = 0.5;
  dy->p0_ = 10.0;
  dy->m_  = 2000.0;
  dy->SetUp();

  // check Hamiltonian
  MatrixXcd H(numI, numI);
  dy->Hamiltonian(&H);
  EXPECT_DOUBLE_EQ(0.0, imag(H(0,0)));
  EXPECT_DOUBLE_EQ(real(H(1,0)), real(H(0,1)));
  EXPECT_DOUBLE_EQ(imag(H(1,0)), -imag(H(0,1)));  
  EXPECT_DOUBLE_EQ(0.0, imag(H(1,1)));
  
}

TEST(TestNac, TestBasic) {
  
  // basic info
  int numA = 1;
  VectorXi ns(numA); ns << 0;
  
  // potential
  int numI = 2;
  OpMat opHeIJ(boost::extents[numI][numI]);
  OpMat opXkIJP(boost::extents[numI][numI]);
  VectorXd xs = VectorXd::LinSpaced(500, -10.0, 10.0);
  double A = 0.01;
  double B = 1.6;
  double C = 0.005;
  double D = 1.0;
  VectorXd y00s(xs.size());
  VectorXd y01s(xs.size());
  VectorXd y10s(xs.size());
  VectorXd y11s(xs.size());
  for(int i = 0; i < xs.size(); i++) {
    auto x = xs(i);
    if(x>0)
      y00s(i) = +A*(1-exp(-B*x));
    else
      y00s(i) = -A*(1-exp(+B*x));
    y01s(i) = C*exp(-D*x*x);
    y10s(i) = y01s(i);
    y11s(i) = -y00s(i);
  }
  opHeIJ[0][0] = new OperatorSpline(xs, y00s);
  opHeIJ[0][1] = new OperatorSpline(xs, y01s);
  opHeIJ[1][0] = new OperatorSpline(xs, y10s);
  opHeIJ[1][1] = new OperatorSpline(xs, y11s);
  opXkIJP[0][0] = nullptr;
  opXkIJP[0][1] = nullptr;
  opXkIJP[1][1] = nullptr;
  opXkIJP[1][0] = nullptr;

  // dynamics
  auto dy = new DyNac(opHeIJ, opXkIJP, ns, "thawed", 2, 0.01);
  dy->type_intenuc_ = "RK4";
  dy->gamma0_ = 1.0;
  dy->q0_ = -1.0;
  dy->p0_ = 15.0;
  dy->m_  = 2000.0;
  dy->SetUp();

  // check matrix
  MatrixXcd H(numI, numI);
  dy->Hamiltonian(dy->id_, &H);
  EXPECT_DOUBLE_EQ(0.0, imag(H(0,0)));
  EXPECT_DOUBLE_EQ(real(H(1,0)), real(H(0,1)));
  EXPECT_DOUBLE_EQ(imag(H(1,0)), -imag(H(0,1)));  
  EXPECT_DOUBLE_EQ(0.0, imag(H(1,1)));

}

TEST(SolveModel, Harmonic) {

  Con& con = Con::getInstance();
  con.set_root("_con_fp");
  auto i = complex<double>(0,1);

  // basic info
  int num = 1;
  VectorXi ns(num); ns << 0;

  // system (Harmonic)
  int numI = 1;
  OpMat opHeIJ(boost::extents[numI][numI]);
  OpMat opXkIJP(boost::extents[numI][numI]);
  auto m = 2.0;
  auto w = 0.8;  
  auto a = m*w/2;
  VectorXd xs = VectorXd::LinSpaced(400, -5.0, 5.0);
  VectorXd ys = (m*w*w/2)*xs.array()*xs.array();
  opHeIJ[0][0]  = new OperatorSpline(xs, ys);
  opXkIJP[0][0] = nullptr;

  // initial value
  auto q0 = -2.0;
  auto p0 = 0.8;
  auto g0 = complex<double>(1.2, 0.0);
  
  // dynamics
  auto dy = new DyNac(opHeIJ, opXkIJP, ns, "thawed", 2, 0.01);
  dy->type_intenuc_ = "RK4";
  dy->gamma0_ = g0;
  dy->q0_ = q0;
  dy->p0_ = p0;
  dy->m_  = m;
  dy->SetUp();

  // time
  auto dt = 0.1;
  auto nt = 100;

  // exact
  double q, p;
  complex<double> g;

  // calculate
  for(int it = 0; it < nt; it++) {
    
    // time
    con.write_f("t", it, it*dt);

    // analytic solution
    auto t = dt*it;
    auto cwt = cos(w*t);
    auto swt = sin(w*t);
    g = a * (g0*cwt+i*a*swt) / (i*g0*swt+a*cwt);
    q = q0*cwt + p0/(m*w)*swt;
    p = p0*cwt - (m*w*q0)*swt;
    con.write_f("exact_q0", it, q);
    con.write_f("exact_p0", it, p);
    con.write_f("exact_gr0", it, real(g));
    con.write_f("exact_gi0", it, imag(g));

    // dump calculation
    dy->DumpCon(it, "");

    // update
    if(it!=nt-1)
      dy->Update(dt);
  }
  
  auto tol = 0.00002;
  EXPECT_NEAR(real(g), real(dy->gamma0_), tol);
  EXPECT_NEAR(imag(g), imag(dy->gamma0_), tol);
  EXPECT_NEAR(q,       dy->q0_,           tol);
  EXPECT_NEAR(p,       dy->p0_,           tol);
  
  delete dy;
}
TEST(SolveModel, TwoState) {
  // Two state model with constant potential.
  // see ${QPBRANCH}/doc/two_state

  Con& con = Con::getInstance();
  con.set_root("_con_2state");
  complex<double> i(0, 1);

  // nuclear basis
  int numA = 1;
  VectorXi ns(numA); ns << 0;

  // potential
  int numI = 2;
  OpMat opHeIJ(boost::extents[numI][numI]);
  OpMat opXkIJP(boost::extents[numI][numI]);
  VectorXd xs = VectorXd::LinSpaced(500, -10.0, 10.0);
  int nx = xs.size();
  auto Ea = 0.0;
  auto Eb = 1.0;
  auto w0 = Eb-Ea;
  auto x = 1.0;
  VectorXd y00s = VectorXd::Constant(nx, Ea);
  VectorXd y11s = VectorXd::Constant(nx, Eb);
  VectorXd y01s = VectorXd::Constant(nx, -x/2);
  VectorXd y10s = VectorXd::Constant(nx, -x/2);
  opHeIJ[0][0] = new OperatorSpline(xs, y00s);
  opHeIJ[0][1] = new OperatorSpline(xs, y01s);
  opHeIJ[1][0] = new OperatorSpline(xs, y10s);
  opHeIJ[1][1] = new OperatorSpline(xs, y11s);
  opXkIJP[0][0] = nullptr;
  opXkIJP[0][1] = nullptr;
  opXkIJP[1][1] = nullptr;
  opXkIJP[1][0] = nullptr;

  // dynamics
  auto dy = new DyNac(opHeIJ, opXkIJP, ns, "frozen", 0, 0.01);
  dy->type_intenuc_ = "RK4";
  dy->gamma0_ = 1.0;
  dy->q0_ = 0.0;
  dy->p0_ = 0.0;
  dy->m_  = 2000.0;
  dy->SetUp();

  // time
  auto dt = 0.1;
  auto nt = 100;

  // exact solution
  VectorXcd cs(2);

  // calculate
  for(int it = 0; it < nt; it++) {
    
    // time
    auto t = dt*it;
    con.write_f("t", it, t);

    // analytic solution
    double Omega = sqrt(w0*w0 + x*x);
    auto arg = Omega*t/2.0;
    cs(0) = exp(-i*Ea*t) * exp(-i*w0*t/2.0) * (cos(arg) + i*w0/Omega*sin(arg));
    cs(1) = exp(-i*Eb*t) * exp(+i*w0*t/2.0) * x/(2*Omega) * 2.0*i*sin(arg);
    con.write_f1("exact_cr", it, cs.real());
    con.write_f1("exact_ci", it, cs.imag());

    // dump calculation results
    dy->DumpCon(it);
    
    if(it!=nt-1)
      dy->Update(dt);
  }

  auto tol = 0.00001;
  EXPECT_NEAR( real(cs[1]/cs[0]), real(dy->cAI_[1]/dy->cAI_[0]), tol);
  EXPECT_NEAR( imag(cs[1]/cs[0]), imag(dy->cAI_[1]/dy->cAI_[0]), tol);
  
}
TEST(SolveModel, TwoConst) {

  // nuclear basis
  int numA = 1;
  VectorXi ns(numA); ns << 0;

  // potential
  int numI = 2;
  OpMat opHeIJ(boost::extents[numI][numI]);
  OpMat opXkIJP(boost::extents[numI][numI]);
  VectorXd xs = VectorXd::LinSpaced(500, -10.0, 10.0);
  int nx = xs.size();
  auto Delta = 1.0;
  auto mu_eps = 1.0;
  VectorXd y00s = VectorXd::Zero(nx);
  VectorXd y11s = VectorXd::Constant(nx, Delta);
  VectorXd y01s = VectorXd::Constant(nx, -mu_eps);
  VectorXd y10s = y01s;
  opHeIJ[0][0] = new OperatorSpline(xs, y00s);
  opHeIJ[0][1] = new OperatorSpline(xs, y01s);
  opHeIJ[1][0] = new OperatorSpline(xs, y10s);
  opHeIJ[1][1] = new OperatorSpline(xs, y11s);
  opXkIJP[0][0] = nullptr;
  opXkIJP[0][1] = nullptr;
  opXkIJP[1][1] = nullptr;
  opXkIJP[1][0] = nullptr;

  // dynamics
  auto dy = new DyNac(opHeIJ, opXkIJP, ns, "frozen", 0, 0.01);
  dy->type_intenuc_ = "RK4";
  dy->gamma0_ = 1.0;
  dy->q0_ = 0.0;
  dy->p0_ = 0.0;
  dy->m_  = 2000.0;
  dy->SetUp();

  // Matrix
  int n = numA*numI;
  MatrixXcd H(n,n), S(n,n);
  dy->Overlap(&S);
  dy->Hamiltonian(dy->id_, &H);

  EXPECT_DOUBLE_EQ(Delta,   real(H(1,1)-H(0,0)));
  EXPECT_DOUBLE_EQ(-mu_eps, real(H(0,1)));
  EXPECT_DOUBLE_EQ(-mu_eps, real(H(1,0)));

  EXPECT_DOUBLE_EQ(1.0, real(S(0,0)));
  EXPECT_DOUBLE_EQ(1.0, real(S(1,1)));
  EXPECT_DOUBLE_EQ(0.0, real(S(0,1)));
  EXPECT_DOUBLE_EQ(0.0, real(S(1,0)));
}
TEST(SolveModel, Tully1) {
  // Calculation for Tully Type I model.
  // See numerically exact resutls at ${WPDY}/example/tully1

  Con& con = Con::getInstance();
  con.set_root("_con_tully1");

  // basic info
  int numA = 1;
  VectorXi ns(numA); ns << 0;
  
  // potential
  int numI = 2;
  OpMat opHeIJ(boost::extents[numI][numI]);
  OpMat opXkIJP(boost::extents[numI][numI]);
  VectorXd xs = VectorXd::LinSpaced(500, -10.0, 10.0);
  double A = 0.01;
  double B = 1.6;
  double C = 0.005;
  double D = 1.0;
  VectorXd y00s(xs.size());
  VectorXd y01s(xs.size());
  VectorXd y10s(xs.size());
  VectorXd y11s(xs.size());
  for(int i = 0; i < xs.size(); i++) {
    auto x = xs(i);
    if(x>0)
      y00s(i) = +A*(1-exp(-B*x));
    else
      y00s(i) = -A*(1-exp(+B*x));
    y01s(i) = C*exp(-D*x*x);
    y10s(i) = y01s(i);
    y11s(i) = -y00s(i);
  }
  opHeIJ[0][0] = new OperatorSpline(xs, y00s);
  opHeIJ[0][1] = new OperatorSpline(xs, y01s);
  opHeIJ[1][0] = new OperatorSpline(xs, y10s);
  opHeIJ[1][1] = new OperatorSpline(xs, y11s);
  opXkIJP[0][0] = nullptr;
  opXkIJP[0][1] = nullptr;
  opXkIJP[1][1] = nullptr;
  opXkIJP[1][0] = nullptr;

  // dynamics
  auto dy = new DyNac(opHeIJ, opXkIJP, ns, "thawed", 1, 0.01);
  dy->type_intenuc_ = "RK4";
  dy->gamma0_ = 1.0;
  dy->q0_ = -7.0;
  dy->p0_ = 15.0;
  dy->m_  = 2000.0;
  dy->SetUp();

  // save
  con.write_f1("xs", 0, xs);
  VectorXcd ys(xs.size());
  ((OperatorSpline*)dy->opHeIJ_[0][0])->At(xs, &ys); con.write_f1("v11", 0, ys.real());
  ((OperatorSpline*)dy->opHeIJ_[0][1])->At(xs, &ys); con.write_f1("v12", 0, ys.real());
  ((OperatorSpline*)dy->opHeIJ_[1][0])->At(xs, &ys); con.write_f1("v21", 0, ys.real());
  ((OperatorSpline*)dy->opHeIJ_[1][1])->At(xs, &ys); con.write_f1("v22", 0, ys.real());

  // time
  auto dt = 20.0;
  auto nt = 100;

  // calculate
  for(int it = 0; it < nt; it++) {
    con.write_f("t", it, it*dt);
    dy->DumpCon(it);
    for(int I = 0; I < 2; I++) {
      VectorXcd ys(xs.size());
      dy->At(I, xs, &ys);
      con.write_f1("psi_"+to_string(I)+"_re", it, ys.real() );
      con.write_f1("psi_"+to_string(I)+"_im", it, ys.imag() );
    }
  
    if(it!=nt-1)
      dy->Update(dt);
  }
  
  delete dy;
}
TEST(SolveModel, Tully2) {
  // Calculation for Tully Type II model.
  // See numerically exact resutls at ${WPDY}/example/tully2_dvr

  Con& con = Con::getInstance();
  con.set_root("_con_tully2");

  // basic info
  int numA = 1;
  VectorXi ns(numA); ns << 0;
  
  // potential
  int numI = 2;
  OpMat opHeIJ(boost::extents[numI][numI]);
  OpMat opXkIJP(boost::extents[numI][numI]);
  VectorXd xs = VectorXd::LinSpaced(500, -10.0, 10.0);
  int nx = xs.size();
  double A = 0.1;
  double B = 0.28;
  double C = 0.015;
  double D = 0.06;
  double E = 0.05;
  VectorXd y00s = VectorXd::Zero(nx);
  VectorXd y01s(nx), y10s(nx), y11s(nx);
  for(int i = 0; i < xs.size(); i++) {
    auto x = xs(i);
    y01s(i) = C*exp(-D*x*x);
    y10s(i) = y01s(i);
    y11s(i) = -A*exp(-B*x*x) + E;
  }
  opHeIJ[0][0] = new OperatorSpline(xs, y00s);
  opHeIJ[0][1] = new OperatorSpline(xs, y01s);
  opHeIJ[1][0] = new OperatorSpline(xs, y10s);
  opHeIJ[1][1] = new OperatorSpline(xs, y11s);
  opXkIJP[0][0] = nullptr;
  opXkIJP[0][1] = nullptr;
  opXkIJP[1][1] = nullptr;
  opXkIJP[1][0] = nullptr;

  // save
  con.write_f1("xs", 0, xs);

  // dynamics
  auto dy = new DyNac(opHeIJ, opXkIJP, ns, "thawed", 1, 0.01);
  dy->type_intenuc_ = "RK4";
  dy->gamma0_ = 1.0;
  dy->q0_ = -7.0;
  dy->p0_ = 20.0;
  dy->m_  = 2000.0;
  dy->SetUp();

  // time
  auto dt = 10.0;
  auto nt = 150;

  // calculate
  for(int it = 0; it < nt; it++) {
    con.write_f("t", it, it*dt);
    dy->DumpCon(it);
    for(int I = 0; I < 2; I++) {
      VectorXcd ys(xs.size());
      dy->At(I, xs, &ys);
      con.write_f1("psi_"+to_string(I)+"_re", it, ys.real() );
      con.write_f1("psi_"+to_string(I)+"_im", it, ys.imag() );
    }
  
    if(it!=nt-1)
      dy->Update(dt);
  }
  
  delete dy;
}

class TestNonCoupled : public ::testing::Test {
protected:
  DyNac *dy;
  int numI;
  void SetUp() {

    // nuclear part
    int numA = 1;
    VectorXi ns(numA); ns << 0;

    // electron part
    int numI = 2;
    OpMat HeIJ(boost::extents[numI][numI]);
    OpMat XkIJP(boost::extents[numI][numI]);

    VectorXd xs = VectorXd::LinSpaced(100, -5.0, 5.0);
    VectorXd ys1 = xs.array()*xs.array();
    VectorXd ys2 = xs.array()*xs.array() + 1.0;
    HeIJ[0][0] = new OperatorSpline(xs, ys1);
    HeIJ[1][1] = new OperatorSpline(xs, ys2);
    HeIJ[0][1] = nullptr;
    HeIJ[1][0] = nullptr;
    XkIJP[0][0] = nullptr;
    XkIJP[0][1] = nullptr;
    XkIJP[1][0] = nullptr;
    XkIJP[1][1] = nullptr;

    dy = new DyNac(HeIJ, XkIJP, ns, "frozen", 2, 0.01);
    dy->m_ = 1.0;
    dy->gamma0_ = 1.0;
    dy->q0_ = 2.0;
    dy->p0_ = 3.0;
    dy->SetUp();
    dy->UpdateBasis();
  }
  void TearDown() {
    delete dy;
  }  
};
TEST_F(TestNonCoupled, Update) {
  VectorXcd c0 = dy->cAI_;
  dy->Update(0.01);
  VectorXcd c1 = dy->cAI_;
  for(int I = 0; I < 2; I++)
    EXPECT_DOUBLE_EQ(abs(c0[I]), abs(c1[I]));
}

TEST(TestCoupled, Matrix) {

  // nuclear part
  int numA = 1;
  VectorXi ns(numA); ns << 0;

  // electron part
  int numI = 2;
  OpMat HeIJ( boost::extents[numI][numI]);
  OpMat XkIJP(boost::extents[numI][numI]);
  VectorXd xs = VectorXd::LinSpaced(200, -10.0, 10.0);
  double A = 0.01;
  double B = 1.6;
  double C = 0.005;
  double D = 1.0;
  VectorXd y00s(xs.size());
  VectorXd y01s(xs.size());  
  VectorXd y11s(xs.size());
  VectorXd y10s(xs.size());
  for(int i = 0; i < xs.size(); i++) {
    double x = xs(i);
    if(x>0) 
      y00s(i) = A*(1-exp(-B*x));
    else
      y00s(i) = -A*(1-exp(+B*x));
    y01s(i) = C*exp(-D*x*x);
    y10s(i) = y01s(i);
    y11s(i) = -y00s(i);
  }
  HeIJ[0][0] = new OperatorSpline(xs, y00s);
  HeIJ[0][1] = new OperatorSpline(xs, y01s);
  HeIJ[1][0] = new OperatorSpline(xs, y10s);
  HeIJ[1][1] = new OperatorSpline(xs, y11s);
  XkIJP[0][0] = nullptr;
  XkIJP[0][1] = nullptr;
  XkIJP[1][1] = nullptr;
  XkIJP[1][0] = nullptr;
  
  DyNac *dy = new DyNac(HeIJ, XkIJP, ns, "frozen", 0, 0.01);
  dy->m_  = 2000.0;
  dy->q0_ = xs(110);
  dy->p0_ = 0.0;
  dy->SetUp();

  MatrixXcd H(numA*numI,numA*numI);
  dy->Hamiltonian(dy->id_, &H);
  MatrixXcd P2(numA,numA);
  dy->basis_->Matrix(dy->id_, dy->p2_, &P2);
  auto e0 = P2(0,0)/(2.0*dy->m_);
  EXPECT_DOUBLE_EQ(real(e0+y00s(110)), real(H(0,0)));  
  EXPECT_DOUBLE_EQ(real(e0+y11s(110)), real(H(1,1)));
  EXPECT_NEAR(real( y01s(110)), real(H(0,1)), pow(10.0,-10));
  EXPECT_NEAR(real( y10s(110)), real(H(1,0)), pow(10.0,-10));
}


