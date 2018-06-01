#include <iostream>
#include <gtest/gtest.h>
#include <qpbranch/gtestplus.hpp>
#include <qpbranch/dy_branch.hpp>
#include <qpbranch/dy_nac.hpp>

using namespace std;
using namespace qpbranch;

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
    dy_old = new DySetPoly(v, ns, "thawed");
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
  dy_nac->DotxQhamilton(&dotx1);
  dy_old->DotxQhamilton(&dotx2);
  double tol = pow(10.0, -10.0);
  EXPECT_VECTORXD_NEAR(dotx1, dotx2, tol);
}

TEST(SolveModel, free_particle) {

  // basic info
  int num = 1;
  VectorXi ns(num); ns << 0;

  // system
  int numI = 1;
  OpMat opHeIJ(boost::extents[numI][numI]);
  OpMat opXkIJP(boost::extents[numI][numI]);
  auto w = 0.8;
  auto m = 2.0;
  auto a = m*w/2;
  VectorXd xs = VectorXd::LinSpaced(100, -5.0, 5.0);
  VectorXd ys = (m*w*w/2)*xs.array()*xs.array();
  opHeIJ[0][0]  = new OperatorSpline(xs, ys);
  opXkIJP[0][0] = nullptr;

  // initial value
  auto q0 = -2.0;
  auto p0 = +1.0;
  auto g0 = complex<double>(1.2, 0.0);
  //auto g0 = a;
  
  // dynamics
  auto dy = new DyNac(opHeIJ, opXkIJP, ns, "thawed", 2, 0.01);
  dy->gamma0_ = g0;
  dy->q0_ = q0;
  dy->p0_ = p0;
  dy->m_  = m;
  dy->SetUp();

  // calculate
  auto dt = 0.001;
  auto nt = 10;
  for(int it = 0; it < nt; it++)
    dy->Update(dt);
  auto t = dt*nt;

  // analytic solution. see Tannor's book p.29
  auto cwt = cos(w*t);
  auto swt = sin(w*t); 
  auto i = complex<double>(0,1);
  auto g = a * (g0*cwt+i*swt) / (a*cwt+i*swt);
  auto q = q0*cwt + p0/(m*w)*swt;
  auto p = p0*cwt - (m*w*q0)*swt;

  auto tol = pow(10.0, -10.0);
  EXPECT_NEAR(real(g), real(dy->gamma0_), tol);
  EXPECT_NEAR(imag(g), imag(dy->gamma0_), tol);
  EXPECT_NEAR(q,       dy->q0_,           tol);
  EXPECT_NEAR(p,       dy->p0_,           tol);
  
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
  cout << H << endl;
  cout << y00s(110) << " " << y01s(110) << endl;
  cout << y10s(110) << " " << y11s(110) << endl;
  MatrixXcd P2(numA,numA);
  dy->basis_->Matrix(dy->id_, dy->p2_, &P2);
  auto e0 = P2(0,0)/(2.0*dy->m_);
  EXPECT_DOUBLE_EQ(real(e0+y00s(110)), real(H(0,0)));  
  EXPECT_DOUBLE_EQ(real(e0+y11s(110)), real(H(1,1)));
  EXPECT_NEAR(real( y01s(110)), real(H(0,1)), pow(10.0,-10));
  EXPECT_NEAR(real( y10s(110)), real(H(1,0)), pow(10.0,-10));
}
