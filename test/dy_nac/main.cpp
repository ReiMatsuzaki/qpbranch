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
    OpMat opXkIJ(boost::extents[numI][numI]);
    opHeIJ[0][0] = v;
    opXkIJ[0][0] = nullptr;
    
    // dynamics for Non Adiabatic Coupled 
    dy_nac = new DyNac(opHeIJ, opXkIJ, ns, "thawed");
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
  
}
TEST_F(TestMonoState, Dotx) {
  VectorXd dotx1(4), dotx2(4);
  dy_nac->DotxQhamilton("tdvp", &dotx1);
  dy_old->DotxQhamilton(&dotx2);
  EXPECT_VECTORXD_NEAR(dotx1, dotx2, pow(10.0,-10));
}

