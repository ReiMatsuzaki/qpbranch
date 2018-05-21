#include <iostream>
#include <gtest/gtest.h>
#include <qpbranch/gtestplus.hpp>
#include <qpbranch/dy_branch.hpp>
#include <qpbranch/dy_aadf.hpp>

using namespace std;
using namespace qpbranch;

class TestMonoGauss : public ::testing::Test {
protected:
  int num;
  DyAadf *dy_aadf;
  DySetPoly *dy_old;
  void SetUp() {
    
    // basic info
    num = 1;
    VectorXi ns(num); ns << 0;
    
    // system
    complex<double> v0(0.225, 0.0);
    complex<double> b  = 1.0;
    complex<double> q0  = 0.0;
    auto v = new OperatorGausspot(v0, b, q0);
    double m = 2000.0;
    
    // AADF
    dy_aadf = new DyAadf(v, ns, "thawed");
    dy_aadf->rho_ = 1.0;
    dy_aadf->lambda_ = 0.0;
    dy_aadf->q0_ = -10.0;
    dy_aadf->p0_ = 21.2312; // = sqrt(real(v0) *2.0*m)
    dy_aadf->m_ = m;
    dy_aadf->SetUp();
    dy_aadf->UpdateBasis();
    
    // Ordinary
    dy_old = new DySetPoly(v, ns, "thawed");
    dy_old->gr0_ = 1.0/(4*pow(dy_aadf->rho_, 2));
    dy_old->gi0_ = -dy_aadf->lambda_/(2*dy_aadf->rho_);
    dy_old->q0_  = dy_aadf->q0_;
    dy_old->p0_  = dy_aadf->p0_;
    dy_old->m_   = dy_aadf->m_;  
    dy_old->SetUp();
    dy_old->UpdateBasis();
  }
  void TearDown() {
    delete dy_aadf;
    delete dy_old;
  }  
};

TEST_F(TestMonoGauss, At) {

  int nx = 2;
  VectorXd xs(nx); xs << 0.3, 0.4;
  VectorXcd cs(num); cs << 1.0;
  VectorXcd y1s(nx), y2s(nx);
  dy_aadf->basis_->At(dy_aadf->id_, cs, xs, &y1s);
  dy_old->basis_->At( dy_aadf->id_, cs, xs, &y2s);
  EXPECT_VECTORXCD_EQ(y1s, y2s);
}
TEST_F(TestMonoGauss, Hamiltonian) {

  MatrixXcd H1(num,num), H2(num,num);
  dy_aadf->Hamiltonian(dy_aadf->id_, &H1);
  dy_old->Hamiltonian( dy_old->id_,  &H2);

  EXPECT_MATRIXXCD_EQ(H1, H2);
  
}
TEST_F(TestMonoGauss, Dotx) {
  VectorXd dotx1(4), dotx2(4);
  dy_aadf->DotxQhamilton(&dotx1);
  dy_old->DotxQhamilton(&dotx2);
  EXPECT_VECTORXD_NEAR(dotx1, dotx2, pow(10.0,-10));
}
