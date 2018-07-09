#include <iostream>
#include <gtest/gtest.h>
#include <qpbranch/gtestplus.hpp>
#include <qpbranch/dy_mono.hpp>
#include <qpbranch/dy_aadf.hpp>

using namespace std;
using namespace qpbranch;

class TestMonoGauss : public ::testing::Test {
protected:
  int num;
  DyAadf *dy_aadf;
  DySetPoly *dy_mono;
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
    //    dy_aadf->p0_ = 21.2312; // = sqrt(real(v0) *2.0*m)
    dy_aadf->m_ = m;
    dy_aadf->SetUp();
    dy_aadf->UpdateBasis();
    
    // Ordinary
    dy_mono = new DySetPoly(v, ns, "thawed", 2, 0.001);
    dy_mono->gr0_ = 1.0/(4*pow(dy_aadf->rho_, 2));
    dy_mono->gi0_ = -dy_aadf->lambda_/(2*dy_aadf->rho_);
    dy_mono->q0_  = dy_aadf->q0_;
    dy_mono->p0_  = dy_aadf->p0_;
    dy_mono->m_   = dy_aadf->m_;  
    dy_mono->SetUp();
    dy_mono->UpdateBasis();
  }
  void TearDown() {
    delete dy_aadf;
    delete dy_mono;
  }  
};

TEST_F(TestMonoGauss, At) {
  int nx = 2;
  VectorXd xs(nx); xs << dy_aadf->q0_+0.3, dy_aadf->q0_;
  VectorXcd y_aadf(nx), y_mono(nx);
  dy_aadf->At(xs, &y_aadf);
  dy_mono->At(xs, &y_mono);
  cerr << y_aadf << endl;
  EXPECT_VECTORXCD_EQ(y_aadf, y_mono);
}
TEST_F(TestMonoGauss, Hamiltonian) {
  MatrixXcd H1(num,num), H2(num,num);
  dy_aadf->Hamiltonian(dy_aadf->id_, &H1);
  dy_mono->Hamiltonian( dy_mono->id_,  &H2);
  EXPECT_MATRIXXCD_EQ(H1, H2);
}
TEST_F(TestMonoGauss, Dotx) {
  VectorXd dotx1(4), dotx2(4);
  dy_aadf->DotxQhamilton(&dotx1);
  dy_mono->DotxQhamilton(&dotx2);
  EXPECT_VECTORXD_NEAR(dotx1, dotx2, pow(10.0,-10));
}
