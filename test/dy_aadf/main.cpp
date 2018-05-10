#include <iostream>
#include <gtest/gtest.h>
#include <qpbranch/dy_branch.hpp>
#include <qpbranch/dy_aadf.hpp>

using namespace std;
using namespace qpbranch;

TEST(utest_dy_aadf, test_first) {

  // basic info
  int num = 2;
  VectorXi ns(num); ns << 0, 2;

  // system
  complex<double> v0(0.225, 0.0);
  complex<double> b  = 1.0;
  complex<double> q0  = 0.0;
  auto v = new OperatorGausspot(v0, b, q0);
  double m = 2000.0;

  // AADF
  auto aadf = new DyAadf(v, ns, "thawed");
  aadf->alpha_ = 1.0;
  aadf->beta_ = 0.0;
  aadf->q0_ = -10.0;
  aadf->p0_ = sqrt(real(v0) *2.0*m); // => 21.2132
  aadf->p0_ = 21.2132;
  aadf->m_ = m;
  aadf->SetUp();

  // Ordinary
  auto dy = new DySetPoly(v, ns, "thawed");
  dy->gr0_ = 1.0/(4*aadf->alpha_);
  dy->gi0_ = aadf->beta_;
  dy->q0_  = aadf->q0_;
  dy->p0_  = aadf->p0_;
  dy->m_   = aadf->m_;  
  dy->SetUp();

  // compare Hamiltonian
  MatrixXcd H1(num,num), H2(num,num);
  aadf->Hamiltonian(aadf->id_, &H1);
  dy->Hamiltonian(  dy->id_,   &H2);

  cerr << H1 << endl;
  cerr << H2 << endl;
  
}
