#include <iostream>
#include <gtest/gtest.h>
#include <qpbranch/mathplus.hpp>

using namespace std;
using namespace qpbranch;
using namespace Eigen;

TEST(utest_math, test_fact) {
  ASSERT_EQ(1, Fact(0));
  ASSERT_EQ(1, Fact(1));
  ASSERT_EQ(2, Fact(2));
  ASSERT_EQ(6, Fact(3));
  ASSERT_EQ(24, Fact(4));
  ASSERT_EQ(120, Fact(5));  
}
TEST(utest_math, test_comb) {
  
  ASSERT_EQ(1, Comb(3,0));
  ASSERT_EQ(1, Comb(3,3));

  ASSERT_EQ(3, Comb(3,1));
  ASSERT_EQ(3, Comb(3,2));

  ASSERT_EQ(4, Comb(4,1));
  ASSERT_EQ(6, Comb(4,2));  // 4.3/2.1 = 6
  ASSERT_EQ(4, Comb(4,3));
  
}
TEST(utest_math, test_gtoint2n) {
  VectorXcd res(5);
  complex<double> a(2.0);
  IntGto2N(1, a, &res);
}
TEST(utest_math, test_gtointn_shift) {
  int maxn = 6;
  complex<double> a(1.1), b(1.2);
  double w(0.3);
  VectorXcd calc(maxn+1);
  VectorXd ref(maxn+1);
  double tol;
  tol = pow(10.0, -13);
  ref << 1.10988699911633,
    -0.173721443339947,
    0.268470964852411,
    -0.117552604646441,
    0.193489297804928,
    -0.132504937609851,
    0.231054357413855;

  IntGtoShift(maxn, a, b, w, &calc);

  for(int n = 0; n <= maxn; n++) {
    ASSERT_NEAR(ref(n), real(calc(n)), tol) << "n: " << n;
  }
  
}
TEST(utest_math, test_dwn_gaussint_shift) {
  // see daily/2018/4/17/qpbranch

  complex<double> a(1.2);
  complex<double> b(2.1);
  complex<double> w(1.3);
  
  int maxN = 3;
  VectorXcd gints(maxN+1);
  DwIntGauss(maxN, a, b, w, &gints);

  double tol = pow(10.0, -14.0);
  ASSERT_NEAR(0.268436270813165,  real(gints(0)), tol);
  ASSERT_NEAR(-0.532968014050865, real(gints(1)), tol);
  ASSERT_NEAR(0.648208370655430,  real(gints(2)), tol);
    
}
