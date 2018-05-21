#include <iostream>
#include <complex>

#include <gtest/gtest.h>

#include <qpbranch/univariate.hpp>

using namespace std;
using namespace qpbranch;

TEST(utest_univariate, test_double) {
  
  Univariate<double> p1({ {1.0, 2}, {1.1, 3} });
  Univariate<double> p2({ {1.2, 0}, {1.3, 3} });
  auto p1p2 = p1 + p2;
  double x = 0.2;
  ASSERT_DOUBLE_EQ(p1p2.Sub(x), p1.Sub(x)+p2.Sub(x));

  double dx = 0.0001;
  auto dp = p1p2.Diff(1);
  double ref = (p1p2.Sub(x+dx)-p1p2.Sub(x-dx))/(2*dx);
  ASSERT_NEAR(ref, dp.Sub(x), dx*dx*10);

  auto dp2(dp);
  dp.Refresh();
  ASSERT_DOUBLE_EQ(dp.Sub(x), dp2.Sub(x));
}
TEST(utest_univariate, test_complex) {
  
  Univariate<complex<double> > p1({ {1.0, 2}, {1.1, 3} });
  Univariate<complex<double> > p2({ {1.2, 0}, {1.3, 3} });
  auto p1p2 = p1 + p2;
  complex<double> x(0.2, 0.3);
  ASSERT_DOUBLE_EQ(p1p2.Sub(x).real(), (p1.Sub(x)+p2.Sub(x)).real());
  ASSERT_DOUBLE_EQ(p1p2.Sub(x).imag(), (p1.Sub(x)+p2.Sub(x)).imag());

  double dx = 0.0001;
  auto dp = p1p2.Diff(1);
  complex<double> ref = (p1p2.Sub(x+dx)-p1p2.Sub(x-dx))/(2*dx);
  ASSERT_NEAR(ref.real(), dp.Sub(x).real(), dx*dx*10);
  ASSERT_NEAR(ref.imag(), dp.Sub(x).imag(), dx*dx*10);

  auto dp2(dp);
  dp.Refresh();
  ASSERT_DOUBLE_EQ(dp.Sub(x).real(), dp2.Sub(x).real());
  ASSERT_DOUBLE_EQ(dp.Sub(x).imag(), dp2.Sub(x).imag());
}
TEST(utest_univariate, test_multi) {
  typedef Univariate<complex<double> > Poly;
  Poly p1({ {1.0, 2}, {1.1, 3} });
  Poly p2({ {1.2, 0}, {1.3, 3} });
  Poly p3({ {1.3, 3}});
  Poly p4({ {1.4, 1}});
  Poly p1p2 = p1 * p2;
  complex<double> c(2.1, 0.1);
  Poly p = c * p1p2 + p3 + (-p4);

  complex<double> x(1.0, 0.1);
  auto ref = c * p1.Sub(x) * p2.Sub(x) + p3.Sub(x) - p4.Sub(x);
  auto calc= p.Sub(x);
  
  ASSERT_DOUBLE_EQ(ref.real(), calc.real());
  ASSERT_DOUBLE_EQ(ref.imag(), calc.imag());
}
