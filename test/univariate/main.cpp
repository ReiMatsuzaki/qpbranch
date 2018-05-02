#include <iostream>

#include <gtest/gtest.h>

#include <qpbranch/univariate.hpp>

using namespace std;
using namespace qpbranch;

TEST(utest_univariate, first) {
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


