#include <iostream>
#include <gtest/gtest.h>
#include "mathplus.hpp"

using namespace std;
using namespace qpbranch;

TEST(utest_math, test_gtoint2n) {
  VectorXcd res(5);
  complex<double> z(2.0);
  gtoint2n(-1, z, &res);
}
