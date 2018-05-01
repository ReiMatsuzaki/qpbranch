#include <iostream>
#include <gtest/gtest.h>

#include <qpbranch/con.hpp>

using namespace std;
using namespace mangan4;

TEST(utest_con, get_instance) {

  Con& con = Con::getInstance();

  VectorXi x(3); x << 1,2,3;
  con.write("x", 11, x);

  VectorXd y(3); y << 1.1, 2.1, -0.1;
  con.write("y", 1, y);

  MatrixXd z(2,3); 
  z << 1.1, 1.2, 1.3,
    2.1, 2.2, 2.3;
  con.write("z", 1, z);
  
}
