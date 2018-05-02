#include <qpbranch/operator.hpp>

namespace qpbranch {
  /*
  void OperatorPot::at0(double x, complex<double> *res) {
    VectorXd xs(1); xs << x;
    VectorXcd ys(1);
    this->at(xs, &ys);
    *res = ys[0];
  }
  */
  void OperatorGausspot::At(const VectorXd& xs, VectorXcd *res) {
    assert(xs.size()==res->size());
    *res = v0_ * (-b_*xs.array().pow(2)).exp();
  }
}



