#include <qpbranch/operator.hpp>

namespace qpbranch {
  
  OperatorRn::OperatorRn(int n) : n_(n) {}
  OperatorPn::OperatorPn(int n) : n_(n) {}
  OperatorDa::OperatorDa(int id) : id_(id) {}
  OperatorGausspot::OperatorGausspot(complex<double> v0, complex<double> b, complex<double> q0):
    v0_(v0), b_(b), q0_(q0) {}
  void OperatorPot::at0(double x, complex<double> *res) {
    VectorXd xs(1); xs << x;
    VectorXcd ys(1);
    this->at(xs, &ys);
    *res = ys[0];
  }
  void OperatorGausspot::at(const VectorXd& xs, VectorXcd *res) {
    assert(xs.size()==res->size());
    *res = v0_ * (-b_*xs.array().pow(2)).exp();
  }
}



