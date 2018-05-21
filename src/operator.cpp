#include <qpbranch/operator.hpp>

using namespace std;

namespace qpbranch {
  string OperatorDa::str() const {
    string buf = "Da[";
    buf += "id = " + to_string(id_) + "]";
    return buf;
  }
  void OperatorGausspot::At(const VectorXd& xs, VectorXcd *res) {
    assert(xs.size()==res->size());
    *res = v0_ * (-b_*xs.array().pow(2)).exp();
  }
  string OperatorGausspot::str() const {
    //string buf = "Gausspot[" + to_string(v0_.real()) + "]";
    string buf = "Gausspot[";
    buf += "v0 = (" + to_string(v0_.real()) + ", " +  to_string(v0_.imag()) + "), ";
    buf += "b  = (" + to_string(b_.real()) + ", " +  to_string(b_.imag()) + "), ";
    buf += "q0 = (" + to_string(q0_.real()) + ", " +  to_string(q0_.imag()) + ") ] ";
    return buf;
  }
}



