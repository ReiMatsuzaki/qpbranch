#include <qpbranch/operator.hpp>
#include <qpbranch/pwgto.hpp>
namespace qpbranch {
  OperatorRn::OperatorRn(int n) : n_(n) {}
  OperatorPn::OperatorPn(int n) : n_(n) {}
  OperatorDa::OperatorDa(int id) : id_(id) {}
  OperatorGausspot::OperatorGausspot(complex<double> v0, complex<double> b, complex<double> q0):
    v0_(v0), b_(b), q0_(q0) {}
}



