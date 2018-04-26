#include "operator.hpp"
#include "pwgto.hpp"
namespace qpbranch {
  void OperatorId::call_new(PlaneWaveGto *basis) { basis->new_op(this); }
  
  OperatorRn::OperatorRn(int n) : n_(n) {}
  void OperatorRn::call_new(PlaneWaveGto *basis) { basis->new_op(this); }
  
  OperatorPn::OperatorPn(int n) : n_(n) {}
  void OperatorPn::call_new(PlaneWaveGto *basis) { basis->new_op(this); }

  OperatorDa::OperatorDa(int id) : id_(id) {}
  void OperatorDa::call_new(PlaneWaveGto *basis) { basis->new_op(this); }
  
  OperatorGausspot::OperatorGausspot(complex<double> v0, complex<double> b, complex<double> q0):
    v0_(v0), b_(b), q0_(q0) {}
  void OperatorGausspot::call_new(PlaneWaveGto *basis) { basis->new_op(this); }
  

}
