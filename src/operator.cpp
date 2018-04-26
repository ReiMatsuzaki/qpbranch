#include "operator.hpp"
#include "pwgto.hpp"
namespace qpbranch {
  void OperatorId::call_new(PlaneWaveGto *basis) { basis->new_op(this); }
  void OperatorId::call_setup(PlaneWaveGto *basis) { basis->setup_op(this); }
  
  OperatorRn::OperatorRn(int n) : n_(n) {}
  void OperatorRn::call_new(PlaneWaveGto *basis) { basis->new_op(this); }
  void OperatorRn::call_setup(PlaneWaveGto *basis) { basis->setup_op(this); }
  
  OperatorPn::OperatorPn(int n) : n_(n) {}
  void OperatorPn::call_new(PlaneWaveGto *basis) { basis->new_op(this); }
  void OperatorPn::call_setup(PlaneWaveGto *basis) { basis->setup_op(this); }

  /*
  PlaneWaveGto::PlaneWaveGto(const VectorXi& ns, vector<Operator*> ops) {}
  void PlaneWaveGto::new_op(OperatorRn *op) {}
  void PlaneWaveGto::new_op(OperatorPn *op) {}
  void PlaneWaveGto::setup_op(OperatorRn *op) {}
  void PlaneWaveGto::setup_op(OperatorPn *op) {}
  */

}
