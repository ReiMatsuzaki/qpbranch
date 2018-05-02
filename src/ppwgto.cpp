#include <iostream>
using namespace std;

namespace qpbranch {

  /*
  OpBuf *MakeOpBuf(const Polys& polys, Operator *op) {
    const type_info& optype = typeid(*op);
    if(optype == typeid(OperatorId)) {
      return new OpBufId(ns);
    } else if(optype == typeid(OperatorRn)) {
      auto opop = static_cast<OperatorRn*>(op);
      return new OpBufRn(ns, opop->n());
    } else if(optype == typeid(OperatorPn)) {
      auto opop = static_cast<OperatorPn*>(op);
      return new OpBufPn(ns, opop->n());
    } else if(optype == typeid(OperatorDa)) {
      auto opop = static_cast<OperatorDa*>(op);
      return new OpBufDa(ns, opop->id());
    } else if(optype == typeid(OperatorGausspot)) {
      auto opop = static_cast<OperatorGausspot*>(op);
      return new OpBufGausspot(ns, opop);
    } else {
      cerr << "Unsupported operator" << endl;
      abort();
    }
  }
  */

  PolyPwgto(const Polys& polys, const vector<Operator*>& ops):
    num_(polys.size()), nops_(ops.size()), polys_(polys), ops(ops_),
    gs_(num_), Rs_(num_), Ps_(num_), is_setup_(false),
    Ns_(num_), maxn_(num_),
    gAB_(num_,num_), eAB_(num_,num_), hAB_(num_,num_), RAB_(num_,num_) {

    // allocate buffer for each operator
    for(auto it = ops.begin(); it!=ops.end(); ++it) {
      Operator *op = *it;
      auto ptr = MakeOpBuf(ns, op);
      buffer_map_[op] = ptr;
    }


  }






  
}
