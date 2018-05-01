#include <qpbranch/pwgto.hpp>
#include <qpbranch/pwgto_buf.hpp>
#include <qpbranch/mathplus.hpp>
using namespace std;
using boost::extents;

namespace qpbranch {
  
  bool is_buf_basic(Operator *op){
    const type_info& op_type = typeid(op);
    if(op_type == typeid(OperatorId) ||
       op_type == typeid(OperatorRn) ||
       op_type == typeid(OperatorPn) ||
       op_type == typeid(OperatorDa) ) {
      return true;
    } else {
      return false;
    }
  }

  PlaneWaveGto::PlaneWaveGto(const VectorXi& ns, const vector<Operator*>& ops):
    num_(ns.size()), nop_(ops.size()), ns_(ns), ops_(ops),
    gs_(num_), Rs_(num_), Ps_(num_), 
    Ns_(num_), maxn_(num_),
    gAB_(num_,num_), eAB_(num_,num_), hAB_(num_,num_), RAB_(num_,num_) {
    
    for(auto it = ops.begin(); it!=ops.end(); ++it) {
      Operator *op = *it;
      buffer_map_[op] = make_op_buff(ns, op);
    }
    
    for(int A = 0; A < num_; A++) {
      int maxn(0);
      for(auto it =  buffer_map_.begin(); it != buffer_map_.end(); ++it) {
	int maxn0 = it->second->maxn(A);
	if(maxn<maxn0)
	  maxn = maxn0;
      }
      maxn_[A] = maxn;
    }
    
    d_ = new multi_array< multi_array<complex<double>,3>*, 2>(extents[num_][num_]);
    for(int A = 0; A < num_; A++) {
      for(int B = 0; B < num_; B++) {
	(*d_)[A][B] = new multi_array<complex<double>,3>(extents[maxn_[A]+1][maxn_[B]+1][maxn_[A]+maxn_[B]+1]);
      }
    }
    
  }
  PlaneWaveGto::~PlaneWaveGto() {}
  void PlaneWaveGto::setup() {

    // normalization term 
    for(int A = 0; A < num_; A++) {
      int n(ns_[A]);
      VectorXcd gg(n+1);
      gtoint2n(n, gs_[A]+conj(gs_[A]), &gg);
      Ns_[A] = real( 1.0/sqrt((gg[n])) );
    }

    // combination values for McCurcie Davidson recursion
    auto ii = complex<double>(0,1);    
    for(int A = 0; A < num_; A++) {
      for(int B = 0; B < num_; B++) {
	auto cgA = conj(gs_(A));
	auto gB = gs_(B);
	auto RA = Rs_[A]; auto RB = Rs_[B]; auto PA = Ps_[A]; auto PB = Ps_[B];
	gAB_(A,B) = cgA + gB;
	RAB_(A,B) = (2.0*cgA*RA + 2.0*gB*RB - ii*PA + ii*PB) / (2.0*gAB_(A,B));
	eAB_(A,B) = exp(-cgA*RA*RA - gB*RB*RB
			+ii*RA*PA - ii*PB*RB
			+gAB_(A,B)*pow(RAB_(A,B), 2));
	hAB_(A,B) = sqrt(M_PI/gAB_(A,B));
	
	hermite_coef_d(gAB_(A,B), RAB_(A,B), Rs_[A], Rs_[B],
		       maxn_[A], maxn_[B], (*d_)[A][B]);
      }
    }    

    // each operator (visitor pattern)
    for(auto it = buffer_map_.begin(); it != buffer_map_.end(); ++it) {
      it->second->setup(this);
    }
    
  }
  void PlaneWaveGto::con() const {
    
  }
  void PlaneWaveGto::matrix(Operator *opbra, Operator *opket, MatrixXcd *res) {
    /* compute overlap matrix. */

    assert(res->rows()>=num_);
    assert(res->cols()>=num_);
    assert(buffer_map_.find(opbra)!=buffer_map_.end());
    assert(buffer_map_.find(opket)!=buffer_map_.end());

    buffer_map_[opket]->matrix(buffer_map_[opbra], this, res);
    
    /*
    if(is_buf_basic(opbra) && is_buf_basic(opket)) {
      OpBufBasic *bra_buf = static_cast<OpBufBasic*>(buffer_map_[opbra]);
      OpBufBasic *ket_buf = static_cast<OpBufBasic*>(buffer_map_[opket]);
      calc_matrix(this, bra_buf, ket_buf, res);
    } else if(is_buf_basic(opbra) && ! is_buf_basic(opket)) {
      OpBufBasic    *bra_buf = static_cast<OpBufBasic*>(   buffer_map_[opbra]);
      OpBufGausspot *ket_buf = static_cast<OpBufGausspot*>(buffer_map_[opket]);
      calc_matrix(this, bra_buf, ket_buf, res);
    } else {
      throw runtime_error("invalid bra type");
    }
    */
    
  }  
  void PlaneWaveGto::at(Operator *op, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res) {

    buffer_map_[op]->at(this, cs, xs, res);
    
    /*
    if(is_buf_basic(op)) {
      OpBufBasic *buf = static_cast<OpBufBasic*>(buffer_map_[op]);
      calc_at(this, buf, cs, xs, res);
    } else {
      throw runtime_error("invalid op");
    }
    */
  }

}

