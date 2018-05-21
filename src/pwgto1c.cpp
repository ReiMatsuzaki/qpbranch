#include <qpbranch/pwgto1c.hpp>
#include <qpbranch/pwgto_buf.hpp>
#include <qpbranch/mathplus.hpp>
using namespace std;

namespace qpbranch {

  Pwgto1c::Pwgto1c(const VectorXi& ns, double R0, double P0,
		   complex<double> g0, const vector<Operator*>& ops):
    num_(ns.size()), nop_(ops.size()), ns_(ns), ops_(ops), R0_(R0), P0_(P0), g0_(g0), Ns_(num_) {

    for(auto it = ops.begin(); it != ops.end(); ++it) {
      auto op = *it;
      buf_map_[op] = MakeOpBuf(this, op);
    }

    maxn_ = 0;
    for(int A = 0; A < num_; A++) {
      for(auto it =  buf_map_.begin(); it != buf_map_.end(); ++it) {
	maxn_ = std::max(maxn_, it->second->Maxn(A));
      }
    }   
    gints2n_ = VectorXcd::Zero(maxn_+1);    
  }
  Pwgto1c::~Pwgto1c() {}
  void Pwgto1c::SetUp() {
    complex<double> gg = conj(g0_) + g0_;

    // integration
    IntGto2N(maxn_, gg, &gints2n_);
    
    // normalization term 
    for(int A = 0; A < num_; A++) {
      int n = ns_[A];
      Ns_[A] = real( 1.0/sqrt((gints2n_[n])) );
    }

    // each operator (visitor pattern)
    for(auto it = buf_map_.begin(); it != buf_map_.end(); ++it) {
      it->second->SetUp(this);
    }
    
  }
  void Pwgto1c::Matrix(Operator *opbra, Operator *opket, MatrixXcd *res) {
    assert(res->rows()>=num_);
    assert(res->cols()>=num_);
    assert(buf_map_.find(opbra)!=buf_map_.end());
    assert(buf_map_.find(opket)!=buf_map_.end());
    
    buf_map_[opket]->Matrix(buf_map_[opbra], this, res);
  }
  void Pwgto1c::At(Operator *op, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res) {
    assert(buf_map_.find(op)!=buf_map_.end());
    assert(cs.size()==num_);
    assert(xs.size()==res->size());
    
    buf_map_[op]->At(this, cs, xs, res);
  }  
}
