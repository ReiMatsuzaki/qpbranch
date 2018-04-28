#include <qpbranch/pwgto.hpp>
#include <qpbranch/pwgto_buf.hpp>
#include <qpbranch/mathplus.hpp>
using namespace std;
using boost::extents;

namespace qpbranch {

  
  void hermite_coef_d(complex<double> gAB, complex<double> wAB, double RA, double RB,
		      int maxnA, int maxnB, multi_array<complex<double>, 3> *res) {
    assert(maxnA>=0);
    assert(maxnB>=0);
    
    int maxNk = maxnA + maxnB;
    
    (*res)[0][0][0] = 1.0;

    for(int nA = 0; nA <= maxnA; nA++) {
      for(int nB = 0; nB <= maxnB; nB++) {
	for(int Nk = nA+nB+1; Nk<=maxNk; Nk++) {
	  (*res)[nA][nB][Nk] = 0.0;
	}
      }
    }

    for(int nAB = 1; nAB <= maxnA+maxnB; nAB++) {
      for(int nA = 0; nA <= min(maxnA, nAB); nA++) {
	int nB = nAB-nA;
	if(nB<=maxnB) {
	  for(int Nk = 0; Nk <= min(nAB, maxNk); Nk++) {
	    complex<double> v(0);
	    if(nA > 0) {
	      v += (wAB-RA) * (*res)[nA-1][nB][Nk];
	      if(Nk-1 >= 0) 
		v += 1.0/(2.0*gAB) * (*res)[nA-1][nB][Nk-1];
	      if(Nk+1 <= maxNk)
		v += (Nk+1.0) * (*res)[nA-1][nB][Nk+1];
	    } else if(nB > 0) {
	      v += (wAB-RB) * (*res)[nA][nB-1][Nk];
	      if(Nk-1>=0)
		v += 1.0/(2.0*gAB) * (*res)[nA][nB-1][Nk-1];
	      if(Nk+1<=maxNk)
		v += (Nk+1.0)   * (*res)[nA][nB-1][Nk+1];
	    }
	    (*res)[nA][nB][Nk] = v;
	  }
	}
      }
    }
  }
  complex<double> hermite_coef_d_0(complex<double> gAB, complex<double> wAB, double RA, double RB,
				   int nA, int nB, int Nk) {
    if(nA==0 && nB==0 && Nk==0)
      return 1.0;
    
    if(Nk<0 || Nk>nA+nB)
      return 0.0;

    if(nA>0) {
      auto res0 = hermite_coef_d_0(gAB, wAB, RA, RB, nA-1, nB, Nk-1);
      auto res1 = hermite_coef_d_0(gAB, wAB, RA, RB, nA-1, nB, Nk  );
      auto res2 = hermite_coef_d_0(gAB, wAB, RA, RB, nA-1, nB, Nk+1);
      return 1.0/(2.0*gAB)*res0 + (wAB-RA)*res1 + (Nk+1.0)*res2;
    } else {
      auto res0 = hermite_coef_d_0(gAB, wAB, RA, RB, nA, nB-1, Nk-1);
      auto res1 = hermite_coef_d_0(gAB, wAB, RA, RB, nA, nB-1, Nk  );
      auto res2 = hermite_coef_d_0(gAB, wAB, RA, RB, nA, nB-1, Nk+1);
      return 1.0/(2.0*gAB)*res0 + (wAB-RB)*res1 + (Nk+1.0)*res2;
    }    
  }

  /*
  void calc_matrix(PlaneWaveGto *basis, OpBufBasic *bra, OpBufBasic *ket, MatrixXcd *res) {
    for(int A = 0; A < basis->size(); A++) {
      for(int B = 0; B < basis->size(); B++) {	
	complex<double> cumsum(0);
	const auto& nAs(bra->ns_[A]); // ncs_bra(bra->ncs_[A]);
	const auto& nBs(ket->ns_[B]); //const auto& ncs_ket(ket->ncs_[B]);
	const auto& cAs(bra->cs_[A]);
	const auto& cBs(ket->cs_[B]);
	for(int i = 0; i < bra->nums_[A]; i++) {
	  for(int j = 0; j < ket->nums_[B]; j++) {	    
	    cumsum += conj(cAs[i]) * cBs[j] * basis->getd(A,B,nAs[i], nBs[j], 0);
	  }
	}
	(*res)(A,B) = cumsum * basis->eAB_(A,B) * basis->hAB_(A,B);
      }
    }

  }
  void calc_matrix(PlaneWaveGto *basis, OpBufBasic *bra, OpBufGausspot *ket, MatrixXcd *res) {
    auto opV = ket->op_;
    for(int A = 0; A < basis->size(); A++) {
      for(int B = 0; B < basis->size(); B++) {
	const auto& nAs(bra->ns_[A]);
	const auto& cAs(bra->cs_[A]);
	int nB = basis->ns_[B];
	complex<double> cB = basis->Ns_[B];

	int max_nA_nB = nAs.maxCoeff() + nB;
	VectorXcd intg(max_nA_nB+1);
	dwn_gaussint_shift(max_nA_nB, opV->b_, basis->gAB_(A,B), basis->RAB_(A,B), &intg);
	
	complex<double> cumsum(0.0);
	for(int i = 0; i < bra->nums_[A]; i++) {
	  int nA = nAs[i];
	  complex<double> cumsum1(0.0);
	  for(int N = 0; N <= max_nA_nB; N++)
	    cumsum1 += basis->getd(A,B,nA,nB,N) * intg(N);
	  cumsum += cumsum1 * conj(cAs[i]) * cB;
	}
	(*res)(A,B) = cumsum * basis->eAB_(A,B) * opV->v0_;
      }
    }
  }
  void calc_at(PlaneWaveGto *basis, OpBufBasic *op, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res) {
    complex<double> ii(0.0, 1.0);
    for (int ix = 0; ix < xs.size(); ix++) {
      complex<double> cumsum(0.0);      
      for (int A = 0; A < basis->size(); A++) {
	double d = xs[ix] - basis->Rs_[A];
	complex<double> cumsum1(0.0);
	for(int i = 0; i < op->nums_[A]; i++) {
	  cumsum1 += op->cs_[A][i] * pow(d, op->ns_[A][i]) * exp(-basis->gs_[A]*d*d - ii*basis->Ps_[A]*d);
	}
	cumsum += cs[A]*cumsum1;
      }
      (*res)(ix) = cumsum;
    }    
  }
  */
  
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


  /*
  class MatrixDispatcher {
  private:
    typedef tuple<type_info, type_info, type_info> TTI;
    map<TTI, void (func*)(OpBuf*, OpBuf*, MatrixXcd*)
  public:
    MatrixDispatcher() {
      
    }
  };
  */
}

