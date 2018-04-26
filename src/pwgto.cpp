#include "pwgto.hpp"
#include "mathplus.hpp"
using namespace std;
  using boost::extents;
  using boost::array;


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

  double calc_nterm(int nd, int nA, double gA) {
    /*
      give difference of normalization term
       (d/d(Re[g]))^(nd)N(t)
     where
       N(t) = 1/sqrt(S)
       S = Int_{-oo}^{+oo} x^{2n} exp[-2Re[g]x^2] dx.
    overlap S can be expressed by Gamma function
       S = (2Re[g])^{-n-1/2} Gamma[n+1/2]
     So, normalization term N becomes
       N = (2Re[g])^{n/2+1/4} Sqrt(1/Gamma[n+1/2])
    */
    assert(false||"not impl");
    return 0.0;
  }

  int PwgtoBufferBasic::maxn(int A) {
    int maxn(0);
    for(auto it = ncs_[A].begin(); it != ncs_[A].end(); ++it) {
      if(maxn < it->first)
	maxn = it->first;
    }
    return maxn;
  }  
  void PwgtoBufferBasic::matrix(PwgtoBuffer *opbra, PlaneWaveGto *basis, MatrixXcd *res) {
    opbra->matrix(this, basis, res);
  }
  void PwgtoBufferBasic::matrix(PwgtoBufferBasic *opket, PlaneWaveGto *basis, MatrixXcd *res) {
    basis->calc_matrix(this, opket, res);
  }
  void PwgtoBufferBasic::matrix(PwgtoBufferGausspot *opket, PlaneWaveGto *basis, MatrixXcd *res) {
    basis->calc_matrix(this, opket, res);
  }

  void PwgtoBufferId::call_setup(PlaneWaveGto *basis) {
    basis->setup_op(this);
  }
  PwgtoBufferRn::PwgtoBufferRn(OperatorRn *op) : op_(op) {}
  void PwgtoBufferRn::call_setup(PlaneWaveGto *basis) {
    basis->setup_op(this);
  }
  PwgtoBufferPn::PwgtoBufferPn(OperatorPn *op) : op_(op) {}
  void PwgtoBufferPn::call_setup(PlaneWaveGto *basis) {
    basis->setup_op(this);
  }
  PwgtoBufferDa::PwgtoBufferDa(OperatorDa *op) : op_(op) {}
  void PwgtoBufferDa::call_setup(PlaneWaveGto *basis) {
    basis->setup_op(this);
  }

  PwgtoBufferGausspot::PwgtoBufferGausspot(OperatorGausspot *op, const VectorXi& ns) {
    op_ = op;
    ns_ = ns;
  }
  int PwgtoBufferGausspot::maxn(int A) { return ns_(A); }
  void PwgtoBufferGausspot::call_setup(PlaneWaveGto *basis) {
    basis->setup_op(this);
  }
  void PwgtoBufferGausspot::matrix(PwgtoBuffer *opbra, PlaneWaveGto *basis, MatrixXcd *res) {
    opbra->matrix(this, basis, res);
  }
  void PwgtoBufferGausspot::matrix(PwgtoBufferBasic *opket, PlaneWaveGto *basis, MatrixXcd *res) {
    assert(false&&"not impl");
  }
  void PwgtoBufferGausspot::matrix(PwgtoBufferGausspot *opket, PlaneWaveGto *basis, MatrixXcd *res) {
    assert(false&&"not impl");
  }
  
  
  
  PlaneWaveGto::PlaneWaveGto(const VectorXi& ns, const vector<Operator*>& ops):
    num_(ns.size()), nop_(ops.size()), ns_(ns), gs_(num_), Rs_(num_), Ps_(num_),
    ops_(ops), Ns_(num_), maxn_(num_),
    gAB_(num_,num_), eAB_(num_,num_), hAB_(num_,num_), RAB_(num_,num_) {

    for(auto it = ops.begin(); it!=ops.end(); ++it) 
      (*it)->call_new(this);

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
      it->second->call_setup(this);
    }
    
  }
  void PlaneWaveGto::matrix(Operator *ibra, Operator *iket, MatrixXcd *res) {
    /* compute overlap matrix. */

    assert(res->rows()>=num_);
    assert(res->cols()>=num_);
    assert(buffer_map_.find(ibra)!=buffer_map_.end());
    assert(buffer_map_.find(iket)!=buffer_map_.end());

    //    buffer_map_[ibra]->matrix(buffer_map_[iket], this, res);
    buffer_map_[iket]->matrix(buffer_map_[ibra], this, res);    
    
  }  
  void PlaneWaveGto::at(Operator *iop, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res) {
    for (int ix = 0; ix < xs.size(); ix++) {
      complex<double> cumsum(0.0);
      complex<double> ii(0.0, 1.0);
      for (int A = 0; A < num_; A++) {
	double d = xs[ix] - Rs_[A];
	cumsum += cs[A]*pow(d, ns_[A]) * exp(-gs_[A]*d*d - ii*Ps_[A]*d);
      }
      (*res)(ix) = cumsum;
    }
  }
  
  void PlaneWaveGto::new_op(OperatorId *op) {
    auto ptr = new PwgtoBufferId();
    for(int A = 0; A < num_; A++) {
      vector<pair<int,complex<double>>> ncs(1);
      ncs[0] = make_pair(ns_[A], 1.0);
      ptr->ncs_.push_back(ncs);
    }
    buffer_map_[op] = ptr;
  }
  void PlaneWaveGto::new_op(OperatorRn *op) {
    auto ptr = new PwgtoBufferRn(op);
    for(int A = 0; A < num_; A++) {
      vector<pair<int,complex<double>>> ncs(1);
      ncs[0] = make_pair(ns_[A]+op->n(), 1.0);
      ptr->ncs_.push_back(ncs);
    }
    buffer_map_[op] = ptr;
  }
  void PlaneWaveGto::new_op(OperatorPn *op) {
    assert(op->n()==1 || op->n()==2);
    auto ptr = new PwgtoBufferPn(op);
    for(int A = 0; A < num_; A++) {
      int nA = ns_(A);
      vector<pair<int,complex<double> > > ncs;
      if(op->n()==1) {
	ncs.push_back(make_pair(nA+1, 1.0));
	ncs.push_back(make_pair(nA,   1.0));
	if(nA>0)
	  ncs.push_back(make_pair(nA-1, 1.0));
      } else if(op->n()==2) {
	ncs.push_back(make_pair(nA+2, 1.0));
	ncs.push_back(make_pair(nA+1, 1.0));
	ncs.push_back(make_pair(nA,   1.0));
	if(nA>0)
	  ncs.push_back(make_pair(nA-1, 1.0));
	if(nA>1)
	  ncs.push_back(make_pair(nA-2, 1.0));
      }
      ptr->ncs_.push_back(ncs);
    }
    buffer_map_[op] = ptr;
  }
  void PlaneWaveGto::new_op(OperatorDa *op) {
    auto ptr = new PwgtoBufferDa(op);
    for(int A = 0; A < num_; A++) {
      int nA = ns_(A);
      vector<pair<int,complex<double>>> ncs;
      switch(op->id_) {
      case kIdDR:
	ncs.push_back(make_pair(nA+1, 1.0));
	ncs.push_back(make_pair(nA,   1.0));
	if(nA>0)
	  ncs.push_back(make_pair(nA-1, 1.0));
	break;
      case kIdDP:
	ncs.push_back(make_pair(nA+1, 1.0));
	break;
      case kIdDgr:
	ncs.push_back(make_pair(nA+2, 1.0));
	ncs.push_back(make_pair(nA+0, 1.0));
	break;
      case kIdDgi:
	ncs.push_back(make_pair(nA+2, 1.0));
	break;
      default:
	assert(false||"invalid id");
      }
      ptr->ncs_.push_back(ncs);
    }
    buffer_map_[op] = ptr;
  }
  void PlaneWaveGto::new_op(OperatorGausspot *op) {
    auto ptr = new PwgtoBufferGausspot(op, ns_);
    buffer_map_[op] = ptr;
  }
  void PlaneWaveGto::setup_op(PwgtoBufferId *op) {
    for(int A = 0; A<num_;A++) {
      op->ncs_[A][0].second = Ns_[A];
    }
  }
  void PlaneWaveGto::setup_op(PwgtoBufferRn *op) {
    for(int A = 0; A<num_;A++) {
      op->ncs_[A][0].second = Ns_[A];
    }
  }
  void PlaneWaveGto::setup_op(PwgtoBufferPn *op) {
    complex<double> ii(0.0, 1.0);
    if(op->op_->n()==1) {
      for(int A = 0; A < num_; A++) {
	int nA = ns_[A];
	complex<double> c = -ii*Ns_[A];
	double pA = Ps_[A];
	complex<double> gA = gs_[A];
	op->ncs_[A][0].second = c * (-2.0*gA);
	op->ncs_[A][1].second = c * (ii*pA);
	if(nA>0)
	  op->ncs_[A][2].second = c * (1.0*nA);
      }
    } else if(op->op_->n()==2) {
      for(int A = 0; A < num_; A++) {
	int nA = ns_[A];
	complex<double> c = -Ns_[A];
	double pA = Ps_[A];
	complex<double> gA = gs_[A];
	op->ncs_[A][0].second = c * 4.0*gA*gA;
	op->ncs_[A][1].second = c * (-4.0*ii*gA*pA);
	op->ncs_[A][2].second = c * (-2.0*gA*(nA+1.0) -2.0*nA*gA -pA*pA);
	if(nA>0)
	  op->ncs_[A][3].second = c * (2.0*nA*ii*pA);
	if(nA>1)
	  op->ncs_[A][4].second = c * (1.0*nA*(nA-1));
      }
    }
  }
  void PlaneWaveGto::setup_op(PwgtoBufferDa *buf) {
    complex<double> ii(0.0, 1.0);
    for(int A = 0; A < num_; A++) {
      int nA = ns_[A];      
      double pA = Ps_[A];
      complex<double> gA = gs_[A];
      double NA = Ns_[A];
      complex<double> c = -Ns_[A];
      switch(buf->op_->id_) {
      case kIdDR: // dR	
	buf->ncs_[A][0].second = c * (-2.0*gA);
	buf->ncs_[A][1].second = c * (ii*pA);
	if(nA>0)
	  buf->ncs_[A][2].second = c * (1.0*nA);
	break;
      case kIdDP: // dP
	buf->ncs_[A][0].second =  NA * ii;
	break;
      case kIdDgr: // dgr
	buf->ncs_[A][0].second =  -NA;
	buf->ncs_[A][1].second =  NA * calc_nterm(1, nA, real(gA));
	break;
      case kIdDgi: // dgi
	buf->ncs_[A][0].second =  -ii*NA;
	break;
      default:
	assert(false||"invalid id");
      }
    }
  }
  void PlaneWaveGto::setup_op(PwgtoBufferGausspot *op) {}
  void PlaneWaveGto::calc_matrix(PwgtoBufferBasic *bra, PwgtoBufferBasic *ket, MatrixXcd *res) {
    
    for(int A = 0; A < num_; A++) {
      for(int B = 0; B < num_; B++) {	
	complex<double> cumsum(0);
	const auto& ncs_bra(bra->ncs_[A]);
	const auto& ncs_ket(ket->ncs_[B]);
	for(auto it = ncs_bra.begin(); it != ncs_bra.end(); ++it)
	  for(auto jt = ncs_ket.begin(); jt != ncs_ket.end(); ++jt) 
	    cumsum += conj(it->second) * jt->second * getd(A,B,it->first, jt->first, 0);
	(*res)(A,B) = cumsum * eAB_(A,B) * hAB_(A,B);
      }
    }

  }
  void PlaneWaveGto::calc_matrix(PwgtoBufferBasic *ibra, PwgtoBufferGausspot *iket, MatrixXcd *res) {
    auto opV = iket->op_;
    for(int A = 0; A < num_; A++) {
      for(int B = 0; B < num_; B++) {
	int nB = ns_[B];
	int nA_nB = 0;
	vector<pair<int,complex<double> > >& ncs(ibra->ncs_[A]);
	for(auto it = ncs.begin(); it != ncs.end(); ++it) 
	  nA_nB = max(nA_nB, it->first);
	nA_nB += nB;

	VectorXcd intg(nA_nB+1);
	dwn_gaussint_shift(nA_nB, opV->b_, gAB_(A,B), RAB_(A,B), &intg);
	complex<double> cumsum(0.0);
	for(auto it = ncs.begin(); it != ncs.end(); ++it) {
	  int nA = it->first;
	  complex<double> cA = it->second;
	  complex<double> cumsum1(0.0);
	  for(int N = 0; N <= nA_nB; N++) {
	    cumsum1 += getd(A,B,nA,nB,N) * intg(N);
	  }
	  cumsum += cumsum1 * conj(cA) * Ns_[B];
	}
	(*res)(A,B) = cumsum * eAB_(A,B) * opV->v0_;
      }
    }
  }
}

