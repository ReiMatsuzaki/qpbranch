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
  
  OpBasis::OpBasis(int num) : num_(num), cs_(num), ns_(num) {}

  /** utils functions */
  /*
  int num_op_basis(Operator op, int n) {
    
    if(op==kOp0 || op==kOp1 || op==kOp2) {
      return 1;
    } else if(op==kOpP1) {
      if(n==0)
	return 2;
      else
	return 3;
    } else if(op==kOpP2) {
      if(n==0)
	return 3;
      else if(n==1)
	return 4;
      else
	return 5;
    } else if(op==kOpdR) {
      if(n==0)
	return 2;
      else
	return 3;
    } else if(op==kOpdP) {
      return 1;
    } else if(op==kOpdgr) {
      return 2;
    } else if(op==kOpdgi) {
      return 1;
    } else {
      assert(false||"unsupported");
    }
  }

  int maxn_op_basis(Operator op, int n) {
    if(op==kOp0) {
      return n;
    } else if(op==kOp1) {
      return n+1;
    } else if(op==kOp2) {
      return n+2;
    } else if(op==kOpP1) {
      return n+1;
    } else if(op==kOpP2) {
      return n+2;
    } else if(op==kOpdR) {
      return n+1;
    } else if(op==kOpdP) {
      return n+1;
    } else if(op==kOpdgr) {
      return n+2;
    } else if(op==kOpdgi) {
      return n+2;
    } else {
      assert(false||"unsupported");
    }
  }
  */
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

  PlaneWaveGto::PlaneWaveGto(const VectorXi& ns, const vector<Operator*>& ops):
    num_(ns.size()), nop_(ops.size()), ns_(ns), gs_(num_), Rs_(num_), Ps_(num_),
    ops_(ops), Ns_(num_), maxn_(num_),
    gAB_(num_,num_), eAB_(num_,num_), hAB_(num_,num_), RAB_(num_,num_) {

    for(auto it = ops.begin(); it!=ops.end(); ++it) {
      (*it)->call_new(this);
      //	int npoly = num_op_basis(*it, ns_(A));
      //	OpBasis *ptr = new OpBasis(npoly);
      //	op_basis_[*it].push_back(ptr);
    }
    for(int A = 0; A < num_; A++) {
      int maxn = 0;
      for(auto it = op_basis_.begin(); it!=op_basis_.end(); ++it) {
	int maxn0 = (it->second)[A]->ns_.maxCoeff();
	if(maxn < maxn0)
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
    for(auto it = ops_.begin(); it != ops_.end(); ++it) {
      (*it)->call_setup(this);
    }
    
  }
  void PlaneWaveGto::matrix(Operator *ibra, Operator *iket, MatrixXcd *res) {
    /* compute overlap matrix. */


    assert(res->rows()>=num_);
    assert(res->cols()>=num_);
    assert(op_basis_.find(ibra)!=op_basis_.end());
    assert(op_basis_.find(iket)!=op_basis_.end());
    
    for(int A = 0; A < num_; A++) {
      for(int B = 0; B < num_; B++) {	
	auto *bra = op_basis_[ibra][A];
	auto *ket = op_basis_[iket][B];
	complex<double> cumsum(0);
	for(int i = 0; i < bra->num_; i++) {
	  for(int j = 0; j < ket->num_; j++) {
	    cumsum += conj(bra->cs_(i))*ket->cs_(j)*getd(A,B,bra->ns_(i),ket->ns_(j),0);
	  }
	}
	(*res)(A,B) = cumsum * eAB_(A,B) * hAB_(A,B);
      }
    }
    
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
    for(int A = 0; A < num_; A++) {
      auto ptr = new OpBasis(1);
      ptr->ns_[0] = ns_[A];
      ptr->cs_[0] = 1.0;
      op_basis_[op].push_back(ptr);      
    }    
  }
  void PlaneWaveGto::new_op(OperatorRn *op) {
    for(int A = 0; A < num_; A++) {
      auto ptr = new OpBasis(1);
      ptr->ns_[0] = ns_[A]+op->n();
      ptr->cs_[0] = 1.0;
      op_basis_[op].push_back(ptr);      
    }    
  }
  void PlaneWaveGto::new_op(OperatorPn *op) {
    for(int A = 0; A < num_; A++) {
      op_basis_[op].push_back(new OpBasis(3));
    }    
  }
  void PlaneWaveGto::setup_op(OperatorId *op) {
    for(int A = 0; A<num_;A++) {
      op_basis_[op][A]->cs_[0] = Ns_[A];
    }
  }
  void PlaneWaveGto::setup_op(OperatorRn *op) {
    for(int A = 0; A<num_;A++) {
      op_basis_[op][A]->cs_[0] = Ns_[A];
    }
  }
  void PlaneWaveGto::setup_op(OperatorPn *op) {
    for(int A = 0; A<num_;A++) {
      op_basis_[op][A]->cs_ = Ns_;
    }    
  }  

  /*
  void PlaneWaveGto::setup_operator() {
    complex<double> ii(0, 1);
    for(int A = 0; A < num_; A++) {
      int nA = ns_[A];
      complex<double> Nt = Ns_[A];
      complex<double> gA = gs_[A];
      double pA = Ps_[A];
      for(auto it = ops_.begin(); it != ops_.end(); ++it) {
	OpBasis *ptr = op_basis_[*it][A];
	if(*it==kOp0) {
	  ptr->cs_[0] = Nt;
	  ptr->ns_[0] = nA;
	} else if(*it==kOp1) {
	  ptr->cs_[0] = Nt;
	  ptr->ns_[0] = nA+1;
	} else if(*it==kOp2) {
	  ptr->cs_[0] = Nt;
	  ptr->ns_[0] = nA+2;
	} else if(*it==kOpP1) {
	  ptr->cs_[0] = Nt;
	  ptr->ns_[0] = nA;
	  ptr->cs_[1] = Nt;
	  ptr->ns_[1] = nA*Ps_[A];
	  if(nA!=0) {
	    ptr->cs_[2] = Nt*ii*(1.0*nA);
	    ptr->ns_[2] = nA;
	  }
	} else if(*it==kOpP2) {
	  ptr->ns_(0) = nA+2;
	  ptr->cs_(0) = Nt*(-4.0*gA*gA);
	  ptr->ns_(1) = nA+1;
	  ptr->cs_(1) = Nt*(4.0*ii*gA*pA);
	  ptr->ns_(2) = nA;
	  ptr->cs_(2) = Nt*(2.0*gA*(nA+1.0) + 2.0*nA*gA + pA*pA);
	  if(nA>0) {
	    ptr->ns_(3) = nA-1;
	    ptr->cs_(3) = Nt*(-2.0*nA*ii*pA);
	  }
	  if(nA>1) {
	    ptr->ns_(4) = nA-2;
	    ptr->cs_(4) = Nt*(-nA*(nA-1.0));
	  }
	} else if(*it==kOpdR) {
	  ptr->cs_(0) = Nt*2.0*gA;
	  ptr->ns_(0) = nA+1;
	  ptr->cs_(1) = Nt*(-ii*pA);
	  ptr->ns_(1) = nA;
	  if(nA>0) {
	    ptr->cs_(2) = Nt*(-nA*1.0);
	    ptr->ns_(2) = nA-1;
	  }
	} else if(*it==kOpdP) {
	  ptr->cs_(0) = Nt*ii;
	  ptr->ns_(0) = nA+1;	  
	} else if(*it==kOpdgr) {
	  ptr->cs_(0) = -Nt;
	  ptr->ns_(0) = nA+2;
	  ptr->cs_(1) = calc_nterm(1, nA, real(gA));
	  ptr->ns_(1) = nA;
	} else if(*it==kOpdgi) {
	  ptr->cs_(0) = Nt*(-ii);
	  ptr->ns_(0) = nA+2;	  
	} else {
	  assert(false&&"unsupported operator");	  
	}
      }
    }    
  }
  */
  /*
  void PlaneWaveGtoMDR::gausspot(Operator ibra, Operator iket, complex<double> b, MatrixXcd *ptr_res) {
     compute matrix element of gauss type potential
       .   res(A,B) = <gA|V|gB>,
       where
       .   V = exp[-bq^2].
       
    
    MatrixXcd& res(*ptr_res);
    assert(res.rows()>=num_);
    assert(res.cols()>=num_);
    
    for(int A = 0; A < num_; A++) {
      for(int B = 0; B < num_; B++) {		
	OpBasis *bra = op_basis_[ibra][A];
	OpBasis *ket = op_basis_[iket][B];

	int nA_nB = bra->ns_.maxCoeff() + ket->ns_.maxCoeff();

	VectorXcd intg(nA_nB+1);
	dwn_gaussint_shift(nA_nB, b, gAB_(A,B), RAB_(A,B), &intg);

	complex<double> cumsum(0);
	for(int i = 0; i < bra->num_; i++) {
	  for(int j = 0; j < ket->num_; j++) {
	    complex<double> cumsum1(0);
	    int nA(bra->ns_[i]);
	    int nB(ket->ns_[j]);
	    for(int N = 0; N <= nA_nB; N++) 
	      cumsum1 += getd(A,B,nA,nB,N) * intg(N);
	    cumsum += cumsum1 * conj(bra->cs_(i)) * ket->cs_(j);
	  }
	}
	res(A,B) = cumsum*eAB_(A,B);
      }
    }    
  }
*/

  
  /** same 
  PlaneWaveGto1Center::PlaneWaveGto1Center(const VectorXi& ns, const vector<Operator>& ops) :
    PlaneWaveGto(ns, ops) {
    maxmaxn_ = maxn_.maxCoeff();
    ints_ = VectorXcd::Zero(maxmaxn_+1);
  }
  PlaneWaveGto1Center::~PlaneWaveGto1Center() {}
  void PlaneWaveGto1Center::setup() {
    setup_normalize();
    setup_operator();    
    gtoint2n(maxmaxn_, conj(gs_(0))+gs_(0), &ints_);
  }
  void PlaneWaveGto1Center::overlap(Operator ibra, Operator iket, MatrixXcd *res){

    assert(res->rows()>=num_);
    assert(res->cols()>=num_);
    assert(op_basis_.find(ibra)!=op_basis_.end());
    assert(op_basis_.find(iket)!=op_basis_.end());

    for(int A = 0; A < num_; A++) {
      for(int B = 0; B < num_; B++) {	
	auto *bra = op_basis_[ibra][A];
	auto *ket = op_basis_[iket][B];
	complex<double> cumsum(0);
	for(int i = 0; i < bra->num_; i++) {
	  for(int j = 0; j < ket->num_; j++) {
	    cumsum += conj(bra->cs_(i))*ket->cs_(j)*ints_(bra->ns_(i)+ket->ns_(j));
	  }
	}
	(*res)(A,B) = cumsum;
      }
    }    
    
  }
  void PlaneWaveGto1Center::at(Operator iop, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res) {
    at_slow(iop, cs, xs, res);
  }
  */  
}

