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

  OpBasis::OpBasis(int np, int nb) : num_(np), cs_(np, nb), ns_(np, nb) {}
    
  GaussBasis::GaussBasis(const VectorXi& ns, const vector<Operator>& ops):
    num_(ns.size()), ns_(ns), gs_(num_), Rs_(num_), Ps_(num_), ops_(ops),
    op_basis_(ops.size()), Ns_(num_),
    gAB_(num_,num_), eAB_(num_,num_), hAB_(num_,num_), RAB_(num_,num_), maxn_(num_) {

    for(int A = 0; A < num_; A++) {
      int maxn = 0;
      for(auto it = ops.begin(); it != ops.end(); ++it) {
	int maxn0 = OpBasis::maxn(*it, ns_[A]);
	if(maxn<maxn0) maxn=maxn0;
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
  GaussBasis::~GaussBasis() {}
  void GaussBasis::setup() {
    setup_normalize();
    setup_combination();
    setup_operator();
  }
  void GaussBasis::overlap(Operator ibra, Operator iket, MatrixXcd *res) {
    /* compute overlap matrix. */
    assert(ibra==kOp0);
    assert(iket==kOp0);
   
    for(int A = 0; A < num_; A++) {
      for(int B = 0; B < num_; B++) {
	complex<double> cumsum(0);
	cumsum = Ns_(A)*Ns_(B)*eAB_(A,B)*hAB_(A,B)*getd(A,B,ns_[A],ns_[B],0);
	(*res)(A,B) = cumsum;
      }
    }
  }
  void GaussBasis::at(Operator iop, const VectorXcd& cs, const VectorXd& xs, MatrixXcd *res) {
    for (int ix = 0; ix < xs.size(); ix++) {
      complex<double> cumsum(0.0);
      complex<double> ii(0.0, 1.0);
      for (int A = 0; A < num_; A++) {
	double d = xs[ix] - Rs_[A];
	cumsum += cs[A]*pow(d, ns_[A]) * exp(-gs_[A]*d*d - ii*Ps_[A]*d);
      }
    }
  }
  void GaussBasis::setup_normalize() {

    for(int A = 0; A < num_; A++) {
      int n(ns_[A]);
      VectorXcd gg(n+1);
      gtoint2n(n, gs_[A]+conj(gs_[A]), &gg);
      Ns_[A] = real( 1.0/sqrt((gg[n])) );
    }
  }
  void GaussBasis::setup_combination() {

    auto ii = complex<double>(0,1);
    
    for(int A = 0; A < num_; A++) {
      for(int B = 0; B < num_; B++) {
	auto cgA = conj(gs_(A));
	auto  gB = gs_(B);
	auto RA = Rs_[A]; auto RB = Rs_[B]; auto PA = Ps_[A]; auto PB = Ps_[B];
	gAB_(A,B) = cgA + gB;
	RAB_(A,B) = (2.0*cgA*RA + 2.0*gB*RB - ii*PA + ii*PB) / (2.0*gAB_(A,B));
	eAB_(A,B) = exp(-cgA*RA*RA - gB*RB*RB
			+ii*RA*PA - ii*PB*RB
			+gAB_(A,B)*pow(RAB_(A,B), 2));
	hAB_(A,B) = sqrt(M_PI/gAB_(A,B));
	/*
	if(A+B==1) {
	  cerr << eAB_(A,B) << RAB_(A,B) << pow(RAB_(A,B), 2) << endl;
	}
	*/
	hermite_coef_d(gAB_(A,B), RAB_(A,B), Rs_[A], Rs_[B],
		       ns_[A], ns_[B], (*d_)[A][B]);
      }
    }    
  }
  void GaussBasis::setup_operator() {
    
    for(auto it = ops_.begin(); it != ops_.end(); ++it) {
      if(*it==kOp0) {
	OpBasis *ptr = new OpBasis(1, num_);
	ptr->cs_.row(0) = Ns_;
	ptr->ns_.row(0) = ns_;
	op_basis_[kOp0] = ptr;
      } else if(*it==kOp1) {
	OpBasis *ptr = new OpBasis(1, num_);
	ptr->cs_.row(0) = Ns_;
	ptr->ns_.row(0) = ns_ + VectorXi::Ones(num_);
	op_basis_[kOp0] = ptr;
      } else if(*it==kOp2) {
	OpBasis *ptr = new OpBasis(1, num_);
	ptr->cs_.row(0) = Ns_;
	ptr->ns_.row(0) = ns_ + 2*VectorXi::Ones(num_);
	op_basis_[kOp0] = ptr;
      } else {
	assert(false&&"unsupported operator");	  
      }
    }    
  }
  
  PlaneWaveGTO::PlaneWaveGTO(const VectorXi& ns, const vector<Operator>& ops):
    GaussBasis(ns, ops) {}
  PlaneWaveGTO::~PlaneWaveGTO() {
  }
  void PlaneWaveGTO::setup() { GaussBasis::setup(); }
  void PlaneWaveGTO::overlap(Operator ibra, Operator iket, MatrixXcd *res) {
    GaussBasis::overlap(ibra, iket, res);
  }
  void PlaneWaveGTO::at(Operator iop, const VectorXcd& cs, const VectorXd& xs, MatrixXcd *res) {
    GaussBasis::at(iop, cs, xs, res);
  }
}

