#include <iostream>

#include <qpbranch/pwgto.hpp>
#include <qpbranch/pwgto1c.hpp>
#include <qpbranch/pwgto_buf.hpp>
#include <qpbranch/mathplus.hpp>
using namespace std;

namespace qpbranch {

  OpBuf* make_op_buff(const VectorXi& ns, Operator *op) {
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
    return 0.0*nd*nA*gA;
  }

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

  

  OpBufBasic::OpBufBasic(int num) : num_(num), nums_(num_), ns_(num_), cs_(num_){}
  int OpBufBasic::maxn(int A) {
    assert(A>=0 && A<num_);
    return ns_[A].maxCoeff();    
  }
  void OpBufBasic::init_zero(int A, int num) {
    nums_[A] = num;
    ns_[A] = VectorXi::Zero(num);
    cs_[A] = VectorXcd::Zero(num);
  }
  void OpBufBasic::matrix(OpBuf *opbra, PlaneWaveGto *basis, MatrixXcd *res) {
    opbra->matrix(this, basis, res);
  }
  void OpBufBasic::matrix(OpBufBasic *ket, PlaneWaveGto *basis, MatrixXcd *res) {
    
    for(int A = 0; A < basis->num_; A++) {
      for(int B = 0; B < basis->num_; B++) {	
	complex<double> cumsum(0);
	const auto& nAs(this->ns_[A]); // ncs_bra(bra->ncs_[A]);
	const auto& nBs(ket->ns_[B]); //const auto& ncs_ket(ket->ncs_[B]);
	const auto& cAs(this->cs_[A]);
	const auto& cBs(ket->cs_[B]);
	for(int i = 0; i < this->nums_[A]; i++) {
	  for(int j = 0; j < ket->nums_[B]; j++) {	    
	    cumsum += conj(cAs[i]) * cBs[j] * basis->getd(A,B,nAs[i], nBs[j], 0);
	  }
	}
	(*res)(A,B) = cumsum * basis->eAB_(A,B) * basis->hAB_(A,B);
      }
    }

  }
  void OpBufBasic::matrix(OpBufGausspot *ket, PlaneWaveGto *basis, MatrixXcd *res) {
    
    auto opV = ket->op_;
    for(int A = 0; A < basis->num_; A++) {
      for(int B = 0; B < basis->num_; B++) {
	const auto& nAs(this->ns_[A]);
	const auto& cAs(this->cs_[A]);
	int nB = basis->ns_[B];
	complex<double> cB = basis->Ns_[B];

	int max_nA_nB = nAs.maxCoeff() + nB;
	VectorXcd intg(max_nA_nB+1);
	dwn_gaussint_shift(max_nA_nB, opV->b_, basis->gAB_(A,B), basis->RAB_(A,B), &intg);
	
	complex<double> cumsum(0.0);
	for(int i = 0; i < this->nums_[A]; i++) {
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
  void OpBufBasic::at(PlaneWaveGto *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res) {
    complex<double> ii(0.0, 1.0);
    for (int ix = 0; ix < xs.size(); ix++) {
      complex<double> cumsum(0.0);      
      for (int A = 0; A < num_; A++) {
	double d = xs[ix] - basis->Rs_[A];
	complex<double> cumsum1(0.0);
	for(int i = 0; i < nums_[A]; i++) {
	  cumsum1 += cs_[A][i] * pow(d, ns_[A][i]) * exp(-basis->gs_[A]*d*d - ii*basis->Ps_[A]*d);
	}
	cumsum += cs[A]*cumsum1;
      }
      (*res)(ix) = cumsum;
    }    
  }
  void OpBufBasic::matrix(OpBuf *opbra, Pwgto1c *basis, MatrixXcd *res) {
    opbra->matrix(this, basis, res);
  }
  void OpBufBasic::matrix(OpBufBasic *ket, Pwgto1c *basis, MatrixXcd *res) {
    
    for(int A = 0; A < basis->num_; A++) {
      for(int B = 0; B < basis->num_; B++) {	
	complex<double> cumsum(0);
	const auto& nAs(this->ns_[A]); // ncs_bra(bra->ncs_[A]);
	const auto& nBs(ket->ns_[B]); //const auto& ncs_ket(ket->ncs_[B]);
	const auto& cAs(this->cs_[A]);
	const auto& cBs(ket->cs_[B]);
	for(int i = 0; i < this->nums_[A]; i++) {
	  for(int j = 0; j < ket->nums_[B]; j++) {
	    int n = nAs[i] + nBs[j];
	    if(n % 2 == 0)
	      cumsum += conj(cAs[i]) * cBs[j] * basis->gints2n_[n/2];
	  }
	}
	(*res)(A,B) = cumsum;
      }
    }

  }
  void OpBufBasic::matrix(OpBufGausspot *ket, Pwgto1c *basis, MatrixXcd *res) {
    auto opV = ket->op_;
    complex<double> gg = conj(basis->g0_) + basis->g0_;
    complex<double> dq = basis->R0_ - opV->q0_;
    for(int A = 0; A < basis->num_; A++) {
      for(int B = 0; B < basis->num_; B++) {
	const auto& nAs(this->ns_[A]);
	const auto& cAs(this->cs_[A]);
	int nB = basis->ns_[B];
	complex<double> cB = basis->Ns_[B];
	int max_nA_nB = nAs.maxCoeff() + nB;
	VectorXcd intg(max_nA_nB+1);
	dwn_gaussint_shift(max_nA_nB, opV->b_, gg, dq, &intg);	
	complex<double> cumsum(0); 
	for(int i = 0; i < this->nums_[A]; i++) {
	  int nA = nAs[i];
	  cumsum += conj(cAs[i]) * cB * intg(nA);
	}
	(*res)(A,B) = cumsum;
      }
    }
  }
  void OpBufBasic::at(Pwgto1c *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res) {
    complex<double> ii(0.0, 1.0);
    for (int ix = 0; ix < xs.size(); ix++) {
      complex<double> cumsum(0.0);      
      for (int A = 0; A < num_; A++) {
	double d = xs[ix] - basis->R0_;
	complex<double> cumsum1(0.0);
	for(int i = 0; i < nums_[A]; i++) {
	  cumsum1 += cs_[A][i] * pow(d, ns_[A][i]) * exp(-basis->g0_*d*d - ii*basis->P0_*d);
	}
	cumsum += cs[A]*cumsum1;
      }
      (*res)(ix) = cumsum;
    }    
  }

  
  OpBufId::OpBufId(const VectorXi& ns) : OpBufBasic(ns.size()) {
    for(int A = 0; A < num_; A++) {
      this->init_zero(A, 1);
      this->ns_[A][0] = ns[A];
    }
  }
  void OpBufId::setup(PlaneWaveGto *basis) {
    for(int A = 0; A < num_; A++) {
      ns_[A][0] = basis->ns_[A];
      cs_[A][0] = basis->Ns_[A];
    }
  }
  void OpBufId::setup(Pwgto1c *basis) {
    for(int A = 0; A < num_; A++) {
      ns_[A][0] = basis->ns_[A];
      cs_[A][0] = basis->Ns_[A];
    }
  }
  OpBufRn::OpBufRn(const VectorXi& ns, int n): OpBufBasic(ns.size()), n_(n) {
    for(int A = 0; A < num_; A++) {
      this->init_zero(A, 1);
      this->ns_[A][0] = ns[A] + n;
    }
  }
  void OpBufRn::setup(PlaneWaveGto *basis) {
    for(int A = 0; A < num_; A++) {
      ns_[A][0] = basis->ns_[A] + n_;
      cs_[A][0] = basis->Ns_[A];
    }
  }
  void OpBufRn::setup(Pwgto1c *basis) {
    for(int A = 0; A < num_; A++) {
      ns_[A][0] = basis->ns_[A] + n_;
      cs_[A][0] = basis->Ns_[A];
    }
  }  
  OpBufPn::OpBufPn(const VectorXi& ns, int n) : OpBufBasic(ns.size()), n_(n) {
    
    assert(n==1 || n==2);

    for(int A = 0; A < num_; A++) {
      int nA = ns(A);
      vector<pair<int,complex<double> > > ncs;
      if(n==1) {
	if(nA == 0) {
	  this->init_zero(A, 2);
	  ns_[A][0] = nA+1;
	  ns_[A][1] = nA;
	} else {
	  this->init_zero(A, 3);
	  ns_[A][2] = nA-1;
	}
      } else if(n==2) {
	if (nA == 0) {
	  this->init_zero(A, 3);
	  ns_[A][0] = nA+2;
	  ns_[A][1] = nA+1;
	  ns_[A][2] = nA;
	} else if(nA == 1) {
	  this->init_zero(A, 4);
	  ns_[A][0] = nA+2;
	  ns_[A][1] = nA+1;
	  ns_[A][2] = nA;
	  ns_[A][3] = nA-1;
	} else {
	  this->init_zero(A, 5);
	  ns_[A][0] = nA+2;
	  ns_[A][1] = nA+1;
	  ns_[A][2] = nA;
	  ns_[A][3] = nA-1;
	  ns_[A][4] = nA-2;
	}

      }
    }
  }
  void OpBufPn::setup(PlaneWaveGto *basis) {
    complex<double> ii(0.0, 1.0);
    for(int A = 0; A < num_; A++) {
      int nA = basis->ns_[A];      
      double pA = basis->Ps_[A];
      complex<double> gA = basis->gs_[A];
      if(n_ == 1) {
	complex<double> c = -ii*basis->Ns_[A];
	ns_[A][0] = nA+1;
	cs_[A][0] = c * (-2.0*gA);
	ns_[A][1] = nA;
	cs_[A][1] = c * (ii*pA);
	if(nA>0) {
	  ns_[A][2] = nA-1;
	  cs_[A][2] = c * (1.0*nA);
	}	
      } else if(n_ == 2) {
	complex<double> c = -basis->Ns_[A];
	ns_[A][0] = nA+2;
	cs_[A][0] = c * 4.0*gA*gA;
	ns_[A][1] = nA+1;
	cs_[A][1] = c * (-4.0*ii*gA*pA);
	ns_[A][2] = nA;
	cs_[A][2] = c * (-2.0*gA*(nA+1.0) -2.0*nA*gA -pA*pA);
	if(nA>0) {
	  ns_[A][3] = nA-1;
	  cs_[A][3] = c * (2.0*nA*ii*pA);
	}
	if(nA>1) {
	  ns_[A][4] = nA-2;
	  cs_[A][4] = c * (1.0*nA*(nA-1));
	}
      }
    }
  }
  void OpBufPn::setup(Pwgto1c *basis) {
    complex<double> ii(0.0, 1.0);
    for(int A = 0; A < num_; A++) {
      int nA = basis->ns_[A];      
      double pA = basis->P0_;
      complex<double> gA = basis->g0_;
      if(n_ == 1) {
	complex<double> c = -ii*basis->Ns_[A];
	ns_[A][0] = nA+1;
	cs_[A][0] = c * (-2.0*gA);
	ns_[A][1] = nA;
	cs_[A][1] = c * (ii*pA);
	if(nA>0) {
	  ns_[A][2] = nA-1;
	  cs_[A][2] = c * (1.0*nA);
	}	
      } else if(n_ == 2) {
	complex<double> c = -basis->Ns_[A];
	ns_[A][0] = nA+2;
	cs_[A][0] = c * 4.0*gA*gA;
	ns_[A][1] = nA+1;
	cs_[A][1] = c * (-4.0*ii*gA*pA);
	ns_[A][2] = nA;
	cs_[A][2] = c * (-2.0*gA*(nA+1.0) -2.0*nA*gA -pA*pA);
	if(nA>0) {
	  ns_[A][3] = nA-1;
	  cs_[A][3] = c * (2.0*nA*ii*pA);
	}
	if(nA>1) {
	  ns_[A][4] = nA-2;
	  cs_[A][4] = c * (1.0*nA*(nA-1));
	}
      }
    }
  }
  OpBufDa::OpBufDa(const VectorXi& ns, int id) : OpBufBasic(ns.size()), id_(id) {
    for(int A = 0; A < ns.size(); A++) {
      int nA = ns[A];
      switch(id) {
      case kIdDR:
	if(nA == 0)  {
	  this->init_zero(A, 2);
	  ns_[A][0] = nA+1;
	  ns_[A][1] = nA;
	} else {
	  this->init_zero(A, 3);
	  ns_[A][2] = nA-1;
	}
	break;
      case kIdDP:
	this->init_zero(A, 1);
	ns_[A][0] = nA+1;
	break;
      case kIdDgr:
	this->init_zero(A, 2);
	ns_[A][0] = nA+2;
	ns_[A][1] = nA;
	break;
      case kIdDgi:
	this->init_zero(A, 1);
	ns_[A][0] = nA+2;
	break;
      default:
	assert(false||"invalid id");
      }
    }
  }
  void OpBufDa::setup(PlaneWaveGto *basis) {
    complex<double> ii(0, 1);
    for(int A = 0; A < basis->num_; A++) {
      int nA = basis->ns_[A];
      double pA = basis->Ps_[A];
      complex<double> gA = basis->gs_[A];
      double NA = basis->Ns_[A];
      switch(id_) {
      case kIdDR:
	ns_[A][0] = nA+1;
	cs_[A][0] = (-NA) * (-2.0*gA);
	ns_[A][1] = nA;
	cs_[A][1] = (-NA) * (ii*pA);
	if(nA>0) {
	  ns_[A][2] = nA-1;
	  cs_[A][1] = (-NA) * (1.0*nA);
	}
	break;
      case kIdDP:
	ns_[A][0] = nA+1;
	cs_[A][0] = NA * ii;
	break;
      case kIdDgr:
	ns_[A][0] = nA+2;
	cs_[A][0] = -NA;
	ns_[A][1] = nA;
	cs_[A][1] = calc_nterm(1, nA, real(gA));
	break;
      case kIdDgi:
	ns_[A][0] = nA+2;
	cs_[A][0] = -ii*NA;
	break;
      default:
	assert(false||"invalid id");
      }
    }
  }
  void OpBufDa::setup(Pwgto1c *basis) {
    complex<double> ii(0, 1);
    for(int A = 0; A < basis->num_; A++) {
      int nA = basis->ns_[A];
      double pA = basis->P0_;
      complex<double> gA = basis->g0_;
      double NA = basis->Ns_[A];
      switch(id_) {
      case kIdDR:
	ns_[A][0] = nA+1;
	cs_[A][0] = (-NA) * (-2.0*gA);
	ns_[A][1] = nA;
	cs_[A][1] = (-NA) * (ii*pA);
	if(nA>0) {
	  ns_[A][2] = nA-1;
	  cs_[A][1] = (-NA) * (1.0*nA);
	}
	break;
      case kIdDP:
	ns_[A][0] = nA+1;
	cs_[A][0] = NA * ii;
	break;
      case kIdDgr:
	ns_[A][0] = nA+2;
	cs_[A][0] = -NA;
	ns_[A][1] = nA;
	cs_[A][1] = calc_nterm(1, nA, real(gA));
	break;
      case kIdDgi:
	ns_[A][0] = nA+2;
	cs_[A][0] = -ii*NA;
	break;
      default:
	assert(false||"invalid id");
      }
    }
  }

  OpBufGausspot::OpBufGausspot(const VectorXi& ns, OperatorGausspot *op) : num_(ns_.size()),
									   ns_(ns), op_(op) {}  
  int OpBufGausspot::maxn(int A) { return ns_(A); }
  void OpBufGausspot::setup(PlaneWaveGto*) {}
  void OpBufGausspot::matrix(OpBuf *opbra, PlaneWaveGto *basis, MatrixXcd *res) {
    opbra->matrix(this, basis, res);
  }
  void OpBufGausspot::matrix(OpBufBasic *, PlaneWaveGto *, MatrixXcd *) {
    assert(false&&"not impl");
  }
  void OpBufGausspot::matrix(OpBufGausspot *, PlaneWaveGto *, MatrixXcd *) {}
  void OpBufGausspot::at(PlaneWaveGto*, const VectorXcd&, const VectorXd&, VectorXcd *) {
    throw runtime_error("not impl");
  }
  void OpBufGausspot::setup(Pwgto1c *) {}
  void OpBufGausspot::matrix(OpBuf *opbra, Pwgto1c *basis, MatrixXcd *res) {
    opbra->matrix(this, basis, res);
  }
  void OpBufGausspot::matrix(OpBufBasic *, Pwgto1c *, MatrixXcd *) {
    assert(false&&"not impl");
  }
  void OpBufGausspot::matrix(OpBufGausspot *, Pwgto1c *, MatrixXcd *) {}
  void OpBufGausspot::at(Pwgto1c*, const VectorXcd&, const VectorXd&, VectorXcd *) {
    throw runtime_error("not impl");
  }
  
}
