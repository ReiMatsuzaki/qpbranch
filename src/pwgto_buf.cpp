#include <iostream>

#include <qpbranch/pwgto.hpp>
#include <qpbranch/pwgto1c.hpp>
#include <qpbranch/pwgto_buf.hpp>
#include <qpbranch/mathplus.hpp>
using namespace std;

namespace qpbranch {

  OpBuf* MakeOpBuf(Pwgto *basis, Operator *op) {
    assert(op!=nullptr);
    const VectorXi& ns = basis->ns();
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
    } else if(optype == typeid(OperatorSpline)) {
      auto opop = static_cast<OperatorSpline*>(op);
      return new OpBufSpline(ns, basis->norder(), opop);
    } else if(optype == typeid(OperatorSplineP1)) {
      auto opop = static_cast<OperatorSplineP1*>(op);
      return new OpBufSplineP1(ns, basis->norder(), opop);      
    } else if(optype == typeid(OperatorGausspot)) {
      auto opop = static_cast<OperatorGausspot*>(op);
      return new OpBufGausspot(ns, opop);
    } else {
      cerr << "Unsupported operator" << endl;
      abort();
    }
  }
  OpBuf* MakeOpBuf(Pwgto1c *basis, Operator *op) {
    const VectorXi& ns = basis->ns();
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
    } else if(optype == typeid(OperatorSpline)) {
      auto opop = static_cast<OperatorSpline*>(op);
      return new OpBufSpline(ns, 0, opop);
    } else if(optype == typeid(OperatorGausspot)) {
      auto opop = static_cast<OperatorGausspot*>(op);
      return new OpBufGausspot(ns, opop);
    } else {
      cerr << "Unsupported operator" << endl;
      abort();
    }
  }
  
  // give difference of normalization term
  //     (d/d(Re[g]))^(nd)N(t)
  // where
  //     N(t) = 1/sqrt(S)
  //     S = Int_{-oo}^{+oo} x^{2n} exp[-2Re[g]x^2] dx.
  // overlap S can be expressed by Gamma function
  //     S = (2Re[g])^{-n-1/2} Gamma[n+1/2]
  // So, normalization term N becomes
  //     N = (2Re[g])^{n/2+1/4} Sqrt(1/Gamma[n+1/2])
  double Nterm(int nd, int nA, double gAr) {

    VectorXd gint(nA+1);
    GammaNhalf(nA, &gint);

    double nn = nA*0.5+0.25;
    double c = pow(2.0, nn) * sqrt(1/gint[nA]);
    switch(nd) {
    case 0:
      return pow(gAr, nn) * c;
      break;
    case 1:
      return pow(gAr, nn-1) * nn * c;
      break;
    default:
      throw runtime_error("unsupported nd");
    }
  }

  int OpBufBasic::Maxn(int A) {
    assert(A>=0 && A<num_);
    return ns_[A].maxCoeff();    
  }
  void OpBufBasic::InitZero(int A, int num) {

    nums_[A] = num;
    ns_[A] = VectorXi::Zero(num);
    cs_[A] = VectorXcd::Zero(num);
  }
  void OpBufBasic::Matrix(OpBuf *opbra, Pwgto *basis, MatrixXcd *res) {
    opbra->Matrix(this, basis, res);
  }
  void OpBufBasic::Matrix(OpBufBasic *ket, Pwgto *basis, MatrixXcd *res) {
    
    for(int A = 0; A < basis->num(); A++) {
      for(int B = 0; B < basis->num(); B++) {	
	complex<double> cumsum(0);
	const auto& nAs(this->ns_[A]); // ncs_bra(bra->ncs_[A]);
	const auto& nBs(ket->ns_[B]); //const auto& ncs_ket(ket->ncs_[B]);
	const auto& cAs(this->cs_[A]);
	const auto& cBs(ket->cs_[B]);
	for(int i = 0; i < this->nums_[A]; i++) {
	  for(int j = 0; j < ket->nums_[B]; j++) {	    
	    cumsum += conj(cAs[i]) * cBs[j] * basis->get_d(A,B,nAs[i], nBs[j], 0);
	  }
	}
	(*res)(A,B) = cumsum * basis->eAB_(A,B) * basis->hAB_(A,B);
      }
    }

  }
  void OpBufBasic::Matrix(OpBufGausspot *ket, Pwgto *basis, MatrixXcd *res) {
    
    auto opV = ket->op();
    for(int A = 0; A < basis->num(); A++) {
      for(int B = 0; B < basis->num(); B++) {
	const auto& nAs(this->ns_[A]);
	const auto& cAs(this->cs_[A]);
	int nB = basis->ns()[B];
	complex<double> cB = basis->Ns()[B];

	int max_nA_nB = nAs.maxCoeff() + nB;
	VectorXcd intg(max_nA_nB+1);
	DwIntGauss(max_nA_nB, opV->b(), basis->gAB_(A,B), basis->RAB_(A,B), &intg);
	
	complex<double> cumsum(0.0);
	for(int i = 0; i < this->nums_[A]; i++) {
	  int nA = nAs[i];
	  complex<double> cumsum1(0.0);
	  for(int N = 0; N <= max_nA_nB; N++)
	    cumsum1 += basis->get_d(A,B,nA,nB,N) * intg(N);
	  cumsum += cumsum1 * conj(cAs[i]) * cB;
	}
	(*res)(A,B) = cumsum * basis->eAB_(A,B) * opV->v0();
      }
    }
  }
  void OpBufBasic::At(Pwgto *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res) {
    complex<double> ii(0.0, 1.0);
    for (int ix = 0; ix < xs.size(); ix++) {
      complex<double> cumsum(0.0);      
      for (int A = 0; A < num_; A++) {
	double d = xs[ix] - basis->Rs()[A];
	complex<double> cumsum1(0.0);
	for(int i = 0; i < nums_[A]; i++) {
	  cumsum1 += cs_[A][i]*pow(d, ns_[A][i]) * exp(-basis->gs()[A]*d*d - ii*basis->Ps()[A]*d);
	}
	cumsum += cs[A]*cumsum1;
      }
      (*res)(ix) = cumsum;
    }
  }
  void OpBufBasic::Matrix(OpBuf *opbra, Pwgto1c *basis, MatrixXcd *res) {
    opbra->Matrix(this, basis, res);
  }
  void OpBufBasic::Matrix(OpBufBasic *ket, Pwgto1c *basis, MatrixXcd *res) {
    
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
  void OpBufBasic::Matrix(OpBufGausspot *ket, Pwgto1c *basis, MatrixXcd *res) {
    auto opV = ket->op();
    complex<double> gg = conj(basis->g0_) + basis->g0_;
    complex<double> dq = basis->R0_ - opV->q0();
    for(int A = 0; A < basis->num_; A++) {
      for(int B = 0; B < basis->num_; B++) {
	const auto& nAs(this->ns_[A]);
	const auto& cAs(this->cs_[A]);
	int nB = basis->ns_[B];
	complex<double> cB = basis->Ns_[B];
	int max_nA_nB = nAs.maxCoeff() + nB;

	VectorXcd intg(max_nA_nB+1);
	IntGtoShift(max_nA_nB, gg, opV->b(), dq, &intg);
	complex<double> cumsum(0); 
	for(int i = 0; i < this->nums_[A]; i++) {
	  int nA = nAs[i];
	  cumsum += conj(cAs[i]) * cB * intg(nA+nB);
	}

	
	/*
	VectorXcd intg(max_nA_nB+1);	
	DwIntGauss(max_nA_nB, gg, opV->b(), dq, &intg);	
	complex<double> cumsum(0); 
	for(int i = 0; i < this->nums_[A]; i++) {
	  int nA = nAs[i];
	  cumsum += conj(cAs[i]) * cB * intg(nA+nB);
	}
	*/
	(*res)(A,B) = cumsum * opV->v0();
      }
    }
  }
  void OpBufBasic::At(Pwgto1c *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res) {
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
  void OpBufBasic::Dump() const {
    using std::cout;
    using std::endl;
    cout << "OpBufBasic" << endl;
    cout << "num: " << num_ << endl;
    for(int A = 0; A < num_; A++) {
      cout << "A=" << A << endl;
      for(int i = 0; i < nums_[A]; i++) {
	cout << ns_[A][i] << " " << cs_[A][i] << endl;
      }
    }
  }
  
  OpBufId::OpBufId(const VectorXi& ns) : OpBufBasic(ns.size()) {
    for(int A = 0; A < num(); A++) {
      this->InitZero(A, 1);
      this->ns_[A][0] = ns[A];
    }
  }
  void OpBufId::SetUp(Pwgto *basis) {
    for(int A = 0; A < num_; A++) {
      ns_[A][0] = basis->ns()[A];
      cs_[A][0] = basis->Ns()[A];
    }
  }
  void OpBufId::SetUp(Pwgto1c *basis) {
    for(int A = 0; A < num_; A++) {
      ns_[A][0] = basis->ns_[A];
      cs_[A][0] = basis->Ns_[A];
    }
  }
  OpBufRn::OpBufRn(const VectorXi& ns, int n): OpBufBasic(ns.size()), n_(n) {
    for(int A = 0; A < num_; A++) {
      this->InitZero(A, 1);
      this->ns_[A][0] = ns[A] + n;
    }
  }
  void OpBufRn::SetUp(Pwgto *basis) {
    for(int A = 0; A < num_; A++) {
      ns_[A][0] = basis->ns()[A] + n_;
      cs_[A][0] = basis->Ns()[A];
    }
  }
  void OpBufRn::SetUp(Pwgto1c *basis) {
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
	  this->InitZero(A, 2);
	  ns_[A][0] = nA+1;
	  ns_[A][1] = nA;
	} else {
	  this->InitZero(A, 3);
	  ns_[A][2] = nA-1;
	}
      } else if(n==2) {
	if (nA == 0) {
	  this->InitZero(A, 3);
	  ns_[A][0] = nA+2;
	  ns_[A][1] = nA+1;
	  ns_[A][2] = nA;
	} else if(nA == 1) {
	  this->InitZero(A, 4);
	  ns_[A][0] = nA+2;
	  ns_[A][1] = nA+1;
	  ns_[A][2] = nA;
	  ns_[A][3] = nA-1;
	} else {
	  this->InitZero(A, 5);
	  ns_[A][0] = nA+2;
	  ns_[A][1] = nA+1;
	  ns_[A][2] = nA;
	  ns_[A][3] = nA-1;
	  ns_[A][4] = nA-2;
	}

      }
    }
  }
  void OpBufPn::SetUp(Pwgto *basis) {
    complex<double> ii(0.0, 1.0);
    for(int A = 0; A < num_; A++) {
      int nA = basis->ns()[A];
      double pA = basis->Ps()[A];
      complex<double> gA = basis->gs()[A];
      if(n_ == 1) {
	complex<double> c = -ii*basis->Ns()[A];
	ns_[A][0] = nA+1;
	cs_[A][0] = c * (-2.0*gA);
	ns_[A][1] = nA;
	cs_[A][1] = c * (ii*pA);
	if(nA>0) {
	  ns_[A][2] = nA-1;
	  cs_[A][2] = c * (1.0*nA);
	}	
      } else if(n_ == 2) {
	complex<double> c = -basis->Ns()[A];
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
  void OpBufPn::SetUp(Pwgto1c *basis) {
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
	  this->InitZero(A, 2);
	  ns_[A][0] = nA+1;
	  ns_[A][1] = nA;
	} else {
	  this->InitZero(A, 3);
	  ns_[A][0] = nA+1;
	  ns_[A][1] = nA;
	  ns_[A][2] = nA-1;
	}
	break;
      case kIdDP:
	this->InitZero(A, 1);
	ns_[A][0] = nA+1;
	break;
      case kIdDgr:
	this->InitZero(A, 2);
	ns_[A][0] = nA+2;
	ns_[A][1] = nA;
	break;
      case kIdDgi:
	this->InitZero(A, 1);
	ns_[A][0] = nA+2;
	break;
      default:
	assert(false||"invalid id");
      }
    }
  }
  void OpBufDa::SetUp(Pwgto *basis) {
    complex<double> ii(0, 1);
    for(int A = 0; A < basis->num(); A++) {
      int nA = basis->ns()[A];
      double pA = basis->Ps()[A];
      complex<double> gA = basis->gs()[A];
      double NA = basis->Ns()[A];
      switch(id_) {
      case kIdDR:
	ns_[A][0] = nA+1;
	cs_[A][0] = (-NA) * (-2.0*gA);
	ns_[A][1] = nA;
	cs_[A][1] = (-NA) * (ii*pA);
	if(nA>0) {
	  ns_[A][2] = nA-1;
	  cs_[A][2] = (-NA) * (1.0*nA);
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
	cs_[A][1] = Nterm(1, nA, real(gA));
	break;
      case kIdDgi:
	ns_[A][0] = nA+2;
	cs_[A][0] = -ii*NA;
	break;
      default:
	throw runtime_error("invalid id");	
      }
    }
  }
  void OpBufDa::SetUp(Pwgto1c *basis) {
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
	  cs_[A][2] = (-NA) * (1.0*nA);
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
	cs_[A][1] = Nterm(1, nA, real(gA));
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
  OpBufSpline::OpBufSpline(const VectorXi& ns, int norder, OperatorSpline *op): OpBufBasic(ns.size()), op_(op)  {

    assert(0<=norder && norder<=2);
    
    for(int A = 0; A < ns.size(); A++) {
      this->InitZero(A, norder+1);
      int nA = ns[A];
      for(int i = 0; i < norder+1; i++) 
	ns_[A][i] = nA + i;      
    }
  }
  void OpBufSpline::SetUp(Pwgto *basis) {
    int norder = basis->norder();
    assert(0<=norder && norder<=2);
    double dx = basis->dx();
    for(int A = 0; A < basis->num(); A++) {

      auto nA = basis->ns()[A];
      auto NA = basis->Ns()[A];

      // principle number
      for(int i = 0; i < basis->norder(); i++) {
	ns_[A][i] = nA+i;
      }

      // numerical derivative
      VectorXd xs(3);
      VectorXcd ys(3);
      xs(1) = basis->Rs()[A];
      xs(0) = xs(1)-dx;
      xs(2) = xs(1)+dx;
      op_->At(xs, &ys);
      cs_[A][0]   = NA * ys(1);
      if(norder>0)
	cs_[A][1] = NA * (ys(2)-ys(0))/(2*dx);
      if(norder>1)
	cs_[A][2] = NA/2 * (ys(2)+ys(0)-2.0*ys(1))/(dx*dx);
    }
  }
  void OpBufSpline::SetUp(Pwgto1c *basis) {

    double dx = 0.01;
    int norder = 2;
    
    for(int A = 0; A < basis->num_; A++) {

      auto nA = basis->ns_[A];

      // principle number
      for(int i = 0; i < norder; i++) {
	ns_[A][i] = nA+i;
      }

      // numerical derivative
      VectorXd xs(3);
      VectorXcd ys(3);
      xs(1) = basis->R0_;
      xs(0) = xs(1)-dx;
      xs(2) = xs(1)+dx;
      op_->At(xs, &ys);
      cs_[A][0] = ys(1);
      cs_[A][1] = (ys(2)-ys(0))/(2*dx);
      cs_[A][2] = (ys(2)+ys(0)-2.0*ys(1))/(dx*dx);
    }
  }
  OpBufSplineP1::OpBufSplineP1(const VectorXi& ns, int norder, OperatorSplineP1 *op):
    OpBufBasic(ns.size()), op_(op) {
    assert(0<=norder && norder<=2);

    for(int A = 0; A < ns.size(); A++) {
      int nA = ns(A);
      if(nA == 0) {
	this->InitZero(A, 2*(norder+1));
	for(int i = 0; i < norder+1; i++) {
	  ns_[A][2*i]   = nA+i+1;
	  ns_[A][2*i+1] = nA+i;
	}
      } else {
	this->InitZero(A, 3*(norder+1));
	for(int i = 0; i < norder+1; i++) {
	  ns_[A][3*i]   = nA+i+1;
	  ns_[A][3*i+1] = nA+i;
	  ns_[A][3*i+2] = nA+i-1;
	}
      }
    }
  }
  void OpBufSplineP1::SetUp(Pwgto *basis) {
    int norder = basis->norder();
    assert(0<=norder && norder<=2);
    double dx = basis->dx();
    complex<double> i(0, 1);
    
    for(int A = 0; A < basis->num(); A++) {

      auto nA = basis->ns()[A];
      auto NA = basis->Ns()[A];
      auto pA = basis->Ps()[A];
      auto gA = basis->gs()[A];

      // numerical derivative
      VectorXd xs(3);
      VectorXcd ys(3);
      VectorXcd ds(3);
      xs(1) = basis->Rs()[A];
      xs(0) = xs(1)-dx;
      xs(2) = xs(1)+dx;
      op_->At(xs, &ys);
      ds(0) = ys(1);
      ds(1) = (ys(2)-ys(0))/(2*dx);
      ds(2) = 0.5 * (ys(2)+ys(0)-2.0*ys(1))/(dx*dx);

      //
      // f = (v0 + v1.x + v2.xx).(-2g.x^(nA+1) + i.pA.x^nA + nA.x^{nA-1})
      if(nA==0) {
	ns_[A][2*0+0]   = nA+1; cs_[A][2*0+0] = -i*NA*ds(0) * (-2.0*gA);
	ns_[A][2*0+1]   = nA+0; cs_[A][2*0+1] = -i*NA*ds(0) * (i*pA);
	if(norder>0) {
	  ns_[A][2*1+0] = nA+2; cs_[A][2*1+0] = -i*NA*ds(1) * (-2.0*gA);
	  ns_[A][2*1+1] = nA+1; cs_[A][2*1+1] = -i*NA*ds(1) * (i*pA);
	}
	if(norder>1) {
	  ns_[A][2*2+0] = nA+3; cs_[A][2*2+0] = -i*NA*ds(2) * (-2.0*gA);
	  ns_[A][2*2+1] = nA+2; cs_[A][2*2+1] = -i*NA*ds(2) * (i*pA);
	}
      } else {
	ns_[A][3*0+0] = nA+1; cs_[A][3*0+0]   = -i*NA*ds(0) * (-2.0*gA);
	ns_[A][3*0+1] = nA+0; cs_[A][3*0+1]   = -i*NA*ds(0) * i*pA;
	ns_[A][3*0+2] = nA-1; cs_[A][3*0+2]   = -i*NA*ds(0) * (1.0*nA);
	if(norder>0) {
	  ns_[A][3*1+0] = nA+2; cs_[A][3*1+0] = -i*NA*ds(1) * (-2.0*gA);
	  ns_[A][3*1+1] = nA+1; cs_[A][3*1+1] = -i*NA*ds(1) * i*pA;
	  ns_[A][3*1+2] = nA+0; cs_[A][3*1+2] = -i*NA*ds(1) * (1.0*nA);
	}
	if(norder>1) {
	  ns_[A][3*2+0] = nA+3; cs_[A][3*2+0] = -i*NA*ds(2) * (-2.0*gA);
	  ns_[A][3*2+1] = nA+2; cs_[A][3*2+1] = -i*NA*ds(2) * (i*pA);
	  ns_[A][3*2+2] = nA+1; cs_[A][3*2+2] = -i*NA*ds(2) * (1.0*nA);
	}
      }

      
    }
  }
  void OpBufSplineP1::SetUp(Pwgto1c *) {
    throw runtime_error("not impl");
  }
  OpBufGausspot::OpBufGausspot(const VectorXi& ns, OperatorGausspot *op) : num_(ns_.size()),
									   ns_(ns), op_(op) {}  
  int OpBufGausspot::Maxn(int A) { return ns_(A); }
  void OpBufGausspot::Dump() const { cout << "gauss" << endl; }
  void OpBufGausspot::SetUp(Pwgto*) {}
  void OpBufGausspot::Matrix(OpBuf *opbra, Pwgto *basis, MatrixXcd *res) {
    opbra->Matrix(this, basis, res);
  }
  void OpBufGausspot::Matrix(OpBufBasic *, Pwgto *, MatrixXcd *) {
    assert(false&&"not impl");
  }
  void OpBufGausspot::Matrix(OpBufGausspot *, Pwgto *, MatrixXcd *) {}
  void OpBufGausspot::At(Pwgto*, const VectorXcd&, const VectorXd&, VectorXcd *) {
    throw runtime_error("not impl");
  }
  void OpBufGausspot::SetUp(Pwgto1c *) {}
  void OpBufGausspot::Matrix(OpBuf *opbra, Pwgto1c *basis, MatrixXcd *res) {
    opbra->Matrix(this, basis, res);
  }
  void OpBufGausspot::Matrix(OpBufBasic *, Pwgto1c *, MatrixXcd *) {
    assert(false&&"not impl");
  }
  void OpBufGausspot::Matrix(OpBufGausspot *, Pwgto1c *, MatrixXcd *) {}
  void OpBufGausspot::At(Pwgto1c*, const VectorXcd&, const VectorXd&, VectorXcd *) {
    throw runtime_error("not impl");
  }  
}
