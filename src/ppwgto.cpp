#include <iostream>
#include <qpbranch/pwgto.hpp>
#include <qpbranch/ppwgto.hpp>
#include <qpbranch/mathplus.hpp>
#include <qpbranch/con.hpp>

using namespace std;

namespace qpbranch {

  using boost::extents;
  using namespace mangan4;

  OpBuff *MakeOpBuff(int num, Operator *op) {
    const type_info& optype = typeid(*op);
    if(optype == typeid(OperatorId)) {
      return new OpBuffId(num);
    } else if(optype == typeid(OperatorRn)) {
      auto opop = static_cast<OperatorRn*>(op);
      return new OpBuffRn(num, opop->n());
    } else if(optype == typeid(OperatorPn)) {
      auto opop = static_cast<OperatorPn*>(op);
      return new OpBuffPn(num, opop->n());
    } else if(optype == typeid(OperatorDa)) {
      auto opop = static_cast<OperatorDa*>(op);
      return new OpBuffDa(num, opop->id());
    } else if(optype == typeid(OperatorGausspot)) {
      auto opop = static_cast<OperatorGausspot*>(op);
      return new OpBuffGausspot(opop);
    } else {
      cerr << "Unsupported operator" << endl;
      abort();
    }
  }  
  void OpBuffBasic::MatrixAsKet(OpBuff *opbra, PolyPwgto *basis, MatrixXcd *res) {
    opbra->MatrixAsBra(this, basis, res);
  }
  void OpBuffBasic::MatrixAsBra(OpBuffBasic *opket, PolyPwgto *basis, MatrixXcd *res) {
    for(int A = 0; A < basis->size(); A++) {
      for(int B = 0; B < basis->size(); B++) {
	complex<double> cumsum(0);
	for(auto it = this->poly_[A].begin(); it != this->poly_[A].end(); ++it) {
	  for(auto jt = opket->poly_[B].begin(); jt != opket->poly_[B].end(); ++jt) {
	    cumsum += conj(it->first) * jt->first * basis->get_d(A,B,it->second,jt->second,0);
	  }
	}
	(*res)(A,B) = cumsum;
      }
    }
  }
  void OpBuffBasic::MatrixAsBra(OpBuffGausspot *opket, PolyPwgto *basis, MatrixXcd *res) {
    auto opV = opket->op();
    for(int A = 0; A < basis->size(); A++) {
      for(int B = 0; B < basis->size(); B++) {
	int max_nA_nB = this->poly_[A].Order() + basis->poly(B).Order();
	VectorXcd intg(max_nA_nB+1);
	DwIntGauss(max_nA_nB, opV->b(), basis->gAB_(A,B), basis->RAB_(A,B), &intg);
	
	complex<double> cumsum(0);
	for(auto it = this->poly_[A].begin(); it != this->poly_[A].end(); ++it) {
	  for(auto jt = basis->poly(B).begin(); jt != this->poly(B).end(); ++jt) {
	    complex<double> cumsum1(0);
	    for(int N = 0; N <= max_nA_nB; N++) 
	      cumsum1 += basis->get_d(A,B,it->second,jt->second,N);
	    cumsum += conj(it->first)*jt->first*cumsum1;
	    }
	  }
	(*res)(A,B) = cumsum * basis->eAB_(A,B) * opV->v0();
	}
      }
    }
  void OpBuffBasic::At(PolyPwgto *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res) {
    complex<double> ii(0.0, 1.0);
    for (int ix = 0; ix < xs.size(); ix++) {
      complex<double> cumsum(0.0);      
      for (int A = 0; A < basis->size(); A++) {
	auto gA = basis->gs()[A];
	auto pA = basis->Ps()[A];
	double d = xs[ix] - basis->Rs()[A];
	complex<double> cumsum1(0.0);
	const Poly& p(basis->poly(A));
	for(auto it = p.begin(); it != p.end(); ++it) {
	  cumsum1 += it->first*pow(d, it->second) * exp(-gA*d*d - ii*pA*d);
	}
	cumsum += cs[A]*cumsum1;
      }
      (*res)[ix] = cumsum;
    }    
  }
  
  void OpBuffId::SetUp(PolyPwgto *basis) {
    for(int A = 0; A < basis->size(); A++) {
      poly_[A] = basis->poly(A);
    }
  }
  void OpBuffRn::SetUp(PolyPwgto *basis) {
    Poly Rn({{1.0, n_}});
    for(int A = 0; A < basis->size(); A++) {
      poly_[A] = basis->poly(A);
      poly_[A] *= Rn;
    }
  }
  void OpBuffPn::SetUp(PolyPwgto *basis) {
    Poly R1({{1.0, 1}});
    Poly R2({{1.0, 2}});
    complex<double> ii(0,1);
    complex<double> c = pow(ii, n_);
    for(int A = 0; A < basis->size(); A++) {      
      auto gA = basis->gs()(A);
      auto pA = basis->Ps()(A);
      const Poly& p = basis->poly(A);
      if(n_==1) {
	auto t1 = p.Diff(n_);
	auto t2 = 2.0*gA*R1*p  - ii*pA*p;
	poly_[A] = c * (t1 + t2);
      } else if(n_==2) {
	auto d1p = p.Diff(1);
	auto d2p = d1p.Diff(1);
	auto t2 = 2.0*gA*R1*d1p  - ii*pA*d1p;
	auto t3 = (4.0*gA*gA)*R2*p  - (2.0*(2.0*gA)*(ii*pA))*d1p;
	poly_[A] = c * (d2p + t2 + t3);
      } else {
	cerr << "not supported n" << endl;
	abort();
      }
    }
  }
  OpBuffDa::OpBuffDa(int num, int id) : OpBuffBasic(num), id_(id) {}  
  int OpBuffDa::Maxn() const {
      if(id_==1)
	return 1;
      if(id_==2)
	return 1;
      cerr << "invalid id" << endl;
      abort();
    }
  void OpBuffDa::SetUp(PolyPwgto *basis) {
    complex<double> ii(0.0, 1.0);
    Poly R1({{1.0,1}});
    for(int A = 0; A < basis->size(); A++) {
      auto gA = basis->gs()(A);
      auto pA = basis->Ps()(A);
      const Poly& p = basis->poly(A);
      if(id_==1) { // DR
	auto t1 = p.Diff(1);
	auto t2 = 2.0*gA*R1*p  - ii*pA*p;
	poly_[A] = t1 + t2;
      } else if(id_==2) {
	poly_[A] = ii * R1 * p;
	break;
      } else {
	cerr << "invalid id" << endl;
	break;
      }
      
    }
  }
  
  PolyPwgto::PolyPwgto(const vector<Poly>& polys, const vector<Operator*>& ops):
    num_(polys.size()), nop_(ops.size()), polys_(polys), ops_(ops),
    gs_(num_), Rs_(num_), Ps_(num_), is_setup_(false),
    Ns_(num_), maxn_(num_),
    gAB_(num_,num_), eAB_(num_,num_), hAB_(num_,num_), RAB_(num_,num_) {

    op_id_ = NULL;
    for(auto it = ops.begin(); it!=ops.end(); ++it) {
      Operator *op = *it;
      buff_map_[op] = MakeOpBuff(num_, op);
      if(typeid(*op) == typeid(OperatorId))
	op_id_ = op;
    }
    if(op_id_==NULL)
      throw runtime_error("ops do not include OperatorId");

    // claculate
    for(int A = 0; A < num_; A++) {
      int maxn(0);
      for(auto it =  buff_map_.begin(); it != buff_map_.end(); ++it) {
	int maxn0 = it->second->Maxn();
	if(maxn<maxn0)
	  maxn = maxn0;
      }
      maxn_[A] = maxn + polys[A].Order();
    }

    // coef d
    d_ = new multi_array< multi_array<complex<double>,3>*, 2>(extents[num_][num_]);
    for(int A = 0; A < num_; A++) {
      for(int B = 0; B < num_; B++) {
	(*d_)[A][B] = new multi_array<complex<double>,3>(extents[maxn_[A]+1][maxn_[B]+1][maxn_[A]+maxn_[B]+1]);
      }
    }

    // initialize variables
    gs_ = VectorXcd::Constant(num_, complex<double>(1.0, 0.0));
    Rs_ = VectorXd::Constant( num_, 0.0);
    Ps_ = VectorXd::Constant( num_, 0.0);

  }
  void PolyPwgto::SetUp() {

    // normalization term 
    for(int A = 0; A < num_; A++) {
      auto& p(poly(A));
      int nA(0);
      for(auto it = p.begin(); it != p.end(); ++it) {
	nA = std::max(nA, it->second);
      }
            
      VectorXcd gg(nA+1);
      IntGto2N(nA, gs_[A]+conj(gs_[A]), &gg);

      complex<double> cumsum(0);
      for(auto it = p.begin(); it != p.end(); ++it) {
	for(auto jt = p.begin(); jt != p.end(); ++jt) {
	  int n = it->second + jt->second;
	  if(n%2==0)
	    cumsum += conj(it->first)*jt->first * gg(n/2);
	}
      }
      Ns_[A] = sqrt(1.0/cumsum.real());
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
	
	HermiteCoefD(gAB_(A,B), RAB_(A,B), Rs_[A], Rs_[B],
		     maxn_[A], maxn_[B], (*d_)[A][B]);
      }
    }    

    // each operator (visitor pattern)
    for(auto it =ops_.begin(); it != ops_.end(); ++it) {
      buff_map_[*it]->SetUp(this);
    }

    is_setup_ = true;
  }
  void PolyPwgto::Matrix(Operator *opbra, Operator *opket, MatrixXcd *res) {
    
    assert(is_setup_);
    assert(res->rows()>=num_);
    assert(res->cols()>=num_);
    assert(buff_map_.find(opbra)!=buff_map_.end());
    assert(buff_map_.find(opket)!=buff_map_.end());

    buff_map_[opket]->MatrixAsKet(buff_map_[opbra], this, res);
    
  }  
  void PolyPwgto::At(Operator *op, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res) {

    assert(is_setup_);
    buff_map_[op]->At(this, cs, xs, res);

  }
  void PolyPwgto::DumpCon(int it) {
    
    assert(is_setup_);
    Con& con = Con::getInstance();
    con.write_f1("R", it, Rs_);
    con.write_f1("P", it, Ps_);
    con.write_f1("gr", it, gs_.real());
    con.write_f1("gi", it, gs_.imag());
  }

}
