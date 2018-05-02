#include <qpbranch/pwgto.hpp>
#include <qpbranch/pwgto_buf.hpp>
#include <qpbranch/mathplus.hpp>
#include <qpbranch/con.hpp>
using namespace std;
using boost::extents;

namespace qpbranch {

  using namespace mangan4;

  void HermiteCoefD(complex<double> gAB, complex<double> wAB, double RA, double RB,
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
  complex<double> HermiteCoefD0(complex<double> gAB, complex<double> wAB, double RA, double RB,
				int nA, int nB, int Nk) {
    if(nA==0 && nB==0 && Nk==0)
      return 1.0;
    
    if(Nk<0 || Nk>nA+nB)
      return 0.0;

    if(nA>0) {
      auto res0 = HermiteCoefD0(gAB, wAB, RA, RB, nA-1, nB, Nk-1);
      auto res1 = HermiteCoefD0(gAB, wAB, RA, RB, nA-1, nB, Nk  );
      auto res2 = HermiteCoefD0(gAB, wAB, RA, RB, nA-1, nB, Nk+1);
      return 1.0/(2.0*gAB)*res0 + (wAB-RA)*res1 + (Nk+1.0)*res2;
    } else {
      auto res0 = HermiteCoefD0(gAB, wAB, RA, RB, nA, nB-1, Nk-1);
      auto res1 = HermiteCoefD0(gAB, wAB, RA, RB, nA, nB-1, Nk  );
      auto res2 = HermiteCoefD0(gAB, wAB, RA, RB, nA, nB-1, Nk+1);
      return 1.0/(2.0*gAB)*res0 + (wAB-RB)*res1 + (Nk+1.0)*res2;
    }    
  }
 
  Pwgto::Pwgto(const VectorXi& ns, const vector<Operator*>& ops):
    num_(ns.size()), nop_(ops.size()), ns_(ns), ops_(ops),
    gs_(num_), Rs_(num_), Ps_(num_), is_setup_(false),
    Ns_(num_), maxn_(num_),
    gAB_(num_,num_), eAB_(num_,num_), hAB_(num_,num_), RAB_(num_,num_) {

    // allocate buffer for each operator
    for(auto it = ops.begin(); it!=ops.end(); ++it) {
      Operator *op = *it;
      auto ptr = MakeOpBuf(ns, op);
      buffer_map_[op] = ptr;
    }

    // calculate
    for(int A = 0; A < num_; A++) {
      int maxn(0);
      for(auto it =  buffer_map_.begin(); it != buffer_map_.end(); ++it) {
	int maxn0 = it->second->Maxn(A);
	if(maxn<maxn0)
	  maxn = maxn0;
      }
      maxn_[A] = maxn;
    }

    // allocate coefficient d
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
  void Pwgto::SetUp() {

    // normalization term 
    for(int A = 0; A < num_; A++) {
      int n(ns_[A]);
      VectorXcd gg(n+1);
      IntGto2N(n, gs_[A]+conj(gs_[A]), &gg);
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
	
	HermiteCoefD(gAB_(A,B), RAB_(A,B), Rs_[A], Rs_[B],
		     maxn_[A], maxn_[B], (*d_)[A][B]);
      }
    }    

    // each operator (visitor pattern)
    for(auto it = buffer_map_.begin(); it != buffer_map_.end(); ++it) {
      it->second->SetUp(this);
    }

    is_setup_ = true;
    
  }
  void Pwgto::Matrix(Operator *opbra, Operator *opket, MatrixXcd *res) {
    /* compute overlap matrix. */

    assert(is_setup_);
    assert(res->rows()>=num_);
    assert(res->cols()>=num_);
    assert(buffer_map_.find(opbra)!=buffer_map_.end());
    assert(buffer_map_.find(opket)!=buffer_map_.end());

    buffer_map_[opket]->Matrix(buffer_map_[opbra], this, res);
    
  }  
  void Pwgto::At(Operator *op, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res) {

    assert(is_setup_);
    buffer_map_[op]->At(this, cs, xs, res);

  }
  void Pwgto::DumpCon(int it) {
    assert(is_setup_);
    Con& con = Con::getInstance();
    con.write_f1("R", it, Rs_);
    con.write_f1("P", it, Ps_);
    con.write_f1("gr", it, gs_.real());
    con.write_f1("gi", it, gs_.imag());
  }

}

