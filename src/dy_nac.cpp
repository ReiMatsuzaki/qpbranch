#include <qpbranch/dy_nac.hpp>
#include <qpbranch/con.hpp>

namespace qpbranch{
  using mangan4::Con;

  DyNac::DyNac(const OpMat& opHeIJ, const OpMat& opXkIJ, const VectorXi& ns, string type_gauss) :
    numA_(ns.size()), numI_(opHeIJ.shape()[0]),
    q0_(0), p0_(0), gamma0_(1.0,0.0), 
    ns_(ns), m_(1.0),
    type_gauss_(type_gauss), is_setup_(false), opHeIJ_(opHeIJ), opXkIJ_(opXkIJ) {

    assert(type_gauss=="frozen"||type_gauss=="thawed");
    assert(numI_ == (int)opHeIJ.shape()[1]);
    assert(numI_ == (int)opXkIJ.shape()[0]);
    assert(numI_ == (int)opXkIJ.shape()[1]);

    id_ = new OperatorId();
    r1_ = new OperatorRn(1);
    r2_ = new OperatorRn(2);
    p1_ = new OperatorPn(1);
    p2_ = new OperatorPn(2);    
    DR_ = new OperatorDa(kIdDR);
    DP_ = new OperatorDa(kIdDP);
    vector<Operator*> ops{id_,r1_,r2_,p1_,p2_,DR_,DP_};
    
    for(int I = 0; I < numI_; I++)
      for(int J = 0; J < numI_; J++) {
	if(I<=J)
	  if(opHeIJ_[I][J] != nullptr)
	    ops.push_back(opHeIJ_[I][J]);
	if(I<J)
	  if(opXkIJ_[I][J] != nullptr)
	    ops.push_back(opXkIJ_[I][J]);
      }

    if(type_gauss == "frozen") {
      Dgr_ = nullptr;
      Dgi_ = nullptr;
    } else if(type_gauss == "thawed") {
      Dgr_ = new OperatorDa(kIdDgr); ops.push_back(Dgr_);
      Dgi_ = new OperatorDa(kIdDgi); ops.push_back(Dgi_);
    }

    basis_ = new Pwgto(ns, ops);
  
  }
  void DyNac::SetUp() {

    const double tol = pow(10.0, -10.0);    

    cAI_ = VectorXcd::Zero(numI_*numA_);
    cAI_(0) = 1;
    
    complex<double> norm2(0);
    MatrixXcd S(numA_,numA_);
    this->UpdateBasis();
    basis_->Matrix(id_, id_, &S);
    for(int A = 0; A < numA_; A++) {
      for(int B = 0; B < numA_; B++) {
	for(int I = 0; I < numI_; I++) {
	  int AI = this->idx(A,I);
	  int BI = this->idx(B,I);
	  norm2 += conj(cAI_[AI]) * cAI_[BI] * S(A,B);
	}
      }
    }
    if(abs(norm2) < tol) {
      cerr << __FILE__ << ":" << __LINE__ << ":Error too small norm for dA" << endl; abort();
    }    
    cAI_ *= 1.0/sqrt(norm2);
    
    is_setup_ = true;
  }
  void DyNac::Update(double dt) {
    assert(is_setup_);
    MatrixXcd H(numA_*numI_,numA_*numI_);
    this->Hamiltonian(id_, &H);
  }
  void DyNac::UpdateBasis() {

    for(int A = 0; A < numA_; A++) {
      basis_->ref_Rs()(A) = q0_;
      basis_->ref_Ps()(A) = p0_;
      basis_->ref_gs()(A) = gamma0_;
    }
    basis_->SetUp();
    
  }
  void DyNac::At(int I, const VectorXd& xs, VectorXcd *ys) const {
    VectorXcd c(numA_);
    for(int A = 0; A < numA_; A++) {
      c(A) = cAI_(this->idx(A,I));
    }
    basis_->At(id_, c, xs, ys);
  }
  void DyNac::Hamiltonian(Operator *op_bra, MatrixXcd *res) {
    assert(is_setup_);
    assert(res->rows()==numA_*numI_);
    assert(res->cols()==numA_*numI_);

    MatrixXcd S(numA_,numA_), P2(numA_,numA_), HeIJ(numA_,numA_), XkIJ(numA_,numA_);

    basis_->Matrix(id_,id_,&S);
    basis_->Matrix(id_,p2_,&P2);

    *res = MatrixXcd::Zero(numA_*numI_, numA_*numI_);
    for(int I = 0; I < numI_; I++) {
      for(int A = 0; A < numA_; A++) {
	for(int B = 0; B < numA_; B++) {
	  (*res)(idx(A,I),idx(B,I)) = 1.0/(2.0*m_) * P2(A,B);
	}
      }
    }

    for(int I = 0; I < numI_; I++) {
      for(int J = I; J < numI_; J++) {
	basis_->Matrix(id_, opHeIJ_[I][J], &HeIJ);
	for(int A = 0; A < numA_; A++) {
	  for(int B = 0; B < numA_; B++) {
	    (*res)(idx(A,I),idx(B,J)) += HeIJ(A,B);
	    if(I!=J)
	      (*res)(idx(A,J),idx(B,I)) += HeIJ(A,B);
	  }
	}
      }
    }

    complex<double> i(0,1);
    for(int I = 0; I < numI_; I++) {
      for(int J = I+1; J < numI_; J++) {
	basis_->Matrix(p1_, opXkIJ_[I][J], &XkIJ);
	for(int A = 0; A < numA_; A++) {
	  for(int B = 0; B < numA_; B++) {
	    (*res)(idx(A,I),idx(B,J)) += -i/m_*XkIJ(A,B);
	    (*res)(idx(A,J),idx(B,I)) += +i/m_*XkIJ(A,B);
	  }
	}
      }
    }

  }
  void DyNac::DotxQhamilton(string vp_name, VectorXd *res) {
  }  
}
