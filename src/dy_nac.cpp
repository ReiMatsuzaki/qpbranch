#include <qpbranch/dy_nac.hpp>
#include <qpbranch/con.hpp>
#include <qpbranch/eigenplus.hpp>

namespace qpbranch{
  using mangan4::Con;

  DyNac::DyNac(const OpMat& opHeIJ, const OpMat& opXkIJP, const VectorXi& ns,
	       string type_gauss, int norder, double dx) :
    numA_(ns.size()), numI_(opHeIJ.shape()[0]),
    q0_(0), p0_(0), gamma0_(1.0,0.0), 
    ns_(ns), m_(1.0),
    type_gauss_(type_gauss), is_setup_(false), opHeIJ_(opHeIJ), opXkIJP_(opXkIJP) {

    assert(type_gauss=="frozen"||type_gauss=="thawed");
    assert(numI_ == (int)opHeIJ.shape()[1]);
    assert(numI_ == (int)opXkIJP.shape()[0]);
    assert(numI_ == (int)opXkIJP.shape()[1]);

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
	  if(opXkIJP_[I][J] != nullptr)
	    ops.push_back(opXkIJP_[I][J]);
      }

    if(type_gauss == "frozen") {
      Dgr_ = nullptr;
      Dgi_ = nullptr;
    } else if(type_gauss == "thawed") {
      Dgr_ = new OperatorDa(kIdDgr); ops.push_back(Dgr_);
      Dgi_ = new OperatorDa(kIdDgi); ops.push_back(Dgi_);
    }

    basis_ = new Pwgto(ns, ops, norder, dx);
  
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

    this->UpdateBasis();

    // non linear parameter
    VectorXd dotx(4);
    this->DotxQhamilton(&dotx);

    // coefficient
    MatrixXcd H(numA_*numI_,numA_*numI_), S(numA_*numI_,numA_*numI_);    
    //    this->EffHamiltonian(dotx, &H);
    this->Hamiltonian(id_, &H);
    this->Overlap(&S);

    // propagate nonlinear
    q0_ += dotx(0) *dt;
    p0_ += dotx(1) *dt;
    if(type_gauss_=="thawed") {
      gamma0_ += complex<double>(dotx(2), dotx(3))*dt;
    }

    // propagate coefficient
    IntetGdiag(H, S, dt, &cAI_);
    
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

    MatrixXcd P2(numA_,numA_), HeIJ(numA_,numA_), XkIJP(numA_,numA_);

    *res = MatrixXcd::Zero(numA_*numI_, numA_*numI_);
    
    basis_->Matrix(op_bra,p2_,&P2);    
    for(int I = 0; I < numI_; I++) {
      for(int A = 0; A < numA_; A++) {
	for(int B = 0; B < numA_; B++) {
	  (*res)(idx(A,I),idx(B,I)) = 1.0/(2.0*m_) * P2(A,B);
	}
      }
    }

    for(int I = 0; I < numI_; I++) {
      for(int J = I; J < numI_; J++) {
	if(opHeIJ_[I][J]==nullptr)
	  continue;
	basis_->Matrix(op_bra, opHeIJ_[I][J], &HeIJ);
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
	if(opXkIJP_[I][J]==nullptr)
	  continue;
	basis_->Matrix(op_bra, opXkIJP_[I][J], &XkIJP);
	for(int A = 0; A < numA_; A++) {
	  for(int B = 0; B < numA_; B++) {
	    (*res)(idx(A,I),idx(B,J)) += -i/m_*XkIJP(A,B);
	    (*res)(idx(A,J),idx(B,I)) += +i/m_*XkIJP(A,B);
	  }
	}
      }
    }

  }
  void DyNac::EffHamiltonian(const VectorXd& dotx, MatrixXcd *ptr_res) {
    assert(is_setup_);

    complex<double> i(0,1);
    MatrixXcd M(numA_*numI_,numA_*numI_);
    MatrixXcd& res(*ptr_res);

    this->Hamiltonian(id_, &res);

    basis_->Matrix(id_, DR_, &M);
    res += -i*M*dotx(0);

    basis_->Matrix(id_, DP_, &M);
    res += -i*M*dotx(1);

    if(type_gauss_=="thawed") {
      basis_->Matrix(id_, Dgr_, &M);
      res += -i*M*dotx(2);
      basis_->Matrix(id_, Dgi_, &M);
      res += -i*M*dotx(3);
    } else if(type_gauss_=="frozen") {
    } else {
      cerr << __FILE__ << ":" << __LINE__ << ":Error invalid type_gauss" << endl; abort();
    }
    
  }
  void DyNac::Overlap(MatrixXcd *ptr_res) {
    assert(ptr_res->rows()==numI_*numA_);
    assert(ptr_res->rows()==ptr_res->cols());

    MatrixXcd& res(*ptr_res);
    res.array() = 0.0;

    MatrixXcd S = MatrixXcd::Zero(numA_,numA_);
    basis_->Matrix(id_, id_, &S);
    
    for(int I = 0; I < numI_; I++) {
      for(int A = 0; A < numA_; A++) {
	for(int B = 0; B < numA_; B++) {
	  res(idx(A,I),idx(A,I)) = S(A,B);
	}
      }      
    }

    
  }
  void DyNac::DotxQhamilton(VectorXd *res) {
    assert(numA_==1);
    assert(res->size()==4);

    MatrixXcd H(numI_,numI_);

    *res = VectorXd::Zero(4);

    this->Hamiltonian(DP_, &H);  // dot{R} = dH/dP
    (*res)[0] = 2*real(cAI_.dot(H*cAI_));  
    this->Hamiltonian(DR_, &H);  // dot{P} = -dH/dR
    (*res)[1] = -2*real(cAI_.dot(H*cAI_));  
    
    if(type_gauss_=="thawed") {
      auto gr2 = pow(real(gamma0_), 2);
      this->Hamiltonian(Dgi_, &H); // dot{gi} = (4gr^2)  dH/dgi
      (*res)[2] = (+4.0*gr2) * 2.0*real(cAI_.dot(H*cAI_));
      this->Hamiltonian(Dgr_, &H); // dot{gr} = (-4gr^2) dH/dgr
      (*res)[3] = (-4.0*gr2) * 2.0*real(cAI_.dot(H*cAI_));
    }
  }
  void DyNac::CalcProb(VectorXd *res) {
    assert(res->size()==numI_);

    for(int I = 0; I < numI_; I++) {
      double cumsum(0);
      for(int A = 0; A < numA_; A++) {
	cumsum += pow(abs(cAI_(idx(A,I))), 2);
      }
      (*res)(I) = cumsum;
    }
    
  }
  void DyNac::DumpCon(int it, string prefix) {

    Con& con = Con::getInstance();

    con.write_f(prefix+"q0", it, q0_);
    con.write_f(prefix+"p0", it, p0_);
    con.write_f(prefix+"gr0", it, gamma0_.real());
    con.write_f(prefix+"gi0", it, gamma0_.imag());
    con.write_f1(prefix+"c_re", it, cAI_.real());
    con.write_f1(prefix+"c_im", it, cAI_.imag());

    VectorXd prob(numI_);
    this->CalcProb(&prob);
    con.write_f1(prefix+"prob", it, prob);

    if(it==0) {
      con.write_i(prefix+"_numA", 0, numA_);
      con.write_i(prefix+"_numI", 0, numI_);
    }
    
  }
}
