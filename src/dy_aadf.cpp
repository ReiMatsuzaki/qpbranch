#include <qpbranch/pwgto.hpp>
#include <qpbranch/dy_aadf.hpp>
#include <qpbranch/con.hpp>
#include <qpbranch/eigenplus.hpp>

namespace qpbranch {
  using mangan4::Con;
  
  DyAadf::DyAadf(OperatorPot *pot, const VectorXi& ns, string type_gauss) {
    assert(type_gauss=="frozen"||type_gauss=="thawed");

    num_ = ns.size();

    q0_ = 0.0;
    p0_ = 0.0;
    rho_ = 1.0;
    lambda_ = 0.0;
    theta_ = 0.0;

    c_ = VectorXcd::Zero(num_);
    c_[0] = 1;

    m_ = 1.0;

    type_gauss_ = type_gauss;

    is_setup_ = false;

    pot_ = pot;
    id_ = new OperatorId();
    p2_ = new OperatorPn(2);
    r1_ = new OperatorRn(1);
    r2_ = new OperatorRn(2);
    DR_ = new OperatorDa(kIdDR);
    DP_ = new OperatorDa(kIdDP);
    vector<Operator*> ops = {id_, pot_, p2_, DR_, DP_, r1_, r2_};
    if(type_gauss == "frozen") {
      Dgr_ = nullptr;
      Dgi_ = nullptr;
    } else if(type_gauss == "thawed") {
      Dgr_ = new OperatorDa(kIdDgr);
      Dgi_ = new OperatorDa(kIdDgi);
      ops.push_back(Dgr_);
      ops.push_back(Dgi_);
    }

    basis_ = new Pwgto(ns, ops);
    
  }
  void DyAadf::SetUp() {
    basis_->SetUp();

    MatrixXcd S(num_,num_);
    basis_->Matrix(id_, id_, &S);

    double norm2 = real(c_.dot(S*c_));
    c_ *= 1/sqrt(norm2);
    
    is_setup_ = true;
  }
  void DyAadf::Update(double dt) {
    assert(is_setup_);

    this->UpdateBasis();

    // time derivative of nonlinear parameters
    VectorXd dotx(4);
    this->DotxQhamilton(&dotx);

    // propagate variables
    q0_ += dotx(0) * dt;
    p0_ += dotx(1) * dt;
    if(type_gauss_=="thawed") {
      rho_ += dotx(2) *dt;
      lambda_  += dotx(3) *dt;
    }
    
  }
  void DyAadf::DotxQhamilton(VectorXd *res) {
    assert(is_setup_);
    assert(res->size()==4);

    const double dx = 0.001;

    int num = 2;
    if(type_gauss_=="frozen") {
      num = 2;
    } else if(type_gauss_=="thawed") {
      num = 4;
    } else {
      cerr << "invalid type gauss" << endl; abort();
    }

    MatrixXcd H(num_,num_);
    VectorXd dH(4);
    double ep, em;
    for(int i = 0; i < num; i++) {
      this->ShiftVar(i, dx);
      this->Hamiltonian(id_, &H); ep = real(c_.dot(H*c_));
      this->ShiftVar(i, -2*dx);
      this->Hamiltonian(id_, &H); em = real(c_.dot(H*c_));
      this->ShiftVar(i, dx);
      dH[i] = (ep-em)/(2*dx);
    }
    (*res)[0] = +dH[1]; // dot{q} = dH/dp
    (*res)[1] = -dH[0]; // dot{p} = -dH/dq

    if(type_gauss_=="frozen") {
      (*res)[2] = 0.0;
      (*res)[3] = 0.0;
    } else if(type_gauss_=="thawed") {
      (*res)[2] = -dH[3]; // dot{rho}    = +dH/d\lambda
      (*res)[3] = +dH[2]; // dot{lambda} = -dH/d\rho
    } else {
      cerr << __FILE__ << ":" << __LINE__ << ":Error invalid type gauss" << endl; abort();
    }
  }
  void DyAadf::HamiltonianForF(Operator *op, MatrixXcd *res) {
    assert(is_setup_);
    MatrixXcd P2(num_,num_), V(num_,num_);
    basis_->Matrix(op, p2_,  &P2);
    basis_->Matrix(op, pot_, &V );
    *res = 1.0/(2*m_) * P2 + V;
  }
  void DyAadf::Hamiltonian(Operator *op, MatrixXcd *res) {
    assert(is_setup_);

    MatrixXcd S(num_, num_), R1(num_,num_), R2(num_,num_), HF(num_,num_);
    basis_->Matrix(op, id_,  &S );
    basis_->Matrix(op, r1_,  &R1 );
    basis_->Matrix(op, r2_,  &R2 );
    this->HamiltonianForF(op, &HF);

    *res = 1.0/(2*m_)*( p0_*p0_*S +
			-2*p0_*lambda_/rho_*R1 +
			pow(lambda_/rho_,2)*R2) + HF;

  }
  void DyAadf::ShiftVar(int var_id, double dx) {
    if(var_id==0) 
      this->q0_ += dx;
    else if(var_id==1)
      this->p0_ += dx;
    else if(var_id==2)
      this->rho_ += dx;
    else if(var_id==3)
      this->lambda_ += dx;
    else {
      cerr << __FILE__ << ":" << __LINE__ << ": Error invalid var_id." << endl;
      cerr << "var_id: " << var_id << endl;
      abort();
    }
    this->UpdateBasis();
  }
  void DyAadf::UpdateBasis() {
    for(int i = 0; i < num_; i++) {
      basis_->ref_Rs()(i) = q0_;
      basis_->ref_Ps()(i) = 0.0;
      basis_->ref_gs()(i) = complex<double>(1/(4*rho_*rho_), -lambda_/(2*rho_));
    }
    basis_->SetUp();
  }
  void DyAadf::DumpCon(int it, string prefix) {
    Con& con = Con::getInstance();

    con.write_f("q0",     it, q0_);
    con.write_f("p0",     it, q0_);
    con.write_f("rho",    it, rho_);
    con.write_f("lambda", it, lambda_);

    if(it==0) {
      con.write_i("_num", 0, num_);
    }

    con.write_f1(prefix+"c_re", it, c_.real());
    con.write_f1(prefix+"c_im", it, c_.imag());

    double norm2 = Norm2();
    con.write_f("norm2", it, norm2);
    
  }
  double DyAadf::Norm2() const {
    MatrixXcd S(num_,num_);
    basis_->Matrix(id_, id_, &S);
    return (c_.dot(S*c_)).real();
  }
}


