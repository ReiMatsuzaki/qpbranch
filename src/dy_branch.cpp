#include <boost/multi_array.hpp>
#include <qpbranch/pwgto.hpp>
#include <qpbranch/dy_branch.hpp>
#include <qpbranch/eigenplus.hpp>
#include <qpbranch/con.hpp>

namespace qpbranch {
  
  using boost::multi_array;
  using boost::extents;
  using namespace Eigen;
  using namespace std;
  using namespace mangan4;

  DySetPoly::DySetPoly(Operator *pot, const VectorXi& ns, string type_gauss) {

    assert(type_gauss=="frozen"||type_gauss=="thawed");

    num_ = ns.size();
    if(type_gauss=="frozen") 
      numopt_ = 2;
    else
      numopt_ = 4;
    
    q0_ = 0.0;
    p0_ = 0.0;
    gr0_ = 1.0;
    gi0_ = 0.0;
    c_ = VectorXcd::Zero(num_);
    c_(0) = 1.0;

    m_ = 1.0;

    type_gauss_ = type_gauss;
    type_eomslow_ = "qhamilton";

    is_setup_ = false;
    
    pot_ = pot;
    id_  = new OperatorId();
    p2_  = new OperatorPn(2);
    DR_ = new OperatorDa(kIdDR);
    DP_ = new OperatorDa(kIdDP);
    vector<Operator*> ops = {id_, pot_, p2_, DR_, DP_};
    ops_opt_.push_back(DR_);
    ops_opt_.push_back(DP_);
    if(type_gauss=="frozen") {
      Dgr_ = nullptr;
      Dgi_ = nullptr;
    } else if(type_gauss=="thawed") {
      Dgr_ = new OperatorDa(kIdDgr);
      Dgi_ = new OperatorDa(kIdDgi);
      ops.push_back(Dgr_);
      ops.push_back(Dgi_);
      ops_opt_.push_back(Dgr_);
      ops_opt_.push_back(Dgi_);
    }

    w_ = VectorXd::Zero(num_);
    U_ = MatrixXcd::Zero(num_,num_);
    Clam_ = VectorXcd::Zero(num_);
    
    basis_ = new Pwgto(ns, ops);
    
  }
  void DySetPoly::SetUp() {

    basis_->SetUp();
    MatrixXcd S(num_,num_);
    basis_->Matrix(id_, id_, &S);

    double norm2 = real(c_.dot(S*c_));
    c_ *= 1/sqrt(norm2);

    is_setup_ = true;
  }
  void DySetPoly::Update(double dt) {
    assert(is_setup_);

    this->UpdateBasis();

    // slow part
    VectorXd dotx(numopt_);
    if(type_eomslow_=="qhamilton") {
      this->DotxQhamilton(&dotx);
    } else if(type_eomslow_=="tdvp") {
      this->DotxQuantum(true, &dotx);					       
    } else if(type_eomslow_=="mp") {
      this->DotxQuantum(false, &dotx);
    }

    // fast part
    MatrixXcd H(num_,num_), S(num_,num_);
    basis_->Matrix(id_, id_, &S);
    this->EffHamiltonian(dotx, &H);    

    // propagate slow part
    q0_ += dotx(0) * dt;
    p0_ += dotx(1) * dt;
    if(type_gauss_=="thawed") {
      gr0_ += dotx(2) * dt;
      gi0_ += dotx(3) * dt;
    }

    // propagate fast part
    IntetGdiag(H, S, dt, &c_);

    // save Eigen functions for internal motion;
    double p0_bak = p0_;
    p0_ = 0.0;
    this->UpdateBasis();
    this->Hamiltonian(id_, &H);
    Zhegv(H, S, &U_, &w_);
    p0_ = p0_bak;
    this->UpdateBasis();
    Clam_ = U_.adjoint() * S * c_;

    this->Hamiltonian(id_, &H);
    cerr  << "E = " << (c_.dot(H*c_)).real() << endl; // => E = 0.11275
    
  }
  //
  // see 2018_aadf/0508_doc for detail.
  // Assume the single gaussian trial function
  //       G(q) = (1/2\pi\alpha)^{1/4} exp[-\gamma(q-q_0)^2 + ip_0(q-q_0)]
  //       \gamma = 1/(4\alpha) + i\beta
  // It can be proven that (q_0,p_0) and (\beta, \alpha) are canonical pair and satisfy
  // Hamilton equation. 
  //       \dot{\alpha} = -dH/d\beta
  //       \dot{\beta}  = +dH/d\alpha
  // In this program, orbital exponent are treated as its real(gr) and imaginary(gi) part.
  // Time derivative of them are
  //       \dot{gr} = d/dt (1/4\alpha)
  //                = -1/(4\alpha^2) \dot{\alpha}
  //                =  1/(4\alpha^2) dH/d\beta
  //                =  4gr^2 dH/d(gi)
  //       \dot{gi} = d\beta/dt
  //                = +dH/d\alpha
  //                = +dH/d(gr) d(gr)/d\alpha
  //                = +dH/d(gr) (-1)/(4\alpha^)
  //                = -4gr^2 dH/d(gr) 
  void DySetPoly::DotxQhamilton(VectorXd *res) {
    assert(is_setup_);
    
    MatrixXcd H(num_,num_);
    
    for(int iop = 0; iop < (int)ops_opt_.size(); iop++) {
      if(ops_opt_[iop] == DR_) {
	this->Hamiltonian(DP_, &H);
	(*res)[iop] = 2.0 * real(c_.dot(H*c_));
      } else if(ops_opt_[iop] == DP_) {
	this->Hamiltonian(DR_, &H);
	(*res)[iop] = -2.0 * real(c_.dot(H*c_));
      } else if(ops_opt_[iop] == Dgr_) {
	this->Hamiltonian(Dgi_, &H);
	(*res)[iop] = (4.0*gr0_*gr0_) * 2.0 * real(c_.dot(H*c_));
      } else if(ops_opt_[iop] == Dgi_) {
	this->Hamiltonian(Dgr_, &H);
	(*res)[iop] = (-4.0*gr0_*gr0_) * 2.0 * real(c_.dot(H*c_));
      } else {
	throw runtime_error("invalid iop");
      }      
    }
    
  }
  void DySetPoly::DotxQuantum(bool is_tdvp, VectorXd *res) {

    assert(is_setup_);

    VectorXcd v(num_), v1(num_);
    MatrixXd C(numopt_,numopt_);
    VectorXd y(numopt_);
    
    multi_array<MatrixXcd, 2> Snm(extents[numopt_+1][numopt_+1]);
    multi_array<MatrixXcd, 1> Hn0(extents[numopt_+1]);
    vector<Operator*> ops(ops_opt_);
    ops.push_back(id_);

    int id = numopt_; // index for operator id_

    for(int n = 0; n < 1+numopt_; n++) {
      Hn0[n] = MatrixXcd::Zero(num_,num_);
      for(int m = 0; m < 1+numopt_; m++) {
	Snm[n][m] = MatrixXcd::Zero(num_,num_);
      }
    }

    for(int n = 0; n < 1+numopt_; n++) {
      this->Hamiltonian(ops[n], &Hn0[n]);
      for(int m = 0; m < 1+numopt_; m++) {	
	basis_->Matrix(ops[n], ops[m], &(Snm[n][m]));
      }
    }

    for(int n = 0; n < numopt_; n++) {
      for(int m = 0; m < numopt_; m++) {
	auto t1 = c_.dot(Snm[n][m]*c_);
	v = Snm[id][m]*c_;
	Zgesv(Snm[id][id], v, &v1);
	auto t2 = c_.dot(v1);
	C(n,m) = is_tdvp ? -imag(t1-t2) : real(t1-t2);
      }
      auto t1 = c_.dot(Hn0[n]*c_);
      v = Hn0[id]*c_;
      Zgesv(Snm[id][id], v, &v1);
      v = Snm[n][0]*v;
      auto t2 = c_.dot(v1);
      y(n) = is_tdvp ? real(t1-t2) : imag(t1-t2);
    }

    Dgesv(C, y, res);
    
  }
  void DySetPoly::Hamiltonian(Operator *op_bra, MatrixXcd *res) {
    assert(is_setup_);
    
    MatrixXcd P2(num_, num_), V(num_,num_);    
    basis_->Matrix(op_bra, p2_,  &P2);
    basis_->Matrix(op_bra, pot_, &V);
    *res = 1.0/(2*m_) * P2 + V;
  }  
  void DySetPoly::EffHamiltonian(const VectorXd& dotx, MatrixXcd *ptr_res) {

    assert(is_setup_);
    
    complex<double> ii(0.0, 1.0);
    MatrixXcd  M(num_,num_);
    MatrixXcd& res(*ptr_res);
    this->Hamiltonian(id_, &res);
    for(int iop = 0; iop < (int)ops_opt_.size(); iop++) {
      basis_->Matrix(id_, ops_opt_[iop], &M);
      res += -ii*M*dotx[iop];
    }
  }
  double DySetPoly::Norm2() const {
    MatrixXcd S(num_,num_);
    basis_->Matrix(id_, id_, &S);
    double norm2 = (c_.dot(S*c_)).real();
    return norm2;
  }
  void DySetPoly::UpdateBasis() {

    assert(is_setup_);

    for(int A = 0; A < basis_->num(); A++) {
      basis_->ref_Rs()(A) = q0_;
      basis_->ref_Ps()(A) = p0_;
      basis_->ref_gs()(A) = complex<double>(gr0_, gi0_);
    }
    basis_->SetUp();
  }
  void DySetPoly::DumpCon(int it, string prefix) {
    Con& con = Con::getInstance();
    
    this->UpdateBasis();
    basis_->DumpCon(it);

    if(it==0) {
      if(xs_.size()!=0) {
	con.write_f1("_x", 0, xs_);
      }
    }
    
    con.write_f1(prefix+"c_re", it, c_.real());
    con.write_f1(prefix+"c_im", it, c_.imag());

    con.write_f1(prefix+"clam_re", it, Clam_.real());
    con.write_f1(prefix+"clam_im", it, Clam_.imag());

    if(xs_.size()>0) {
      VectorXcd ys(xs_.size());;
      basis_->At(id_, c_, xs_, &ys);
      con.write_f1(prefix+"psi_re", it, ys.real());
      con.write_f1(prefix+"psi_im", it, ys.imag());
      for(int lam = 0; lam < num_; lam++) {
	VectorXcd c = U_.col(lam);
	basis_->At(id_, c, xs_, &ys);
	string lbl = prefix + "phi" + to_string(lam);
	con.write_f1(lbl + "_re", it, ys.real());
	con.write_f1(lbl + "_im", it, ys.imag());
      }      
    }

    double norm2 = Norm2();
    con.write_f("norm2", it, norm2);
  }

  /*
  DyPsa::DyPsa(int max_path, Operator *pot, const VectorXi& ns, string type_gauss):
    max_path_(max_path), ns_(ns)
    type_gauss_(type_gauss), type_eomslow_("qhamilton"), m_(1.0), pot_(NULL),
    dys_(max_path_), xs_(0) {

    assert(ns.size()>0);
    
    for(int K = 0; K < max_apth_; K++) {
      auto ptr = new DySetPoly(pot, ns, type_gauss);
      dys_[K] = ptr;
      psa_data_[ptr] = new PsaData(ns.size());
    }

  }
  void DyPsa::SetUp() { }
  void DyPsa::Update(double dt) {    
    for(auto ptr_dy : dys_) {
      PsaData *data = psa_data_[ptr_dy];
      if(data->mode_ == kNone) {
      } else if(data->mode_ == kStr) {
	ptr_dy->Update(dt);
      } else if(data->mode_ == kPsa) {
	throw runtime_error("not impl");
      } else {
	throw runtime_error("not impl");
      }
    }
  }
  void DyPsa::DumpCon(int it) {
    for(int K = 0; K < max_path_; K++) {
      string prefix = to_string(K) + "_";
      dys_[K]->DumpCon(it, prefix);
    }
  }
  void DyPsa::AddPath(DySetPoly *dy, double q0, double p0, complex<double> g0, const VectorXcd& c) {
    assert(psa_datas_.find(dy) != psa_data_.end());
    assert(psa_datas_[dy]->mode_ == kNone);
    
    dy->m = m_;
    dy->q0_ = q0;
    dy->p0_ = p0;
    dy->gr0_ = g0.real();
    dy->gi0_ = g0.imag();
    dy->c_ = c;
    dy->type_eomslow_ = type_emoslow_;
    dy->SetUp();

    psa_datas_[dy]->mode_ = kSet;

  }
  void DyPsa::Branch(DySetPoly* path) {
    assert(dys_.find(path)!=dys_.end());

    auto q0 = path->q0_;
    auto p0 = path->q0_;
    complex<double> g0(path->gr0_, gi0= path->gi0_);
    
    int n(path->num_);
    VectorXd w(n);
    MatrixXcd H(n,n), S(n,n);

    MatrixXcd& U = path->U_;
    
    path->c_ = U.col(0);
    
    for(int lam = 1; lam < n; lam++) {
      this->AddPath(q0, p0, g0, U.col(lam));
    }
    
  }
  void DyPsa::set_xs(int nx, double x0, double x1) {
    xs_ = VectorXd::LinSpaced(nx, x0, x1);
  }  
  */
}
