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

  // Main
  DySetPoly::DySetPoly(Operator *pot, const VectorXi ns, string type_gauss) {

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
    
    basis_ = new Pwgto(ns, ops);
    
  }
  void DySetPoly::setup() {

    basis_->setup();
    MatrixXcd S(num_,num_);
    basis_->matrix(id_, id_, &S);

    double norm2 = real(c_.dot(S*c_));
    c_ *= 1/sqrt(norm2);

    is_setup_ = true;
  }
  void DySetPoly::update(double dt) {
    assert(is_setup_);

    this->update_basis();

    // slow part
    VectorXd dotx(numopt_);
    if(type_eomslow_=="qhamilton") {
      this->calc_dotx_qhamilton(&dotx);
    } else if(type_eomslow_=="tdvp") {
      this->calc_dotx_quantum(true, &dotx);					       
    } else if(type_eomslow_=="mp") {
      this->calc_dotx_quantum(false, &dotx);
    }

    // fast part
    MatrixXcd H(num_,num_), S(num_,num_);
    basis_->matrix(id_, id_, &S);
    this->calc_eff_H(dotx, &H);

    // propagate slow part
    q0_ += dotx(0) * dt;
    p0_ += dotx(1) * dt;
    if(type_gauss_=="thawed") {
      gr0_ += dotx(2) * dt;
      gi0_ += dotx(3) * dt;
    }

    // propagate fast part
    intet_gdiag(H, S, dt, &c_);
    
  }
  // Calc
  void DySetPoly::calc_dotx_qhamilton(VectorXd *res) {
    assert(is_setup_);
    
    MatrixXcd H(num_,num_);
    
    for(int iop = 0; iop < (int)ops_opt_.size(); iop++) {
      double c;
      Operator *op;
      if(ops_opt_[iop] == DR_) {
	op = DP_;
	c = 1;
      } else if(ops_opt_[iop] == DP_) {
	op = DR_;
	c = 1;
      } else if(ops_opt_[iop] == Dgr_) {
	op = Dgi_;
	c = 4.0*gr0_*gr0_;
      } else if(ops_opt_[iop] == Dgi_) {
	op = Dgr_;
	c = 4.0*gr0_*gr0_;
      } else {
	assert(false||"invalid op");
      }      
      this->calc_H(op, &H);
      (*res)[iop] = 2.0 * c * real(c_.dot(H*c_));
    }
    
  }
  void DySetPoly::calc_dotx_quantum(bool is_tdvp, VectorXd *res) {

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
      this->calc_H(ops[n], &Hn0[n]);
      for(int m = 0; m < 1+numopt_; m++) {	
	basis_->matrix(ops[n], ops[m], &(Snm[n][m]));
      }
    }

    for(int n = 0; n < numopt_; n++) {
      for(int m = 0; m < numopt_; m++) {
	auto t1 = c_.dot(Snm[n][m]*c_);
	v = Snm[id][m]*c_;
	zgesv(Snm[id][id], v, &v1);
	auto t2 = c_.dot(v1);
	C(n,m) = is_tdvp ? -imag(t1-t2) : real(t1-t2);
      }
      auto t1 = c_.dot(Hn0[n]*c_);
      v = Hn0[id]*c_;
      zgesv(Snm[id][id], v, &v1);
      v = Snm[n][0]*v;
      auto t2 = c_.dot(v1);
      y(n) = is_tdvp ? real(t1-t2) : imag(t1-t2);
    }

    dgesv(C, y, res);
    
  }
  void DySetPoly::calc_H(Operator *op_bra, MatrixXcd *res) {
    assert(is_setup_);
    
    MatrixXcd P2(num_, num_), V(num_,num_);    
    basis_->matrix(op_bra, p2_,  &P2);
    basis_->matrix(op_bra, pot_, &V);
    *res = 1.0/(2*m_) * P2 + V;
  }  
  void DySetPoly::calc_eff_H(const VectorXd& dotx, MatrixXcd *ptr_res) {

    assert(is_setup_);
    
    complex<double> ii(0.0, 1.0);
    MatrixXcd  M(num_,num_);
    MatrixXcd& res(*ptr_res);
    this->calc_H(id_, &res);
    for(int iop = 0; iop < (int)ops_opt_.size(); iop++) {
      basis_->matrix(id_, ops_opt_[iop], &M);
      res += -ii*M*dotx[iop];
    }
  }
  void DySetPoly::update_basis() {

    assert(is_setup_);

    for(int A = 0; A < basis_->num_; A++) {
      basis_->Rs_(A) = q0_;
      basis_->Ps_(A) = p0_;
      basis_->gs_(A) = complex<double>(gr0_, gi0_);
    }
    basis_->setup();
  }
  void DySetPoly::con(int it) {
    Con& con = Con::getInstance();
    
    this->update_basis();
    basis_->con(it);
    con.write_f1("c_re", it, c_.real());
    con.write_f1("c_im", it, c_.imag());
    
  }
}
