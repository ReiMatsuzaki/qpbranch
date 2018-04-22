#include "pwgto.hpp"
#include "mathplus.hpp"
using namespace std;
  using boost::extents;
  using boost::array;


namespace qpbranch {

  OpBasis::OpBasis(int np, int nb) : num_(np), cs_(np, nb), ns_(np, nb) {}
  //  Polynomial::atp
  
  GaussBasis::GaussBasis(const VectorXi& ns, const vector<Operator>& ops):
    num_(ns.size()), ns_(ns), gs_(num_), Rs_(num_), Ps_(num_), ops_(ops),
    op_basis_(ops.size()), Ns_(num_),
    gAB(num_,num_), eAB(num_,num_), hAB(num_,num_), RAB(num_,num_) {
    int maxnd = 2;
    d_ = new multi_array<double,5>(extents[num_][num_][maxnd][maxnd][2*maxnd]);
  }
  GaussBasis::~GaussBasis() {}
  void GaussBasis::setup() {

    for(int A = 0; A < num_; A++) {
      gtoint2n(2*)
      Ns_[A] = 1.0;
    }
    
    for(auto it = ops_.begin(); it != ops_.end(); ++it) {
      if(*it==kOp0) {
	OpBasis *ptr = new OpBasis(1, num_);
	ptr->cs_.row(0) = Ns_;
	ptr->ns_.row(0) = ns_;
	op_basis_[kOp0] = ptr;
      } else {
	throw std::runtime_error("aa");
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
  
  PlaneWaveGTO::PlaneWaveGTO(const VectorXi& ns, const vector<Operator>& ops):
    GaussBasis(ns, ops) {}
  PlaneWaveGTO::~PlaneWaveGTO() {
  }
  void PlaneWaveGTO::setup() {}  
  void PlaneWaveGTO::overlap(Operator ibra, Operator iket, MatrixXcd *res) {
    *res = MatrixXcd::Zero(num_, num_);
  }
  void PlaneWaveGTO::at(Operator iop, const VectorXcd& cs, const VectorXd& xs, MatrixXcd *res) {
  }
}

