#ifndef PWGTO1C_HPP_
#define PWGTO1C_HPP_

#include <vector>
#include <map>
#include <Eigen/Core>
#include "basis_set.hpp"
#include "operator.hpp"

namespace qpbranch {

  using std::vector;
  using std::map;
  using Eigen::VectorXi;
  using Eigen::VectorXd;
  using Eigen::VectorXcd;
  using Eigen::MatrixXcd;

  // forward declaration
  class OpBuf;

  // Plane Wave Gauss Type Orbital basis set. Each basis share the same gauss parameters
  class Pwgto1c : public BasisSet {
  public:
    // data size
    int num_, nop_;
    // cont
    VectorXi ns_;
    vector<Operator*> ops_;
    // variable
    double R0_, P0_;
    complex<double> g0_;
    // intermediate
    VectorXd Ns_;
    int maxn_;
    map<Operator*, OpBuf*> buf_map_;
    VectorXcd gints2n_; // gints[n] = Int[q^{2n} Exp[-2Re[q0]q^2]]
    // Main
    Pwgto1c(const VectorXi& ns, double R0, double P0, complex<double> g0, const vector<Operator*>& ops);
    ~Pwgto1c();
    void SetUp();
    int size() const { return num_; }
    void Matrix(Operator *opbra, Operator *opket, MatrixXcd *res);
    void At(Operator *op, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
  };
  
} 

#endif
