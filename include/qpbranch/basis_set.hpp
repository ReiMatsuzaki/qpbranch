#ifndef BASIS_SET_HPP_
#define BASIS_SET_HPP_

#include <Eigen/Core>

namespace qpbranch {

  // forward declaration
  class Operator;
  using Eigen::VectorXd;
  using Eigen::VectorXcd;
  using Eigen::MatrixXcd;
  
  // interface
  class BasisSet {
  public:
    virtual ~BasisSet() {}
    virtual void setup() = 0;    
    virtual int size() const = 0;
    virtual void matrix(Operator *ibra, Operator *iket, MatrixXcd *res) = 0;
    virtual void at(Operator *iop, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res) = 0;    
  };
  
}

#endif

