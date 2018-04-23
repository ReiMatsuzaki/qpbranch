#ifndef EIGENPLUS_HPP_
#define EIGENPLUS_HPP_
#include <Eigen/Core>
namespace qpbranch {
  using Eigen::VectorXd;
  using Eigen::MatrixXcd;
  void zhegv(const MatrixXcd& H, const MatrixXcd& S, MatrixXcd *U, VectorXd *w);
}
#endif
