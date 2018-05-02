#ifndef EIGENPLUS_HPP_
#define EIGENPLUS_HPP_
#include <Eigen/Core>
namespace qpbranch {
  using Eigen::VectorXd;
  using Eigen::MatrixXd;
  using Eigen::VectorXcd;
  using Eigen::MatrixXcd;
  // Solve generalized eigen value problem for hermitian complex matrix.
  void Zhegv(const MatrixXcd& H, const MatrixXcd& S, MatrixXcd *U, VectorXd *w);
  // Solve linear problem for complex matrix.
  void Zgesv(const MatrixXcd& M, const VectorXcd& b, VectorXcd *v);
  // Solve linear problem for real matrix
  void Dgesv(const MatrixXd&  M, const VectorXd&  b, VectorXd  *v);
  // From time derivative problem
  //        i S.dot{c} = H.c
  //        c(0) = c0
  // integrate c(dt).
  void IntetGdiag(const MatrixXcd& H, const MatrixXcd& S, double dt, VectorXcd *c);
}
#endif
