#ifndef EIGENPLUS_HPP_
#define EIGENPLUS_HPP_
#include <Eigen/Core>
namespace qpbranch {
  using Eigen::VectorXd;
  using Eigen::MatrixXcd;
  void zhegv(const MatrixXcd& H, const MatrixXcd& S, MatrixXcd *U, VectorXd *w);
  void zgesv(const MatrixXcd& M, VectorXcd *v);
  void dgesv(const MatrixXd&  M, VectorXd  *v);
  void intet_gdiag(const MatrixXcd& H, const MatrixXcd& S, double dt, VectorXcd *c);
}
#endif
