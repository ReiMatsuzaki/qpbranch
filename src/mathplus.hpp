#ifndef MATHPLUS_H_
#define MATHPLUS_H_
#include <complex>
#include <Eigen/Core>

namespace qpbranch {
  using Eigen::VectorXcd;
  using std::complex;
  void gtoint2n(int maxn, complex<double> z, VectorXcd *res);
  void gammaint_nhalf(int maxn, VectorXcd *res);
}

#endif
