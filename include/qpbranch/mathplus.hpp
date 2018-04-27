#ifndef MATHPLUS_H_
#define MATHPLUS_H_
#include <complex>
#include <Eigen/Core>

namespace qpbranch {
  using Eigen::VectorXcd;
  using std::complex;
  void gtoint2n(int maxn, complex<double> z, VectorXcd *res);
  void gammaint_nhalf(int maxn, VectorXcd *res);
  void gtointn_shift(int maxn, complex<double> a, complex<double> b, complex<double> q0, VectorXcd *res);
  void dwn_gaussint_shift(int maxN, complex<double> b, complex<double> gAB, complex<double> wAB,
			  VectorXcd *ptr_res);
}

#endif
