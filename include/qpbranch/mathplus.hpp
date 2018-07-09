#ifndef MATHPLUS_H_
#define MATHPLUS_H_
#include <complex>
#include <Eigen/Core>

namespace qpbranch {
  using Eigen::VectorXd;
  using Eigen::VectorXcd;  
  using std::complex;

  // give a combination
  int Fact(int n);
  // give number of combination nCk
  int Comb(int n, int k);
  // Give an array of integration
  //       { Int_{-oo}^{+oo} x^n Exp[-zx^2] dx   |   n = 0,...,maxn }.
  void IntGto2N(int maxN, complex<double> z, VectorXcd *res);
  // Give an array of gamma function with half integer
  //       { Gamma(n+1/2)  |  n = 0,...,maxn}
  void GammaNhalf(int maxN, VectorXd *res);
  // Give an array of integration
  //       { Int_{-oo}^{+oo} (q-q0)^n Exp[-z(q-w)^2 - b q^2] dq  |  n = 0,...,nmax }
  void IntGtoShift(int maxN, complex<double> a, complex<double> b, complex<double> w,
		   VectorXcd *res);
  // Give an array of derivative of gauss function integration
  //       { (d/dw)^N Int_{-oo}^{+oo} Exp[-a(q-w)^2 - b q^2]  | N = 0,...,maxN }
  void DwIntGauss(int maxN, complex<double> a, complex<double> b, complex<double> w,
		  VectorXcd *ptr_res);
}

#endif
