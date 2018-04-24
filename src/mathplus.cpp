#include <iostream>
#include <math.h>
#include "mathplus.hpp"
#include <Eigen/Core>

namespace qpbranch {
  using namespace std;
  using namespace Eigen;
  using std::complex;
  using std::runtime_error;
  using std::invalid_argument;
  void gtoint2n(int maxn, complex<double> z, VectorXcd *res) {
    /* 
       gives the integrations : { Int_{-oo}^{+oo} x^{2n}Exp[-zx^2] dx | n = 0,...,maxn}
    */

    assert(maxn >= 0);
    assert(res->size() >= maxn+1);

    (*res)(0) = sqrt(M_PI/z);
    if(maxn == 0)
      return;

    for(int n = 1; n <= maxn; n++) {
      (*res)[n] = (2*n-1.0)/(2.0*z) * (*res)[n-1];
    }
    
  }
  void gammaint_nhalf(int maxn, VectorXcd *res) {
    /*
      gives {Gamma(n+1/2) | n = 0,...,maxn}
    */

    assert(maxn >= 0);
    assert(res->size() >= maxn+1);
    
    (*res)(0) = sqrt(M_PI);
    for(int n = 0; n < maxn; n++) {
      (*res)(n+1) = (n+0.5) * (*res)(n);
    }
  }  
  void gtointn_shift(int maxn, complex<double> a, complex<double> b, complex<double> q0, VectorXcd *res) {
    /*
      compute itegrations
      .      J_n = Int_{-oo}^{+oo} dq (q-q0)^n Exp[-a (q-q0)^2 - b q^2]
      by the recursion formula
      .      (1+b/a) J_n = (n-1)/(2a) J_{n-2} - b.q_0/a. J_{n-1}.
      See 2018/4/16/qpbranch for detail
     */
    assert(maxn>=0);

    VectorXcd& vec(*res);

    vec(0) = sqrt(M_PI) * exp(-(a*b*q0*q0)/(a+b)) / sqrt(a+b);
    if(maxn==0) return;

    vec(1) = sqrt(M_PI) * b*(-q0) * exp(a*q0*q0*(a/(a+b)-1.0))/pow(a+b,1.5);
    if(maxn==1) return;

    for(int n = 2; n <= maxn; n++) 
      vec(n) = ((n-1.0)/(2.0*a)*vec(n-2) - b*q0/a*vec(n-1)) / (1.0+b/a);
  }
  void dwn_gaussint_shift(int maxN, complex<double> b, complex<double> gAB, complex<double> wAB,
			  VectorXcd *ptr_res) {
    /** calculation for integration of Hermite polynomial and gauss.
	.   J(N)   = (d/dw)^N Int[ Exp[-g(q-w)^2 - bq^2 ] ]
	To calculate this seriese of integration, we define following integrations
	.   I(N,n) = (d/dw)^N Int[(q-w)^n     Exp[-g(q-w)^2 - bq^2 ]
	I(N,n) satisfie the recursion relation below
	.   I(N,n) = (d/dw)^{N-1} (-n). Int[(q-w)^{n-1} Exp[-g(q-w)^2 - bq^2 ]
	.           +(d/dw)^{N-1} 2g  . Int[(q-w)^{n+1} Exp[-g(q-w)^2 - bq^2 ]
	.          = -n.I(N-1,n-1) + 2g.I(N-1,n+1)
    */
    VectorXcd& res(*ptr_res);
    assert(maxN>=0);
    assert(res.size()>=maxN+1);

    MatrixXcd I(maxN+1, maxN+1);
    gtointn_shift(maxN, gAB, b, wAB, &res);
    I.row(0) = res;
    for(int N = 1; N<=maxN; N++) {
      for(int n = 0; n<maxN; n++) {
	complex<double> cumsum(0);
	if(n!=maxN)
	  cumsum += 2.0*gAB*I(N-1,n+1);
	if(n!=0)
	  cumsum += -double(n)*I(N-1,n-1);
	I(N,n) = cumsum;
      }
    }
    res = I.col(0);
    
  }  
}
