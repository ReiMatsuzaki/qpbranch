#include <math.h>

#include <qpbranch/mathplus.hpp>

namespace qpbranch {
  using namespace std;
  using namespace Eigen;
  using std::complex;
  using std::runtime_error;
  using std::invalid_argument;
  void IntGto2N(int maxN, complex<double> a, VectorXcd *res) {
    assert(maxN >= 0);
    assert(a.real() > 0);
    assert(res->size() >= maxN+1);

    (*res)(0) = sqrt(M_PI/a);
    if(maxN == 0)
      return;

    for(int n = 1; n <= maxN; n++) {
      (*res)[n] = (2*n-1.0)/(2.0*a) * (*res)[n-1];
    }
    
  }
  void GammaNhalf(int maxN, VectorXd *res) {
    assert(maxN >= 0);
    assert(res->size() >= maxN+1);
    
    (*res)(0) = sqrt(M_PI);
    for(int n = 0; n < maxN; n++) {
      (*res)(n+1) = (n+0.5) * (*res)(n);
    }
  }  
  void IntGtoShift(int maxN, complex<double> a, complex<double> b, complex<double> w,
		   VectorXcd *res) {
    // The integrations
    //      J_n = Int_{-oo}^{+oo} dq (q-q0)^n Exp[-a (q-q0)^2 - b q^2]
    // satisfy following recursion relation
    //      (1+b/a) J_n = (n-1)/(2a) J_{n-2} - b.q_0/a. J_{n-1}.
    // See 2018/4/16/qpbranch for detail
    assert(maxN>=0);
    assert(a.real()>0);
    assert(b.real()>0);
    assert(res->size() >= maxN+1);

    VectorXcd& vec(*res);

    vec(0) = sqrt(M_PI) * exp(-(a*b*w*w)/(a+b)) / sqrt(a+b);
    if(maxN==0) return;

    vec(1) = sqrt(M_PI) * b*(-w) * exp(a*w*w*(a/(a+b)-1.0))/pow(a+b,1.5);
    if(maxN==1) return;

    for(int n = 2; n <= maxN; n++) 
      vec(n) = ((n-1.0)/(2.0*a)*vec(n-2) - b*w/a*vec(n-1)) / (1.0+b/a);
    
  }
  void DwIntGauss(int maxN, complex<double> a, complex<double> b, complex<double> w,
		  VectorXcd *ptr_res) {
    // Consider seriese of function
    //      I(N,n) = (d/dw)^N Int[(q-w)^n Exp[-g(q-w)^2 - bq^2 ].
    // I(N,n) satisfy following recursion relation
    //      I(N,n) = (d/dw)^{N-1} (-n). Int[(q-w)^{n-1} Exp[-g(q-w)^2 - bq^2 ]
    //              +(d/dw)^{N-1} 2g  . Int[(q-w)^{n+1} Exp[-g(q-w)^2 - bq^2 ]
    //             = -n.I(N-1,n-1) + 2g.I(N-1,n+1)
    VectorXcd& res(*ptr_res);
    assert(maxN>=0);
    assert(res.size()>=maxN+1);

    MatrixXcd I(maxN+1, maxN+1);
    IntGtoShift(maxN, a, b, w, &res);
    I.row(0) = res;
    for(int N = 1; N<=maxN; N++) {
      for(int n = 0; n<maxN; n++) {
	complex<double> cumsum(0);
	if(n!=maxN)
	  cumsum += 2.0*a*I(N-1,n+1);
	if(n!=0)
	  cumsum += -double(n)*I(N-1,n-1);
	I(N,n) = cumsum;
      }
    }
    res = I.col(0);
  }  
}
