#include <iostream>
#include <math.h>
#include "mathplus.hpp"

namespace qpbranch {
  using namespace std;
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
}
