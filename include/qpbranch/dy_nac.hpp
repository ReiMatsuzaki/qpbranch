#ifndef DY_NAC_HPP_
#define DY_NAC_HPP_

#include "operator.hpp"
#include "pwgto.hpp"

namespace qpbranch {

  using std::string;

  // Nuclear electron wave function propagation.
  // Wave function is assumed to
  //     \Psi(r,Q,t) = \chi(Q,t)\Phi(r,t;Q)
  //     \chi(Q,t)   = \sum_A D_A G_A(Q,t)
  //     \Phi(r,t;Q) = \sum_I C_I \Phi_I(r;Q)
  class DyNac {
    // size
    int numA_, numI_;
    // variable
    double q0_, p0_;
    std::complex<double> gamma0_;
    // const
    double m_;
    // options
    std::string type_gauss_;
    // intermediate
    bool is_setup_;
    std::vector<std::vector<OperatorPot*> > potHIJ_;
    std::vector<std::vector<OperatorPot*> > potXIJ_;
    Operator *id_, *p2_, *DR_, *DP_, *Dgr_, *Dgi_, *r1_, *r2_;
    pwgto *basis_;
    // Main
    DyNac(Operator)
  };
  
}
