#ifndef DY_NAC_HPP_
#define DY_NAC_HPP_

#include <boost/multi_array.hpp>

#include "operator.hpp"
#include "pwgto.hpp"

namespace qpbranch {

  using std::string;
  typedef boost::multi_array<Operator*, 2> OpMat;

  // Nuclear electron wave function propagation.
  // Wave function is assumed to
  //     \Psi(r,Q,t) = sum_AI C_AI \chi_A(Q,t)\Phi_I(r,t;Q)
  class DyNac {
  public:    
    // size
    int numA_, numI_;
    // variable
    double q0_, p0_;
    std::complex<double> gamma0_;
    VectorXcd cAI_;
    // const
    VectorXi ns_;
    double m_;
    // basis functions
    Operator *id_, *r1_, *r2_, *p1_, *p2_, *DR_, *DP_, *Dgr_, *Dgi_;
    Pwgto *basis_;
  private:
    // options
    std::string type_gauss_;
    // intermediate
    bool is_setup_;
    OpMat opHeIJ_;
    OpMat opXkIJP_;    
  public:
    // Main
    DyNac(const OpMat& HeIJ, const OpMat& XkIJP, const VectorXi& ns,
	  string type_gauss, int norder, double dx);	  
    void SetUp();
    void Update(double dt);
    void UpdateBasis();
    // Calc
    void At(int I, const VectorXd& xs, VectorXcd *ys) const;
    // calculate electron-nuclear Hamiltonian
    //     H_{AI,BJ} := (2m)^{-1}(GA,P2.GB) d_{IJ} - i.m^{-1}(GA,XIJ.P.GB) + (GA,HIJ.GB)
    void Hamiltonian(Operator *op_bra, MatrixXcd *res);
    void DotxQhamilton(VectorXd *res);
    // private
    int idx(int A, int I) const { return A+numA_*I; }
  };
  
}
#endif
