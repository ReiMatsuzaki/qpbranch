#ifndef DY_AADF_HPP_
#define DY_AADF_HPP_

#include "operator.hpp"
#include "pwgto.hpp"

namespace qpbranch {

  using std::string;
  
  // Averaged Action Decomposed Function implementation.
  // Wave function is  represented as
  //        \Psi(q,t) = F(q,t)Exp[ iS(q,t)/\hbar], 
  // where the action S describe slow motion and the amplitude F does fast motion.
  // In this class, F are approximated as sum of gauss function
  //        F(q,t) = \sum_i ci Gi
  //        Gi   = N[i] (q-q0)^ni Exp[ -(1/4\alpha) (q-q0)^2]
  // and S are approximated 2nd order taylor expansion
  //        S(q,t) = p0(q-q0) - \beta(q-q0)^2
  class DyAadf {
  public:
    int num_;
    // variable
    double q0_, p0_, alpha_, beta_, theta_;
    VectorXcd c_;
    // const
    double m_;
    // options
    string type_gauss_;
    /*VectorXd xs_;*/
    // intermediate    
    bool is_setup_;
    OperatorPot *pot_;
    Operator *id_, *p2_, *DR_, *DP_, *Dgr_, *Dgi_, *r1_, *r2_;
    Pwgto *basis_;
    // Main
    DyAadf(OperatorPot *pot, const VectorXi& ns, string type_gauss);
    void SetUp();
    void Update(double dt);
    // Accessor
    /*void set_xs(int nx, double x0, double x1) { xs_=VectorXd::LinSpaced(nx,x0,x1); }*/
    // calc
    void DotxQhamilton(VectorXd *res);
    void DotxQuantum(bool is_tdvp, VectorXd *res);
    // calculate AADF hamiltonian for fast motion F.
    //      H^F_{ij} = <Gi| (1/2m) p^2 + V | Gj>
    void HamiltonianForF(Operator *op, MatrixXcd *res);
    // calculate AADF hamiltonian matrix
    //     H_{ij} = (1/2)<Gi|(dS/dq)^2|Gj> + <Gi|H^F|Gj>
    //            = (1/2)<Gi|(p0-2\beta(q-q0))^2|Gj> + <Gi|H^F|Gj>
    void Hamiltonian(Operator *op_bra, MatrixXcd *res);
    void EffHamiltonian(const VectorXd& dotx, MatrixXcd *res);
    double Norm2() const;
    // shift i th variable with dx and update basis. Integer i and variable are connected with
    // following table
    // | i | var   |
    // +---+-------+
    // | 0 | q0    |
    // | 1 | p0    |
    // | 2 | alpha |
    // | 3 | beta  |
    void ShiftVar(int i, double dx);
    void UpdateBasis();
    void DumpCon(int it, string prefix);
  };
}

#endif
