#ifndef DY_BRANCH_HPP_
#define DY_BRANCH_HPP_

#include "operator.hpp"
#include "pwgto.hpp"

namespace qpbranch {
  using std::string;

  // Nuclear wave packet dynamics with polynomial gaussian basis set.
  // all gaussian parameters are same.
  class DySetPoly {
  public:
    // size
    int num_, numopt_;
    // variable     
    double q0_, p0_, gr0_, gi0_;
    VectorXcd c_;
    // const
    double m_;
    // options
    string type_gauss_, type_eomslow_;
    // intermediate
    bool is_setup_;
    Operator *pot_, *id_, *p2_, *DR_, *DP_, *Dgr_, *Dgi_;
    vector<Operator*> ops_opt_;
    Pwgto *basis_;
    // Main
    DySetPoly(Operator *pot, const VectorXi ns, string type_gauss);
    void SetUp();
    void Update(double dt);
    // Calc
    void DotxQhamilton(VectorXd *res);
    void DotxQuantum(bool is_tdvp, VectorXd *res);
    void Hamiltonian(Operator *op_bra, MatrixXcd *res);
    void EffHamiltonian(const VectorXd& dotx, MatrixXcd *res);
    double Norm2() const;
    void UpdateBasis();
    void DumpCon(int it);
  };
}

#endif
