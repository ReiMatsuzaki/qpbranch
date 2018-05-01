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
    Operator *pot_, *id_, *p2_, *DR_, *DP_, *Dgr_, *Dgi_;
    vector<Operator*> ops_opt_;
    Pwgto *basis_;
    // Main
    DySetPoly(Operator *pot, const VectorXi ns, string type_gauss);
    void setup();
    void update(double dt);
    // Calc
    void calc_dotx_qhamilton(VectorXd *res);
    void calc_dotx_quantum(bool is_tdvp, VectorXd *res);
    void calc_H(Operator *op_bra, MatrixXcd *res);
    void calc_eff_H(const VectorXd& dotx, MatrixXcd *res);
    void update_basis();
  };

  /*
  class DyBranch {
  public:
    // size
    int maxpath_;
    // - variable -
    AadfBasis *basis_;
    VectorXd q0_, p0_, gr0_, gi0_;
    VectorXcd c_tot_;
    MatrixXcd c_;
    // - const -
    int nt_, n1t_;
    double m_, pot_v0_, pot_b_, dt_, dal_;
    // - options -
    string type_gauss_, type_hamiltonian_, type_eomslow_;
    int verbose_;
    // - intermediate -
    int nppath_;
    double dydt_;
    DyBranch(AadfBasis *basis, string type_gauss_, int maxpath);
    ~DyBranch();
    void setup();
    void update();
  };
  */
}

#endif
