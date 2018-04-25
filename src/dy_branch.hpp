#ifndef DY_BRANCH_HPP_
#define DY_BRANCH_HPP_

namespace {
  using std::string;
  class DyBranch {
  public:
    // - variable -
    GaussBasis *basis_;
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
    DyBranch
  };
}

#endif
