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
    int num_;
    // variable     
    double q0_, p0_, gr0_, gi0_;
    VectorXcd c_;
    // const
    double m_;
    // basis functions
    Operator *pot_, *id_, *p2_, *DR_, *DP_, *Dgr_, *Dgi_;
    Pwgto *basis_;
    // options
    string type_gauss_, type_dotx_, type_intenuc_;
    // intermediate
    bool is_setup_;
    VectorXd w_;    
    MatrixXcd U_;
    VectorXcd Clam_;
    VectorXd xs_; // grid for dumping wave functions.
    // Main
    DySetPoly(Operator *pot, const VectorXi& ns, string type_gauss, int norder, double dx);
    void SetUp();
    void Update(double dt);
    // Accessor
    void set_xs(int nx, double x0, double x1) { xs_=VectorXd::LinSpaced(nx,x0,x1); }
    // Calc
    // calculate time derivative of non linear real parameters with Hamilton equation for
    // quantum Hamiltonian. The results are the same when single gaussian is applied.  
    void Dotx(VectorXd *res);
    void DotxQhamilton(VectorXd *res);
    void DotxQuantum(bool is_tdvp, VectorXd *res);
    void Hamiltonian(Operator *op_bra, MatrixXcd *res);
    void EffHamiltonian(const VectorXd& dotx, MatrixXcd *res);
    double Norm2() const;
    void UpdateBasis();
    void At(const VectorXd& xs, VectorXcd *ys) const;
    void DumpCon(int it, string prefix);
  };
  /*
  // Dynamics mode for each path
  enum DyMode { kNone, kSet, kPsa };

  // Path Space Averaging data
  class PsaData {
    DyMode mode_;
    VectorXd qlam_, plam_;
    VectorXcd clam_;
    PsaData(int n) : mode_(kNone), qlam_(n), plam_(n), clam_(n) {}
  };
  
  // Support branching dynamics
  class DyPsa {
  public:
    // data
    int max_path_;
    VectorXi ns_;
    string type_gauss_, type_eomslow_;
    double m_;
    Operator *pot_;
    // intermediate 
    vector<DySetPoly*> dys_;
    map<DySetPoly*, PsaData*> psa_datas_;
    VectorXd xs_; // grid for dumping wave functions.
    
    DyPsa(int max_path, Operator *pot, const VectorXi& ns, string type_gauss);
    void SetUp();
    void Update();
    void DumpCon();
    void AddPath(i, double q0, double p0, complex<double> g0);
    void Branch(DySetPoly* path);
    void set_xs(int nx, double x0, double x1);
  }
  */
}

#endif
