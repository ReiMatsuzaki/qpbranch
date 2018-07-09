#ifndef DY_NAC_HPP_
#define DY_NAC_HPP_

#include <boost/multi_array.hpp>

#include "operator.hpp"
#include "pwgto.hpp"

namespace qpbranch {

  using std::string;
  typedef boost::multi_array<Operator*, 2> OpMat;
  typedef boost::multi_array<OperatorSpline*, 2> OpSpMat;

  // Nuclear electron wave function propagation.
  // Wave function is assumed to
  //     \Psi(r,Q,t) = \chi(Q,t) sum_AI C_AI \Phi_I(r,t;Q)
  //     \chi = delta(Q-Q')
  class DyNacDelta {
  public:    
    // size
    int numI_;
    // variable
    double q0_, p0_;
    VectorXcd cI_;
    // const
    double m_;
    // options
    std::string type_dotx_, type_intenuc_;
    double dx_;
  private:
    // intermediate
    bool is_setup_;
    OpSpMat opHeIJ_;
  public:
    // Main
    DyNacDelta(const OpSpMat& HeIJ);
    void SetUp();
    void Update(double dt);
    // Calc
    void Hamiltonian(MatrixXcd *res);
    void Dotx(VectorXd *res);
    void CalcProb(VectorXd *res);
    // IO
    void DumpCon(int it, string prefix="");
  };  
  
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
    // options
    std::string type_gauss_, type_dotx_, type_intenuc_;
    //  private:
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
    void EffHamiltonian(const VectorXd& dotx, MatrixXcd *res);
    void Overlap(MatrixXcd *res);
    void Dotx(VectorXd *res);
    void CalcProb(VectorXd *res);
    // IO
    void DumpCon(int it, string prefix="");
    // private
    int idx(int A, int I) const { return A+numA_*I; }
  };

  enum PathState {
    kSingle, // single path
    kBranch  // branching path
  };
  
  // Nuclear electron wave function propagation with branching
  // Each path is expressed by DyNac object.
  class DyBranch {
  public:
    int numA_, numI_, nlam_;
    double dx_;
    vector<DyNac*> paths_;
    map<DyNac*,PathState> modes_;
    map<DyNac*,VectorXd>  qlam_, plam_;
    //    map<DyNac*,VectorXcd> glam_;
    map<DyNac*,MatrixXcd> U_;
    map<DyNac*,VectorXcd> Clam_;
    map<DyNac*,complex<double> > Ctot_;
    // Main
    DyBranch(const OpMat& HeIJ, const OpMat& XkIJP, const VectorXi& ns, string type_gauss, int norder, double dx);    
    void SetUp();
    void Update(double dt);
    void DumpCon(int it, string prefix="");
    // calc
    void BeginBranch(int K);
    void EndBranch(int K);
  };
}
#endif
