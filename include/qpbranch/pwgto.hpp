#ifndef PWGTO_HPP_
#define PWGTO_HPP_

#include <iostream>
#include <Eigen/Core>
#include <vector>
#include <map>
#include <boost/multi_array.hpp>

#include "basis_set.hpp"
#include "operator.hpp"

namespace qpbranch {
  
  using namespace std;
  using namespace Eigen;
  using boost::multi_array;

  // gauss functions processes
  void HermiteCoefD(complex<double> gP, complex<double> wP, double RA, double RB,
		    int maxnA, int maxnB, multi_array<complex<double>, 3> *res);
  complex<double> HermiteCoefD0(complex<double> gP, complex<double> wP, double RA, double RB,
				int nA, int nB, int Nk);

  // id used in derivative of PWGTO variables.
  const int kIdDR  = 1;
  const int kIdDP  = 2;
  const int kIdDgr = 3;
  const int kIdDgi = 4;

  // forward declaration
  class OpBuf;

  // Plane Wave Gauss Type Orbital set.
  // A th basis GA are defined as follow
  //       GA(q) = NA (q-RA)^nA Exp[-gA(q-RA)^2 + iPA(q-RA)]
  // Matrix elements are evaluated with extended McMurcy-Davidson recursion formula.
  // See Tachikawa et al.
  class Pwgto : public BasisSet {
  protected:
    // data size
    int num_, nop_;
    // const
    VectorXi ns_;
    vector<Operator*> ops_;
    // variable    
    VectorXcd gs_;
    VectorXd Rs_, Ps_;
    // intermediate
    bool is_setup_;
    map<Operator*, OpBuf*> buffer_map_;
    VectorXd Ns_;
    VectorXi maxn_;
    multi_array<multi_array<complex<double>,3>*,2> *d_;
  public:
    MatrixXcd gAB_, eAB_, hAB_,  RAB_;    
  public:
    // Main
    Pwgto(const VectorXi& ns, const vector<Operator*>& ops);
    ~Pwgto() {}
    int num() const { return num_; }
    int nop() const { return nop_; }
    const vector<Operator*>& ops() const { return ops_; }
    const VectorXi& ns() const { return ns_; }
    const VectorXcd& gs() const { return gs_; }
    const VectorXd& Rs() const { return Rs_; }
    const VectorXd& Ps() const { return Ps_; }
    const VectorXd& Ns() const { return Ns_; }
    VectorXcd& ref_gs() { return gs_; }
    VectorXd&  ref_Rs() { return Rs_; }
    VectorXd&  ref_Ps() { return Ps_; }
    void SetUp();
    int size() const { return num_; }
    void Matrix(Operator *ibra, Operator *iket, MatrixXcd *res);
    void At(Operator *iop, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);    
    void DumpCon(int it);
    inline complex<double> get_d(int A, int B,int na,int nb,int Nk) {
      return (*(*d_)[A][B])[na][nb][Nk];
    }
  };
  
}

#endif  // FOO_BAR_BAZ_H_
