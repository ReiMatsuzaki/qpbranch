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
  void hermite_coef_d(complex<double> gP, complex<double> wP, double RA, double RB,
		      int maxnA, int maxnB, multi_array<complex<double>, 3> *res);
  complex<double> hermite_coef_d_0(complex<double> gP, complex<double> wP, double RA, double RB,
				   int nA, int nB, int Nk);

  // id used in derivative of PWGTO variables.
  const int kIdDR  = 1;
  const int kIdDP  = 2;
  const int kIdDgr = 3;
  const int kIdDgi = 4;

  // forward declaration
  class OpBuf;
  
  class Pwgto : public BasisSet {
  public:    
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
    MatrixXcd gAB_, eAB_, hAB_,  RAB_;
    multi_array<multi_array<complex<double>,3>*,2> *d_;
    // Main
    Pwgto(const VectorXi& ns, const vector<Operator*>& ops);
    ~Pwgto();
    void setup();
    int size() const { return num_; }
    void matrix(Operator *ibra, Operator *iket, MatrixXcd *res);
    void at(Operator *iop, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
    inline complex<double> getd(int A, int B,int na,int nb,int Nk) {
      return (*(*d_)[A][B])[na][nb][Nk];
    }
    void con(int it);
  };
  
}

#endif  // FOO_BAR_BAZ_H_
