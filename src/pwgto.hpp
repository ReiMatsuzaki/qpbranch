#ifndef PWGTO_HPP_
#define PWGTO_HPP_

#include <iostream>
#include <Eigen/Core>
#include <vector>
#include <map>
#include <boost/multi_array.hpp>

#include "operator.hpp"

namespace qpbranch {
  using namespace std;
  using namespace Eigen;
  using boost::multi_array;

  // gauss functions processes
  void prod_gauss(complex<double> gA, double RA, double PA, double tA,
		  complex<double> gB, double RB, double PB, double tB,
		  complex<double> *gAB, complex<double> *RAB, complex<double> *eP);
  void hermite_coef_d(complex<double> gP, complex<double> wP, double RA, double RB,
		      int maxnA, int maxnB, multi_array<complex<double>, 3> *res);
  complex<double> hermite_coef_d_0(complex<double> gP, complex<double> wP, double RA, double RB,
				   int nA, int nB, int Nk);

  /*
  class PlaneWaveGtoData {
  public:
    int num_;
    VectorXi ns_;
    VectorXcd gs_;
    VectorXd Rs_, Ps_;
    PlaneWaveGtoData(int num);
    ~PlaneWaveGtoData();
  };

  class OpPwgto {
  protected:
    Operator *op_;
    PlaneWaveGtoData *pwgto_;
  public:
    OpPwgto(Operator *op, PlaneWaveGtoData *pwgto_);
  };
  class OpPwgtoBasic : public OpPwgto {
    vector<pair<int, complex<double>>> ncs_;
  };
  
  enum Operator {
    kNone, kOp0, kOp1, kOp2, kOpP1, kOpP2, kOpdR, kOpdP, kOpdgr, kOpdgi
  };
  */

  class OpBasis {
  public:
    int num_;    
    VectorXcd cs_;
    VectorXi  ns_;
    OpBasis(int num);
  };
 
  class PlaneWaveGto {
  public:    
    // - data size -
    int num_, nop_;
    // - variable -
    VectorXi ns_;
    VectorXcd gs_;
    VectorXd Rs_, Ps_;
    vector<Operator*> ops_;
    // - intermediate -
    map<Operator*, vector<OpBasis*> > op_basis_;
    VectorXd Ns_;
    VectorXi maxn_;
    MatrixXcd gAB_, eAB_, hAB_,  RAB_;
    multi_array<multi_array<complex<double>,3>*,2> *d_;
    // - method -
    PlaneWaveGto(const VectorXi& ns, const vector<Operator*>& ops);
    virtual ~PlaneWaveGto();
    void setup();
    void matrix(Operator *ibra, Operator *iket, MatrixXcd *res);
    void at(Operator *iop, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
    // - operator specific operation -
    void new_op(OperatorId *op);
    void new_op(OperatorRn *op);
    void new_op(OperatorPn *op);
    void setup_op(OperatorId *op);
    void setup_op(OperatorRn *op);
    void setup_op(OperatorPn *op);
  protected:
    inline complex<double> getd(int A, int B,int na,int nb,int Nk) {
      return (*(*d_)[A][B])[na][nb][Nk];
    }
  };


  
  /**
     Plane Wave Gauss Type Orbitals by the McMurchie-Davidson Recursion formula.

  class PlaneWaveGtoMDR : public PlaneWaveGto {
  private:
    // - intermediate -
    
  public:
    // - method -
    PlaneWaveGtoMDR(const VectorXi& ns, const vector<Operator>& ops);
    ~PlaneWaveGtoMDR();
    void setup();
    void overlap(Operator ibra, Operator iket, MatrixXcd *res);
    void gausspot(Operator ibra, Operator iket, complex<double> b, MatrixXcd *res);
    void at(Operator iop, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
    inline complex<double> getd(int A, int B,int na,int nb,int Nk) {
      return (*(*d_)[A][B])[na][nb][Nk];
    }
  protected:    
  };
  */
  
  /**
     all gauss parameters are same

  class PlaneWaveGto1Center : public PlaneWaveGto {
  public:
    int maxmaxn_;
    VectorXcd ints_; // ints_[n] : Int[q^n Exp[-2Re[g]q^2]]
    PlaneWaveGto1Center(const VectorXi& ns, const vector<Operator>& ops);
    ~PlaneWaveGto1Center();
    void setup();
    void overlap(Operator ibra, Operator iket, MatrixXcd *res);
    void at(Operator iop, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
  };
   */

}

#endif  // FOO_BAR_BAZ_H_
