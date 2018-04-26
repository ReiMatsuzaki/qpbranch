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

  // id used in derivative of PWGTO variables.
  const int kIdDR = 1;
  const int kIdDP = 2;
  const int kIdDgr = 3;
  const int kIdDgi = 4;

  class PlaneWaveGto;
  class PwgtoBufferBasic;
  class PwgtoBufferGausspot;

  // Buffer for each operator applying PlaneWaveGto.
  class PwgtoBuffer {
  public:
    virtual int maxn(int A) = 0;
    virtual void call_setup(PlaneWaveGto *basis) = 0;
    virtual void matrix(PwgtoBuffer *opbra, PlaneWaveGto *basis, MatrixXcd *res) = 0;
    virtual void matrix(PwgtoBufferBasic *opket, PlaneWaveGto *basis, MatrixXcd *res) = 0;
    virtual void matrix(PwgtoBufferGausspot *opket, PlaneWaveGto *basis, MatrixXcd *res) = 0;
  };
  // Buffer for basic case. op.phi_A can be represented by polynomial times gauss function
  //     op.phi_A = sum_i c[A][i] (q-qA)^n[A][i] Exp(-gA(q-qA)^2 + ipA(q-qA))
  class PwgtoBufferBasic : public PwgtoBuffer {
  public:
    int num_;
    vector<vector<pair<int,complex<double> > > > ncs_; // ncs_[A][i]
    int maxn(int A);
    virtual void call_setup(PlaneWaveGto *basis) = 0;
    void matrix(PwgtoBuffer *opbra, PlaneWaveGto *basis, MatrixXcd *res);
    void matrix(PwgtoBufferBasic *opket, PlaneWaveGto *basis, MatrixXcd *res);
    void matrix(PwgtoBufferGausspot *opket, PlaneWaveGto *basis, MatrixXcd *res);
  };
  class PwgtoBufferId : public PwgtoBufferBasic {
  public:
    void call_setup(PlaneWaveGto *basis);
  };
  class PwgtoBufferRn : public PwgtoBufferBasic {
  public:
    OperatorRn *op_;
    PwgtoBufferRn(OperatorRn *op);
    void call_setup(PlaneWaveGto *basis);
  };
  class PwgtoBufferPn : public PwgtoBufferBasic {
  public:
    OperatorPn *op_;
    PwgtoBufferPn(OperatorPn *op);
    void call_setup(PlaneWaveGto *basis);
  };
  class PwgtoBufferDa : public PwgtoBufferBasic {
  public:
    OperatorDa *op_;
    PwgtoBufferDa(OperatorDa *op);
    void call_setup(PlaneWaveGto *basis);
  };
  class PwgtoBufferGausspot : public PwgtoBuffer {
  public:
    OperatorGausspot *op_;
    VectorXi ns_;
    PwgtoBufferGausspot(OperatorGausspot *op, const VectorXi& ns);
    int maxn(int A);
    void call_setup(PlaneWaveGto *basis);
    void matrix(PwgtoBuffer *opbrat, PlaneWaveGto *basis, MatrixXcd *res);
    void matrix(PwgtoBufferBasic *opket, PlaneWaveGto *basis, MatrixXcd *res);
    void matrix(PwgtoBufferGausspot *opket, PlaneWaveGto *basis, MatrixXcd *res);
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
    map<Operator*, PwgtoBuffer*> buffer_map_;
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
    void new_op(OperatorDa *op);
    void new_op(OperatorGausspot *op);
    void setup_op(PwgtoBufferId *op);
    void setup_op(PwgtoBufferRn *op);
    void setup_op(PwgtoBufferPn *op);
    void setup_op(PwgtoBufferDa *op);
    void setup_op(PwgtoBufferGausspot *op);
    void calc_matrix(PwgtoBufferBasic *ibra, PwgtoBufferBasic *iket, MatrixXcd *res);
    void calc_matrix(PwgtoBufferBasic *ibra, PwgtoBufferGausspot *iket, MatrixXcd *res);
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
