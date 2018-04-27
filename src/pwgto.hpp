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
  class OpBufBasic;
  class OpBufGausspot;

  // Buffer for each operator applying PlaneWaveGto.
  class OpBuf {
  public:
    virtual int maxn(int A) = 0;
    virtual void setup(PlaneWaveGto *basis) = 0;
    virtual void matrix(OpBuf *opbra, PlaneWaveGto *basis, MatrixXcd *res) = 0;
    virtual void matrix(OpBufBasic *opket, PlaneWaveGto *basis, MatrixXcd *res) = 0;
    virtual void matrix(OpBufGausspot *opket, PlaneWaveGto *basis, MatrixXcd *res) = 0;
  };
  // Buffer for basic case. op.phi_A can be represented by polynomial times gauss function
  //     op.phi_A = sum_i c[A][i] (q-qA)^n[A][i] Exp(-gA(q-qA)^2 + ipA(q-qA))
  class OpBufBasic : public OpBuf {
  public:
    int num_;
    vector<int> nums_;
    vector<VectorXi> ns_;   // ns_[A][i]
    vector<VectorXcd> cs_;  // cs_[A][i]
    /*    vector<vector<pair<int,complex<double> > > > ncs_; // ncs_[A][i] */
    OpBufBasic(int num);
    int maxn(int A);
    virtual void setup(PlaneWaveGto *basis) = 0;    
    void matrix(OpBuf *opbra, PlaneWaveGto *basis, MatrixXcd *res);
    void matrix(OpBufBasic *opket, PlaneWaveGto *basis, MatrixXcd *res);
    void matrix(OpBufGausspot *opket, PlaneWaveGto *basis, MatrixXcd *res);
    void init_zero(int A, int num);
  };
  class OpBufId : public OpBufBasic {
  public:
    OpBufId(int num);
    void setup(PlaneWaveGto *basis);
  };
  class OpBufRn : public OpBufBasic {
  public:
    int n_;
    OpBufRn(int num, int n);
    void setup(PlaneWaveGto *basis);
  };
  class OpBufPn : public OpBufBasic {
  public:
    int n_;
    OpBufPn(const VectorXi& ns, int n);
    void setup(PlaneWaveGto *basis);
  };
  class OpBufDa : public OpBufBasic {
  public:
    int id_;
    OpBufDa(const VectorXi& ns, int id);
    void setup(PlaneWaveGto *basis);
  };
  class OpBufGausspot : public OpBuf {
  public:
    VectorXi ns_;
    OperatorGausspot *op_;    
    OpBufGausspot(const VectorXi& ns, OperatorGausspot *op);
    
    int maxn(int A);
    void setup(PlaneWaveGto *basis);
    void matrix(OpBuf *opbrat, PlaneWaveGto *basis, MatrixXcd *res);
    void matrix(OpBufBasic *opket, PlaneWaveGto *basis, MatrixXcd *res);
    void matrix(OpBufGausspot *opket, PlaneWaveGto *basis, MatrixXcd *res);
  };

  OpBuf* make_op_buff(PlaneWaveGto *basis, Operator *op);
  
  class PlaneWaveGto {
  public:    
    // data size
    int num_, nop_;
    // variable
    VectorXi ns_;
    VectorXcd gs_;
    VectorXd Rs_, Ps_;
    vector<Operator*> ops_;
    // intermediate
    map<Operator*, OpBuf*> buffer_map_;
    VectorXd Ns_;
    VectorXi maxn_;
    MatrixXcd gAB_, eAB_, hAB_,  RAB_;
    multi_array<multi_array<complex<double>,3>*,2> *d_;
    // Main
    PlaneWaveGto(const VectorXi& ns, const vector<Operator*>& ops);
    virtual ~PlaneWaveGto();
    void setup();
    void matrix(Operator *ibra, Operator *iket, MatrixXcd *res);
    void at(Operator *iop, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
    // Calculation
    void calc_matrix(OpBufBasic *ibra, OpBufBasic *iket, MatrixXcd *res);
    void calc_matrix(OpBufBasic *ibra, OpBufGausspot *iket, MatrixXcd *res);
  protected:
    inline complex<double> getd(int A, int B,int na,int nb,int Nk) {
      return (*(*d_)[A][B])[na][nb][Nk];
    }
  };
  
}

#endif  // FOO_BAR_BAZ_H_
