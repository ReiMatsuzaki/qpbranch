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

  // interface
  class BasisSet {
  public:
    virtual ~BasisSet() {}
    virtual void setup() = 0;    
    virtual int size() const = 0;
    virtual void matrix(Operator *ibra, Operator *iket, MatrixXcd *res) = 0;
    virtual void at(Operator *iop, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res) = 0;
    
  };

  class PlaneWaveGto;
  class Pwgto1c;
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
    virtual void at(PlaneWaveGto *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res)=0;
    virtual void setup(Pwgto1c *basis) = 0;
    virtual void matrix(OpBuf *opbra, Pwgto1c *basis, MatrixXcd *res) = 0;
    virtual void matrix(OpBufBasic *opket, Pwgto1c *basis, MatrixXcd *res) = 0;
    virtual void matrix(OpBufGausspot *opket, Pwgto1c *basis, MatrixXcd *res) = 0;
    virtual void at(Pwgto1c *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res)=0;
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
    void init_zero(int A, int num);
    virtual void setup(PlaneWaveGto *basis) = 0;    
    void matrix(OpBuf *opbra, PlaneWaveGto *basis, MatrixXcd *res);
    void matrix(OpBufBasic *opket, PlaneWaveGto *basis, MatrixXcd *res);
    void matrix(OpBufGausspot *opket, PlaneWaveGto *basis, MatrixXcd *res);
    void at(PlaneWaveGto *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
    virtual void setup(Pwgto1c *basis) = 0;    
    void matrix(OpBuf *opbra, Pwgto1c *basis, MatrixXcd *res);
    void matrix(OpBufBasic *opket, Pwgto1c *basis, MatrixXcd *res);
    void matrix(OpBufGausspot *opket, Pwgto1c *basis, MatrixXcd *res);
    void at(Pwgto1c *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
  };
  class OpBufId : public OpBufBasic {
  public:
    OpBufId(const VectorXi& ns);
    void setup(PlaneWaveGto *basis);
    void setup(Pwgto1c *basis);
  };
  class OpBufRn : public OpBufBasic {
  public:
    int n_;
    OpBufRn(const VectorXi& ns, int n);
    void setup(PlaneWaveGto *basis);
    void setup(Pwgto1c *basis);
  };
  class OpBufPn : public OpBufBasic {
  public:
    int n_;
    OpBufPn(const VectorXi& ns, int n);
    void setup(PlaneWaveGto *basis);
    void setup(Pwgto1c *basis);
  };
  class OpBufDa : public OpBufBasic {
  public:
    int id_;
    OpBufDa(const VectorXi& ns, int id);
    void setup(PlaneWaveGto *basis);
    void setup(Pwgto1c *basis);
  };
  class OpBufGausspot : public OpBuf {
  public:
    int num_;
    VectorXi ns_;    
    OperatorGausspot *op_;
    OpBufGausspot(const VectorXi& ns, OperatorGausspot *op);
    
    int maxn(int A);
    void setup(PlaneWaveGto *basis);
    void matrix(OpBuf *opbrat, PlaneWaveGto *basis, MatrixXcd *res);
    void matrix(OpBufBasic *opket, PlaneWaveGto *basis, MatrixXcd *res);
    void matrix(OpBufGausspot *opket, PlaneWaveGto *basis, MatrixXcd *res);
    void at(PlaneWaveGto *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
    void setup(Pwgto1c *basis);
    void matrix(OpBuf *opbra, Pwgto1c *basis, MatrixXcd *res);
    void matrix(OpBufBasic *opket, Pwgto1c *basis, MatrixXcd *res);
    void matrix(OpBufGausspot *opket, Pwgto1c *basis, MatrixXcd *res);
    void at(Pwgto1c *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
  };

  OpBuf* make_op_buff(PlaneWaveGto *basis, Operator *op);
  
  class PlaneWaveGto : public BasisSet {
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
    map<Operator*, OpBuf*> buffer_map_;
    VectorXd Ns_;
    VectorXi maxn_;
    MatrixXcd gAB_, eAB_, hAB_,  RAB_;
    multi_array<multi_array<complex<double>,3>*,2> *d_;
    // Main
    PlaneWaveGto(const VectorXi& ns, const vector<Operator*>& ops);
    ~PlaneWaveGto();
    void setup();
    int size() const { return num_; }
    void matrix(Operator *ibra, Operator *iket, MatrixXcd *res);
    void at(Operator *iop, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
    inline complex<double> getd(int A, int B,int na,int nb,int Nk) {
      return (*(*d_)[A][B])[na][nb][Nk];
    }
  };

  class Pwgto1c : public BasisSet {
  public:
    // data size
    int num_, nop_;
    // cont
    VectorXi ns_;
    vector<Operator*> ops_;
    // variable
    double R0_, P0_;
    complex<double> g0_;
    // intermediate
    VectorXd Ns_;
    int maxn_;
    map<Operator*, OpBuf*> buf_map_;
    VectorXcd gints2n_; // gints[n] = Int[q^{2n} Exp[-2Re[q0]q^2]]
    // Main
    Pwgto1c(const VectorXi& ns, double R0, double P0, complex<double> g0, const vector<Operator*>& ops);
    ~Pwgto1c();
    void setup();
    int size() const { return num_; }
    void matrix(Operator *opbra, Operator *opket, MatrixXcd *res);
    void at(Operator *op, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
  };
  
}

#endif  // FOO_BAR_BAZ_H_
