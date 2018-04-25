#ifndef PWGTO_HPP_
#define PWGTO_HPP_

#include <iostream>
#include <Eigen/Core>
#include <vector>
#include <map>
#include <boost/multi_array.hpp>

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
  class OperatorForGauss {
  public:

    static Operator* Op0();
    static Operator* Op1();
    static Operator* Op2();
    static Operator* OpP1();
    static Operator* OpP2();
    static Operator* OpdR();
    static Operator* OpdP();
    static Operator* Opdgr();
    static Operator* Opdgr();

    virtual int num_poly(int n);
    virtual int maxn(int n);
    virtual void calc_
  };
    */
  
  enum Operator {
    kNone, kOp0, kOp1, kOp2, kOpP1, kOpP2, kOpdR, kOpdP, kOpdgr, kOpdgi
  };
  
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
    vector<Operator> ops_;
    // - intermediate -
    map<Operator, vector<OpBasis*> > op_basis_;
    VectorXd Ns_;
    VectorXi maxn_;
    // - method -
    PlaneWaveGto(const VectorXi& ns, const vector<Operator>& ops);
    virtual ~PlaneWaveGto();
    virtual void setup() = 0;
    virtual void overlap(Operator ibra, Operator iket, MatrixXcd *res) = 0;
    virtual void gausspot(Operator ibra, Operator iket, complex<double> b, MatrixXcd *res) = 0;
    virtual void at(Operator iop, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res) = 0;
  protected:
    void at_slow(Operator iop, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);    
    void setup_normalize();
    void setup_operator();
  };

  /**
     Plane Wave Gauss Type Orbitals by the McMurchie-Davidson Recursion formula.
   */
  class PlaneWaveGtoMDR : public PlaneWaveGto {
  public:
    // - intermediate -
    MatrixXcd gAB_, eAB_, hAB_,  RAB_;
    multi_array<multi_array<complex<double>,3>*,2> *d_;
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
    void setup_combination();    
  };

  /**
     all gauss parameters are same
   */
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

}

#endif  // FOO_BAR_BAZ_H_
