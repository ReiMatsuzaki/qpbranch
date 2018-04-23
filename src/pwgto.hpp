#include <iostream>
#include <Eigen>
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
		
  enum Operator {
    kNone, kOp0, kOp1, kOp2, kOpR1, kOpR2
  };
  
  class OpBasis {
  public:
    int num_;
    MatrixXcd cs_; // cs[iop,A]
    MatrixXi  ns_; // ns[iop,A]    
    OpBasis(int np, int nb);
    static int maxn(Operator op, int n) {
      if(op==kOp0)
	return n;
      else if(op==kOp1)
	return n+1;
      else if(op==kOp2)
	return n+2;
      else
	assert("unsupported operator");
    }
  };
  
  class GaussBasis {
  public:
    GaussBasis(const VectorXi& ns, const vector<Operator>& ops);
    // - data size -
    int num_;
    // - variable -
    VectorXi ns_;
    VectorXcd gs_;
    VectorXd Rs_, Ps_;
    vector<Operator> ops_;
    // - operatorator -
    map<Operator, OpBasis*> op_basis_;
    VectorXd Ns_;
    MatrixXcd gAB_, eAB_, hAB_,  RAB_;
    //    multi_array<complex<double>,5> *d_; // d(A,B,nA,nB,N)
    VectorXi maxn_;
    multi_array<multi_array<complex<double>,3>*,2> *d_;
    virtual ~GaussBasis();
    virtual void setup();
    virtual void overlap(Operator ibra, Operator iket, MatrixXcd *res);
    virtual void at(Operator iop, const VectorXcd& cs, const VectorXd& xs, MatrixXcd *res);
    inline complex<double> getd(int A, int B,int na,int nb,int Nk) {
      return (*(*d_)[A][B])[na][nb][Nk];
    }
  protected:
    void setup_normalize();
    void setup_combination();
    void setup_operator();
    
  };
  
  class PlaneWaveGTO : public GaussBasis {
  public:    
    PlaneWaveGTO(const VectorXi& ns, const vector<Operator>& ops);
    ~PlaneWaveGTO();
    void setup();
    void overlap(Operator ibra, Operator iket, MatrixXcd *res);
    void at(Operator iop, const VectorXcd& cs, const VectorXd& xs, MatrixXcd *res);
  };
}
