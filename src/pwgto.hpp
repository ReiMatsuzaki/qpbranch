#include <iostream>
#include <Eigen>
#include <vector>
#include <boost/multi_array.hpp>

namespace qpbranch {
  using namespace std;
  using namespace Eigen;
  using boost::multi_array;

  enum Operator {
    kNone, kOp0, kOp1, kOpR1, kOpR2
  };

  class OpBasis {
  public:
    int num_;
    MatrixXcd cs_; // cs[iop,A]
    MatrixXi  ns_; // ns[iop,A]    
    OpBasis(int np, int nb);
    //    void at(const VectorXd& xs, VectorXcd *ys);
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
    vector<OpBasis*> op_basis_;
    VectorXd Ns_;
    MatrixXcd gAB, eAB, hAB;
    MatrixXd  RAB;
    multi_array<double,5> *d_; // d(A,B,nA,nB,N)
    virtual ~GaussBasis();
    virtual void setup();
    virtual void overlap(Operator ibra, Operator iket, MatrixXcd *res) = 0;
    virtual void at(Operator iop, const VectorXcd& cs, const VectorXd& xs, MatrixXcd *res);
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
