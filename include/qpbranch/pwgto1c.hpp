#ifndef PWGTO1C_HPP_
#define PWGTO1C_HPP_

namespace qpbranch {

  class Pwgto1c;
  
namespace pwgto_1c {

  typedef Pwgto1c Basis;

  // Buffer for each operator applying PlaneWaveGto.
  class OpBuf {
  public:
    virtual int maxn(int A) = 0;
    virtual void setup(Basis *basis) = 0;
    virtual void matrix(OpBuf *opbra, Basis *basis, MatrixXcd *res) = 0;
    virtual void matrix(OpBufBasic *opket, Basis *basis, MatrixXcd *res) = 0;
    virtual void matrix(OpBufGausspot *opket, Basis *basis, MatrixXcd *res) = 0;
    virtual void at(Basis *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res)=0;
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
    virtual void setup(Basis *basis) = 0;    
    void matrix(OpBuf *opbra, Basis *basis, MatrixXcd *res);
    void matrix(OpBufBasic *opket, Basis *basis, MatrixXcd *res);
    void matrix(OpBufGausspot *opket, Basis *basis, MatrixXcd *res);
    void at(Basis *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
    void init_zero(int A, int num);
  };
  class OpBufId : public OpBufBasic {
  public:
    OpBufId(const VectorXi& ns);
    void setup(Basis *basis);
  };
  class OpBufRn : public OpBufBasic {
  public:
    int n_;
    OpBufRn(const VectorXi& ns, int n);
    void setup(Basis *basis);
  };
  class OpBufPn : public OpBufBasic {
  public:
    int n_;
    OpBufPn(const VectorXi& ns, int n);
    void setup(Basis *basis);
  };
  class OpBufDa : public OpBufBasic {
  public:
    int id_;
    OpBufDa(const VectorXi& ns, int id);
    void setup(Basis *basis);
  };
  class OpBufGausspot : public OpBuf {
  public:
    int num_;
    VectorXi ns_;    
    OperatorGausspot *op_;
    OpBufGausspot(const VectorXi& ns, OperatorGausspot *op);
    
    int maxn(int A);
    void setup(Basis *basis);
    void matrix(OpBuf *opbrat, Basis *basis, MatrixXcd *res);
    void matrix(OpBufBasic *opket, Basis *basis, MatrixXcd *res);
    void matrix(OpBufGausspot *opket, Basis *basis, MatrixXcd *res);
    void at(Basis *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
  };

  OpBuf* make_op_buff(Basis *basis, Operator *op);
  
} /* end namespace pwgto_1c */
  
  class Pwgto1c {
  public:
    // data size
    int num_, nop_;
    // cont
    VectorXi ns_;
    vector<Operator*> ops_;
    // variable
    VectorXi ns_;
    double R0_, P0_;
    complex<double> g0_;
    // intermediate
    map<Operator*, OpBuf*> buf_map_;
    VectorXcd gints_; // gints[n] = Int[q^n Exp[-2Re[q0]q^2]]
    // Main
    Pwgto1c(const VectorXi& ns, double R0, double P0, complex<double> g0, const vector<Operator*>& ops);
    ~Pwgto1c();
    void setup();
    void matrix(Operator *opbra, Operator *opket, MatrixXcd *res);
    void at(Operator *op, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
  };


  
} /* end namespace qpbranch */

#endif
