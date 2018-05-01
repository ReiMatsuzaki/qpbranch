#ifndef PWGTO_BUF_HPP_
#define PWGTO_BUF_HPP_

#include <vector>
#include <Eigen/Core>
#include "operator.hpp"

namespace qpbranch {
  using std::vector;
  using Eigen::MatrixXcd;
  using Eigen::VectorXcd;
  using Eigen::VectorXd;
  using Eigen::VectorXi;
    
  class OpBufBasic;

  class Pwgto;
  class Pwgto1c;
  
  class OpBufGausspot;

  // Buffer for each operator applying PlaneWaveGto.
  class OpBuf {
  public:
    virtual int maxn(int A) = 0;
    virtual void setup(Pwgto *basis) = 0;
    virtual void matrix(OpBuf *opbra, Pwgto *basis, MatrixXcd *res) = 0;
    virtual void matrix(OpBufBasic *opket, Pwgto *basis, MatrixXcd *res) = 0;
    virtual void matrix(OpBufGausspot *opket, Pwgto *basis, MatrixXcd *res) = 0;
    virtual void at(Pwgto *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res)=0;
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
    OpBufBasic(int num);
    int maxn(int A);
    void init_zero(int A, int num);    
    virtual void setup(Pwgto *basis) = 0;
    void matrix(OpBuf *opbra, Pwgto *basis, MatrixXcd *res);
    void matrix(OpBufBasic *opket, Pwgto *basis, MatrixXcd *res);
    void matrix(OpBufGausspot *opket, Pwgto *basis, MatrixXcd *res);
    void at(Pwgto *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
    virtual void setup(Pwgto1c *basis) = 0;
    void matrix(OpBuf *opbra, Pwgto1c *basis, MatrixXcd *res);
    void matrix(OpBufBasic *opket, Pwgto1c *basis, MatrixXcd *res);
    void matrix(OpBufGausspot *opket, Pwgto1c *basis, MatrixXcd *res);
    void at(Pwgto1c *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
  };
  class OpBufId : public OpBufBasic {
  public:
    OpBufId(const VectorXi& ns);
    void setup(Pwgto *basis);
    void setup(Pwgto1c *basis);
  };
  class OpBufRn : public OpBufBasic {
  public:
    int n_;
    OpBufRn(const VectorXi& ns, int n);
    void setup(Pwgto *basis);
    void setup(Pwgto1c *basis);
  };
  class OpBufPn : public OpBufBasic {
  public:
    int n_;
    OpBufPn(const VectorXi& ns, int n);
    void setup(Pwgto *basis);
    void setup(Pwgto1c *basis);
  };
  class OpBufDa : public OpBufBasic {
  public:
    int id_;
    OpBufDa(const VectorXi& ns, int id);
    void setup(Pwgto *basis);
    void setup(Pwgto1c *basis);
  };
  class OpBufGausspot : public OpBuf {
  public:
    int num_;
    VectorXi ns_;    
    OperatorGausspot *op_;
    OpBufGausspot(const VectorXi& ns, OperatorGausspot *op);
    
    int maxn(int A);
    void setup(Pwgto *basis);
    void matrix(OpBuf *opbrat, Pwgto *basis, MatrixXcd *res);
    void matrix(OpBufBasic *opket, Pwgto *basis, MatrixXcd *res);
    void matrix(OpBufGausspot *opket, Pwgto *basis, MatrixXcd *res);
    void at(Pwgto *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
    void setup(Pwgto1c *basis);
    void matrix(OpBuf *opbra, Pwgto1c *basis, MatrixXcd *res);
    void matrix(OpBufBasic *opket, Pwgto1c *basis, MatrixXcd *res);
    void matrix(OpBufGausspot *opket, Pwgto1c *basis, MatrixXcd *res);
    void at(Pwgto1c *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
  };

  OpBuf* make_op_buff(const VectorXi& ns, Operator *op);
  
  // factory function for OpBuf
  //  OpBuf* make_op_buff(PlaneWaveGto *basis, Operator *op);
  
}

#endif


