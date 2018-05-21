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

  // factory function for buffer
  OpBuf* MakeOpBuf(Pwgto *basis, Operator *op);
  OpBuf* MakeOpBuf(Pwgto1c *basis, Operator *op);

  // Buffer for each operator applying Pwgto or Pwgto1c.
  class OpBuf {
  public:
    virtual int Maxn(int A) = 0;
    virtual void SetUp(Pwgto *basis) = 0;
    virtual void Matrix(OpBuf *opbra, Pwgto *basis, MatrixXcd *res) = 0;
    virtual void Matrix(OpBufBasic *opket, Pwgto *basis, MatrixXcd *res) = 0;
    virtual void Matrix(OpBufGausspot *opket, Pwgto *basis, MatrixXcd *res) = 0;
    virtual void At(Pwgto *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res)=0;
    virtual void SetUp(Pwgto1c *basis) = 0;
    virtual void Matrix(OpBuf *opbra, Pwgto1c *basis, MatrixXcd *res) = 0;
    virtual void Matrix(OpBufBasic *opket, Pwgto1c *basis, MatrixXcd *res) = 0;
    virtual void Matrix(OpBufGausspot *opket, Pwgto1c *basis, MatrixXcd *res) = 0;
    virtual void At(Pwgto1c *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res)=0;
  };
  // Buffer for basic case. op.phi_A can be represented by polynomial times gauss function
  //     op.phi_A = sum_i c[A][i] (q-qA)^n[A][i] Exp(-gA(q-qA)^2 + ipA(q-qA)) 
  class OpBufBasic : public OpBuf {
  protected:
    int num_;
    vector<int> nums_;
    vector<VectorXi> ns_;   // ns_[A][i]
    vector<VectorXcd> cs_;  // cs_[A][i]
  public:    
    OpBufBasic(int num) : num_(num), nums_(num), ns_(num_), cs_(num_){}
    int num() const { return num_; }
    const vector<VectorXi>& ns() const { return ns_; }
    const vector<VectorXcd>& cs() const { return cs_; }
    int Maxn(int A);
    void InitZero(int A, int num);    
    virtual void SetUp(Pwgto *basis) = 0;
    void Matrix(OpBuf *opbra, Pwgto *basis, MatrixXcd *res);
    void Matrix(OpBufBasic *opket, Pwgto *basis, MatrixXcd *res);
    void Matrix(OpBufGausspot *opket, Pwgto *basis, MatrixXcd *res);
    void At(Pwgto *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
    virtual void SetUp(Pwgto1c *basis) = 0;
    void Matrix(OpBuf *opbra, Pwgto1c *basis, MatrixXcd *res);
    void Matrix(OpBufBasic *opket, Pwgto1c *basis, MatrixXcd *res);
    void Matrix(OpBufGausspot *opket, Pwgto1c *basis, MatrixXcd *res);
    void At(Pwgto1c *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
  };
  // Buffer for id operator
  class OpBufId : public OpBufBasic {
  public:
    OpBufId(const VectorXi& ns);
    //    OpBufId(const Polys& poly);
    void SetUp(Pwgto *basis);
    void SetUp(Pwgto1c *basis);
  };
  // Buffer for Rn operator
  class OpBufRn : public OpBufBasic {
    int n_;
  public:    
    OpBufRn(const VectorXi& ns, int n);
    //    OpBufRn(const Polys& poly, int n);
    int n() const { return n_; }
    void SetUp(Pwgto *basis);
    void SetUp(Pwgto1c *basis);
  };
  // Buffer for Pn operator
  class OpBufPn : public OpBufBasic {
    int n_;
  public:    
    OpBufPn(const VectorXi& ns, int n);
    //    OpBufPn(const Polys& poly, int n);
    int n() const { return n_; }
    void SetUp(Pwgto *basis);
    void SetUp(Pwgto1c *basis);
  };
  // Buffer for Da operator
  class OpBufDa : public OpBufBasic {
    int id_;
  public:    
    OpBufDa(const VectorXi& ns, int id);
    //    OpBufDa(const Polys& poly, int n);
    int id() const { return id_; }
    void SetUp(Pwgto *basis);
    void SetUp(Pwgto1c *basis);
  };
  // Buffer for Gaussian potential operator
  class OpBufGausspot : public OpBuf {
    int num_;
    VectorXi ns_;    
    OperatorGausspot *op_;
  public:    
    OpBufGausspot(const VectorXi& ns, OperatorGausspot *op);
    //    OpBufGausspot(const Polys& poly,  OperatorGausspot *op);
    int num() const { return num_; }
    const VectorXi& ns() const { return ns_; }
    OperatorGausspot *op() const { return op_; }
    int Maxn(int A);
    void SetUp(Pwgto *basis);
    void Matrix(OpBuf *opbrat, Pwgto *basis, MatrixXcd *res);
    void Matrix(OpBufBasic *opket, Pwgto *basis, MatrixXcd *res);
    void Matrix(OpBufGausspot *opket, Pwgto *basis, MatrixXcd *res);
    void At(Pwgto *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
    void SetUp(Pwgto1c *basis);
    void Matrix(OpBuf *opbra, Pwgto1c *basis, MatrixXcd *res);
    void Matrix(OpBufBasic *opket, Pwgto1c *basis, MatrixXcd *res);
    void Matrix(OpBufGausspot *opket, Pwgto1c *basis, MatrixXcd *res);
    void At(Pwgto1c *basis, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
  };
  // Buffer for numerical potential
  class OpBufSpline : public OpBufBasic {    
    OperatorSpline *op_;
  public:
    OpBufSpline(const VectorXi& ns, int norder, OperatorSpline *op);
    OperatorSpline *op() const { return op_; }
    void SetUp(Pwgto *basis);
    void SetUp(Pwgto1c *basis);
  };
}

#endif


