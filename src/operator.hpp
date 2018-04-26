#ifndef OPERATOR_HPP_
#define OPERATOR_HPP_

namespace qpbranch {

  class PlaneWaveGto;
  
  class Operator {
  public:
    //    Operator();
    //    virtual ~Operator();
    virtual void call_new(PlaneWaveGto *ptr) =0;
    virtual void call_setup(PlaneWaveGto *ptr) =0;
  };
  class OperatorId : public Operator {
  public:
    OperatorId() {}
    ~OperatorId() {}
    void call_new(PlaneWaveGto *ptr);
    void call_setup(PlaneWaveGto *ptr);
  };
  class OperatorRn : public Operator {
    int n_;
  public:
    OperatorRn(int n);
    inline int n() const { return n_; }
    void call_new(PlaneWaveGto *ptr);
    void call_setup(PlaneWaveGto *ptr);
  };
  class OperatorPn : public Operator {
    int n_;
  public:
    OperatorPn(int n);
    inline int n() const { return n_; }
    void call_new(PlaneWaveGto *ptr);
    void call_setup(PlaneWaveGto *ptr);
  };
  class OperatorDR : public Operator {
  public:
    void call_new(PlaneWaveGto *ptr);
    void call_setup(PlaneWaveGto *ptr);
  };
  class OperatorDP : public Operator {
  public:
    void call_new(PlaneWaveGto *ptr);
    void call_setup(PlaneWaveGto *ptr);
  };
  class OperatorDgr : public Operator {
  public:
    void call_new(PlaneWaveGto *ptr);
    void call_setup(PlaneWaveGto *ptr);
  };
  class OperatorDgi : public Operator {
  public:
    void call_new(PlaneWaveGto *ptr);
    void call_setup(PlaneWaveGto *ptr);
  };

  /*
  class BasisSet {
  public:
    virtual void NewOpRn(OperatorRn *op);
    virtual void NewOpPn(OperatorPn *op);
  };
  class PlaneWaveGto : public BasisSet {
  public:
    PlaneWaveGto(const VectorXi& ns, vector<Operator*> ops);
    void new_op(OperatorRn *op);
    void new_op(OperatorPn *op);
    void setup_op(OperatorRn *op);
    void setup_op(OperatorPn *op);
  };
  */
}
#endif
