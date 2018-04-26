#ifndef OPERATOR_HPP_
#define OPERATOR_HPP_

#include <complex>

namespace qpbranch {

  using std::complex;

  class PlaneWaveGto;
  
  class Operator {
  public:
    //    Operator();
    //    virtual ~Operator();
    virtual void call_new(PlaneWaveGto *ptr) =0;
  };
  class OperatorId : public Operator {
  public:
    OperatorId() {}
    ~OperatorId() {}
    void call_new(PlaneWaveGto *ptr);
  };
  class OperatorRn : public Operator {
    int n_;
  public:
    OperatorRn(int n);
    inline int n() const { return n_; }
    void call_new(PlaneWaveGto *ptr);
  };
  class OperatorPn : public Operator {
    int n_;
  public:
    OperatorPn(int n);
    inline int n() const { return n_; }
    void call_new(PlaneWaveGto *ptr);
  };
  class OperatorDa : public Operator {
  public:
    int id_; // id for variable.
    OperatorDa(int id);
    void call_new(PlaneWaveGto *ptr);
  };
  class OperatorGausspot : public Operator {
  public:
    complex<double> v0_, b_, q0_;
    OperatorGausspot(complex<double> v0, complex<double> b, complex<double> q0);
    void call_new(PlaneWaveGto *basis);
  };

}
#endif
