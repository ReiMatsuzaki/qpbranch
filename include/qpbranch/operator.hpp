#ifndef OPERATOR_HPP_
#define OPERATOR_HPP_

#include <complex>

namespace qpbranch {

  using std::complex;
  
  class Operator {
  public:
    virtual std::string str() const = 0;
  };
  class OperatorId : public Operator {
  public:
    std::string str() const { return "id"; }
  };
  class OperatorRn : public Operator {
    int n_;
  public:
    OperatorRn(int n);    
    inline int n() const { return n_; }
    std::string str() const { return "Rn"; }
  };
  class OperatorPn : public Operator {
    int n_;
  public:
    OperatorPn(int n);
    inline int n() const { return n_; }
    std::string str() const { return "Pn"; }
  };
  class OperatorDa : public Operator {
  public:
    int id_; // id for variable.
    OperatorDa(int id);
    inline int id() const { return id_; }
    std::string str() const { return "Da"; }
  };
  class OperatorGausspot : public Operator {
  public:
    complex<double> v0_, b_, q0_;
    OperatorGausspot(complex<double> v0, complex<double> b, complex<double> q0);
    std::string str() const { return "Gausspot"; }
  };

}
#endif
