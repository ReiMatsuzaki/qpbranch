#ifndef OPERATOR_HPP_
#define OPERATOR_HPP_

#include <complex>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <Eigen/Core>

namespace qpbranch {

  using std::string;
  using std::complex;
  using Eigen::VectorXd;
  using Eigen::VectorXcd;
  
  class Operator {
  public:
    virtual string str() const = 0;
  };
  class OperatorId : public Operator {
  public:
    string str() const { return "id"; }
  };
  class OperatorRn : public Operator {
    int n_;
  public:
    OperatorRn(int n) : n_(n) {}
    inline int n() const { return n_; }
    string str() const { return "Rn"; }
  };
  class OperatorPn : public Operator {
    int n_;
  public:
    OperatorPn(int n) : n_(n) {}
    inline int n() const { return n_; }
    string str() const { return "Pn"; }
  };
  class OperatorDa : public Operator {
    int id_; // id for variable.
  public:
    OperatorDa(int id) : id_(id) {}
    inline int id() const { return id_; }
    string str() const;
  };
  class OperatorPot : public Operator {
  public:
    virtual void At(const VectorXd& xs, VectorXcd *res) const = 0;
    virtual string str() const = 0;
  };
  class OperatorGausspot : public OperatorPot {
    complex<double> v0_, b_, q0_;
  public:
    complex<double> v0() const { return v0_; }
    complex<double> b() const { return b_; }
    complex<double> q0() const { return q0_; }
    OperatorGausspot(complex<double> v0, complex<double> b, complex<double> q0):
      v0_(v0), b_(b), q0_(q0) {}
    void At(const VectorXd& xs, VectorXcd *res) const;
    string str() const;
  };
  class OperatorSpline : public OperatorPot {
    boost::math::cubic_b_spline<double> *spline_;
  public:
    OperatorSpline(const VectorXd& xs, const VectorXd& ys);
    void At(const VectorXd& xs, VectorXcd *res) const;
    string str() const;
  };
}
#endif
