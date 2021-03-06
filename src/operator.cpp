#include <qpbranch/operator.hpp>

using namespace std;

namespace qpbranch {
  string OperatorDa::str() const {
    string buf = "Da[";
    buf += "id = " + to_string(id_) + "]";
    return buf;
  }
  void OperatorGausspot::At(const VectorXd& xs, VectorXcd *res) const {
    assert(xs.size()==res->size());
    *res = v0_ * (-b_*xs.array().pow(2)).exp();
  }
  string OperatorGausspot::str() const {
    //string buf = "Gausspot[" + to_string(v0_.real()) + "]";
    string buf = "Gausspot[";
    buf += "v0 = (" + to_string(v0_.real()) + ", " +  to_string(v0_.imag()) + "), ";
    buf += "b  = (" + to_string(b_.real()) + ", " +  to_string(b_.imag()) + "), ";
    buf += "q0 = (" + to_string(q0_.real()) + ", " +  to_string(q0_.imag()) + ") ] ";
    return buf;
  }
  OperatorPoly::OperatorPoly(int maxn, const VectorXd& cs) : maxn_(maxn), cs_(cs) {
    assert(cs_.size()==maxn+1);
  }
  void OperatorPoly::At(const VectorXd& xs, VectorXcd *res) const {
    assert(xs.size()==res->size());
    for(auto i = 0; i < xs.size(); i++) {
      auto cumsum = 0.0;
      for(int n = 0; n <= maxn_; n++) {
	cumsum += cs_[n] * pow(xs[i], n);
      }
      (*res)(i) = cumsum;
    }
  }
  OperatorSpline::OperatorSpline(const VectorXd& xs, const VectorXd& ys) {
    assert(xs.size()==ys.size());

    // convert Eigen::VectorXcd to STL vector.
    int n = ys.size();
    std::vector<double> yys(n);
    for(int i = 0; i < n; i++)
      yys[i] = ys[i];

    // build spline object
    spline_ = new boost::math::cubic_b_spline<double>
      (yys.begin(), yys.end(), xs[0], xs[1]-xs[0]);
  }
  void OperatorSpline::At(const VectorXd& xs, VectorXcd* res) const {
    assert(xs.size()==res->size());
    for(auto i = 0; i < xs.size(); i++) {
      (*res)(i) = (*spline_)(xs[i]);
    }
  }
  OperatorSplineP1::OperatorSplineP1(const VectorXd& xs, const VectorXd& ys) {
    assert(xs.size()==ys.size());

    // convert Eigen::VectorXcd to STL vector.
    int n = ys.size();
    std::vector<double> yys(n);
    for(int i = 0; i < n; i++)
      yys[i] = ys[i];
    
    // build spline object
    spline_ = new boost::math::cubic_b_spline<double>
      (yys.begin(), yys.end(), xs[0], xs[1]-xs[0]);

  }
  void OperatorSplineP1::At(const VectorXd& xs, VectorXcd *res) const {
    assert(xs.size()==res->size());
    for(int i = 0; i < xs.size(); i++) {
      (*res)(i) = (*spline_)(xs[i]);
    }
  }
}



