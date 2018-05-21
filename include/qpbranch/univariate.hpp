#include <vector>
#include <algorithm>
#include <string>
#include <ostream>
#include <cmath>
#include <cassert>

namespace qpbranch {

  using std::string;
  using std::vector;
  using std::pair;
  using std::make_pair;
  using std::to_string;
  using std::ostream;

  // This class express Univariate polynomial
  //      c0 x^n0 + c1 x^n1 + ...
  // This class supports some arithmetic.
  template<class Field>
  class Univariate {
  private:
    typedef vector<pair<Field, int> > Data;
    Data ncs_;
  public:
    typedef typename Data::iterator iterator;
    typedef typename Data::const_iterator const_iterator;
    Univariate() {}
    Univariate(int num) : ncs_(num) {}
    Univariate(const Data& data) : ncs_(data) {}
    Univariate(const Univariate<Field>& o) : ncs_(o.ncs_) {}
    ~Univariate() {}
    // Gives numerical value with substitution.
    const Field Sub(const Field& x) {
      Field cumsum(0);
      for(auto nc : ncs_) {
	cumsum += nc.first * pow(x, nc.second);
      }
      return cumsum;
    }
    // Give short string expression
    string str() const {
      string buf;
      for(const_iterator it = begin(); it != end(); ++it) {	
	Field c(it->first);
	int n = it->second;
	
	if(it != begin()) {
	  buf += " + ";
	}
	buf += to_string(c) + "X^" + to_string(n);
	
      }
      return buf;
    }
    // Below methods are vector methods.
    int size() const { return ncs_.size(); }
    void resize(int n) { ncs_.resize(n); }
    const_iterator begin() const { return ncs_.begin(); }
    const_iterator end() const { return ncs_.end(); }
    iterator begin() { return ncs_.begin(); }
    iterator end() { return ncs_.end(); }
    // Arithmetric
    Univariate<Field>& operator+=(const Univariate<Field>& o) {
      for(auto it = o.begin(); it != o.end(); ++it) {
	ncs_.push_back(make_pair(it->first, it->second));
      }
      return *this;
    }
    Univariate<Field>& operator*=(const Field& c) {
      for(auto it = begin(); it != end(); ++it) {
	it->first *= c;
      }
      return *this;
    }
    Univariate<Field>& operator*=(const Univariate<Field>& o) {
      Data orig;
      this->ncs_.swap(orig);
      for(auto it = orig.begin(); it != orig.end(); ++it) {
	for(auto jt = o.begin(); jt != o.end(); ++jt) {
	  ncs_.push_back(make_pair(it->first*jt->first, it->second+jt->second));
	}
      }
      return *this;
    }
    Univariate<Field>  operator-() {
      Univariate<Field> a(*this);
      a.Negative();
      return a;
    }
    Univariate<Field>  operator-(const  Univariate<Field>& o) {
      Univariate<Field> a(o);
      a.Negative();
      a += *this;
      return a;
    }
    void Negative() {
      for(iterator it = ncs_.begin(); it != ncs_.end(); ++it) {
	it->first *= -Field(1);
      }
    }
    // calc
    Univariate<Field> Diff(int n) const {
      assert(n>0);

      Univariate<Field> res(*this);

      if(n == 1) {
	for(int i = 0; i < size(); i++) {
	  Field ci = ncs_[i].first;
	  int ni = ncs_[i].second;
	  if(ncs_[i].second>n) {
	    res.ncs_[i].first  = ci * Field(ni);
	    res.ncs_[i].second = ni - 1;
	  } else {
	    res.ncs_[i].first  = 0.0;
	    res.ncs_[i].second = 0;
	  }
	}
      } else {
	assert(false||"not impl");
      }
      return res;
    }
    int Order() const {
      int res(0);
      for(auto nc : ncs_) {
	res = std::max(res, nc.second);
      }
      return res;
    }
    // Remove unnecessary data.
    void Refresh() {

      // reduce coefficients for same power
      for(auto it = begin(); it != end(); ++it) {
	for(auto jt = it+1; jt != end(); ++jt) {
	  if(it->second == jt->second) {
	    it->first += jt->first;
	    jt->first = Field(0);
	  }
	}
      }

      // remove zero coefficient
      typedef pair<Field,int> FI;
      double tol = std::pow(10.0, -14.0);
      auto is_zero = [&](FI cn) { return std::abs(cn.first)<tol; } ;
      ncs_.erase(std::remove_if(ncs_.begin(), ncs_.end(), is_zero), ncs_.end());
      
    }
  };

  template<class Field>
  ostream& operator<<(ostream& os, const Univariate<Field>& uni) {
    os << uni.str();
    return os;
  }
  template<class Field>
  Univariate<Field> operator+(const Univariate<Field>& a, const Univariate<Field>& b) {
    Univariate<Field> ab(a);
    ab += b;
    return ab;
  }
  template<class Field>
  Univariate<Field> operator*(const Field& c, const Univariate<Field>& a) {
    Univariate<Field> ca(a);
    ca *= c;
    return ca;
  }
  template<class Field>
  Univariate<Field> operator*(const Univariate<Field>& a, const Univariate<Field>& b) {
    Univariate<Field> ab(a);
    ab *= b;
    return ab;
  }
}
