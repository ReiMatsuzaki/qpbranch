#include <qpbranch/gtestplus.hpp>

using namespace std;
using namespace Eigen;

static double eps(pow(10.0, -12.0));

bool ComplexIsNear(std::complex<double> a, std::complex<double> b, double eps) {
  return abs(a-b)<eps;
}

/*
bool ComplexIsEq(std::complex<double> a, std::complex<double> b) {
static eps(pow(10.0, -12.0));
  return ComplexIsNear(a, b, eps);
}
*/

::testing::AssertionResult AssertComplexNear(const char *a_expr,
					     const char *b_expr,
					     const char *eps_expr,
					     std::complex<double> a,
					     std::complex<double> b,
					     double eps) {
  if(ComplexIsNear(a,b,eps))
    return ::testing::AssertionSuccess();

  return ::testing::AssertionFailure()
    << a_expr << " and "  << b_expr << " are not near." << endl
    << a_expr << " : " << a << endl
    << b_expr << " : " << b << endl    
    << "|" << a_expr << "-" << b_expr << "| : " << abs(a-b) << endl
    << eps_expr << " : " << eps;
}
::testing::AssertionResult AssertComplexEq(const char *a_expr,
					   const char *b_expr,
					   std::complex<double> a,
					   std::complex<double> b) {
  return AssertComplexNear(a_expr, b_expr, "eps", a, b, eps);
}

/*
::testing::AssertionResult AssertMatrixXcdNear(const char *a_expr,
					       const char *b_expr,
					       const char *eps_expr,
					       const MatrixXcd& a,
					       const MatrixXcd& b,
					       double eps) {
  if(a.cols() != b.cols() || a.rows() != b.rows()) {
    return ::testing::AssertionFailure()
      << a_expr << " and "  << b_expr << " are different size." << endl
      << a_expr << " : " << a.cols() << ", " << a.rows() << endl
      << b_expr << " : " << b.cols() << ", " << b.rows() << endl;
  }

  for(int i = 0; i < a.rows(); i++) {
    for(int j = 0; j < a.cols(); j++) {
      if(!ComplexIsNear(a(i,j),b(i,j),eps)) {
      string str_a(a_expr);
	str_a += "(" + to_string(i) + "," + to_string(j)+ ")";
	string str_b(a_expr);
	str_b += "(" + to_string(i) + "," + to_string(j)+ ")";
	auto resij = AssertComplexNear(str_a.c_str(), str_b.c_str(), eps_expr,
				       a(i,j), b(i,j), eps);
	if(!resij) {
	  return resij;
	}
      }
    }
  }
  return ::testing::AssertionSuccess();
}
::testing::AssertionResult AssertMatrixXcdEq(const char *a_expr,
					     const char *b_expr,
					     const MatrixXcd& a,
					     const MatrixXcd& b) {
  return AssertMatrixXcdNear(a_expr, b_expr, "eps", a, b, eps);
}
::testing::AssertionResult AssertVectorXdNear(const char *a_expr,
					      const char *b_expr,
					      const char *eps_expr,
					      const VectorXd& a,
					      const VectorXd& b,
					      double eps) {
  if(a.size() != b.size() ) {
    return ::testing::AssertionFailure()
      << a_expr << " and "  << b_expr << " are different size." << endl
      << a_expr << " : " << a.size() << endl
      << b_expr << " : " << b.size() << endl;
  }

  for(int i = 0; i < a.rows(); i++) {
    if(!ComplexIsNear(a(i),b(i),eps)) {
      string str_a(a_expr);
      str_a += "(" + to_string(i) + ")";
      string str_b(b_expr);
      str_b += "(" + to_string(i) + ")";
      auto resij = AssertComplexNear(str_a.c_str(), str_b.c_str(), eps_expr,
				     a(i), b(i), eps);
      if(!resij) {
	return resij;
      }
    }
  }
  return ::testing::AssertionSuccess();
}
::testing::AssertionResult AssertVectorXdEq(const char *a_expr,
					     const char *b_expr,
					     const VectorXd& a,
					     const VectorXd& b) {
  return AssertVectorXdNear(a_expr, b_expr, "eps", a, b, eps);
}

*/


