#ifndef GTEST_PLUS_HPP_
#define GTEST_PLUS_HPP_

#include <complex>
#include <gtest/gtest.h>
#include <Eigen/Core>

::testing::AssertionResult AssertComplexNear(const char *a_expr,
					     const char *b_expr,
					     const char *eps_expr,
					     std::complex<double> a,
					     std::complex<double> b,
					     double tol);
::testing::AssertionResult AssertComplexEq(const char *a_expr,
					   const char *b_expr, 
					   std::complex<double> a,
					   std::complex<double> b);
#define EXPECT_C_EQ(a, b) EXPECT_PRED_FORMAT2(  AssertComplexEq, a, b);
#define EXPECT_C_NEAR(a, b, eps) EXPECT_PRED_FORMAT3(AssertComplexNear, a, b, eps);

template<typename Scalar, int Rows, int Cols>
::testing::AssertionResult
AssertMatrixNear(const char *a_expr,
		 const char *b_expr,
		 const char *eps_expr,
		 const Eigen::Matrix<Scalar, Rows, Cols>& a,
		 const Eigen::Matrix<Scalar, Rows, Cols>& b,
		 double eps) {

  using namespace std;
  
  if(a.cols() != b.cols() || a.rows() != b.rows()) {
    return ::testing::AssertionFailure()
      << a_expr << " and "  << b_expr << " are different size." << endl
      << a_expr << " : " << a.cols() << ", " << a.rows() << endl
      << b_expr << " : " << b.cols() << ", " << b.rows() << endl;
  }
  
  for(int i = 0; i < a.rows(); i++) {
    for(int j = 0; j < a.cols(); j++) {
      double diff = std::abs(a(i,j)-b(i,j));
      if(diff > eps) {
	string str_a(a_expr);
	string str_b(b_expr);
	if(Cols==1) {
	  str_a += "(" + to_string(i) + ")";	
	  str_b += "(" + to_string(i) + ")";
	} else if(Rows==1) {
	  str_a += "(" + to_string(j)+ ")";	
	  str_b += "(" + to_string(j)+ ")";
	} else {
	  str_a += "(" + to_string(i) + "," + to_string(j)+ ")";	
	  str_b += "(" + to_string(i) + "," + to_string(j)+ ")";
	}
	
	return ::testing::AssertionFailure()
	  << str_a << " and " << str_b << " are not near." << endl
	  << str_a << " : "   << a(i,j) << endl
	  << str_b << " : "   << b(i,j) << endl    
	  << "|" << a_expr << "-" << b_expr << "| : " << diff << endl
	  << eps_expr << " : " << eps;
      }
    }
  }
  
  return ::testing::AssertionSuccess();
}
template<typename Scalar, int Rows, int Cols>
::testing::AssertionResult
AssertMatrixEqual(const char *a_expr,
		  const char *b_expr,
		  const Eigen::Matrix<Scalar, Rows, Cols>& a,
		  const Eigen::Matrix<Scalar, Rows, Cols>& b) {
  static double eps(pow(10.0,-13.0));
  return AssertMatrixNear(a_expr, b_expr, "eps", a, b, eps);
}



// MatrixXcd
::testing::AssertionResult AssertMatrixXcdEqual(const char *a_expr, const char *b_expr,						
						const Eigen::MatrixXcd& a, const Eigen::MatrixXcd& b) {
  return AssertMatrixEqual<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic>(a_expr, b_expr, a, b);
}

::testing::AssertionResult AssertMatrixXcdNear(const char *a_expr, const char *b_expr, const char *eps_expr,
					       const Eigen::MatrixXcd& a, const Eigen::MatrixXcd& b, double eps) {
  return AssertMatrixNear<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic>(a_expr, b_expr, eps_expr, a, b, eps);
}
#define EXPECT_MATRIXXCD_EQ(a, b) EXPECT_PRED_FORMAT2(       AssertMatrixXcdEqual, a, b);
#define EXPECT_MATRIXXCD_NEAR(a, b, eps) EXPECT_PRED_FORMAT3(AssertMatrixXcdNear, a, b, eps);


// VectorXd
::testing::AssertionResult AssertVectorXdEqual(const char *a_expr, const char *b_expr,						
						const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
  return AssertMatrixEqual<double,Eigen::Dynamic,1>(a_expr, b_expr, a, b);
}

::testing::AssertionResult AssertVectorXdNear(const char *a_expr, const char *b_expr, const char *eps_expr,
					      const Eigen::VectorXd& a, const Eigen::VectorXd& b, double eps) {
  return AssertMatrixNear<double,Eigen::Dynamic,1>(a_expr, b_expr, eps_expr, a, b, eps);
}
#define EXPECT_VECTORXD_EQ(a, b) EXPECT_PRED_FORMAT2(       AssertVectorXdEqual, a, b);
#define EXPECT_VECTORXD_NEAR(a, b, eps) EXPECT_PRED_FORMAT3(AssertVectorXdNear, a, b, eps);

// VectorXcd
::testing::AssertionResult AssertVectorXcdEqual(const char *a_expr, const char *b_expr,						
						const Eigen::VectorXcd& a, const Eigen::VectorXcd& b) {
  return AssertMatrixEqual<std::complex<double>,Eigen::Dynamic,1>(a_expr, b_expr, a, b);
}
::testing::AssertionResult AssertVectorXcdNear(const char *a_expr, const char *b_expr, const char *eps_expr,
					       const Eigen::VectorXcd& a, const Eigen::VectorXcd& b, double eps) {
  return AssertMatrixNear<std::complex<double>,Eigen::Dynamic,1>(a_expr, b_expr, eps_expr, a, b, eps);
}
#define EXPECT_VECTORXCD_EQ(a, b) EXPECT_PRED_FORMAT2(       AssertVectorXcdEqual, a, b);
#define EXPECT_VECTORXCD_NEAR(a, b, eps) EXPECT_PRED_FORMAT3(AssertVectorXcdNear, a, b, eps);


#endif
