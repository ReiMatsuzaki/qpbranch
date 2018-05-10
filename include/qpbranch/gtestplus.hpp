#ifndef GTEST_PLUS_HPP_
#define GTEST_PLUS_HPP_

#include <complex>
#include <gtest/gtest.h>
#include <Eigen/Core>

::testing::AssertionResult AssertComplexNear(const char *a_expr,
					     const char *b_expr,
					     std::complex<double> a,
					     std::complex<double> b,
					     double tol);
::testing::AssertionResult AssertComplexEq(const char *a_expr,
					   const char *b_expr,
					   std::complex<double> a,
					   std::complex<double> b);
::testing::AssertionResult AssertMatrixXcdNear(const char *a_expr,
					       const char *b_expr,
					       const Eigen::MatrixXcd& a,
					       const Eigen::MatrixXcd& b,
					       double tol);
::testing::AssertionResult AssertMatrixXcdEq(const char *a_expr,
					     const char *b_expr,
					     const Eigen::MatrixXcd& a,
					     const Eigen::MatrixXcd& b);
#define EXPECT_C_EQ(a, b) EXPECT_PRED_FORMAT2(  AssertComplexEq, a, b);
#define EXPECT_C_NEAR(a, b) EXPECT_PRED_FORMAT2(AssertComplexEq, a, b);
#define EXPECT_MATRIXXCD_EQ(a, b) EXPECT_PRED_FORMAT2(  AssertMatrixXcdEq, a, b);
#define EXPECT_MATRIXXCD_NEAR(a, b) EXPECT_PRED_FORMAT2(AssertMatrixXcdEq, a, b);

#endif
