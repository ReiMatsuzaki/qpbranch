#include <iostream>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

namespace qpbranch {
  using Eigen::MatrixXcd;
  using Eigen::MatrixXd;
  using Eigen::VectorXcd;
  using Eigen::VectorXd;  
  using Eigen::Success;
  using std::complex;
  using std::cerr;
  using std::endl;
  
  void Zheev(const MatrixXcd& H, MatrixXcd *U, VectorXd *w) {
    int n = H.cols();
    assert(n==H.rows());
    assert(n==U->cols());
    assert(n==U->rows());
    assert(n==w->size());

    Eigen::SelfAdjointEigenSolver<MatrixXcd> es(H);
    if(es.info() != Success) {
      cerr << __FILE__ << ":" << __LINE__ << ":Error on Zheev" << endl; abort();
    }

    *U = es.eigenvectors();
    *w = es.eigenvalues();
    
  }
  void Zhegv(const MatrixXcd& H, const MatrixXcd& S, MatrixXcd *U, VectorXd *w) {
    int n = H.cols();
    assert(n==H.rows());
    assert(n==S.rows());
    assert(n==S.cols());
    
    Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXcd> es;
    es.compute(H, S);
    *w = es.eigenvalues();
    *U = es.eigenvectors();
    
  }
  void Zgesv(const MatrixXcd& M, const VectorXcd& b, VectorXcd *x) {
    Eigen::FullPivLU<MatrixXcd> lu(M);
    *x = lu.solve(b);
  }
  void Dgesv(const MatrixXd&  M, const VectorXd& b, VectorXd *x) {
    Eigen::FullPivLU<MatrixXd> lu(M);
    *x = lu.solve(b);    
  }
  void IntetDiag(const MatrixXcd& H, double dt, VectorXcd *ptr_res) {
    complex<double> ii(0.0, 1.0);

    int n = H.rows();
    assert(H.cols()==n);
    assert(ptr_res->size()==n);

    VectorXd w(n);
    MatrixXcd U(n,n);

    VectorXcd& c(*ptr_res);

    Zheev(H, &U, &w);

    c = U.adjoint() * c;
    c = (-ii*w*dt).array().exp() * c.array();
    c = U*c;
  }
  void IntetGdiag(const MatrixXcd& H, const MatrixXcd& S, double dt, VectorXcd *ptr_res) {

    complex<double> ii(0.0, 1.0);    
    
    int n = H.rows();
    assert(H.cols()==n);
    assert(S.rows()==n);
    assert(S.cols()==n);
    assert(ptr_res->size()==n);

    VectorXd w(n);
    MatrixXcd U(n,n);

    VectorXcd& c(*ptr_res);
    
    Zhegv(H, S, &U, &w);

    c = U.adjoint() * S * c;
    c = (-ii*w*dt).array().exp() * c.array();
    c = U*c;
    
  }
}
