#include <Eigen/Eigenvalues>

namespace qpbranch {
  using Eigen::MatrixXcd;
  using Eigen::VectorXd;
  using Eigen::GeneralizedSelfAdjointEigenSolver;
  void zhegv(const MatrixXcd& H, const MatrixXcd& S, MatrixXcd *U, VectorXd *w) {
    GeneralizedSelfAdjointEigenSolver<MatrixXcd> es;
    es.compute(H, S);
    *w = es.eigenvalues();
    *U = es.eigenvectors();
  }
}
