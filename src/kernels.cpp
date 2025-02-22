#include <RcppEigen.h>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace Eigen;

//--------------------------------------------------------------
// Internal worker function: computes the learner estimate - thread safe for OMP
List learner_worker(const Eigen::MatrixXd &Y_source,
                               const Eigen::MatrixXd &Y_target,
                               int r, double lambda1, double lambda2,
                               double step_size, int max_iter, double threshold,
                               double max_value) {
  int p = Y_source.rows();
  int q = Y_source.cols();

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Y_source, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd U_full = svd.matrixU();
  Eigen::MatrixXd V_full = svd.matrixV();
  int r_use = std::min(r, static_cast<int>(U_full.cols()));
  Eigen::MatrixXd U_trunc = U_full.leftCols(r_use);
  Eigen::MatrixXd V_trunc = V_full.leftCols(r_use);
  Eigen::VectorXd singular_vals = svd.singularValues().head(r_use);

  // Initialize as in R: U_init = U_trunc * sqrt(D), V_init = V_trunc * sqrt(D)
  Eigen::MatrixXd sqrtD = singular_vals.array().sqrt().matrix().asDiagonal();
  Eigen::MatrixXd U = U_trunc * sqrtD;
  Eigen::MatrixXd V = V_trunc * sqrtD;

  // Projection matrices (using truncated basis)
  Eigen::MatrixXd P_U = Eigen::MatrixXd::Identity(p, p) - U_trunc * U_trunc.transpose();
  Eigen::MatrixXd P_V = Eigen::MatrixXd::Identity(q, q) - V_trunc * V_trunc.transpose();

  // For now, we assume no missing data (perc_nonmissing = 1).
  // double perc_nonmissing = (Y_target.array().isNaN().count() > 0) ? ( (p*q - num_missing) / double(p*q)) : 1.0;
  double perc_nonmissing = 1.0;
  bool missing = Y_target.hasNaN();

  double obj_init = 0.0;
  {
    Eigen::MatrixXd diff = U * V.transpose() - Y_target;
    if (missing) {
      diff = diff.array().isNaN().select(Eigen::MatrixXd::Zero(diff.rows(), diff.cols()), diff);
    }
    obj_init = diff.squaredNorm() / perc_nonmissing +
      lambda1 * (P_U * U).squaredNorm() +
      lambda1 * (P_V * V).squaredNorm() +
      lambda2 * (U.transpose() * U - V.transpose() * V).squaredNorm();
  }

  double obj_best = obj_init;
  Eigen::MatrixXd U_best = U;
  Eigen::MatrixXd V_best = V;
  double U_norm = U.norm();
  double V_norm = V.norm();

  int convergence_criterion = 2;
  std::vector<double> obj_values;
  obj_values.reserve(max_iter);
  for (int iter = 0; iter < max_iter; ++iter) {
    Eigen::MatrixXd U_tilde = U.transpose() * U;
    Eigen::MatrixXd V_tilde = V.transpose() * V;

    Eigen::MatrixXd adjusted_theta = Y_target;
    if (missing) {
      Eigen::MatrixXd temp = U * V.transpose();
      adjusted_theta = adjusted_theta.array().isNaN().select(temp, adjusted_theta);
    }

    Eigen::MatrixXd grad_U = (2.0 / perc_nonmissing) * (U * V_tilde - adjusted_theta * V)
      + lambda1 * 2 * P_U * U
    + lambda2 * 4 * U * (U_tilde - V_tilde);
    Eigen::MatrixXd grad_V = (2.0 / perc_nonmissing) * (V * U_tilde - adjusted_theta.transpose() * U)
      + lambda1 * 2 * P_V * V
    + lambda2 * 4 * V * (V_tilde - U_tilde);

    double grad_U_norm = grad_U.norm();
    double grad_V_norm = grad_V.norm();

    U = U - (step_size * U_norm / (grad_U_norm + 1e-12)) * grad_U;
    V = V - (step_size * V_norm / (grad_V_norm + 1e-12)) * grad_V;
    U_norm = U.norm();
    V_norm = V.norm();

    double obj = 0.0;
    {
      Eigen::MatrixXd diff = U * V.transpose() - Y_target;
      if (missing) {
        diff = diff.array().isNaN().select(Eigen::MatrixXd::Zero(diff.rows(), diff.cols()), diff);
      }
      obj = diff.squaredNorm() / perc_nonmissing +
        lambda1 * (P_U * U).squaredNorm() +
        lambda1 * (P_V * V).squaredNorm() +
        lambda2 * (U.transpose() * U - V.transpose() * V).squaredNorm();
    }
    obj_values.push_back(obj);

    if (obj < obj_best) {
      obj_best = obj;
      U_best = U;
      V_best = V;
    }

    // Simple convergence check (you might refine this).
    if (iter > 0 && std::abs(obj - obj_values[iter - 1]) < threshold) {
      convergence_criterion = 1;
      break;
    }
    if (iter > 0 && obj > max_value * obj_init) {
      convergence_criterion = 3;
      break;
    }
    obj_init = obj;
  }

  return List::create(
    Named("learner_estimate") = U_best * V_best.transpose(),
    Named("objective_values") = obj_values,
    Named("convergence_criterion") = convergence_criterion
  );
}

//--------------------------------------------------------------
// Exported learner function.
// [[Rcpp::export]]
List learner_cpp(const Eigen::MatrixXd &Y_source, const Eigen::MatrixXd &Y_target,
                 int r, double lambda1, double lambda2, double step_size, int max_iter, double threshold, double max_value) {
  List out = learner_worker(Y_source, Y_target, r, lambda1, lambda2, step_size, max_iter, threshold, max_value);
  return List::create(
    Named("learner_estimate") = out["learner_estimate"],
    Named("objective_values") = out["objective_values"],
    Named("convergence_criterion") = out["convergence_criterion"],
    Named("r") = r
  );
}

//--------------------------------------------------------------
// [[Rcpp::export]]
List cv_learner_cpp(const Eigen::MatrixXd &Y_source, const Eigen::MatrixXd &Y_target,
                    const std::vector<double> &lambda1_all, const std::vector<double> &lambda2_all,
                    double step_size, int n_folds, int max_iter, double threshold,
                    int n_cores, int r, double max_value) {
  int p = Y_source.rows();
  int q = Y_source.cols();
  int n_lambda1 = lambda1_all.size();
  int n_lambda2 = lambda2_all.size();

  std::vector<int> indices(p * q);
  std::iota(indices.begin(), indices.end(), 0);
  std::shuffle(indices.begin(), indices.end(), std::mt19937{std::random_device{}()});

  std::vector<std::vector<int>> index_set(n_folds);
  for (int fold = 0; fold < n_folds; ++fold) {
    int start = fold * (indices.size() / n_folds);
    int end = (fold + 1) * (indices.size() / n_folds);
    index_set[fold] = std::vector<int>(indices.begin() + start, indices.begin() + end);
  }

  Eigen::MatrixXd mse_all = Eigen::MatrixXd::Zero(n_lambda1, n_lambda2);

#pragma omp parallel for collapse(2) schedule(dynamic) num_threads(n_cores)
  for (int i = 0; i < n_lambda1; ++i) {
    for (int j = 0; j < n_lambda2; ++j) {
      double mse = 0.0;

      for (int fold = 0; fold < n_folds; ++fold) {
        Eigen::MatrixXd Y_train = Y_target;
        for (int idx : index_set[fold]) {
          int row = idx / q;
          int col = idx % q;
          Y_train(row, col) = std::numeric_limits<double>::quiet_NaN();
        }
        List temp = learner_worker(Y_source, Y_train, r,
                                                          lambda1_all[i], lambda2_all[j],
                                                                                     step_size, max_iter, threshold, max_value);
        Eigen::MatrixXd learner_estimate = temp["learner_estimate"];
        double fold_mse = 0.0;
        for (int idx : index_set[fold]) {
          int row = idx / q;
          int col = idx % q;
          double diff = learner_estimate(row, col) - Y_target(row, col);
          fold_mse += diff * diff;
        }
        mse += fold_mse;
      }
      mse_all(i, j) = mse;
    }
  }

  Eigen::Index minRow, minCol;
  mse_all.minCoeff(&minRow, &minCol);

  return List::create(
    Named("lambda1_min") = lambda1_all[minRow],
                                      Named("lambda2_min") = lambda2_all[minCol],
                                                                        Named("mse_all") = mse_all,
                                                                        Named("r") = r
  );
}

// --------------------------------------------------------------
// [[Rcpp::export]]
int omp_max_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

