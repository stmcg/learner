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

// Source decompositions are prepared once on the R thread and shared read-only.
struct SourceSpaces {
  std::vector<Eigen::MatrixXd> u, v;
  std::vector<double> weights;
  Eigen::MatrixXd initial_u, initial_v;
};

SourceSpaces prepare_spaces(const Rcpp::List &sources,
                            const std::vector<double> &weights, int r) {
  SourceSpaces spaces;
  spaces.weights = weights;
  for (R_xlen_t k = 0; k < sources.size(); ++k) {
    Eigen::MatrixXd source = Rcpp::as<Eigen::MatrixXd>(sources[k]);
    Eigen::BDCSVD<Eigen::MatrixXd> svd(source, Eigen::ComputeThinU | Eigen::ComputeThinV);
    spaces.u.push_back(svd.matrixU().leftCols(r));
    spaces.v.push_back(svd.matrixV().leftCols(r));
    if (k == 0) {
      Eigen::MatrixXd root = svd.singularValues().head(r).array().sqrt().matrix().asDiagonal();
      spaces.initial_u = spaces.u[0] * root;
      spaces.initial_v = spaces.v[0] * root;
    }
  }
  return spaces;
}

// Apply I - sum(w_k B_k B_k') without allocating a p-by-p projector.
Eigen::MatrixXd space_residual(const std::vector<Eigen::MatrixXd> &bases,
                               const std::vector<double> &weights,
                               const Eigen::MatrixXd &value) {
  Eigen::MatrixXd projected = Eigen::MatrixXd::Zero(value.rows(), value.cols());
  for (size_t k = 0; k < bases.size(); ++k) {
    projected.noalias() += weights[k] * (bases[k] * (bases[k].transpose() * value));
  }
  return value - projected;
}

// Workers use no R objects or allocations, including inside the OpenMP loop.
struct FitResult {
  Eigen::MatrixXd estimate;
  std::vector<double> objective;
  int criterion;
};

FitResult learner_worker(const SourceSpaces &spaces,
                         const Eigen::MatrixXd &Y_target,
                         double lambda1_row, double lambda1_col, double lambda2,
                         double step_size, int max_iter, double threshold,
                         double max_value) {
  Eigen::MatrixXd U = spaces.initial_u;
  Eigen::MatrixXd V = spaces.initial_v;
  double perc_nonmissing = 1.0 - (static_cast<double>((Y_target.array().isNaN()).count()) / Y_target.size());
  bool missing = Y_target.hasNaN();

  double obj_init = 0.0;
  {
    Eigen::MatrixXd diff = U * V.transpose() - Y_target;
    if (missing) {
      diff = diff.array().isNaN().select(Eigen::MatrixXd::Zero(diff.rows(), diff.cols()), diff);
    }
    obj_init = diff.squaredNorm() / perc_nonmissing +
      lambda1_row * space_residual(spaces.u, spaces.weights, U).squaredNorm() +
      lambda1_col * space_residual(spaces.v, spaces.weights, V).squaredNorm() +
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

    Eigen::MatrixXd residual_u = space_residual(spaces.u, spaces.weights, U);
    Eigen::MatrixXd residual_v = space_residual(spaces.v, spaces.weights, V);
    // A weighted average of projectors is generally not idempotent.
    // The derivative of ||A U||^2 is 2 A' A U; here A is symmetric.
    Eigen::MatrixXd penalty_u = spaces.u.size() == 1 ? residual_u :
      space_residual(spaces.u, spaces.weights, residual_u);
    Eigen::MatrixXd penalty_v = spaces.v.size() == 1 ? residual_v :
      space_residual(spaces.v, spaces.weights, residual_v);
    Eigen::MatrixXd grad_U = (2.0 / perc_nonmissing) * (U * V_tilde - adjusted_theta * V)
      + lambda1_row * 2 * penalty_u
      + lambda2 * 4 * U * (U_tilde - V_tilde);
    Eigen::MatrixXd grad_V = (2.0 / perc_nonmissing) * (V * U_tilde - adjusted_theta.transpose() * U)
      + lambda1_col * 2 * penalty_v
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
        lambda1_row * space_residual(spaces.u, spaces.weights, U).squaredNorm() +
        lambda1_col * space_residual(spaces.v, spaces.weights, V).squaredNorm() +
        lambda2 * (U.transpose() * U - V.transpose() * V).squaredNorm();
    }
    obj_values.push_back(obj);
    if (!std::isfinite(obj) || !U.allFinite() || !V.allFinite()) {
      return FitResult{U_best * V_best.transpose(), obj_values, 3};
    }

    if (obj < obj_best) {
      obj_best = obj;
      U_best = U;
      V_best = V;
    }

    // Checking for convergence conditions
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

  return FitResult{U_best * V_best.transpose(), obj_values, convergence_criterion};
}

// [[Rcpp::export]]
List learner_cpp(const Rcpp::List &Y_source, const Eigen::MatrixXd &Y_target,
                 const std::vector<double> &source_weights, int r,
                 double lambda1_row, double lambda1_col, double lambda2,
                 double step_size, int max_iter, double threshold, double max_value) {
  SourceSpaces spaces = prepare_spaces(Y_source, source_weights, r);
  FitResult out = learner_worker(spaces, Y_target, lambda1_row, lambda1_col,
                                 lambda2, step_size, max_iter, threshold, max_value);
  return List::create(Named("learner_estimate") = out.estimate,
                      Named("objective_values") = out.objective,
                      Named("convergence_criterion") = out.criterion,
                      Named("r") = r);
}

// All three vectors contain one entry per candidate combination, in R array order.
// [[Rcpp::export]]
std::vector<double> cv_learner_cpp(const Rcpp::List &Y_source,
                                  const Eigen::MatrixXd &Y_target,
                                  const std::vector<double> &source_weights,
                                  const std::vector<double> &row_grid,
                                  const std::vector<double> &col_grid,
                                  const std::vector<double> &balance_grid,
                                  double step_size, int max_iter, double threshold,
                                  int n_cores, int r, double max_value,
                                  const std::vector<std::vector<int>> &index_set) {
  SourceSpaces spaces = prepare_spaces(Y_source, source_weights, r);
  int p = Y_target.rows();
  int n_grid = row_grid.size();
  std::vector<double> scores(n_grid, 0.0);

#pragma omp parallel for schedule(dynamic) num_threads(n_cores)
  for (int g = 0; g < n_grid; ++g) {
    double score = 0.0;
    for (const auto &fold : index_set) {
      Eigen::MatrixXd training = Y_target;
      for (int idx : fold) training(idx % p, idx / p) = std::numeric_limits<double>::quiet_NaN();
      FitResult fit = learner_worker(spaces, training, row_grid[g], col_grid[g],
                                      balance_grid[g], step_size, max_iter, threshold, max_value);
      if (!fit.estimate.allFinite() ||
          !std::all_of(fit.objective.begin(), fit.objective.end(),
                       [](double x) { return std::isfinite(x); })) {
        score = std::numeric_limits<double>::infinity();
        break;
      }
      double fold_score = 0.0;
      for (int idx : fold) {
        double difference = fit.estimate(idx % p, idx / p) - Y_target(idx % p, idx / p);
        fold_score += difference * difference;
      }
      score += fold_score;
    }
    scores[g] = score;
  }
  return scores;
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
