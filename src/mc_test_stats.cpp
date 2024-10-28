#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;


// Function to calculate the modified Frobenius norm
// [[Rcpp::export]]
double mc_indpd_stats_cpp(const Eigen::MatrixXd &dat) {
  Eigen::MatrixXd centered = dat.rowwise() - dat.colwise().mean();
  Eigen::MatrixXd cov_mat = (centered.adjoint() * centered) / double(dat.rows() - 1);
  cov_mat.diagonal().array() -= 1.0;
  double TT = cov_mat.norm() / sqrt(dat.cols() * (dat.cols() - 0.5));
  return TT;
}


// Function to permute the columns of the dataset independently
// [[Rcpp::export]]
Eigen::MatrixXd colwise_perm_cpp(const Eigen::MatrixXd &dat) {
  Eigen::MatrixXd dat_perm = dat;
  for (int j = 0; j < dat.cols(); ++j) {
    std::random_shuffle(dat_perm.col(j).data(), dat_perm.col(j).data() + dat_perm.rows());
  }
  return dat_perm;
}


// Function to permute the rows of the dataset independently
// [[Rcpp::export]]
Eigen::MatrixXd rowwise_perm_cpp(const Eigen::MatrixXd &dat) {
  Eigen::MatrixXd dat_perm = dat.transpose();
  for (int i = 0; i < dat.rows(); ++i) {
    std::random_shuffle(dat_perm.col(i).data(), dat_perm.col(i).data() + dat_perm.rows());
  }
  return dat_perm.transpose();
}


// Function to scale a matrix (similar to R's scale function)
// [[Rcpp::export]]
Eigen::MatrixXd scale_cpp(const Eigen::MatrixXd& M) {
  Eigen::MatrixXd centered = M.rowwise() - M.colwise().mean();
  Eigen::VectorXd sd = ((centered.array().square().colwise().sum()) / (M.rows() - 1)).sqrt();
  for (int j = 0; j < M.cols(); ++j) {
    if (sd(j) > 0) {
      centered.col(j) = centered.col(j).array() / sd(j);
    }
  }
  return centered;
}


// Main function replicating mc_test_stats
// [[Rcpp::export]]
Eigen::Vector3d mc_test_stats(const Eigen::MatrixXd & counts, 
                              const double gene_select_thre = 0.1) {

  // Filter rows based on the condition
  std::vector<int> selectedRows;
  for (int i = 0; i < counts.rows(); ++i) {
    int nonZeroCount = (counts.row(i).array() > 0).count();
    if (nonZeroCount > counts.cols() * gene_select_thre) {
      selectedRows.push_back(i);
    }
  }
  
  Eigen::MatrixXd counts_var(selectedRows.size(), counts.cols());
  for (size_t i = 0; i < selectedRows.size(); ++i) {
    counts_var.row(i) = counts.row(selectedRows[i]);
  }
  
  // Transpose and scale
  const Eigen::MatrixXd dat = scale_cpp(counts_var.transpose());
  
  // Calculate original and permuted stats
  Eigen::Vector3d res;
  res(0) = mc_indpd_stats_cpp(dat);  // T_org
  res(1) = mc_indpd_stats_cpp(colwise_perm_cpp(dat));  // T_perm
  res(2) = res(0) / res(1);  // T_org / T_perm
  
  // Return the results
  return res;
}
 
