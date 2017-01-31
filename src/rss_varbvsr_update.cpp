// For a description of this C code, see rss_varbvsr_update.m.
#include <RcppEigen.h>
#include "rssvarbvsr.hpp"

                     

// [[Rcpp::export]]
Eigen::MatrixXd rss_varbvsr (Eigen::SparseMatrix<double> SiRiS,
                  const Eigen::ArrayXd sigma_beta,
                  const Eigen::ArrayXd logodds,
                  const Eigen::ArrayXd betahat,
                  const Eigen::ArrayXd se,
                  const Eigen::ArrayXd &alpha0,
                  const Eigen::ArrayXd &mu0,
                  const Eigen::ArrayXd &SiRiSr0
                  ){


  // Get the number of SNPs (p) and coordinate ascent updates (m).
  const size_t p = betahat.size();

  // Initialize outputs.

  // Store a single column of matrix inv(S)*R*inv(S).

  Eigen::ArrayXd alpha=alpha0;
  Eigen::ArrayXd mu=mu0;
  Eigen::ArrayXd SiRiSr=SiRiSr0;
  
  Eigen::ArrayXd  SiRiS_snp(p);
  Eigen::VectorXd  SiRiS_snp_v(p);
  Eigen::MatrixXd retmat(p,3);
  // Run coordinate ascent updates.
  // Repeat for each coordinate ascent update.
  Eigen::ArrayXd se_square = se * se;
  Eigen::ArrayXd sigma_beta_square = sigma_beta * sigma_beta;
  
  // Compute the variational estimate of the posterior variance.
  Eigen::ArrayXd  sigma_square = (se_square * sigma_beta_square) / (se_square + sigma_beta_square);
  for (size_t j = 0; j < p; j++) {
    // Copy the kth column of matrix inv(S)*R*inv(S).
    // copySparseColumn(SiRiS_snp, k, SiRiS.elems, Ir, Jc, p);
    
    // Copy the kth element of vector inv(S)*R*inv(S)*r.
    double SiRiSr_snp = SiRiSr(j);
    // Perform the mean-field variational update.
    rss_varbvsr_update(betahat(j), sigma_square(j), se_square(j), sigma_beta_square(j),
		       SiRiS.col(j), SiRiSr, SiRiSr_snp, 
		       logodds(j), alpha(j), mu(j));
  }
  retmat << alpha,mu,SiRiSr;
  return(retmat);
  // Free the dynamically allocated memory.
}
