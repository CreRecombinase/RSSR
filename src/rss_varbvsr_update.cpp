// For a description of this C code, see rss_varbvsr_update.m.
#include <RcppEigen.h>
#include "rssvarbvsr.hpp"

                     

// [[Rcpp::export]]
Eigen::MatrixXd rss_varbvsr (Eigen::SparseMatrix<double> SiRiS,
                  Rcpp::NumericVector sigma_beta,
                  const Eigen::VectorXd logodds,
                  const Eigen::VectorXd betahat,
                  const Eigen::VectorXd se,
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
  for (size_t j = 0; j < p; j++) {
    SiRiS_snp_v=(SiRiS.col(j));
    SiRiS_snp=SiRiS_snp_v.array();
    // Copy the kth column of matrix inv(S)*R*inv(S).
    // copySparseColumn(SiRiS_snp, k, SiRiS.elems, Ir, Jc, p);
    
    // Copy the kth element of vector inv(S)*R*inv(S)*r.
    double SiRiSr_snp = SiRiSr(j);

    // Perform the mean-field variational update.
    rss_varbvsr_update(betahat(j), se(j), sigma_beta[0], 
		       SiRiS_snp, SiRiSr, SiRiSr_snp, 
		       logodds(0), alpha(j), mu(j));
  }
  retmat << alpha,mu,SiRiSr;
  return(retmat);
  // Free the dynamically allocated memory.
}
