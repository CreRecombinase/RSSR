# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rcppeigen_hello_world <- function() {
    .Call('RSSR_rcppeigen_hello_world', PACKAGE = 'RSSR')
}

rcppeigen_outerproduct <- function(x) {
    .Call('RSSR_rcppeigen_outerproduct', PACKAGE = 'RSSR', x)
}

rcppeigen_innerproduct <- function(x) {
    .Call('RSSR_rcppeigen_innerproduct', PACKAGE = 'RSSR', x)
}

rcppeigen_bothproducts <- function(x) {
    .Call('RSSR_rcppeigen_bothproducts', PACKAGE = 'RSSR', x)
}

rss_varbvsr <- function(SiRiS, sigma_beta, logodds, betahat, se, alpha0, mu0, SiRiSr0) {
    .Call('RSSR_rss_varbvsr', PACKAGE = 'RSSR', SiRiS, sigma_beta, logodds, betahat, se, alpha0, mu0, SiRiSr0)
}

rss_varbvsr_update <- function(betahat, se, sigma_beta, SiRiS_snp, SiRiSr, SiRiSr_snp, logodds, alpha, mu) {
    invisible(.Call('RSSR_rss_varbvsr_update', PACKAGE = 'RSSR', betahat, se, sigma_beta, SiRiS_snp, SiRiSr, SiRiSr_snp, logodds, alpha, mu))
}

SiRSi <- function(R, Si) {
    .Call('RSSR_SiRSi', PACKAGE = 'RSSR', R, Si)
}

