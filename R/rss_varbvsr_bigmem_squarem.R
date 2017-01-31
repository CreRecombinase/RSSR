rss_varbvsr_bigmem_squarem <- function(datafile,sigb=0.058,logodds=-2.9/log(10),options=list()){
  require(h5)
  tolerance <- 1e-4
  sumstat <- h5file(datafile,'r')
  betahat <- sumstat['betahat'][]
  p <- length(betahat)
  se <- sumstat['se'][]
  h5close(sumstat)
  alpha <- options[["alpha"]]
  mu <- options[["mu"]]
  if(is.null(alpha)){
    alpha <- runif(p)
    alpha <- alpha/sum(alpha)
  }
  if(is.null(mu)){
    mu <- rnorm(p)
  }
  new()
  oSiRiS_tmp <- as(SiRiS_tmp,"dgTMatrix")
  SiRiS_tmp <- gen_SiRSi(datafile)
  SiRiSr_cell <- SiRiS_tmp%*%(alpha*mu)
  sigbsquare <- sigb*sigb
  sesquare <- se*se
  q <- betahat/sesquare
  s <- (sesquare*sigbsquare)/(sesquare+sigbsquare)

  lnZ <- -Inf
  iter <- 0
  loglik <- 0
  cont <- T
  
  while(cont){
    iter <- iter +1
    lnZ0 <- lnZ
    alpha0 <- alpha
    mu0 <- mu
    SiRiSr0 <- SiRiSr_cell@x
    system.time(res <- rss_varbvsr(SiRiS  =SiRiS_tmp,
                                   sigma_beta = sigb,
                                   logodds = logodds,
                                   betahat = betahat,
                                   se = se,
                                   alpha = alpha0,
                                   mu = mu0,
                                   SiRiSr = SiRiSr0))
#    rss_va(betahat = )
    maxerr <- 0
    absr <- 0
    asum <- 0
  
  }
}