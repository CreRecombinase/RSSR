library(testthat)
library(RSSR)
library(Matrix)
alpha_input <- scan("~/Downloads/alpha0.txt",what=numeric())
mu_input <- scan("~/Downloads/mu0.txt",what=numeric())
logodds_input <- scan("~/Downloads/logodds.txt",what=numeric())
betahat_input <- scan("~/Downloads/betahat.txt",what=numeric())
se_input <- scan("~/Downloads/se.txt",what=numeric())
SiRiSr_input <- scan("~/Downloads/SiRiSr.txt",what=numeric())
sigbeta_input <- scan("~/Downloads/sigma_beta.txt",what=numeric())
SiRiS_input <- gen_SiRSi("/media/nwknoblauch/Data/1kg/1kg_19.mat")


alpha_update_matlab <- scan("~/Dropbox/RSSR/analyses/alpha_update.txt",what=numeric())
mu_update_matlab <- scan("~/Dropbox/RSSR/analyses/mu_update.txt",what=numeric())
SiRiSr_update_matlab <- scan("~/Dropbox/RSSR/analyses/SiRiSr_update.txt",what=numeric())


system.time(ret <- rss_varbvsr(SiRiS  =SiRiS_input,
                               sigma_beta = sigbeta_input,
                               logodds = logodds_input,
                               betahat = betahat_input,
                               se = se_input,
                               alpha = alpha_input,
                               mu = mu_input,
                               SiRiSr = SiRiSr_input))

alpha_update_R <- ret[,1]
mu_update_R <- ret[,2]
SiRiSr_update_R <- ret[,3]

expect_equal(alpha_update_R,alpha_update_matlab,tolerance=1e-6)
expect_equal(mu_update_matlab,mu_update_R,tolerance=1e-6)
expect_equal(SiRiSr_update_matlab,SiRiSr_update_R)
