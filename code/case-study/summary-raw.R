# summary-raw.R: this script summarizes the full posterior chains coming 
#                from all candidate models fit in the case study. No final
#                figures are produced from this script. Rather, the goal
#                of this script is to assess convergence, extract samples 
#                used for posterior summaries, and compare estimates
#                across the different models
rm(list = ls())
library(coda)
library(spOccupancy)
library(tidyverse)

# Load results from all models --------------------------------------------
# lfJSDM ---------------------------
load("results/bbs-lfJSDM-1-chain-2022-03-27.R")
out.1.lfJSDM <- out
load("results/bbs-lfJSDM-2-chain-2022-03-27.R")
out.2.lfJSDM <- out
load("results/bbs-lfJSDM-3-chain-2022-03-27.R")
out.3.lfJSDM <- out
# sfJSDM ---------------------------
load("results/bbs-sfJSDM-1-chain-2022-03-28.R")
out.1.sfJSDM <- out
load("results/bbs-sfJSDM-2-chain-2022-03-28.R")
out.2.sfJSDM <- out
load("results/bbs-sfJSDM-3-chain-2022-03-28.R")
out.3.sfJSDM <- out
# msPGOcc ---------------------------
load("results/bbs-msPGOcc-1-chain-2022-03-29.R")
out.1.msPGOcc <- out
load("results/bbs-msPGOcc-2-chain-2022-03-29.R")
out.2.msPGOcc <- out
load("results/bbs-msPGOcc-3-chain-2022-03-28.R")
out.3.msPGOcc <- out
# lfMsPGOcc ---------------------------
load("results/bbs-lfMsPGOcc-1-chain-2022-03-28.R")
out.1.lfMsPGOcc <- out
load("results/bbs-lfMsPGOcc-2-chain-2022-03-28.R")
out.2.lfMsPGOcc <- out
load("results/bbs-lfMsPGOcc-3-chain-2022-03-28.R")
out.3.lfMsPGOcc <- out
# sfMsPGOcc ---------------------------
load("results/bbs-sfMsPGOcc-1-chain-2022-03-29.R")
out.1.sfMsPGOcc <- out
load("results/bbs-sfMsPGOcc-2-chain-2022-03-29.R")
out.2.sfMsPGOcc <- out
load("results/bbs-sfMsPGOcc-3-chain-2022-03-29.R")
out.3.sfMsPGOcc <- out

# Get initial values based on small model runs ----------------------------
# lfJSDM ------------------------------
# beta.comm.inits <- apply(out.1.lfJSDM$beta.comm.samples, 2, mean)
# tau.sq.beta.inits <- apply(out.1.lfJSDM$tau.sq.beta.samples, 2, mean)
# p.occ <- length(beta.comm.inits)
# N <- nrow(out.1.lfJSDM$y)
# beta.inits <- matrix(apply(out.1.lfJSDM$beta.samples, 2, mean), N, p.occ)
# q <- 5
# sigma.sq.psi.inits <- apply(out.1.lfJSDM$sigma.sq.psi.samples, 2, mean)
# lambda.inits <- matrix(apply(out.1.lfJSDM$lambda.samples, 2, mean), N, q)
# inits.lfJSDM <- list(beta.comm = beta.comm.inits, 
# 		     tau.sq.beta = tau.sq.beta.inits, 
# 		     beta = beta.inits, 
# 		     sigma.sq.psi = sigma.sq.psi.inits, 
# 		     lambda = lambda.inits)
# save(inits.lfJSDM, file = 'data/inits-lfJSDM.rda')
# # sfJSDM ------------------------------
# beta.comm.inits <- apply(out.1.sfJSDM$beta.comm.samples, 2, mean)
# tau.sq.beta.inits <- apply(out.1.sfJSDM$tau.sq.beta.samples, 2, mean)
# p.occ <- length(beta.comm.inits)
# N <- nrow(out.1.sfJSDM$y)
# beta.inits <- matrix(apply(out.1.sfJSDM$beta.samples, 2, mean), N, p.occ)
# phi.inits <- apply(out.1.sfJSDM$theta.samples, 2, mean)
# sigma.sq.psi.inits <- apply(out.1.sfJSDM$sigma.sq.psi.samples, 2, mean)
# lambda.inits <- matrix(apply(out.1.sfJSDM$lambda.samples, 2, mean), N, q)
# inits.sfJSDM <- list(beta.comm = beta.comm.inits, 
# 		     tau.sq.beta = tau.sq.beta.inits, 
# 		     beta = beta.inits, 
# 		     phi = phi.inits, 
# 		     sigma.sq.psi = sigma.sq.psi.inits,
# 		     lambda = lambda.inits)
# save(inits.sfJSDM, file = 'data/inits-sfJSDM.rda')
# # msPGOcc ------------------------------
# beta.comm.inits <- apply(out.1.msPGOcc$beta.comm.samples, 2, mean)
# tau.sq.beta.inits <- apply(out.1.msPGOcc$tau.sq.beta.samples, 2, mean)
# alpha.comm.inits <- apply(out.1.msPGOcc$alpha.comm.samples, 2, mean)
# tau.sq.alpha.inits <- apply(out.1.msPGOcc$tau.sq.alpha.samples, 2, mean)
# p.occ <- length(beta.comm.inits)
# p.det <- length(alpha.comm.inits)
# N <- nrow(out.1.msPGOcc$y)
# beta.inits <- matrix(apply(out.1.msPGOcc$beta.samples, 2, mean), N, p.occ)
# alpha.inits <- matrix(apply(out.1.msPGOcc$alpha.samples, 2, mean), N, p.det)
# sigma.sq.p.inits <- apply(out.1.msPGOcc$sigma.sq.p.samples, 2, mean)
# z.inits <- apply(out.1.msPGOcc$y, c(1, 2), max, na.rm = TRUE)
# inits.msPGOcc <- list(beta.comm = beta.comm.inits, 
# 		      tau.sq.beta = tau.sq.beta.inits, 
# 		      beta = beta.inits, 
# 		      alpha.comm = alpha.comm.inits, 
# 		      tau.sq.alpha = tau.sq.alpha.inits,
# 		      alpha = alpha.inits,
# 		      sigma.sq.p = sigma.sq.p.inits, 
# 		      z = z.inits)
# save(inits.msPGOcc, file = 'data/inits-msPGOcc.rda')
# # lfMsPGOcc ------------------------------
# beta.comm.inits <- apply(out.1.lfMsPGOcc$beta.comm.samples, 2, mean)
# tau.sq.beta.inits <- apply(out.1.lfMsPGOcc$tau.sq.beta.samples, 2, mean)
# alpha.comm.inits <- apply(out.1.lfMsPGOcc$alpha.comm.samples, 2, mean)
# tau.sq.alpha.inits <- apply(out.1.lfMsPGOcc$tau.sq.alpha.samples, 2, mean)
# p.occ <- length(beta.comm.inits)
# p.det <- length(alpha.comm.inits)
# N <- nrow(out.1.lfMsPGOcc$y)
# beta.inits <- matrix(apply(out.1.lfMsPGOcc$beta.samples, 2, mean), N, p.occ)
# alpha.inits <- matrix(apply(out.1.lfMsPGOcc$alpha.samples, 2, mean), N, p.det)
# sigma.sq.p.inits <- apply(out.1.lfMsPGOcc$sigma.sq.p.samples, 2, mean)
# lambda.inits <- matrix(apply(out.1.lfMsPGOcc$lambda.samples, 2, mean), N, q)
# z.inits <- apply(out.1.lfMsPGOcc$y, c(1, 2), max, na.rm = TRUE)
# inits.lfMsPGOcc <- list(beta.comm = beta.comm.inits, 
# 		        tau.sq.beta = tau.sq.beta.inits, 
# 		        beta = beta.inits, 
# 		        alpha.comm = alpha.comm.inits, 
# 		        tau.sq.alpha = tau.sq.alpha.inits,
# 		        alpha = alpha.inits,
# 		        sigma.sq.p = sigma.sq.p.inits, 
# 		        lambda = lambda.inits,
# 		        z = z.inits)
# save(inits.lfMsPGOcc, file = 'data/inits-lfMsPGOcc.rda')
# # sfMsPGOcc ------------------------------
# beta.comm.inits <- apply(out.1.sfMsPGOcc$beta.comm.samples, 2, mean)
# tau.sq.beta.inits <- apply(out.1.sfMsPGOcc$tau.sq.beta.samples, 2, mean)
# alpha.comm.inits <- apply(out.1.sfMsPGOcc$alpha.comm.samples, 2, mean)
# tau.sq.alpha.inits <- apply(out.1.sfMsPGOcc$tau.sq.alpha.samples, 2, mean)
# p.occ <- length(beta.comm.inits)
# p.det <- length(alpha.comm.inits)
# N <- nrow(out.1.sfMsPGOcc$y)
# beta.inits <- matrix(apply(out.1.sfMsPGOcc$beta.samples, 2, mean), N, p.occ)
# alpha.inits <- matrix(apply(out.1.sfMsPGOcc$alpha.samples, 2, mean), N, p.det)
# sigma.sq.p.inits <- apply(out.1.sfMsPGOcc$sigma.sq.p.samples, 2, mean)
# lambda.inits <- matrix(apply(out.1.sfMsPGOcc$lambda.samples, 2, mean), N, q)
# phi.inits <- apply(out.1.sfMsPGOcc$theta.samples, 2, mean)
# z.inits <- apply(out.1.sfMsPGOcc$y, c(1, 2), max, na.rm = TRUE)
# inits.sfMsPGOcc <- list(beta.comm = beta.comm.inits, 
# 		        tau.sq.beta = tau.sq.beta.inits, 
# 		        beta = beta.inits, 
# 		        alpha.comm = alpha.comm.inits, 
# 		        tau.sq.alpha = tau.sq.alpha.inits,
# 		        alpha = alpha.inits,
# 		        sigma.sq.p = sigma.sq.p.inits, 
# 			phi = phi.inits,
# 		        lambda = lambda.inits,
# 		        z = z.inits)
# save(inits.sfMsPGOcc, file = 'data/inits-sfMsPGOcc.rda')

# Assess Convergence ------------------------------------------------------
# beta --------------------------------
out.lfJSDM.beta.samples <- mcmc.list(out.1.lfJSDM$beta.samples, 
				     out.2.lfJSDM$beta.samples, 
				     out.3.lfJSDM$beta.samples)
out.sfJSDM.beta.samples <- mcmc.list(out.1.sfJSDM$beta.samples, 
				     out.2.sfJSDM$beta.samples, 
				     out.3.sfJSDM$beta.samples)
out.msPGOcc.beta.samples <- mcmc.list(out.1.msPGOcc$beta.samples, 
				     out.2.msPGOcc$beta.samples, 
				     out.3.msPGOcc$beta.samples)
out.lfMsPGOcc.beta.samples <- mcmc.list(out.1.lfMsPGOcc$beta.samples, 
				     out.2.lfMsPGOcc$beta.samples, 
				     out.3.lfMsPGOcc$beta.samples)
out.sfMsPGOcc.beta.samples <- mcmc.list(out.2.sfMsPGOcc$beta.samples, 
				        out.3.sfMsPGOcc$beta.samples)
# Rhat
r.hat.beta.samples.lfJSDM <- gelman.diag(out.lfJSDM.beta.samples)
sum(r.hat.beta.samples.lfJSDM$psrf[, 2] > 1.1)
r.hat.beta.samples.sfJSDM <- gelman.diag(out.sfJSDM.beta.samples)
sum(r.hat.beta.samples.sfJSDM$psrf[, 2] > 1.1)
r.hat.beta.samples.msPGOcc <- gelman.diag(out.msPGOcc.beta.samples)
sum(r.hat.beta.samples.msPGOcc$psrf[, 2] > 1.1)
r.hat.beta.samples.lfMsPGOcc <- gelman.diag(out.lfMsPGOcc.beta.samples)
sum(r.hat.beta.samples.lfMsPGOcc$psrf[, 2] > 1.1)
r.hat.beta.samples.sfMsPGOcc <- gelman.diag(out.sfMsPGOcc.beta.samples)
sum(r.hat.beta.samples.sfMsPGOcc$psrf[, 2] > 1.1)
# ESS
ess.beta.lfJSDM <- effectiveSize(out.lfJSDM.beta.samples)
ess.beta.sfJSDM <- effectiveSize(out.sfJSDM.beta.samples)
ess.beta.msPGOcc <- effectiveSize(out.msPGOcc.beta.samples)
ess.beta.lfMsPGOcc <- effectiveSize(out.lfMsPGOcc.beta.samples)
ess.beta.sfMsPGOcc <- effectiveSize(out.sfMsPGOcc.beta.samples)
# beta.comm --------------------------------
out.lfJSDM.beta.comm.samples <- mcmc.list(out.1.lfJSDM$beta.comm.samples, 
				     out.2.lfJSDM$beta.comm.samples, 
				     out.3.lfJSDM$beta.comm.samples)
out.sfJSDM.beta.comm.samples <- mcmc.list(out.1.sfJSDM$beta.comm.samples, 
				     out.2.sfJSDM$beta.comm.samples, 
				     out.3.sfJSDM$beta.comm.samples)
out.msPGOcc.beta.comm.samples <- mcmc.list(out.1.msPGOcc$beta.comm.samples, 
				     out.2.msPGOcc$beta.comm.samples, 
				     out.3.msPGOcc$beta.comm.samples)
out.lfMsPGOcc.beta.comm.samples <- mcmc.list(out.1.lfMsPGOcc$beta.comm.samples, 
				     out.2.lfMsPGOcc$beta.comm.samples, 
				     out.3.lfMsPGOcc$beta.comm.samples)
out.sfMsPGOcc.beta.comm.samples <- mcmc.list(out.2.sfMsPGOcc$beta.comm.samples, 
				     out.3.sfMsPGOcc$beta.comm.samples)
# Rhat
r.hat.beta.comm.samples.lfJSDM <- gelman.diag(out.lfJSDM.beta.comm.samples)
sum(r.hat.beta.comm.samples.lfJSDM$psrf[, 2] > 1.1)
r.hat.beta.comm.samples.sfJSDM <- gelman.diag(out.sfJSDM.beta.comm.samples)
sum(r.hat.beta.comm.samples.sfJSDM$psrf[, 2] > 1.1)
r.hat.beta.comm.samples.msPGOcc <- gelman.diag(out.msPGOcc.beta.comm.samples)
sum(r.hat.beta.comm.samples.msPGOcc$psrf[, 2] > 1.1)
r.hat.beta.comm.samples.lfMsPGOcc <- gelman.diag(out.lfMsPGOcc.beta.comm.samples)
sum(r.hat.beta.comm.samples.lfMsPGOcc$psrf[, 2] > 1.1)
r.hat.beta.comm.samples.sfMsPGOcc <- gelman.diag(out.sfMsPGOcc.beta.comm.samples)
sum(r.hat.beta.comm.samples.sfMsPGOcc$psrf[, 2] > 1.1)
# ESS
ess.beta.comm.lfJSDM <- effectiveSize(out.lfJSDM.beta.comm.samples)
ess.beta.comm.sfJSDM <- effectiveSize(out.sfJSDM.beta.comm.samples)
ess.beta.comm.msPGOcc <- effectiveSize(out.msPGOcc.beta.comm.samples)
ess.beta.comm.lfMsPGOcc <- effectiveSize(out.lfMsPGOcc.beta.comm.samples)
ess.beta.comm.sfMsPGOcc <- effectiveSize(out.sfMsPGOcc.beta.comm.samples)
# theta --------------------------------
out.sfJSDM.theta.samples <- mcmc.list(out.1.sfJSDM$theta.samples, 
				     out.2.sfJSDM$theta.samples, 
				     out.3.sfJSDM$theta.samples)
out.sfMsPGOcc.theta.samples <- mcmc.list(out.2.sfMsPGOcc$theta.samples, 
				     out.3.sfMsPGOcc$theta.samples)
# Rhat
r.hat.theta.samples.sfJSDM <- gelman.diag(out.sfJSDM.theta.samples)
sum(r.hat.theta.samples.sfJSDM$psrf[, 2] > 1.1)
r.hat.theta.samples.sfMsPGOcc <- gelman.diag(out.sfMsPGOcc.theta.samples)
sum(r.hat.theta.samples.sfMsPGOcc$psrf[, 2] > 1.1)
# ESS
ess.theta.sfJSDM <- effectiveSize(out.sfJSDM.theta.samples)
ess.theta.sfMsPGOcc <- effectiveSize(out.sfMsPGOcc.theta.samples)
geweke.diag(mcmc(out.1.lfJSDM$w.samples[, , 1200]))
pdf("traceplots.pdf")
plot(out.sfMsPGOcc.theta.samples, density = FALSE)
dev.off()
# lambda --------------------------------
tmp <- function(a) {
  ifelse(var(a) > 0, 1, 0)
}
# Get indices of the lower triangular elements
lower.triag.indx <- which(apply(out.1.sfJSDM$lambda.samples, 2, tmp) == 1)
out.lfJSDM.lambda.samples <- mcmc.list(out.1.lfJSDM$lambda.samples[, lower.triag.indx], 
				     out.2.lfJSDM$lambda.samples[, lower.triag.indx], 
				     out.3.lfJSDM$lambda.samples[, lower.triag.indx])
out.sfJSDM.lambda.samples <- mcmc.list(out.1.sfJSDM$lambda.samples[, lower.triag.indx], 
				     out.2.sfJSDM$lambda.samples[, lower.triag.indx], 
				     out.3.sfJSDM$lambda.samples[, lower.triag.indx])
out.lfMsPGOcc.lambda.samples <- mcmc.list(out.1.lfMsPGOcc$lambda.samples[, lower.triag.indx], 
				     out.2.lfMsPGOcc$lambda.samples[, lower.triag.indx], 
				     out.3.lfMsPGOcc$lambda.samples[, lower.triag.indx])
out.sfMsPGOcc.lambda.samples <- mcmc.list(out.2.sfMsPGOcc$lambda.samples[, lower.triag.indx], 
				     out.3.sfMsPGOcc$lambda.samples[, lower.triag.indx])
# save(lambda.means, file = 'data/lambda-starting.rda')
# Rhat
r.hat.lambda.samples.lfJSDM <- gelman.diag(out.lfJSDM.lambda.samples)
sum(r.hat.lambda.samples.lfJSDM$psrf[, 2] > 1.1)
r.hat.lambda.samples.sfJSDM <- gelman.diag(out.sfJSDM.lambda.samples)
sum(r.hat.lambda.samples.sfJSDM$psrf[, 2] > 1.1)
r.hat.lambda.samples.lfMsPGOcc <- gelman.diag(out.lfMsPGOcc.lambda.samples)
sum(r.hat.lambda.samples.lfMsPGOcc$psrf[, 2] > 1.1)
r.hat.lambda.samples.sfMsPGOcc <- gelman.diag(out.sfMsPGOcc.lambda.samples)
sum(r.hat.lambda.samples.sfMsPGOcc$psrf[, 2] > 1.1)
# ESS
ess.lambda.lfJSDM <- effectiveSize(out.lfJSDM.lambda.samples)
ess.lambda.sfJSDM <- effectiveSize(out.sfJSDM.lambda.samples)
ess.lambda.lfMsPGOcc <- effectiveSize(out.lfMsPGOcc.lambda.samples)
ess.lambda.sfMsPGOcc <- effectiveSize(out.sfMsPGOcc.lambda.samples)
pdf("traceplots.pdf")
plot(out.sfMsPGOcc.lambda.samples[, 310:315], density = FALSE)
dev.off()

# Compare beta estimates --------------------------------------------------
beta.mean.msPGOcc <- apply(do.call('rbind', out.msPGOcc.beta.samples), 2, mean)
beta.mean.lfJSDM <- apply(do.call('rbind', 
				  out.lfJSDM.beta.samples[, 1:length(beta.mean.msPGOcc)]), 2, mean)
beta.mean.lfMsPGOcc <- apply(do.call('rbind', out.lfMsPGOcc.beta.samples), 2, mean)
beta.mean.sfMsPGOcc <- apply(do.call('rbind', out.sfMsPGOcc.beta.samples), 2, mean)
beta.mean.sfJSDM <- apply(do.call('rbind', 
				  out.sfJSDM.beta.samples[, 1:length(beta.mean.sfMsPGOcc)]), 2, mean)
pdf('traceplots.pdf')
plot(beta.mean.sfMsPGOcc[1:98], beta.mean.sfJSDM[1:98], pch = 19)
abline(0, 1)
dev.off()

pdf('traceplots.pdf')
plot(mcmc(out$w.samples[, 1, 101:110]), density = FALSE)
# plot(out$beta.comm.samples, density = FALSE)
dev.off()

# Compute WAIC ------------------------------------------------------------
# lfJSDM
waicOcc(out.2.lfJSDM)
# sfJSDM
waicOcc(out.2.sfJSDM)
# msPGOcc
waicOcc(out.2.msPGOcc)
# lfMsPGOcc
waicOcc(out.2.lfMsPGOcc)
# sfMsPGOcc
waicOcc(out.2.sfMsPGOcc)

