# Model selection

rm(list=ls())
library(asreml)
library(tidyverse)
source('Functions_MTME.R')

# Load model ----
load('Data/MTME_NFA5.RData')

# Get gebvs ----
mod <- MTME_NFA5.asr
k <- 5
data <- ILYT_Pheno
#pev <- MTME_NFA5.pev

# variance parameters
(vparams <- mod$vparameters)
# latent environmental covariates (loadings)
(Lam <- matrix(vparams[grep("^rr.*fa", names(vparams))], ncol = k))
rownames(Lam) <- levels(data$TraitEnv)
# specific variances
(Psi <- diag(vparams[grep("^TraitEnv.*vm", names(vparams))]))
# in this example we dont have much/any specific variance.
colnames(Psi) <- rownames(Psi) <- levels(data$TraitEnv)

# coefficients
coefs <- mod$coefficients$random
# genotype scores (slopes)
f <- coefs[grep('Comp', rownames(coefs))]
Lamf <- c(matrix(f, ncol = k) %*% t(Lam)) # common GET effects
# genetic regressions residuals
delta <- coefs[grep(paste0("^TraitEnv.*vm"), rownames(coefs))] # specific GET effects
# total GET effects
Lamfdelta <- Lamf + delta
Lamfdelta

# get gebv for each genotype, by trait-environment combination

# vector of GET
effects_all <- rownames(coefs)
head(effects_all)
effects_GET <- effects_all[grepl(paste0("^TraitEnv.*vm"), effects_all)]
head(effects_GET)
GET <- sub(paste0('TraitEnv_'), '', effects_GET)
GET <- sub(':vm\\(.*?\\)', '', GET)
head(GET)

split_GET <- do.call(rbind, strsplit(GET, '_', fixed = TRUE))

data |> group_by(TraitEnv) |>
  summarise(Pheno_mean=mean(Pheno_mean),
            Pheno_sd=mean(Pheno_sd),
            Trait=unique(Trait),
            Env=unique(Env)
            ) |>
  glimpse()

gebv <- data.frame(TraitEnv=split_GET[,1],
                   Gkeep=split_GET[,2],
                   blup=Lamfdelta) |>
  left_join(
    data |> group_by(TraitEnv) |>
      summarise(Pheno_mean=mean(Pheno_mean),
                Pheno_sd=mean(Pheno_sd),
                Trait=unique(Trait),
                Env=unique(Env)
      )
  ) |>
  mutate(gebv=blup*Pheno_sd) |>
  select(Gkeep, Trait, Env, TraitEnv, gebv) |>
  mutate_if(is.character, ~as.factor(.)) |>
  glimpse()

# Accuracy ----
# additive effects
# cve effects
(BLUPs_cve <- mod$coef$random[grep("rr.*vm", rownames(mod$coef$random)),])
(PEV_cve_diag <- mod$vcoef$random[grep("rr.*vm", rownames(mod$coef$random))])
(PEV_cve_diag <- PEV_cve_diag[grep("Comp", names(BLUPs_cve), invert = T)])
(BLUPs_cve <- BLUPs_cve[grep("Comp", names(BLUPs_cve), invert = T)])
# specific effects
(BLUPs_sve <- mod$coef$random[grep("^TraitEnv.*vm", rownames(mod$coef$random)),])
(PEV_sve_diag <- mod$vcoef$random[grep("^TraitEnv.*vm", rownames(mod$coef$random))])
# total (common + specific) effects
(BLUPs_tve <- BLUPs_cve + BLUPs_sve)
plot(BLUPs_tve, BLUPs_cve)
(PEV_tve_diag_approx <- PEV_cve_diag + PEV_sve_diag)
# non-additive effects
(BLUPs_ide <- mod$coef$random[grep("^TraitEnv.*ide", rownames(mod$coef$random)),])

# Common + specific GET effects
head(pev$pvals)
PEV_tve <- as.matrix(tve_pev$vcov)

# Quick check of BLUPs and PEVs
# BLUPs first
plot(BLUPs_tve, tve_pev$pvals$predicted.value); abline(a=0, b=1)
range(BLUPs_tve - tve_pev$pvals$predicted.value)
# PEVs next
plot(diag(PEV_tve), tve_pev$pvals$std.error^2); abline(a=0, b=1)
plot(diag(PEV_tve), PEV_tve_diag_approx); abline(a=0, b=1)
range(diag(PEV_tve) - PEV_tve_diag_approx)

# obtain accuracies
VAR_tve <- kronecker(Lam %*% t(Lam) + Psi, Ginv)
ACC_tve <- sqrt(1 - diag(PEV_tve)/diag(VAR_tve))
plot(ACC_tve, ACC_cve) # difference arising for that one environment where 
# Psi is non-zero, as expected.
# note that this assumes the variance is given by g_ii (lamlam' + Psi),
# where g_ii is the ith diagonal element of G
plot(ACC_tve, BLUPs_tve)
