# Functions

# Update ASReml-R ----
# This function updates ASReml-R until it converges
update_asreml <- function(mod, max_updates = 5, save_path) {
  count <- 0
  
  while (!mod$converge && count < max_updates) {
    count <- count + 1
    mod <- update(mod)
    
    # Print model update information
    print(paste('Update', count))
    print(paste('Convergence =', mod$converge))
    
    # Print LogLik value
    loglik <- mod$trace |>
      as.data.frame() |>
      rownames_to_column('Iteration') |>
      filter(Iteration == 'LogLik')
    
    print(loglik)
    
    # Save model state after each update
    save.image(save_path)
  }
  
  if (mod$converge) {
    print('Model successfully converged!')
  } else {
    print('Maximum updates reached. Model did not converge.')
  }
  
  return(mod) # Return the final model
}

# % Va----
# Percentage of additive genetic variance explained by the model

VaPct <- function(mod,k,data,TE_fct) {
  # variance parameters
  vparams <- mod$vparameters
  # latent environmental covariates (loadings)
  Lam <- matrix(vparams[grep('^rr.*vm.*fa', names(vparams))], ncol = k)
  # specific variances
  Psi <- diag(vparams[grep(paste0('^',TE_fct,'.*vm'), names(vparams))])
  
  # % var explained:
  # mean of ratios
  meanVaPct <- mean(diag(Lam %*% t(Lam))/diag(Lam %*% t(Lam)+Psi)) * 100
  
  # for each trait x env combo
  TE <- names(vparams)[grep(paste0('^',TE_fct,'.*vm'), names(vparams))] |>
    sub(paste0('.*',TE_fct,'_'), '', x=_)
  
  TraitEnv_VaPct <- data.frame(TE_fct = TE,
                               VaPct = diag(Lam %*% t(Lam))/diag(Lam %*% t(Lam)+Psi) * 100)
  
  return(list(meanVaPct = meanVaPct, TraitEnv_VaPct = TraitEnv_VaPct))
}


# GEBVs ----
# Get genetic estimated breeding values (GEBVs)

gebvs_asreml <- function(mod,k,data,TE_fct) {
  # variance parameters
  vparams <- mod$vparameters
  # latent environmental covariates (loadings)
  Lam <- matrix(vparams[grep('^rr.*vm.*fa', names(vparams))], ncol = k)
  # specific variances
  Psi <- diag(vparams[grep(paste0('^',TE_fct,'.*vm'), names(vparams))])
  
  # coefficients
  coefs <- mod$coefficients$random
  # genotype scores (slopes)
  f <- coefs[grep('Comp', rownames(coefs))]
  Lamf <- c(matrix(f, ncol = k) %*% t(Lam)) # common GET effects
  # genetic regressions residuals
  delta <- coefs[grep(paste0('^',TE_fct,'.*vm'), rownames(coefs))] # specific GET effects
  # total GET effects
  Lamfdelta <- Lamf + delta
  
  # get gebv for each genotype, by trait-environment combination
  
  # vector of GET
  effects_all <- rownames(coefs)
  head(effects_all)
  effects_GET <- effects_all[grepl(paste0('^',TE_fct,'.*vm'), effects_all)]
  head(effects_GET)
  GET <- sub(paste0(TE_fct,'_'), '', effects_GET)
  GET <- sub(':vm\\(.*?\\)', '', GET)
  head(GET)
  
  split_GET <- do.call(rbind, strsplit(GET, '_', fixed = TRUE))
  
  blup <- data.frame(TE_fct=split_GET[,1],
                      G_fct=split_GET[,2],
                      blup=Lamfdelta)
  rownames(blup) <- GET
  return(blup)
}

gebvs_common_asreml <- function(mod,k,data,TE_fct) {
  # variance parameters
  vparams <- mod$vparameters
  # latent environmental covariates (loadings)
  Lam <- matrix(vparams[grep('^rr.*vm.*fa', names(vparams))], ncol = k)
  # specific variances
  Psi <- diag(vparams[grep(paste0('^',TE_fct,'.*vm'), names(vparams))])
  
  # coefficients
  coefs <- mod$coefficients$random
  # genotype scores (slopes)
  f <- coefs[grep('Comp', rownames(coefs))]
  Lamf <- c(matrix(f, ncol = k) %*% t(Lam)) # common GET effects

  # get gebv for each genotype, by trait-environment combination
  
  # vector of GET
  effects_all <- rownames(coefs)
  head(effects_all)
  effects_GET <- effects_all[grepl(paste0('^',TE_fct,'.*vm'), effects_all)]
  head(effects_GET)
  GET <- sub(paste0(TE_fct,'_'), '', effects_GET)
  GET <- sub(':vm\\(.*?\\)', '', GET)
  head(GET)
  
  split_GET <- do.call(rbind, strsplit(GET, '_', fixed = TRUE))
  
  blup <- data.frame(TE_fct=split_GET[,1],
                     G_fct=split_GET[,2],
                     blup=Lamf)
  rownames(blup) <- GET
  return(blup)
}

# Gen Corr ----
gcorr_asreml <- function(mod,k,data,TE_fct) {
  # variance parameters
  vparams <- mod$vparameters
  # latent environmental covariates (loadings)
  Lam <- matrix(vparams[grep('^rr.*vm.*fa', names(vparams))], ncol = k)
  # Rotate
  # Perform Singular Value Decomposition (SVD) for rotation
  svd <- svd(Lam) # Perform SVD on the loadings matrix
  # Compute rotated estimated loadings
  LamStar <- Lam %*% svd$v
  # specific variances
  Psi <- diag(vparams[grep(paste0('^',TE_fct,'.*vm'), names(vparams))])
  # Variance covariance matrix
  Gvar <- LamStar%*%t(LamStar)+Psi
  rownames(Gvar) <- levels(data$TraitEnv)
  colnames(Gvar) <- levels(data$TraitEnv)
  # Genetic correlation matrix
  Cmat <- cov2cor(Gvar)
  return(list(Gvar = Gvar, gcorr = Cmat))
}

# Accuracy ----