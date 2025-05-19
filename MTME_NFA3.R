# Multi-Trait Multi-Environment models ----
# This script fits Multi-Trait Multi-Environment models

# Clean workspace
rm(list = objects())  # Removes all objects from the environment.

# Packages ----
library(tidyverse) # R packages for data science.
library(asreml) # ASReml-R package.
source('Functions_MTME.R')  # Load functions

# Use for HPC only
setwd('~/mtme-gsi4wheat/')

# Load data ----
## Pheno & Ginv
ILYT_Pheno <- readRDS('Data/ILYT_Pheno.rds') # Load the phenotypic data.
Ginv <- readRDS('Data/Ginv.rds') # Load the relationship matrix

# Set k
k <- 3

# Fit NFA3 model ----
## Run model ----
MTME_NFA3.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,k):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) + 
    diag(TraitEnv):Block, 
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno,
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '72gb'
)

# Print model info
print(paste('convergence =', MTME_NFA3.asr$converge))
MTME_NFA3.asr$trace |>
  as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Save
save.image('Data/MTME_NFA3.RData')

# Update model ----
MTME_NFA3.asr <- update_asreml(MTME_NFA3.asr, 
                               max_updates = 20,
                               save_path = "Data/MTME_NFA3.RData")

save.image('Data/MTME_NFA3.RData')

# Save
save.image('Data/MTME_NFA3.RData')