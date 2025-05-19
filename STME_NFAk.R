# Single-Trait Multi-TraitEnvironment models ----
# This script fits Single-Trait Multi-TraitEnvironment models

# Clean workspace
rm(list = objects())  # Removes all objects from the TraitEnvironment.

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

# Fit models ----

## GY ----
### NFA1----
# Run model
GY_STME_NFA1.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,1):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) + 
    diag(TraitEnv):Block, 
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'GY') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '16gb'
)
# Print model info
print('GY_STME_NFA1')
print(paste('convergence =',GY_STME_NFA1.asr$converge))
GY_STME_NFA1.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
GY_STME_NFA1.asr <- update_asreml(GY_STME_NFA1.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA2----
# Run model
GY_STME_NFA2.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,2):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'GY') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '16gb'
)
# Print model info
print('GY_STME_NFA2')
print(paste('convergence =',GY_STME_NFA2.asr$converge))
GY_STME_NFA2.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
GY_STME_NFA2.asr <- update_asreml(GY_STME_NFA2.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA3----
# Run model
GY_STME_NFA3.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,3):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'GY') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '20gb'
)
# Print model info
print('GY_STME_NFA3')
print(paste('convergence =',GY_STME_NFA3.asr$converge))
GY_STME_NFA3.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
GY_STME_NFA3.asr <- update_asreml(GY_STME_NFA3.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA4----
# Run model
GY_STME_NFA4.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,4):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'GY') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '20gb'
)
# Print model info
print('GY_STME_NFA4')
print(paste('convergence =',GY_STME_NFA4.asr$converge))
GY_STME_NFA4.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
GY_STME_NFA4.asr <- update_asreml(GY_STME_NFA4.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA5----
# Run model
GY_STME_NFA5.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,5):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'GY') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '20gb'
)
# Print model info
print('GY_STME_NFA5')
print(paste('convergence =',GY_STME_NFA5.asr$converge))
GY_STME_NFA5.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
GY_STME_NFA5.asr <- update_asreml(GY_STME_NFA5.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA6----
# Run model
GY_STME_NFA6.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,6):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'GY') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '20gb'
)
# Print model info
print('GY_STME_NFA6')
print(paste('convergence =',GY_STME_NFA6.asr$converge))
GY_STME_NFA6.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
GY_STME_NFA6.asr <- update_asreml(GY_STME_NFA6.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

## TW ----
### NFA1----
# Run model
TW_STME_NFA1.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,1):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'TW') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '16gb'
)
# Print model info
print('TW_STME_NFA1')
print(paste('convergence =',TW_STME_NFA1.asr$converge))
TW_STME_NFA1.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
TW_STME_NFA1.asr <- update_asreml(TW_STME_NFA1.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA2----
# Run model
TW_STME_NFA2.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,2):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'TW') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '16gb'
)
# Print model info
print('TW_STME_NFA2')
print(paste('convergence =',TW_STME_NFA2.asr$converge))
TW_STME_NFA2.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
TW_STME_NFA2.asr <- update_asreml(TW_STME_NFA2.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA3----
# Run model
TW_STME_NFA3.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,3):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'TW') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '20gb'
)
# Print model info
print('TW_STME_NFA3')
print(paste('convergence =',TW_STME_NFA3.asr$converge))
TW_STME_NFA3.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
TW_STME_NFA3.asr <- update_asreml(TW_STME_NFA3.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA4----
# Run model
TW_STME_NFA4.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,4):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'TW') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '20gb'
)
# Print model info
print('TW_STME_NFA4')
print(paste('convergence =',TW_STME_NFA4.asr$converge))
TW_STME_NFA4.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
TW_STME_NFA4.asr <- update_asreml(TW_STME_NFA4.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA5----
# Run model
TW_STME_NFA5.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,5):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'TW') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '20gb'
)
# Print model info
print('TW_STME_NFA5')
print(paste('convergence =',TW_STME_NFA5.asr$converge))
TW_STME_NFA5.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
TW_STME_NFA5.asr <- update_asreml(TW_STME_NFA5.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA6----
# Run model
TW_STME_NFA6.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,6):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'TW') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '20gb'
)
# Print model info
print('TW_STME_NFA6')
print(paste('convergence =',TW_STME_NFA6.asr$converge))
TW_STME_NFA6.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
TW_STME_NFA6.asr <- update_asreml(TW_STME_NFA6.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

## HD ----
### NFA1----
# Run model
HD_STME_NFA1.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,1):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'HD') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '16gb'
)
# Print model info
print('HD_STME_NFA1')
print(paste('convergence =',HD_STME_NFA1.asr$converge))
HD_STME_NFA1.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
HD_STME_NFA1.asr <- update_asreml(HD_STME_NFA1.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA2----
# Run model
HD_STME_NFA2.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,2):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'HD') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '16gb'
)
# Print model info
print('HD_STME_NFA2')
print(paste('convergence =',HD_STME_NFA2.asr$converge))
HD_STME_NFA2.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
HD_STME_NFA2.asr <- update_asreml(HD_STME_NFA2.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA3----
# Run model
HD_STME_NFA3.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,3):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'HD') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '20gb'
)
# Print model info
print('HD_STME_NFA3')
print(paste('convergence =',HD_STME_NFA3.asr$converge))
HD_STME_NFA3.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
HD_STME_NFA3.asr <- update_asreml(HD_STME_NFA3.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

## HT ----
### NFA1----
# Run model
HT_STME_NFA1.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,1):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'HT') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '16gb'
)
# Print model info
print('HT_STME_NFA1')
print(paste('convergence =',HT_STME_NFA1.asr$converge))
HT_STME_NFA1.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
HT_STME_NFA1.asr <- update_asreml(HT_STME_NFA1.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA2----
# Run model
HT_STME_NFA2.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,2):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'HT') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '16gb'
)
# Print model info
print('HT_STME_NFA2')
print(paste('convergence =',HT_STME_NFA2.asr$converge))
HT_STME_NFA2.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
HT_STME_NFA2.asr <- update_asreml(HT_STME_NFA2.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA3----
# Run model
HT_STME_NFA3.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,3):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'HT') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '20gb'
)
# Print model info
print('HT_STME_NFA3')
print(paste('convergence =',HT_STME_NFA3.asr$converge))
HT_STME_NFA3.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
HT_STME_NFA3.asr <- update_asreml(HT_STME_NFA3.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA4----
# Run model
HT_STME_NFA4.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,4):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'HT') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '20gb'
)
# Print model info
print('HT_STME_NFA4')
print(paste('convergence =',HT_STME_NFA4.asr$converge))
HT_STME_NFA4.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
HT_STME_NFA4.asr <- update_asreml(HT_STME_NFA4.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

### NFA5----
# Run model
HT_STME_NFA5.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,5):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'HT') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '20gb'
)
# Print model info
print('HT_STME_NFA5')
print(paste('convergence =',HT_STME_NFA5.asr$converge))
HT_STME_NFA5.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
HT_STME_NFA5.asr <- update_asreml(HT_STME_NFA5.asr, 
                                     max_updates = 10,
                                     save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

## MAT ----
### NFA1----
# Run model
MAT_STME_NFA1.asr <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~ rr(TraitEnv,1):vm(Gkeep, Ginv) + diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
  sparse = ~ TraitEnv:Gdrop,
  data = ILYT_Pheno |> filter(Trait == 'MAT') |> droplevels(),
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '16gb'
)
# Print model info
print('MAT_STME_NFA1')
print(paste('convergence =',MAT_STME_NFA1.asr$converge))
MAT_STME_NFA1.asr$trace |> as.data.frame() |> rownames_to_column('Iteration') |>
  filter(Iteration=='LogLik') |> print()

# Update model
MAT_STME_NFA1.asr <- update_asreml(MAT_STME_NFA1.asr, 
                                      max_updates = 10,
                                      save_path = 'Data/STME_NFAk.RData')

save.image('Data/STME_NFAk.RData')

# End ----