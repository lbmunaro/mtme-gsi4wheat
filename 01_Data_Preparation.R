# Data wrangling ----
# This script organizes and checks phenotypic data for further analysis.

# Clean workspace
rm(list = objects())  # Removes all objects from the environment.

# Packages ----
library(tidyverse) # R packages for data science.
library(janitor) # Simple Tools for Examining and Cleaning Dirty Data.

# Phenotypic data ----
## Phenotypic data in the wide format ----
ILYT_Pheno_w <- 
  read.csv('Data/ILYT_22-24-2024.07.19.csv') |> # Load raw data in csv format.
  remove_empty(which = c('cols')) |> # Remove empty columns.
  clean_names() |> # Clean names.
  # Replace location name.
  mutate(
    study_name = str_replace(study_name,'YT_Stj_22','YT_Addie_22'),
    location_name = str_replace(location_name,'St. Jacob Township, IL','Addieville, IL')
  ) |>
  # Select relevant variables.
  dplyr::select(
    study_name, study_year, location_name, germplasm_name,
    block_number, col_number, row_number, plot_number,
    grain_yield_kg_ha_co_321_0001218, grain_test_weight_g_l_co_321_0001210,
    heading_time_julian_date_jd_co_321_0001233, plant_height_cm_co_321_0001301,
    maturity_time_spike_estimation_julian_date_jd_co_321_0501101
  ) |>
  # Simplify study name and convert it to a factor.
  mutate(study_name=as.factor(gsub('^YT_', '', study_name)),
         study_name=as.factor(gsub('^Addie_', 'Adv_', study_name)),
         study_name=str_replace(study_name, '(\\w+)_(\\d+)', '\\2-\\1'),
         ) |>
  # Rename columns for easier reference.
  rename(
    Env=study_name,
    Year=study_year,
    Loc=location_name,
    Gen=germplasm_name,
    Block=block_number,
    Col=col_number,
    Row=row_number,
    Plot=plot_number,
    GY=grain_yield_kg_ha_co_321_0001218,
    TW=grain_test_weight_g_l_co_321_0001210,
    HD=heading_time_julian_date_jd_co_321_0001233,
    HT=plant_height_cm_co_321_0001301,
    MAT=maturity_time_spike_estimation_julian_date_jd_co_321_0501101
  ) |>
  # Arrange data by Env, Col, and Row.
  arrange(Env,Col,Row) |>
  # Convert specific columns to factors.
  mutate_at(vars(Env:Plot),as.factor) |>
  # Convert trait data to numeric format.
  mutate_at(vars(GY:MAT), as.numeric) |>
  # Replace zero values with NA in trait columns.
  mutate_at(vars(GY:MAT), ~ifelse(.<=0,NA,.)) |>
  mutate_at(vars(GY:MAT), ~ifelse(Gen=='FILL',NA,.)) |>
  mutate(Gen = as.factor(ifelse(Gen=='FILL',NA,as.character(Gen)))) |>
  mutate(Gen = as.factor(ifelse(Plot=='958'&Env=='24-Urb',NA,as.character(Gen)))) |>
  filter(Loc!='Belleville, IL') |>
  droplevels() |>
  # Display the dataset.
  glimpse()

## Phenotypic data in the long format ----
# Reshape the wide format data into long format for easier analysis.
ILYT_Pheno <- ILYT_Pheno_w |>
  pivot_longer(
    cols = GY:MAT,  # Convert trait columns into a single column.
    names_to = 'Trait',  # Name of the new column for trait types.
    values_to = 'Pheno',  # Name of the column holding trait values.
    values_drop_na = FALSE  # Keep NA values in the data.
  ) |>
  # Convert the trait column to a factor.
  mutate(Trait=as.factor(Trait)) |>
  # Arrange data by trait, environment, column, and row.
  arrange(Trait, Env, as.numeric(Col), as.numeric(Row)) |>
  # Group data by environment and trait, and remove groups with all NA values.
  group_by(Env, Trait) |>
  filter(!is.na(mean(Pheno, na.rm = TRUE))) |>
  ungroup() |>
  # Assign NA to outliers
  # GY
  mutate(Pheno= ifelse(Env=='22-Adv' & Col=='3' & Row=='17', NA, Pheno)) |>
  #  mutate(Pheno= ifelse(Env=='Neo_24' & Col=='16' & Row=='24', NA, Pheno)) |>
  # TW
  mutate(Pheno= ifelse(Trait=='TW' & Env=='24-Adv' & Col=='2' & Row=='24', NA, Pheno)) |>
  mutate(Pheno= ifelse(Trait=='TW' & Env=='24-Stp' & Col=='12' & Row=='16', NA, Pheno)) |>
  mutate(Pheno= ifelse(Trait=='TW' & Env=='24-Neo' & Col=='8' & Row=='6', NA, Pheno)) |>
  # Standardize phenotypic values (Z-scores) for each trait.
  group_by(Trait) |>
  mutate(Pheno_SI = Pheno, # standard units
         Pheno = ifelse(Trait=='GY', Pheno/c(60 * 0.453592 * 2.47105), Pheno), # GY in bu/ac and TW in lbs/bu
         Pheno = ifelse(Trait=='TW', Pheno/1000 *2.2046 *35.2391, Pheno),
         Pheno_z = as.vector(scale(Pheno_SI, center = T)),
         Pheno_std = as.vector(scale(Pheno_SI, center = F)),
         Pheno_mean = mean(Pheno_SI, na.rm = TRUE),
         Pheno_sd = sd(Pheno_SI, na.rm = TRUE),
         IDEU = as.factor(paste(Env,Plot, sep='_'))
         ) |>
  ungroup() |>
  # Create a new factor combining trait and environment.
  mutate(TraitEnv=as.factor(paste(Trait,Env,sep = '-'))) |>
  # Select relevant columns.
  dplyr::select(Pheno, Pheno_SI, Pheno_z,Pheno_std, Pheno_mean, Pheno_sd,
                TraitEnv,Trait,Env,Year,Loc,Block,Col,Row,Plot,IDEU,Gen) |>
  # Display the dataset.
  glimpse()

# Save data ----
#save(ILYT_Pheno, file='Data/ILYT_Pheno.RData')
saveRDS(ILYT_Pheno, file = 'Data/ILYT_Pheno.rds')

# End ----