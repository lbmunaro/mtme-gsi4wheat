# Genotypic data ----
# This script subsets the phenotypic and genotypic data, performs quality control, 
# and creates a genomic relationship matrix (G matrix) for further analysis.

# Clean workspace
rm(list = objects())  # Removes all objects from the environment.

# Packages ----
library(tidyverse) # R packages for data science.
library(gaston) # Genetic Data Handling.
library(ASRgenomics) # R package for Genomic Selection Analysis in R.

# Load data ----
## Phenotypic data ----
ILYT_Pheno <- readRDS('Data/ILYT_Pheno.rds') # Load the phenotypic data.

# Vector of genotypes with phenotypic data
genotypes <- ILYT_Pheno |>
  filter(!is.na(Pheno)) |>
  ungroup()  |>
  droplevels() |>
  select(Gen) |>
  reframe(Gen=unique(Gen)) |>
  as_vector() |>
  glimpse()
genotypes <- unique(genotypes)

## Marker data ----
Markers <- read.vcf('Data/imputed_24IL_all_regions_filtV3_filtered.vcf.gz') # Load VCF file containing genotype markers.

# Subset marker data to match the genotypes with phenotypic data
Markers <- select.inds(Markers, id %in% genotypes)

genoM <- as.matrix(Markers) # Convert marker data to a matrix.
genoM[1:5,1:5]
dim(genoM)
# Filter the molecular matrix based on quality control thresholds
genoM_filter <- qc.filtering(
  M = genoM,  # Input molecular matrix.
  maf = 0.05,  # Minimum allele frequency threshold.
  marker.callrate = 0.2,  # Minimum call rate per marker.
  ind.callrate = 0.2,  # Minimum call rate per individual.
  impute = F,  # No imputation of missing data
  plots = T  # Generate plots for quality control.
)
dim(genoM_filter$M.clean)
genoM_filter$M.clean[1:5,1:5]

# Calculate the G matrix (genomic relationship matrix) using the VanRaden method
Gg <- G.matrix(
  M = genoM_filter$M.clean,  # Use cleaned marker data.
  method = 'VanRaden',  # VanRaden method for calculating G matrix.
  na.string = NA  # Handle missing data as NA.
)$G
Gg[1:5, 1:5]

check_G <- kinship.diagnostics(K = Gg, duplicate.thr = 0.95)

check_G$list.diagonal
check_G$list.duplicate

check_G$plot.diag
check_G$plot.offdiag

# Find the row and column indices of off-diagonal values > 2
off_diag_indices <- which(Gg > 1.75 & row(Gg) != col(Gg), arr.ind = TRUE)

# Extract the row and column names for the identified indices
off_diag_combinations <- data.frame(
  RowName = rownames(Gg)[off_diag_indices[, 1]],
  ColName = colnames(Gg)[off_diag_indices[, 2]],
  Value = Gg[off_diag_indices]
)
# Display the combinations with values > 2
print(off_diag_combinations)

to_remove <- c('IL2022-10404','IL2022-10405','IL2022-10413')
Gg <- Gg[setdiff(rownames(Gg), to_remove), setdiff(colnames(Gg), to_remove)]

check_G <- kinship.diagnostics(K = Gg, duplicate.thr = 0.95)

#check_G$list.duplicate

#check_G$plot.diag
#check_G$plot.offdiag

Gg.blend <- G.tuneup(G = Gg, blend = T, pblend=0.05)$Gb
Gg.blend[1:5,1:5]

check_Gg.blend <- kinship.diagnostics(K = Gg.blend, duplicate.thr = 0.95)

Ginv <- G.inverse(G = Gg.blend)$Ginv  # Compute the inverse of the bent G matrix.

Ginv.sparse <- G.inverse(G = Gg.blend, sparseform = T)$Ginv.sparse  # Compute the inverse of the bent G matrix.

# Add Gkeep and Gdrop factors for genotypes with and without marker data
ILYT_Pheno <- ILYT_Pheno |>
  mutate(
    Gkeep = factor(if_else(Gen %in% colnames(Ginv), Gen, NA)),  # Genotypes with marker data.
    Gdrop = factor(if_else(Gen %in% colnames(Ginv), NA, Gen))  # Genotypes without marker data.
  ) |>
  glimpse()

saveRDS(ILYT_Pheno, file = 'Data/ILYT_Pheno.rds')
saveRDS(Ginv, file = 'Data/Ginv.rds')
saveRDS(Ginv.sparse, file = 'Data/Ginv_sparse.rds')

#save(Ginv.sparse, ILYT_Pheno, file='Data/ILYT_Pheno-Gmatrix.RData')

# # heatmaps
# heatmapGg <- kinship.heatmap(K = Gg, dendrogram = TRUE,
#                              row.label = FALSE, col.label = FALSE)
# heatmapGg_blend <- kinship.heatmap(K = Gg.blend, dendrogram = TRUE, 
#                                    row.label = FALSE, col.label = FALSE)
# heatmapGinv <- kinship.heatmap(K = Ginv, dendrogram = TRUE, 
#                                row.label = FALSE, col.label = FALSE)

# End ----