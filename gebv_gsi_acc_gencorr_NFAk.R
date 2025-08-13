# Model selection

rm(list=ls())
library(asreml)
library(tidyverse)
source('Functions_MTME.R')

library(ggpubr)
library(ggforce)

# Load model ----
load('Data/MTME_NFA5.RData')

mod <- MTME_NFA5.asr

# Accuracy ----
effects_all <- rownames(mod$coefficients$random)
effects_GET <- effects_all[grepl(paste0('^TraitEnv.*vm'), effects_all)]
GET <- sub(paste0('TraitEnv_'), '', effects_GET)
GET <- sub(':vm\\(.*?\\)', '', GET)
split_GET <- do.call(rbind, strsplit(GET, '_', fixed = TRUE))

ACC_MTME <- data.frame(TraitEnv=split_GET[,1],
                       Gkeep=split_GET[,2],
                       ACC_cve=readRDS('Data/ACC_cve_diag_NFAk.rds'),
                       ACC_tve=readRDS('Data/ACC_tve_diag_NFAk.rds')) |>
  separate(TraitEnv, into = c("Trait","Env"), sep = "-", extra = "merge", remove = F) |>
  mutate(Trait = factor(Trait,
                        levels = c('GSI','GY','TW','HD','HT','MAT'))) |>
  glimpse()

## mean accuracy per trait
ACC_MTME |>
  group_by(Trait) |>
  summarise(ACC_cve=mean(ACC_cve),
            ACC_tve=mean(ACC_tve))

ggplot(ACC_MTME, aes(x=Trait, y=ACC_tve)) +
  geom_boxplot( width = 0.5, alpha=0.5, size = 0.5) + 
  scale_y_continuous(name = expression(italic(r[vps])),
                     limits = c(0:1),
                     breaks = seq(0, 1, by = .2)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    legend.title = element_text(size = 12, family = 'Times New Roman'),
    legend.text = element_text(size = 10, family = 'Times New Roman'),
    axis.text = element_text(size = 10, family = 'Times New Roman'),
    axis.title = element_text(size = 10, family = 'Times New Roman')
  )
ggsave('Figures/Figure6.tif', width = 3.5, height = 3.5, units = 'in', dpi = 300)

# Get gebvs ----
MTME_gebvs <- gebvs_asreml(mod, k, ILYT_Pheno, 'TraitEnv') |>
  rename(TraitEnv=TE_fct,
         Gkeep=G_fct) |>
  left_join(
    ILYT_Pheno |>
      group_by(TraitEnv) |>
      summarise(Pheno_mean=mean(Pheno_mean),
                Pheno_sd=mean(Pheno_sd),
                Trait=unique(Trait),
                Env=unique(Env))
  ) |>
  mutate(MTME_gebv=blup*Pheno_sd) |>
  select(Gkeep, Trait, Env, MTME_gebv) |>
  glimpse()

 # MTME vs. STME

## gebvs ----
STME_gebvs <- readRDS('Data/STME_NFAk_gebvs.rds')

gebvs <- MTME_gebvs |>
  left_join(STME_gebvs)

gebvs_long <- gebvs |>
  rename(
    MTME = MTME_gebv,
    STME = STME_gebv
  ) |>
  select(Gkeep, Trait, MTME, STME)

## gsi ----
gsi <- gebvs |>
  group_by(Gkeep, Trait) |>
  summarise(MTME_gebv=mean(MTME_gebv),
            STME_gebv=mean(STME_gebv)) |>
  ungroup() |>
  pivot_wider(names_from = Trait, values_from = c(MTME_gebv, STME_gebv)) |>
  mutate(MTME_GSI=((213.75/1000*MTME_gebv_GY)+(6.5*0.41*MTME_gebv_TW)+(-13.32*MTME_gebv_HD)),
         STME_GSI=((213.75/1000*STME_gebv_GY)+(6.5*0.41*STME_gebv_TW)+(-13.32*STME_gebv_HD))
  ) |>
  arrange(desc(MTME_GSI)) |>
  glimpse()

gsi_long <- gsi |> 
  rename(
    MTME = MTME_GSI,
    STME = STME_GSI
  ) |> 
  select(Gkeep, MTME, STME) |> 
  mutate(Trait = 'GSI')

gebv_gsi <- bind_rows(gebvs_long, gsi_long) |>
  mutate(Trait = factor(Trait,
                        levels = c('GSI','GY','TW','HD','HT','MAT'))) |>
  glimpse()

## Plot ----
ggplot(gebv_gsi, aes(x = MTME, y = STME)) +
  geom_hex(alpha = 1) +
  geom_smooth(method = 'lm', se = FALSE, linewidth = 0.25, color = 'gray', alpha = 0.5) +
  stat_cor(
    method = "pearson",
    label.x.npc = "left",
    label.y.npc = "top",
    aes(label = after_stat(r.label))
  ) +
  scale_fill_gradient(low = '#13294B', high = '#FF5F05') +
  facet_wrap(~Trait, scales = 'free', ncol = 3) +
  xlab('MTME') + ylab('STME') +
  theme_bw() +
  theme(
    text = element_text(family = 'Times New Roman'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'grey90'),
    strip.text = element_text(size = 12, face = 'bold'),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave('Figures/Figure7.tif',
       width = 7, height = 4, units = 'in', dpi = 300)

# RespSel ----

# Selection candidates
g <- ILYT_Pheno |>
  filter(str_detect(Env,'24-')) |>
  filter(!is.na(Gkeep)) |>
  droplevels() |>
  group_by(Gkeep) |>
  reframe() |>
  pull() |> as.vector()

# Selection intensity 0.1
s <- round(0.1*length(g),0)
s

# Averaged across Trait and Gkeep
gebv_gsi_summary <- gebv_gsi |>
  filter(Gkeep %in% g) |>
  droplevels() |>
  group_by(Trait, Gkeep) |>
  summarise(MTME=mean(MTME),
            STME=mean(STME)) |>
  ungroup()

## Selected individuals ----
### MTME-GSI
rank_MT_GSI <- gebv_gsi_summary |>
  filter(Trait=='GSI') |>
  arrange(desc(MTME)) |>
  mutate(rank=rank(-MTME)) |>
  filter(rank<=s) |>
  pull(Gkeep)
### MTME-GY
rank_MT_GY <- gebv_gsi_summary |>
  filter(Trait=='GY') |>
  arrange(desc(MTME)) |>
  mutate(rank=rank(-MTME)) |>
  filter(rank<=s) |>
  pull(Gkeep)
### STME-GSI
rank_ST_GSI <- gebv_gsi_summary |>
  filter(Trait=='GSI') |>
  arrange(desc(STME)) |>
  mutate(rank=rank(-STME)) |>
  filter(rank<=s) |>
  pull(Gkeep) 
### STME-GY
rank_ST_GY <- gebv_gsi_summary |>
  filter(Trait=='GY') |>
  arrange(desc(STME)) |>
  mutate(rank=rank(-STME)) |>
  filter(rank<=s) |>
  pull(Gkeep) 

## Sel vs All
none <- gebv_gsi_summary |>
  mutate(Sel='NA') 

MT_GSI <- gebv_gsi_summary |>
  filter(Gkeep %in% rank_MT_GSI) |>
  mutate(Sel='GSI-MTME') 

MT_GY <- gebv_gsi_summary |>
  filter(Gkeep %in% rank_MT_GY) |>
  mutate(Sel='GY-MTME') 

ST_GSI <- gebv_gsi_summary |>
  filter(Gkeep %in% rank_ST_GSI) |>
  mutate(Sel='GSI-STME') 

ST_GY <- gebv_gsi_summary |>
  filter(Gkeep %in% rank_ST_GY) |>
  mutate(Sel='GY-STME') 

MTcomb <- bind_rows(none, MT_GSI, MT_GY) |>
  mutate(Sel = factor(Sel,
                      levels = c('NA', 'GY-MTME', 'GSI-MTME'))) |>
  glimpse()

MTSTcomb <- bind_rows(none,MT_GSI, MT_GY, ST_GSI, ST_GY) |>
  mutate(Sel = factor(Sel,
                      levels = c('NA', 'GY-STME',  'GY-MTME', 'GSI-STME', 'GSI-MTME'))) |>
  glimpse()

means_MTcomb <- MTcomb %>%
  group_by(Trait, Sel) %>%
  summarize(mean_MTME = mean(MTME), .groups = "drop") |>
  glimpse()

means_MTSTcomb <- MTSTcomb %>%
  group_by(Trait, Sel) %>%
  summarize(mean_MTME = mean(MTME), .groups = "drop") |>
  glimpse()

## Plot ----
### NA vs MTME ----
ggplot(MTcomb, aes(x = MTME, fill = Sel)) +
  geom_histogram(
    aes(color = Sel),
    position = "identity",
    bins     = 20,
    alpha    = 0.4,
    size     = 0
  ) +
  geom_vline(
    data = means_MTcomb,
    aes(xintercept = mean_MTME, color = Sel),
    linetype = "dashed",
    size = 0.6
  ) +
  scale_fill_manual(name   = "Selection",
                    labels = c('NA', 'GY', 'GSI'),
                    values = c("gray50", "#13294B", "#FF5F05")) +
  scale_color_manual(name   = "Selection",
                     labels = c('NA', 'GY', 'GSI'),
                     values = c("gray50", "#13294B", "#FF5F05")) +
  facet_wrap(
    ~Trait, scales = 'free_x', 
    labeller = as_labeller(
      c(
        GSI = "GSI~(USD~ha^{-1})",
        GY  = "GY~(kg~ha^{-1})",
        HD  = "HD~(day)",
        MAT = "MAT~(day)",
        HT  = "HT~(cm)",
        TW  = "TW~(g~L^{-1})"
      ), 
      default = label_parsed
    )
  ) +
  xlab(NULL) + ylab('No. Genotypes') +
  theme_bw() +
  theme(
    text = element_text(family = 'Times New Roman'),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = 'grey90'),
    strip.text = element_text(size = 12, face = 'bold', family = 'Times New Roman'),  # Ensure facet labels use Times New Roman
    legend.title = element_text(size = 12, family = 'Times New Roman'),
    legend.text = element_text(size = 10, family = 'Times New Roman'),
    legend.position = 'bottom',
    axis.text = element_text(size = 10, family = 'Times New Roman'),
    axis.title = element_text(size = 12, family = 'Times New Roman'),
    axis.text.x = element_text(angle = 0, hjust = 0.5, family = 'Times New Roman')  # Keep x-axis labels centered
  )
ggsave('Figures/Figure8.tiff', width = 7, height = 5, units = 'in', dpi = 300)

# MTSTcomb2 <- MTSTcomb |>
#   separate(Sel,into = c('Crit', 'Mod'), remove = F) |>
#   mutate(Mod = factor(Mod,
#                       levels = c('NA', 'STME', 'MTME'))) |>
#   mutate(Crit = factor(Crit,
#                       levels = c('NA', 'GY', 'GSI'))) |>
#   glimpse()
# 
# ### MTME vs. STME ----
# ggplot(MTSTcomb2, aes(x = MTME, y=Crit, fill = Mod, color=Mod)) +
#   ggridges::stat_density_ridges(
#     scale=1,
#     quantile_lines = TRUE, 
#                                 quantiles = 0.5, # median line
#                                 alpha=0.5) + 
#   facet_wrap(
#     ~Trait, scales = 'free_x', 
#     labeller = as_labeller(
#       c(
#         GSI = "GSI~(USD~ha^{-1})",
#         GY  = "GY~(kg~ha^{-1})",
#         HD  = "HD~(day)",
#         MAT = "MAT~(day)",
#         HT  = "HT~(cm)",
#         TW  = "TW~(g~L^{-1})"
#       ), 
#       default = label_parsed
#     )
#   ) +
#   scale_fill_manual(name   = "Model",
#                     values = c( "#13294B", "#FF5F05",NULL)) +
#   scale_color_manual(name   = "Model",
#                      values = c( "#13294B", "#FF5F05",NULL)) +
#   xlab(NULL) + ylab('Selection') +
#   theme_bw() +
#   theme(
#     text = element_text(family = 'Times New Roman'),
#     panel.grid = element_blank(),
#     strip.background = element_rect(fill = 'grey90'),
#     strip.text = element_text(size = 12, face = 'bold', family = 'Times New Roman'),  # Ensure facet labels use Times New Roman
#     legend.title = element_text(size = 12, family = 'Times New Roman'),
#     legend.text = element_text(size = 10, family = 'Times New Roman'),
#     legend.position = 'bottom',
#     axis.text = element_text(size = 10, family = 'Times New Roman'),
#     axis.title = element_text(size = 12, family = 'Times New Roman'),
#     axis.text.x = element_text(angle = 0, hjust = 0.5, family = 'Times New Roman')  # Keep x-axis labels centered
#   )
# ggsave('Figures/Figure9.tiff', width = 7, height = 5, units = 'in', dpi = 300)

RespSel <- gebv_gsi_summary |>
  group_by(Trait) |>
  mutate(MTmeanAll = mean(MTME),
         STmeanAll = mean(STME)) |>

  summarise(
    GSI_MTME = mean(MTME[Gkeep %in% rank_MT_GSI]) - unique(MTmeanAll),
    GSI_STME = mean(MTME[Gkeep %in% rank_ST_GSI]) - unique(MTmeanAll),
    GY_MTME = mean(MTME[Gkeep %in% rank_MT_GY]) - unique(MTmeanAll),
    GY_STME = mean(MTME[Gkeep %in% rank_ST_GY]) - unique(MTmeanAll),
  ) |>
  pivot_longer(cols = GSI_MTME:GY_STME, names_to = 'ModCrit', values_to = 'RespSel') |>
  separate(ModCrit, into = c('Crit','Model'), sep = '_', remove = F) |>
  mutate(RespSel=round(RespSel,1)) |>
  glimpse()
  
RespSel |>
  select(Trait,ModCrit,RespSel) |>
  pivot_wider(names_from = ModCrit, values_from = RespSel)

# GenCorr ----

## Data ----
gcorr <- gcorr_asreml(mod = MTME_NFA5.asr, k=5, data = ILYT_Pheno, TE_fct = 'TraitEnv')$gcorr
gcorr

library(corrplot)

## TraitEnv ----
# Define color palette
my_colors <- colorRampPalette(c('#13294b', 'white', '#FF5F0F'))(20)

# Create correlation plot with specified colors
png('Figures/Figure10a.tiff', width = 7, height = 7, units = 'in', res = 300)

# Set font to Times New Roman
par(family = "Times")

corrplot(gcorr, 
         method = 'pie', 
         type = 'lower', 
         col = my_colors,
         diag = T,
         tl.col = 'black', 
         tl.cex = 0.75,
         tl.srt = 90)
dev.off()

## Trait ----

corr_df <- as.data.frame(as.table(as.matrix(gcorr))) |>
  rename(Row = Var1, Col = Var2, Correlation = Freq) |>
  mutate(
    Row = as.character(Row),
    Col = as.character(Col),
    Trait_row = sapply(strsplit(Row, "-"), `[`, 1),
    Trait_col = sapply(strsplit(Col, "-"), `[`, 1),
    Year_row = sapply(strsplit(Row, "-"), `[`, 2),
    Year_col = sapply(strsplit(Col, "-"), `[`, 2),
    Loc_row = sapply(strsplit(Row, "-"), `[`, 3),
    Loc_col = sapply(strsplit(Col, "-"), `[`, 3),
    Env_row = sapply(strsplit(Row, "-"), function(x) paste(x[-1], collapse = "-")),
    Env_col = sapply(strsplit(Col, "-"), function(x) paste(x[-1], collapse = "-"))
  ) |>
  glimpse()


# Step 3: Create Trait × Trait average correlation matrix
trait_corr_mat <- corr_df |>
  group_by(Trait_row, Trait_col) |>
  summarise(mean_corr = mean(Correlation, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = Trait_col, values_from = mean_corr) |>
  column_to_rownames("Trait_row") |>
  as.matrix()
trait_corr_mat


png('Figures/Figure10b.tiff', width = 2, height = 2, units = 'in', res = 300)

# Set font to Times New Roman
par(family = "Times")

corrplot(trait_corr_mat, 
         method = 'pie', 
         type = 'upper', 
         col = my_colors,
         number.cex = .75, 
         number.digits = 2,
         diag = T,
         tl.col = 'black', 
         tl.cex = 0.75,
         tl.srt = 90,
         cl.pos = 'n',
         addCoef.col = "black")
dev.off()

### Grain yield ----
GY_corr_df <- corr_df |>
  filter(Trait_row=='GY') |>
  filter(Trait_col=='GY') |>
  glimpse()

loc_corr_mat <- GY_corr_df |>
  group_by(Loc_row, Loc_col) |>
  summarise(mean_corr = mean(Correlation, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = Loc_col, values_from = mean_corr) |>
  column_to_rownames("Loc_row") |>
  as.matrix()
loc_corr_mat

png('Figures/Figure10c.tiff', width = 2/5*4, height = 2/5*4, units = 'in', res = 300)

# Set font to Times New Roman
par(family = "Times")
corrplot(loc_corr_mat, 
         method = 'pie', 
         type = 'upper', 
         col = my_colors,
         number.cex = .75, 
         number.digits = 2,
         diag = T,
         tl.col = 'black', 
         tl.cex = 0.75,
         tl.srt = 90,
         cl.pos = 'n',
         addCoef.col = "black")
dev.off()

