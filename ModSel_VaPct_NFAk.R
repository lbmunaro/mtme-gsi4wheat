# Model selection

rm(list=ls())
library(asreml)
library(tidyverse)
source('Functions_MTME.R')

# Load models ----

files <- sprintf("Data/MTME_NFA%d.RData", 1:9)
invisible(lapply(files, load, envir = .GlobalEnv))

# Models ----
models <- mget(paste0("MTME_NFA", 1:9, ".asr"))

# LRT ----
exec(lrt.asreml, !!!models)

# Summary ----
#nPar, LogLik, AIC
ModSummary <- models %>%
  enframe(name = "model", value = "fit") |>
  mutate(
    npar   = map_dbl(fit, ~ attr(summary(.x)$aic, "parameters")),
    AIC    = map_dbl(fit, ~ summary(.x)$aic),
    logLik = map_dbl(fit, "loglik")
  ) |>
  select(model, npar, AIC, logLik)
ModSummary


# VaPct ----
VaPct_list <- map2(
  .x = models,
  .y = seq_along(models),
  ~ VaPct(mod    = .x,
          k      = .y,
          data   = ILYT_Pheno,
          TE_fct = "TraitEnv")$TraitEnv_VaPct
)

names(VaPct_list) <- paste0("NFA", 1:9)

VaPct_df <- VaPct_list |> as.data.frame() |> pivot_longer(
  cols       = matches("^NFA\\d+\\."),
  names_to   = c("NFAk", ".value"),           # extract the k‐level and the measurement name
  names_pattern = "NFA(\\d+)\\.(.*)"          # “NFA” + (digits) + “.” + (TE_fct or VaPct)
) |>
  rename(TraitEnv = TE_fct) |>                # rename for clarity
  mutate(model = paste0("NFA", NFAk)) |>
  mutate(
    Trait = str_extract(TraitEnv, "^[^-]+"),        # text up to first hyphen
    Env   = str_remove(TraitEnv, "^[^-]+-")         # drop “trait-” prefix
  ) |>
  glimpse()

## VaPct mean ----
VaPct_df |> group_by(model,Trait) |> summarise(mVaPct = mean(VaPct)) |> 
  pivot_wider(names_from = Trait, values_from = mVaPct) |>
  left_join(VaPct_df |> group_by(model) |> summarise(mVaPct = mean(VaPct)))

map2(
  .x = models,
  .y = seq_along(models),
  ~ VaPct(mod    = .x,
          k      = .y,
          data   = ILYT_Pheno,
          TE_fct = "TraitEnv")$meanVaPct
)

## Plot ----

ggplot(VaPct_df, aes(x=Trait, y=VaPct, fill=model)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.5, alpha=0.5) +  # Dodge boxplots
  geom_jitter(aes(color=model),position = position_dodge(width = 0.8),
              alpha=0.75, shape=3, size=1) +  # Dodge boxplots
  
  scale_y_continuous(name = bquote(V[a] ~ 'explained (%)'),
                     breaks = seq(0, 100, by = 20)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    legend.title = element_text(size = 12, family = 'Times New Roman'),
    legend.text = element_text(size = 10, family = 'Times New Roman'),
    axis.text = element_text(size = 10, family = 'Times New Roman'),
    axis.title = element_text(size = 10, family = 'Times New Roman')
  )
#ggsave('Figures/Figure5.png', width = 7, height = 4, units = 'in', dpi = 300)
ggsave('Figures/Figure5.tif', width = 7, height = 4, units = 'in', dpi = 300)
