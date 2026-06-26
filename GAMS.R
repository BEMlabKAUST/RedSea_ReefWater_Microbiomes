library(mgcv)
library(tidyverse)
library(gratia)
library(patchwork)
library(performance)
library(marginaleffects)
library(cowplot)
library(ggeffects)
library(ragg)
library(broom)

# ----------------------------
# Read data
# ----------------------------
df <- read_csv("~/Documents/R/EffortsGloria/DataAnalysis/February_2026/RF/model_ave_relab_01032026.csv")

# ----------------------------
# Output path: Submission folder
# ----------------------------
path_output <- "~/Documents/R/EffortsGloria/DataAnalysis/April_2026/GAMs/Submission/"
dir.create(path_output, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Function to run and save GAM outputs
# ----------------------------
run_and_save_gam <- function(model_obj, response_name, file_prefix, path_output) {
  
  model_summary <- summary(model_obj)
  
  # 1. Partial effects plot
  p <- draw(model_obj, parametric = TRUE, rug = TRUE) +
    plot_annotation(title = response_name)
  
  ggsave(
    filename = file.path(path_output, paste0(file_prefix, "_partial_effects.png")),
    plot = p,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(path_output, paste0(file_prefix, "_partial_effects.pdf")),
    plot = p,
    width = 10,
    height = 8
  )
  
  # 2. Parametric terms
  parametric_tbl <- as.data.frame(model_summary$p.table) %>%
    rownames_to_column("term")
  
  write_csv(
    parametric_tbl,
    file.path(path_output, paste0(file_prefix, "_parametric_terms.csv"))
  )
  
  # 3. Smooth terms
  smooth_tbl <- as.data.frame(model_summary$s.table) %>%
    rownames_to_column("term")
  
  write_csv(
    smooth_tbl,
    file.path(path_output, paste0(file_prefix, "_smooth_terms.csv"))
  )
  
  # 4. Collapse smooth p-values into one column
  smooth_pvalues <- smooth_tbl %>%
    select(term, `p-value`) %>%
    mutate(
      p_value_formatted = case_when(
        `p-value` < 0.001 ~ "<0.001",
        TRUE ~ as.character(signif(`p-value`, 3))
      ),
      term_p = paste0(term, " p=", p_value_formatted)
    )
  
  # 5. Model stats
  model_stats <- tibble(
    response = response_name,
    formula = paste(deparse(formula(model_obj)), collapse = " "),
    family = family(model_obj)$family,
    link = family(model_obj)$link,
    n = nobs(model_obj),
    r_sq_adj = model_summary$r.sq,
    deviance_explained_percent = model_summary$dev.expl * 100,
    scale_est = model_summary$scale,
    reml_score = model_obj$gcv.ubre,
    smooth_p_values = paste(smooth_pvalues$term_p, collapse = "; ")
  )
  
  write_csv(
    model_stats,
    file.path(path_output, paste0(file_prefix, "_model_stats.csv"))
  )
  
  # 6. K-check
  kcheck_tbl <- as.data.frame(k.check(model_obj)) %>%
    rownames_to_column("smooth_term")
  
  write_csv(
    kcheck_tbl,
    file.path(path_output, paste0(file_prefix, "_kcheck.csv"))
  )
  
  # 7. Collinearity
  collinearity_tbl <- check_collinearity(model_obj) %>%
    as.data.frame()
  
  write_csv(
    collinearity_tbl,
    file.path(path_output, paste0(file_prefix, "_collinearity.csv"))
  )
  
  # 8. Concurvity
  conc <- concurvity(model_obj, full = FALSE)
  
  conc_worst <- as.data.frame(conc$worst) %>%
    rownames_to_column("term")
  
  conc_observed <- as.data.frame(conc$observed) %>%
    rownames_to_column("term")
  
  conc_estimate <- as.data.frame(conc$estimate) %>%
    rownames_to_column("term")
  
  write_csv(
    conc_worst,
    file.path(path_output, paste0(file_prefix, "_concurvity_worst.csv"))
  )
  
  write_csv(
    conc_observed,
    file.path(path_output, paste0(file_prefix, "_concurvity_observed.csv"))
  )
  
  write_csv(
    conc_estimate,
    file.path(path_output, paste0(file_prefix, "_concurvity_estimate.csv"))
  )
  
  # 9. GAM check plot
  png(
    filename = file.path(path_output, paste0(file_prefix, "_gamcheck.png")),
    width = 2400,
    height = 2000,
    res = 300
  )
  
  gam.check(model_obj)
  dev.off()
  
  return(list(
    model = model_obj,
    plot = p,
    parametric_terms = parametric_tbl,
    smooth_terms = smooth_tbl,
    stats = model_stats,
    kcheck = kcheck_tbl,
    collinearity = collinearity_tbl
  ))
}

# ----------------------------
# 1. Shannon index taxonomy
# ----------------------------
gam_shannon <- gam(
  WAT_Shannon_Index_ln ~ 
    s(log2(Chla), k = 3, bs = "tp") +
    s(SST, k = 5, bs = "tp") +
    s(Coral_richness, k = 3, bs = "tp"),
  data = df,
  family = gaussian(),
  method = "REML",
  select = TRUE
)

# ----------------------------
# 2. Function richness
# ----------------------------
gam_function <- gam(
  WAT_Function_richness ~ 
    s(HC_bleached, k = 3, bs = "tp") +
    s(LOAP, k = 3, bs = "tp") +
    s(SST, k = 6, bs = "tp"),
  data = df,
  family = gaussian(),
  method = "REML",
  select = TRUE
)

# ----------------------------
# 3. Rhodobacterales
# ----------------------------
gam_rhodo <- gam(
  WAT_Rhodobacterales ~ 
    s(log2(DTS_km), k = 6, bs = "tp") +
    s(Phosphate, k = 6, bs = "tp") +
    s(Mn, k = 3, bs = "tp"),
  data = df,
  family = gaussian(),
  method = "REML",
  select = TRUE
)

# ----------------------------
# 4. Prochlorococcus marinus
# ----------------------------
gam_pro <- gam(
  WAT_Prochlorococcus_marinus ~ 
    s(log2(DTS_km), k = 3, bs = "tp") +
    s(SST, k = 6, bs = "tp") +
    s(Ni, k = 6, bs = "tp"),
  data = df,
  family = betar(link = "logit"),
  method = "REML",
  select = TRUE
)

# ----------------------------
# 5. Bacillota
# ----------------------------
gam_bacillota <- gam(
  WAT_Bacillota ~ 
    s(LOAP, k = 9, bs = "tp") +
    s(HC_nonbleached, k = 3, bs = "tp") +
    s(Mn, k = 4, bs = "tp"),
  data = df,
  family = gaussian(),
  method = "REML",
  select = TRUE
)

# ----------------------------
# 6. Vibrio spp.
# ----------------------------
gam_vibrio <- gam(
  WAT_Vibrio_spp ~ 
    s(HC_bleached, k = 3, bs = "tp") +
    s(DTS_km, k = 5, bs = "tp") +
    s(CCA, k = 3, bs = "tp"),
  data = df,
  family = gaussian(),
  method = "REML",
  select = TRUE
)

# ----------------------------
# 7. Flavobacteriales
# ----------------------------
gam_flavo <- gam(
  WAT_Flavobacteriales ~ 
    s(Algae, k = 6, bs = "tp") +
    s(Mn, k = 6, bs = "tp") +
    s(SST, k = 6, bs = "tp"),
  data = df,
  family = gaussian(),
  method = "REML",
  select = TRUE
)

# ----------------------------
# 8. Synechococcus
# ----------------------------
gam_syne <- gam(
  WAT_Synechococcus ~ 
    s(DTS_km, k = 3, bs = "tp") +
    s(SST, k = 9, bs = "tp") +
    s(log2(Chla), k = 5, bs = "tp"),
  data = df,
  family = gaussian(),
  method = "REML",
  select = TRUE
)

# ----------------------------
# 9. Hyphomicrobiales
# ----------------------------
gam_hypho <- gam(
  WAT_Hyphomicrobiales ~ 
    s(Phosphate, k = 3, bs = "tp") +
    s(log2(DTS_km), k = 6, bs = "tp") +
    s(log2(Chla), k = 9, bs = "tp"),
  data = df,
  family = gaussian(),
  method = "REML",
  select = TRUE
)

# ----------------------------
# 10. Nisaea spp.
# ----------------------------
gam_nisaea <- gam(
  WAT_Nisaea_spp ~ 
    s(Turf, k = 3, bs = "tp") +
    s(SST, k = 6, bs = "tp") +
    s(TN, k = 3, bs = "tp"),
  data = df,
  family = gaussian(),
  method = "REML",
  select = TRUE
)

# ----------------------------
# Run and save all GAMs
# ----------------------------
res_shannon <- run_and_save_gam(
  model_obj = gam_shannon,
  response_name = "Shannon index (taxonomy)",
  file_prefix = "gam_shannon_taxonomy",
  path_output = path_output
)

res_function <- run_and_save_gam(
  model_obj = gam_function,
  response_name = "Function richness",
  file_prefix = "gam_function_richness",
  path_output = path_output
)

res_rhodo <- run_and_save_gam(
  model_obj = gam_rhodo,
  response_name = "Rhodobacterales",
  file_prefix = "gam_rhodobacterales",
  path_output = path_output
)

res_pro <- run_and_save_gam(
  model_obj = gam_pro,
  response_name = "Prochlorococcus marinus",
  file_prefix = "gam_prochlorococcus_marinus",
  path_output = path_output
)

res_bacillota <- run_and_save_gam(
  model_obj = gam_bacillota,
  response_name = "Bacillota",
  file_prefix = "gam_bacillota",
  path_output = path_output
)

res_vibrio <- run_and_save_gam(
  model_obj = gam_vibrio,
  response_name = "Vibrio spp.",
  file_prefix = "gam_vibrio_spp",
  path_output = path_output
)

res_flavo <- run_and_save_gam(
  model_obj = gam_flavo,
  response_name = "Flavobacteriales",
  file_prefix = "gam_flavobacteriales",
  path_output = path_output
)

res_syne <- run_and_save_gam(
  model_obj = gam_syne,
  response_name = "Synechococcus",
  file_prefix = "gam_synechococcus",
  path_output = path_output
)

res_hypho <- run_and_save_gam(
  model_obj = gam_hypho,
  response_name = "Hyphomicrobiales",
  file_prefix = "gam_hyphomicrobiales",
  path_output = path_output
)

res_nisaea <- run_and_save_gam(
  model_obj = gam_nisaea,
  response_name = "Nisaea spp.",
  file_prefix = "gam_nisaea_spp",
  path_output = path_output
)

# ----------------------------
# Unified model statistics table
# ----------------------------
all_model_stats <- bind_rows(
  res_shannon$stats,
  res_function$stats,
  res_rhodo$stats,
  res_pro$stats,
  res_bacillota$stats,
  res_vibrio$stats,
  res_flavo$stats,
  res_syne$stats,
  res_hypho$stats,
  res_nisaea$stats
)

write_csv(
  all_model_stats,
  file.path(path_output, "GAM_all_model_stats_submission.csv")
)

# ----------------------------
# Unified smooth terms table
# ----------------------------
all_smooth_terms <- bind_rows(
  res_shannon$smooth_terms %>% mutate(response = "Shannon index (taxonomy)"),
  res_function$smooth_terms %>% mutate(response = "Function richness"),
  res_rhodo$smooth_terms %>% mutate(response = "Rhodobacterales"),
  res_pro$smooth_terms %>% mutate(response = "Prochlorococcus marinus"),
  res_bacillota$smooth_terms %>% mutate(response = "Bacillota"),
  res_vibrio$smooth_terms %>% mutate(response = "Vibrio spp."),
  res_flavo$smooth_terms %>% mutate(response = "Flavobacteriales"),
  res_syne$smooth_terms %>% mutate(response = "Synechococcus"),
  res_hypho$smooth_terms %>% mutate(response = "Hyphomicrobiales"),
  res_nisaea$smooth_terms %>% mutate(response = "Nisaea spp.")
) %>%
  relocate(response, .before = term)

write_csv(
  all_smooth_terms,
  file.path(path_output, "GAM_all_smooth_terms_submission.csv")
)

# ----------------------------
# Unified parametric terms table
# ----------------------------
all_parametric_terms <- bind_rows(
  res_shannon$parametric_terms %>% mutate(response = "Shannon index (taxonomy)"),
  res_function$parametric_terms %>% mutate(response = "Function richness"),
  res_rhodo$parametric_terms %>% mutate(response = "Rhodobacterales"),
  res_pro$parametric_terms %>% mutate(response = "Prochlorococcus marinus"),
  res_bacillota$parametric_terms %>% mutate(response = "Bacillota"),
  res_vibrio$parametric_terms %>% mutate(response = "Vibrio spp."),
  res_flavo$parametric_terms %>% mutate(response = "Flavobacteriales"),
  res_syne$parametric_terms %>% mutate(response = "Synechococcus"),
  res_hypho$parametric_terms %>% mutate(response = "Hyphomicrobiales"),
  res_nisaea$parametric_terms %>% mutate(response = "Nisaea spp.")
) %>%
  relocate(response, .before = term)

write_csv(
  all_parametric_terms,
  file.path(path_output, "GAM_all_parametric_terms_submission.csv")
)