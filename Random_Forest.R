##########################################
# Random Forest Analysis for Water Data 
# #######################################


# 1. Setup ---

# Load necessary libraries
library(here)
library(tidyverse)
library(randomForest)
library(purrr)

# Set working directory
setwd(here("/Users/gilramgl/Documents/R/EffortsGloria/DataAnalysis/February_2026/RF"))

# 2. Load Data ---

water_file <- "model_ave_relab_01032026.csv"
if (!file.exists(water_file)) {
  stop("Data file not found! Please check the file path.")
}

water_data <- read_csv(water_file, show_col_types = FALSE)


#3. Define Variables ---

# Predictor variables for water samples
predictor_variables <- c(
  "Region", "LOAP",   "LOAP_class", "SST", "DTS_km", "TN", "TOC",
  "Cr", "Mn", "Ni", "Cu", "Zn", "As", "Cd", "Pb",
  "Chla", "Nitrite", "Nitrate", "Phosphate",
  "Algae", "HC_dead", "HC_bleached", "HC_nonbleached",
  "Turf", "CCA", "Live_cover", "Coral_richness", "Shannon_corals"
)

# Response variables (WATER groups)
response_variables <- c(
  "WAT_Species_Richness", "WAT_Shannon_Index_ln", "WAT_Shannon_functions","WAT_Function_richness",
  "WAT_Vibrio_spp", "WAT_Pseudoalteromonas", "WAT_Nisaea_spp",
  "WAT_Alteromonadales", "WAT_Rhodobacterales", "WAT_Bacillota",
  "WAT_Bacteroidota", "WAT_Synechococcus", "WAT_Prochlorococcus_marinus",
  "WAT_Cyanobacteriota", "WAT_Hyphomicrobiales",
  "WAT_Flavobacteriales", "WAT_Desulfovibrionales"
)

# Define variable categories for coloring
physical_vars <- c("SST", "DTS_km")
chemical_vars <- c("TN", "TOC", "Cr", "Mn", "Ni", "Cu", "Zn", 
                   "As", "Cd", "Pb", "Chla", "Nitrite", "Nitrate", "Phosphate")
biological_vars <- c("Algae", "HC_dead", "HC_bleached", "HC_nonbleached", "Turf", 
                     "CCA", "Live_cover", "Coral_richness", "Shannon_corals")
anthropogenic_vars <- c("Region", "LOAP", "LOAP_class")

category_levels <- c("Physical", "Chemical", "Biological", "Anthropogenic")
category_colors <- c(
  "Physical" ="#c5bd96",    # golden
  "Chemical" = "#173f5f",   # blue
  "Biological" ="#93b0ac", #Green
  "Anthropogenic" = "gray" # Gray
)

# Function to assign variable category
get_category <- function(variable) {
  if (variable %in% physical_vars) return("Physical")
  if (variable %in% chemical_vars) return("Chemical")
  if (variable %in% biological_vars) return("Biological")
  if (variable %in% anthropogenic_vars) return("Anthropogenic")
  return(NA_character_)
}

#4. Output Folder ---

# Use today's date automatically (e.g., 13OCT2025)
.today_tag <- toupper(format(Sys.Date(), "%d%b%Y"))
plot_dir <- paste0("RF_Importance_Plots_Wat_", .today_tag)
dir.create(plot_dir, showWarnings = FALSE)

#5. Set Seed ---

set.seed(2025)

#6. Run Random Forests and Save Results ---

rf_summary_water <- map_df(response_variables, function(response_var) {
  # Check variable exists
  if (!response_var %in% names(water_data)) {
    warning(paste("Missing response variable:", response_var))
    return(NULL)
  }
  
  # Prepare model data
  rf_data <- water_data %>%
    select(all_of(c(predictor_variables, response_var))) %>%
    drop_na()
  
  if (nrow(rf_data) < 10) {
    warning(paste("Too few observations for:", response_var))
    return(NULL)
  }
  
  # Ensure sensible types for predictors (RF can handle numeric + factors)
  # If LOAP_class is character, coerce to factor
  if ("LOAP_class" %in% names(rf_data) && !is.factor(rf_data$LOAP_class)) {
    rf_data <- rf_data %>% mutate(LOAP_class = as.factor(LOAP_class))
  }
  
  # Train RF model
  rf_model <- randomForest(
    formula = as.formula(paste(response_var, "~ .")),
    data = rf_data,
    ntree = 1000,
    importance = TRUE
  )
  
  # Extract %IncMSE and sort
  importance_df <- importance(rf_model, type = 1) %>%
    as.data.frame() %>%
    rownames_to_column("Variable") %>%
    rename(Importance = `%IncMSE`) %>%
    arrange(desc(Importance)) %>%
    mutate(Category = factor(sapply(Variable, get_category), levels = category_levels))
  
  # Plot
  p <- ggplot(importance_df, aes(x = Importance, y = reorder(Variable, Importance), fill = Category)) +
    geom_col() +
    scale_fill_manual(values = category_colors, na.value = "grey50") +
    labs(
      title = paste("Random Forest Variable Importance for", response_var),
      subtitle = paste0("Variance Explained: ", round(100 * rf_model$rsq[length(rf_model$rsq)], 1), "%"),
      x = "% Increase in MSE", y = NULL, fill = "Variable Category"
    ) +
    theme_bw(base_size = 14)
  
  ggsave(
    filename = file.path(plot_dir, paste0("RF_Importance_Wat_", response_var, "_", .today_tag, ".jpg")),
    plot = p, width = 9, height = 6, dpi = 300
  )
  
  saveRDS(importance_df, file = file.path(plot_dir, paste0("RF_Immdf_", response_var, ".RDS")))
  saveRDS(p, file = file.path(plot_dir, paste0("RF_RFplot_", response_var, ".RDS")))
  
  # Top 10 predictors
  top_predictors <- importance_df %>% slice(1:10)
  
  tibble(
    Response = response_var,
    VarianceExplained = round(100 * rf_model$rsq[length(rf_model$rsq)], 1),
    TopPredictor1 = top_predictors$Variable[1], IncMSE1 = round(top_predictors$Importance[1], 1),
    TopPredictor2 = top_predictors$Variable[2], IncMSE2 = round(top_predictors$Importance[2], 1),
    TopPredictor3 = top_predictors$Variable[3], IncMSE3 = round(top_predictors$Importance[3], 1),
    TopPredictor4 = top_predictors$Variable[4], IncMSE4 = round(top_predictors$Importance[4], 1),
    TopPredictor5 = top_predictors$Variable[5], IncMSE5 = round(top_predictors$Importance[5], 1),
    TopPredictor6 = top_predictors$Variable[6], IncMSE6 = round(top_predictors$Importance[6], 1),
    TopPredictor7 = top_predictors$Variable[7], IncMSE7 = round(top_predictors$Importance[7], 1),
    TopPredictor8 = top_predictors$Variable[8], IncMSE8 = round(top_predictors$Importance[8], 1),
    TopPredictor9 = top_predictors$Variable[9], IncMSE9 = round(top_predictors$Importance[9], 1),
    TopPredictor10 = top_predictors$Variable[10], IncMSE10 = round(top_predictors$Importance[10], 1)
  )
})

# --- 7. Save Summary Table ---

summary_path <- file.path(plot_dir, paste0("RF_Summary", .today_tag, ".csv"))
write_csv(rf_summary_water, summary_path)

# --- 8. Final Message ---

message("Random Forest for WATER complete! Plots saved to '", plot_dir, "' and summary to '", basename(summary_path), "'.")

multi_taxa_plot <- wrap_plots(
  res_syne$plot,
  res_flavo$plot,
  res_hypho$plot,
  res_fourth$plot,
  ncol = 2
) +
  plot_annotation(
    title = "GAM partial effects across focal bacterial groups"
  )

ggsave(
  filename = file.path(path_output, "gam_multitaxa_partial_effects.png"),
  plot = multi_taxa_plot,
  width = 18,
  height = 16,
  dpi = 300
)

ggsave(
  filename = file.path(path_output, "gam_multitaxa_partial_effects.pdf"),
  plot = multi_taxa_plot,
  width = 18,
  height = 16
)

all_model_stats <- bind_rows(
  res_syne$stats,
  res_flavo$stats,
  res_hypho$stats,
  res_fourth$stats
)

write_csv(
  all_model_stats,
  file.path(path_output, "gam_all_model_stats.csv")
)


