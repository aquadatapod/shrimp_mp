# ============================================================================
# HETEROGENEITY ANALYSIS FOR MICROPLASTIC META-ANALYSIS
# Includes: Season, Method, Environment, Catch Method, Country, Year
# ============================================================================

# Load required packages
library(tidyverse)
library(metafor)
library(ggplot2)
library(patchwork)
library(cowplot)
library(openxlsx)

# Set consistent theme for all plots
theme_hetero_base <- theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(face = "bold", size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )


# ============================================================================
# MULTILEVEL META-ANALYSIS WITH SENSITIVITY CHECKS
# ============================================================================

df <- read_excel("shrimp_data001.xlsx", sheet = "shrimp")

# 1. CREATE UNIQUE IDs
df <- df %>%
  mutate(
    StudyID = as.factor(Reference),  # cluster by study
    ObsID = 1:nrow(df)  # unique observation ID
  )

# 2. CHECK DEPENDENCY STRUCTURE
dependency_summary <- df %>%
  group_by(StudyID) %>%
  summarise(
    n_effects = n(),
    n_species = n_distinct(Species),
    n_sites = n_distinct(Site_number),
    n_years = n_distinct(Years)
  )

print(dependency_summary)
write.xlsx(dependency_summary, "Table_Dependency_Summary.xlsx")



# ----------------------------------------------------------------------------
# PART 1: DATA PREPARATION
# ----------------------------------------------------------------------------

cat("\n📋 STEP 1: Preparing data...\n")

# Calculate effect sizes if not already done
if(!"yi" %in% names(df)) {
  df$yi <- log(df$MPperind + 0.1)
  df$vi <- (0.5^2) / df$SpeciesNo
  df$sei <- sqrt(df$vi)
}

# Create unique study ID
df <- df %>%
  mutate(
    StudyID = as.factor(Reference),
    ObsID = 1:nrow(df)
  )
#-------------------------------------------------------------------------------

# 3. MAIN ANALYSIS: 3-LEVEL META-ANALYSIS
res_mlma <- rma.mv(
  yi = yi,
  V = vi,
  random = list(~ 1 | StudyID, ~ 1 | ObsID),
  data = df,
  method = "REML",
  test = "t",  # t-tests for small samples
  dfs = "contain"  # containment method for degrees of freedom
)

summary(res_mlma)

# 4. COMPARE WITH 2-LEVEL MODEL (ignoring dependency)
res_2level <- rma(yi = yi, vi = vi, data = df, method = "REML")

# 5. TEST IF WITHIN-STUDY VARIANCE IS SIGNIFICANT
# Likelihood ratio test
res_null <- rma.mv(
  yi = yi,
  V = vi,
  random = list(~ 1 | StudyID),
  data = df,
  method = "REML"
)

# Compare models
anova(res_null, res_mlma)  # Significant? Then within-study variance matters!

# 6. SENSITIVITY ANALYSIS 1: AGGREGATE BY STUDY
df_agg <- df %>%
  group_by(StudyID, Reference, Years, Country, quality_score) %>%
  summarise(
    yi_agg = weighted.mean(yi, 1/vi),
    vi_agg = 1 / sum(1/vi),
    n_original = n(),
    .groups = "drop"
  )

res_agg <- rma(yi = yi_agg, vi = vi_agg, data = df_agg, method = "REML")

library(clubSandwich)
# 7. SENSITIVITY ANALYSIS 2: ROBUST VARIANCE ESTIMATION
res_rve <- robust(res_mlma, cluster = StudyID, clubSandwich = TRUE)

# 8. COMPARE ALL RESULTS
comparison <- data.frame(
  Method = c("3-level MLMA", "2-level (ignoring dependency)", 
             "Aggregated", "Robust Variance"),
  Estimate = c(coef(res_mlma), coef(res_2level), coef(res_agg), coef(res_rve)),
  SE = c(res_mlma$se, res_2level$se, res_agg$se, res_rve$se),
  p_value = c(res_mlma$pval, res_2level$pval, res_agg$pval, res_rve$pval),
  CI_lower = c(res_mlma$ci.lb, res_2level$ci.lb, res_agg$ci.lb, res_rve$ci.lb),
  CI_upper = c(res_mlma$ci.ub, res_2level$ci.ub, res_agg$ci.ub, res_rve$ci.ub)
)

write.xlsx(comparison, "Table_Method_Comparison.xlsx")


# ============================================================================
# FOREST PLOT FOR SPECIES (Panel A)
# ============================================================================

cat("\n📊 Creating Species Forest Plot (Panel A)...\n")

# Run meta-analysis for each species
species_list <- unique(df$Species[!is.na(df$Species)])
species_results <- data.frame()

for(sp in species_list) {
  df_sub <- subset(df, Species == sp)
  if(nrow(df_sub) >= 3) {  # Need at least 3 observations
    res_sub <- tryCatch(
      rma(yi = yi, vi = vi, data = df_sub, method = "REML"),
      error = function(e) NULL
    )
    if(!is.null(res_sub)) {
      species_results <- rbind(species_results,
                               data.frame(
                                 Species = sp,
                                 k = res_sub$k,
                                 Estimate = res_sub$b,
                                 ci.lb = res_sub$ci.lb,
                                 ci.ub = res_sub$ci.ub,
                                 stringsAsFactors = FALSE
                               ))
    }
  }
}

# Create effect size groups for coloring (based on tertiles)
species_results <- species_results %>%
  mutate(
    # Create clean label with species name and k
    label = paste0(Species, " (k=", k, ")"),
    
    # Create effect size groups for coloring
    effect_group = case_when(
      Estimate > quantile(Estimate, 0.66) ~ "High",
      Estimate > quantile(Estimate, 0.33) ~ "Medium",
      TRUE ~ "Low"
    ),
    effect_group = factor(effect_group, levels = c("High", "Medium", "Low"))
  ) %>%
  arrange(desc(Estimate))

# Create species forest plot
p_species_forest <- ggplot(species_results, 
                           aes(x = Estimate, y = reorder(label, Estimate), 
                               color = effect_group)) +
  # Vertical reference line at 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.3) +
  
  # Add points and error bars
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci.lb, xmax = ci.ub), 
                 height = 0.2, size = 0.6) +
  
  # Color by effect size group
  scale_color_manual(values = c("High" = "#DB5829", "Medium" = "#5F6D7A", "Low" = "#E9DC6D")) +
  
  # Labels and title
  labs(x = "Log Response Ratio (95% CI)", 
       y = "",
       title = "b. Species") +
  
  # Theme 
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.text.y = element_text(size = 8, face = "italic"),  # Italic for species names
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(face = "bold", size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", size = 0.3),
    legend.position = "none"
  )

# Save the plot
#ggsave("Panel_A_Species_Forest.png", p_species_forest, 
#width = 6, height = nrow(species_results)*0.25 + 2, dpi = 300)

print(p_species_forest)
cat("✅ Saved: Panel_A_Species_Forest.png\n")


# ============================================================================
# FOREST PLOT FOR STUDIES (Panel B) - Aggregated by Study
# ============================================================================

cat("\n📊 Creating Study Forest Plot (Panel B)...\n")

# Aggregate by study (one effect per study)
study_results <- df %>%
  group_by(Reference, Years, Country, quality_score) %>%
  summarise(
    yi_agg = weighted.mean(yi, 1/vi, na.rm = TRUE),
    vi_agg = 1 / sum(1/vi, na.rm = TRUE),
    n_effects = n(),
    ci_lower = yi_agg - 1.96 * sqrt(vi_agg),
    ci_upper = yi_agg + 1.96 * sqrt(vi_agg),
    .groups = "drop"
  ) %>%
  # Create effect size groups for coloring (based on tertiles)
  mutate(
    # Truncate long study names for readability
    short_ref = ifelse(nchar(Reference) > 30,
                       paste0(substr(Reference, 1, 27), "..."),
                       Reference),
    label = paste0(short_ref, " (", n_effects, " effects)"),
    
    # Create effect size groups for coloring
    effect_group = case_when(
      yi_agg > quantile(yi_agg, 0.66) ~ "High",
      yi_agg > quantile(yi_agg, 0.33) ~ "Medium",
      TRUE ~ "Low"
    ),
    effect_group = factor(effect_group, levels = c("High", "Medium", "Low"))
  ) %>%
  arrange(desc(yi_agg)) %>%
  slice(1:min(20, n()))  # Show top 20 studies for readability

# Create study forest plot
p_study_forest <- ggplot(study_results, 
                         aes(x = yi_agg, y = reorder(label, yi_agg), 
                             color = effect_group)) +
  # Vertical reference line at 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.3) +
  
  # Add overall meta-analysis estimate as dotted line
  geom_vline(xintercept = coef(res_mlma), 
             linetype = "dotted", color = "red", size = 0.5) +
  
  # Add points and error bars
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), 
                 height = 0.2, size = 0.6) +
  
  # Color by effect size group
  scale_color_manual(values = c("High" = "#DB5829", "Medium" = "#5F6D7A", "Low" = "#E9DC6D")) +
  
  # Labels and title
  labs(x = "Log Response Ratio (95% CI)", 
       y = "",
       title = "a. Studies") +
  
  # Theme 
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(face = "bold", size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", size = 0.3),
    legend.position = "none"
  )



# ----------------------------------------------------------------------------
# PART 2: SEASON ANALYSIS
# ----------------------------------------------------------------------------

cat("\n📊 STEP 2: Analyzing Season...\n")

# Clean and categorize season data
df <- df %>%
  mutate(
    season_simple = case_when(
      # Single seasons
      season == "Summer" ~ "Summer",
      season == "Winter" ~ "Winter",
      season == "Autumn" ~ "Autumn",
      season == "Spring" ~ "Spring",
      season == "Monsoon" ~ "Monsoon",
      
      # Two-season combinations
      season == "Spring/Summer" ~ "Spring/Summer",
      season == "Autumn/Winter" ~ "Autumn/Winter",
      season == "Winter/Monsoon" ~ "Winter/Monsoon",
      season == "Spring/Autumn" ~ "Spring/Autumn",
      
      # Three-season combinations
      season == "Summer/Autumn/Winter" ~ "Summer/Autumn/Winter",
      season == "Summer/Monsoon/Winter" ~ "Summer/Monsoon/Winter",
      
      # Four-season (annual)
      season == "Autumn/Winter/Spring/Summer" ~ "Annual (All Seasons)",
      
      # Typos/variations
      season == "Monsoom/Winter" ~ "Winter/Monsoon",
      season == "Visula sorting" ~ NA_character_,
      
      # NA or missing
      is.na(season) | season == "NA" ~ "Not specified",
      
      # Default
      TRUE ~ "Other"
    )
  )

# Season groups for coloring
df <- df %>%
  mutate(
    season_group = case_when(
      season_simple %in% c("Summer", "Spring/Autumn") ~ "High_Season",
      season_simple %in% c("Autumn", "Winter", "Spring/Summer", "Annual (All Seasons)") ~ "Medium_Season",
      season_simple %in% c("Monsoon", "Winter/Monsoon") ~ "Low_Season",
      season_simple %in% c("Autumn/Winter", "Summer/Monsoon/Winter") ~ "Minimal_Season",
      TRUE ~ "Other"
    ),
    season_group = factor(season_group, 
                          levels = c("High_Season", "Medium_Season", "Low_Season", "Minimal_Season", "Other"))
  )

# Season boxplot
season_plot_data <- df %>%
  filter(!season_simple %in% c("Not specified", "Other"), !is.na(season_simple)) %>%
  group_by(season_simple) %>%
  filter(n() >= 3)

n_seasons <- length(unique(season_plot_data$season_simple))
season_colors <- colorRampPalette(c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
                                    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"))(n_seasons)

p_season_box <- ggplot(season_plot_data, aes(x = season_simple, y = yi, fill = season_simple)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = season_simple), width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_manual(values = season_colors) +
  scale_color_manual(values = season_colors) +
  labs(x = "Season", y = "Effect Size (Log Response Ratio)") +
  theme_hetero_base +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_blank()) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  stat_summary(fun = mean, geom = "text", aes(label = paste0("mean=", round(after_stat(y), 2))),
               vjust = -0.5, color = "black", size = 3)

#ggsave("Heterogeneity_Season_Boxplot.png", p_season_box, width = 12, height = 6, dpi = 300)

# Season forest plot (subgroup meta-analysis)
season_categories <- unique(df$season_simple[!is.na(df$season_simple) & 
                                               !df$season_simple %in% c("Not specified", "Other")])

season_results <- data.frame()

for(s in season_categories) {
  df_sub <- subset(df, season_simple == s)
  if(nrow(df_sub) >= 3) {
    res_sub <- tryCatch(
      rma(yi = yi, vi = vi, data = df_sub, method = "REML"),
      error = function(e) NULL
    )
    if(!is.null(res_sub)) {
      season_results <- rbind(season_results,
                              data.frame(
                                Season = s,
                                k = res_sub$k,
                                Estimate = res_sub$b,
                                ci.lb = res_sub$ci.lb,
                                ci.ub = res_sub$ci.ub,
                                stringsAsFactors = FALSE
                              ))
    }
  }
}


# First, ensure season_results has proper grouping
season_results <- season_results %>%
  mutate(
    label = paste0(Season, " (k=", k, ")"),
    group = case_when(
      Season %in% c("Summer", "Spring/Autumn") ~ "High_Season",
      Season %in% c("Autumn", "Winter", "Spring/Summer", "Annual (All Seasons)") ~ "Medium_Season",
      Season %in% c("Monsoon", "Winter/Monsoon") ~ "Low_Season",
      Season %in% c("Autumn/Winter", "Summer/Monsoon/Winter") ~ "Minimal_Season",
      TRUE ~ "Other"
    ),
    # Create clean group names for facet labels
    group_label = case_when(
      group == "High_Season" ~ "High Season",
      group == "Medium_Season" ~ "Medium Season", 
      group == "Low_Season" ~ "Low Season",
      group == "Minimal_Season" ~ "Minimal Season",
      TRUE ~ "Other"
    ),
    # Create ordering within groups
    group_order = case_when(
      group == "High_Season" ~ 1,
      group == "Medium_Season" ~ 2,
      group == "Low_Season" ~ 3,
      group == "Minimal_Season" ~ 4,
      TRUE ~ 5
    )
  ) %>%
  # Arrange by group, then by Estimate within group
  arrange(group_order, desc(Estimate)) %>%
  mutate(label = factor(label, levels = label))

# Define colors 
season_group_colors <- c(
  "High_Season" = "#DB5829",    # Burnt orange
  "Medium_Season" = "#5F6D7A",  # Slate gray
  "Low_Season" = "#E9DC6D",     # Soft yellow
  "Minimal_Season" = "#800080", # Purple
  "Other" = "#B0C4DE"           # Light steel blue
)

# Create faceted season plot
p_season_facet <- ggplot(season_results, 
                         aes(x = Estimate, y = reorder(label, Estimate), 
                             color = group)) +
  # Reference lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
  geom_vline(xintercept = mean(season_results$Estimate), 
             linetype = "dotted", color = "red", size = 0.5) +
  
  # Data points
  geom_point(size = 3.5) +
  geom_errorbarh(aes(xmin = ci.lb, xmax = ci.ub), height = 0.2, size = 0.8) +
  
  # Color
  scale_color_manual(values = season_group_colors) +
  
  # Facet by group - THIS IS THE KEY!
  facet_grid(group_label ~ ., scales = "free_y", space = "free") +
  
  # Labels
  labs(x = "Log Response Ratio (95% CI)", 
       y = "", 
       title = "e. Season") +
  
  # Theme
  theme_hetero_base +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95", color = NA),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    panel.spacing = unit(0.5, "lines")  # Add space between facets
  )

# Display
print(p_season_facet)

# Save
ggsave("Season_Forest_Plot_Faceted.png", p_season_facet, 
       width = 8, height = 8, dpi = 300)

# Also save a version with white background for publication
ggsave("Season_Forest_Plot_Faceted_publication.png", p_season_facet, 
       width = 8, height = 8, dpi = 600, bg = "white")


# ----------------------------------------------------------------------------
# PART 3: EXTRACTION METHOD ANALYSIS
# ----------------------------------------------------------------------------

cat("\n📊 STEP 3: Analyzing Extraction Method...\n")

# Clean and categorize method data
df <- df %>%
  mutate(
    method_simple = case_when(
      # KOH-based methods
      extraction_method == "KOH" ~ "KOH only",
      extraction_method == "KOH/H2O2" ~ "KOH + H2O2",
      extraction_method == "KOH/HNO3" ~ "KOH + HNO3",
      extraction_method == "KOH+ASE" ~ "KOH + ASE",
      
      # NaOH-based methods
      extraction_method == "NaOH" ~ "NaOH only",
      extraction_method == "H2O2+NAOH" ~ "NaOH + H2O2",
      
      # H2O2-based methods
      extraction_method == "H2O2" ~ "H2O2 only",
      
      # Acid-based methods
      extraction_method == "HNO3" ~ "HNO3 only",
      extraction_method == "HNO3/HCIO4" ~ "HNO3 + HClO4",
      extraction_method == "H2O2/HNO3" ~ "HNO3 + H2O2",
      extraction_method == "H2O2/HNO3/HCIO4" ~ "HNO3 + HClO4 + H2O2",
      
      # Multi-reagent combinations
      extraction_method == "H2O2/KOH/HCIO4/HNO3" ~ "Multiple (H2O2+KOH+HClO4+HNO3)",
      
      # No digestion
      extraction_method %in% c("Visual sorting", "Visula sorting") ~ "Visual sorting only",
      
      # NA or missing
      is.na(extraction_method) | extraction_method == "NA" ~ "Not specified",
      
      # Default
      TRUE ~ "Other"
    )
  )

# Method groups for coloring
df <- df %>%
  mutate(
    method_group = case_when(
      method_simple %in% c("Multiple (H2O2+KOH+HClO4+HNO3)", "NaOH + H2O2", "HNO3 only") ~ "High_Recovery",
      method_simple %in% c("HNO3 + HClO4", "H2O2 only", "KOH only") ~ "Medium_Recovery",
      method_simple %in% c("KOH + HNO3", "HNO3 + HClO4 + H2O2", "Visual sorting only") ~ "Low_Recovery",
      method_simple == "NaOH only" ~ "Problematic",
      TRUE ~ "Other"
    ),
    method_group = factor(method_group, 
                          levels = c("High_Recovery", "Medium_Recovery", "Low_Recovery", "Problematic", "Other"))
  )

# Method boxplot
method_plot_data <- df %>%
  filter(!method_simple %in% c("Not specified", "Other"), !is.na(method_simple)) %>%
  group_by(method_simple) %>%
  filter(n() >= 3)



n_methods <- length(unique(method_plot_data$method_simple))
method_colors <- colorRampPalette(c("#DB5829", "#5F6D7A", "#E9DC6D", "#4682B4", 
                                    "#9370DB", "#CD5C5C", "#20B2AA"))(n_methods)

p_method_box <- ggplot(method_plot_data, aes(x = reorder(method_simple, yi, median), y = yi, fill = method_simple)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = method_simple), width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  labs(x = "Extraction Method", y = "Effect Size (Log Response Ratio)") +
  theme_hetero_base +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_blank()) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  stat_summary(fun = mean, geom = "text", aes(label = paste0("mean=", round(after_stat(y), 2))),
               vjust = -0.5, color = "black", size = 3)

#ggsave("Heterogeneity_Method_Boxplot.png", p_method_box, width = 14, height = 6, dpi = 300)

# Method forest plot
method_categories <- unique(df$method_simple[!is.na(df$method_simple) & 
                                               !df$method_simple %in% c("Not specified", "Other")])

method_results <- data.frame()

for(m in method_categories) {
  df_sub <- subset(df, method_simple == m)
  if(nrow(df_sub) >= 3) {
    res_sub <- tryCatch(
      rma(yi = yi, vi = vi, data = df_sub, method = "REML"),
      error = function(e) NULL
    )
    if(!is.null(res_sub)) {
      method_results <- rbind(method_results,
                              data.frame(
                                Method = m,
                                k = res_sub$k,
                                Estimate = res_sub$b,
                                ci.lb = res_sub$ci.lb,
                                ci.ub = res_sub$ci.ub,
                                stringsAsFactors = FALSE
                              ))
    }
  }
}

method_results <- method_results %>%
  mutate(
    label = paste0(Method, " (k=", k, ")"),
    group = case_when(
      Method %in% c("Multiple (H2O2+KOH+HClO4+HNO3)", "NaOH + H2O2", "HNO3 only") ~ "High_Recovery",
      Method %in% c("HNO3 + HClO4", "H2O2 only", "KOH only") ~ "Medium_Recovery",
      Method %in% c("KOH + HNO3", "HNO3 + HClO4 + H2O2", "Visual sorting only") ~ "Low_Recovery",
      Method == "NaOH only" ~ "Problematic",
      TRUE ~ "Other"
    )
  )

method_order <- c("Multiple (H2O2+KOH+HClO4+HNO3)", "NaOH + H2O2", "HNO3 only",
                  "HNO3 + HClO4", "H2O2 only", "KOH only",
                  "KOH + HNO3", "HNO3 + HClO4 + H2O2", "KOH + H2O2", "KOH + ASE",
                  "Visual sorting only", "NaOH only")

method_results <- method_results %>%
  mutate(Method = factor(Method, levels = method_order[method_order %in% Method])) %>%
  arrange(Method) %>%
  mutate(label = factor(label, levels = label))


method_group_colors <- c("High_Recovery" = "#DB5829", "Medium_Recovery" = "#5F6D7A",
                         "Low_Recovery" = "#E9DC6D", "Problematic" = "#800080", "Other" = "#B0C4DE")

p_method_facet <- ggplot(method_results, aes(x = Estimate, y = reorder(label, Estimate), color = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
  geom_vline(xintercept = mean(method_results$Estimate), 
             linetype = "dotted", color = "red", size = 0.5) +
  geom_point(size = 3.5) +
  geom_errorbarh(aes(xmin = ci.lb, xmax = ci.ub), height = 0.2, size = 0.8) +
  scale_color_manual(values = method_group_colors) +
  
  # Facet by group
  facet_grid(group ~ ., scales = "free_y", space = "free") +
  
  labs(x = "Log Response Ratio (95% CI)", y = "", title = "d. Extraction Method") +
  
  theme_hetero_base +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "gray95", color = NA),
    axis.text.x = element_text(size = 8),
    legend.position = "none"
  )

ggsave("Method_Forest_Plot_Faceted.png", p_method_facet, 
       width = 9, height = 8, dpi = 300)


# ============================================================================
# PART 4: ENVIRONMENT ANALYSIS 
# ============================================================================

cat("\n📊 STEP 4: Analyzing Environment with color coding...\n")

# Decode Environment codes
df <- df %>%
  mutate(
    Environment_full = case_when(
      Environment == "ES" ~ "Estuarine",
      Environment == "MA" ~ "Marine",
      Environment == "FW" ~ "Freshwater",
      Environment == "PO" ~ "Pond",
      Environment == "LN" ~ "Lagoon",
      Environment == "RF" ~ "Reef",
      TRUE ~ Environment
    )
  )

env_categories <- unique(df$Environment_full[!is.na(df$Environment_full)])
env_results <- data.frame()

for(e in env_categories) {
  df_sub <- subset(df, Environment_full == e)
  if(nrow(df_sub) >= 3) {
    res_sub <- tryCatch(
      rma(yi = yi, vi = vi, data = df_sub, method = "REML"),
      error = function(e) NULL
    )
    if(!is.null(res_sub)) {
      env_results <- rbind(env_results,
                           data.frame(
                             Environment = e,
                             k = res_sub$k,
                             Estimate = res_sub$b,
                             ci.lb = res_sub$ci.lb,
                             ci.ub = res_sub$ci.ub,
                             stringsAsFactors = FALSE
                           ))
    }
  }
}

# Add effect size groups for coloring (based on tertiles)
env_results <- env_results %>%
  mutate(
    label = paste0(Environment, " (k=", k, ")"),
    # Create effect size groups for consistent coloring
    effect_group = case_when(
      Estimate > quantile(Estimate, 0.66) ~ "High",
      Estimate > quantile(Estimate, 0.33) ~ "Medium",
      TRUE ~ "Low"
    ),
    effect_group = factor(effect_group, levels = c("High", "Medium", "Low"))
  ) %>%
  arrange(desc(Estimate))




# Create forest plot with proper color scheme
p_env_forest <- ggplot(env_results, 
                       aes(x = Estimate, y = reorder(label, Estimate), 
                           color = effect_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
  geom_vline(xintercept = mean(env_results$Estimate), 
             linetype = "dotted", color = "red", size = 0.5) +
  geom_point(size = 3.5) +
  geom_errorbarh(aes(xmin = ci.lb, xmax = ci.ub), 
                 height = 0.2, size = 0.8) +
  scale_color_manual(values = c("High" = "#DB5829", 
                                "Medium" = "#5F6D7A", 
                                "Low" = "#E9DC6D")) +
  labs(x = "Log Response Ratio (95% CI)", 
       y = "", 
       title = "f. Environment") +
  theme_hetero_base

# Save the plot
#ggsave("Forest_Plot_Environment.png", p_env_forest, width = 7, height = 4.5, dpi = 300)
print(p_env_forest)
cat("✅ Saved: Forest_Plot_Environment.png\n")

# ============================================================================
# PART 5: CATCH METHOD ANALYSIS 
# ============================================================================

cat("\n📊 STEP 5: Analyzing Catch Method with color coding...\n")

# Decode Catch_method codes
df <- df %>%
  mutate(
    Catch_method_full = case_when(
      Catch_method == "TN" ~ "Trawl Net",
      Catch_method == "FM" ~ "Fish Market",
      Catch_method == "DD" ~ "Dredge",
      Catch_method == "BT" ~ "Beam Trawl",
      Catch_method == "CN" ~ "Cast Net",
      Catch_method == "ZN" ~ "Zooplankton Net",
      Catch_method == "AP" ~ "Aquaculture Pond",
      Catch_method == "MN" ~ "Macroinvertebrate Net",
      Catch_method == "PN" ~ "Pond Net",
      Catch_method == "AT" ~ "Agassiz Trawl",
      Catch_method == "HN" ~ "Hoop Net",
      TRUE ~ Catch_method
    )
  )

catch_categories <- unique(df$Catch_method_full[!is.na(df$Catch_method_full)])
catch_results <- data.frame()

for(c in catch_categories) {
  df_sub <- subset(df, Catch_method_full == c)
  if(nrow(df_sub) >= 3) {
    res_sub <- tryCatch(
      rma(yi = yi, vi = vi, data = df_sub, method = "REML"),
      error = function(e) NULL
    )
    if(!is.null(res_sub)) {
      catch_results <- rbind(catch_results,
                             data.frame(
                               Catch_Method = c,
                               k = res_sub$k,
                               Estimate = res_sub$b,
                               ci.lb = res_sub$ci.lb,
                               ci.ub = res_sub$ci.ub,
                               stringsAsFactors = FALSE
                             ))
    }
  }
}

# Add effect size groups for coloring (based on tertiles)
catch_results <- catch_results %>%
  mutate(
    label = paste0(Catch_Method, " (k=", k, ")"),
    # Create effect size groups for consistent coloring
    effect_group = case_when(
      Estimate > quantile(Estimate, 0.66) ~ "High",
      Estimate > quantile(Estimate, 0.33) ~ "Medium",
      TRUE ~ "Low"
    ),
    effect_group = factor(effect_group, levels = c("High", "Medium", "Low"))
  ) %>%
  arrange(desc(Estimate))


# Create forest plot with proper color scheme
p_catch_forest <- ggplot(catch_results, 
                         aes(x = Estimate, y = reorder(label, Estimate), 
                             color = effect_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
  geom_vline(xintercept = mean(catch_results$Estimate), 
             linetype = "dotted", color = "red", size = 0.5) +
  geom_point(size = 3.5) +
  geom_errorbarh(aes(xmin = ci.lb, xmax = ci.ub), 
                 height = 0.2, size = 0.8) +
  scale_color_manual(values = c("High" = "#DB5829", 
                                "Medium" = "#5F6D7A", 
                                "Low" = "#E9DC6D")) +
  labs(x = "Log Response Ratio (95% CI)", 
       y = "", 
       title = "g. Catch Method") +
  theme_hetero_base

# Save the plot
plot_height_catch <- max(4, nrow(catch_results) * 0.3)
#ggsave("Forest_Plot_Catch_Method.png", p_catch_forest, 
#      width = 8, height = plot_height_catch, dpi = 300)
print(p_catch_forest)
cat("✅ Saved: Forest_Plot_Catch_Method.png\n")


# ============================================================================
# PART 6: COUNTRY ANALYSIS
# ============================================================================

cat("\n📊 STEP 6: Analyzing Country with color coding...\n")

country_categories <- unique(df$Country[!is.na(df$Country) & df$Country != "NA"])
country_results <- data.frame()

for(c in country_categories) {
  df_sub <- subset(df, Country == c)
  if(nrow(df_sub) >= 3) {
    res_sub <- tryCatch(
      rma(yi = yi, vi = vi, data = df_sub, method = "REML"),
      error = function(e) NULL
    )
    if(!is.null(res_sub)) {
      country_results <- rbind(country_results,
                               data.frame(
                                 Country = c,
                                 k = res_sub$k,
                                 Estimate = res_sub$b,
                                 ci.lb = res_sub$ci.lb,
                                 ci.ub = res_sub$ci.ub,
                                 stringsAsFactors = FALSE
                               ))
    }
  }
}

# Add effect size groups for coloring (based on tertiles)
country_results <- country_results %>%
  mutate(
    label = paste0(Country, " (k=", k, ")"),
    # Create effect size groups for consistent coloring
    effect_group = case_when(
      Estimate > quantile(Estimate, 0.66) ~ "High",
      Estimate > quantile(Estimate, 0.33) ~ "Medium",
      TRUE ~ "Low"
    ),
    effect_group = factor(effect_group, levels = c("High", "Medium", "Low"))
  ) %>%
  arrange(desc(Estimate)) %>%
  slice(1:min(8, n()))  # Top 8 countries


# Create forest plot with proper color scheme
p_country_forest <- ggplot(country_results, 
                           aes(x = Estimate, y = reorder(label, Estimate), 
                               color = effect_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
  geom_vline(xintercept = mean(country_results$Estimate), 
             linetype = "dotted", color = "red", size = 0.5) +
  geom_point(size = 3.5) +
  geom_errorbarh(aes(xmin = ci.lb, xmax = ci.ub), 
                 height = 0.2, size = 0.8) +
  scale_color_manual(values = c("High" = "#DB5829", 
                                "Medium" = "#5F6D7A", 
                                "Low" = "#E9DC6D")) +
  labs(x = "Log Response Ratio (95% CI)", 
       y = "", 
       title = "h. Country (Top 8)") +
  theme_hetero_base

# Save the plot
#ggsave("Forest_Plot_Country.png", p_country_forest, 
#      width = 7, height = 4, dpi = 300)
print(p_country_forest)
cat("✅ Saved: Forest_Plot_Country.png\n")

# ----------------------------------------------------------------------------
# PART 7: QUALITY SCORE ANALYSIS
# ----------------------------------------------------------------------------

cat("\n📊 STEP 7: Analyzing Quality Score...\n")

quality_groups <- c("High", "Medium", "Low")
quality_results <- data.frame()

for(q in quality_groups) {
  df_sub <- subset(df, quality_score == q)
  if(nrow(df_sub) >= 3) {
    res_sub <- tryCatch(
      rma(yi = yi, vi = vi, data = df_sub, method = "REML"),
      error = function(e) NULL
    )
    if(!is.null(res_sub)) {
      quality_results <- rbind(quality_results,
                               data.frame(
                                 Quality = q,
                                 k = res_sub$k,
                                 Estimate = res_sub$b,
                                 ci.lb = res_sub$ci.lb,
                                 ci.ub = res_sub$ci.ub,
                                 stringsAsFactors = FALSE
                               ))
    }
  }
}

quality_results <- quality_results %>%
  mutate(
    label = paste0(Quality, " (k=", k, ")"),
    Quality = factor(Quality, levels = c("High", "Medium", "Low"))
  ) %>%
  arrange(Quality)


p_quality_forest <- ggplot(quality_results, aes(x = Estimate, y = label, color = Quality)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
  geom_point(size = 3.5) +
  geom_errorbarh(aes(xmin = ci.lb, xmax = ci.ub), height = 0.2, size = 0.8) +
  scale_color_manual(values = c("High" = "#DB5829", "Medium" = "#5F6D7A", "Low" = "#E9DC6D")) +
  labs(x = "Log Response Ratio (95% CI)", y = "", title = "c. Quality Score") +
  theme_hetero_base +
  geom_vline(xintercept = mean(quality_results$Estimate), 
             linetype = "dotted", color = "red", size = 0.5)

#ggsave("Forest_Plot_Quality.png", p_quality_forest, width = 6, height = 3.5, dpi = 300)



# ============================================================================
# MASTER FIGURE WITH PROPER ALIGNMENT
# ============================================================================

library(patchwork)

master_figure <- 
  # Row 1: Two panels - Studies and Species
  (p_study_forest | p_species_forest) /
  
  # Row 2: Three panels - Quality, Method, Season
  (p_quality_forest | p_method_facet | p_season_facet) /
  
  # Row 3: Three panels - Environment, Catch, Country
  (p_env_forest | p_catch_forest | p_country_forest) +
  
  # Add layout specifications to ensure proper alignment
  plot_layout(heights = c(1.2, 1, 1)) +  # Row 1 slightly taller for better balance
  
  # Add annotation
  plot_annotation(
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30")
    )
  )

# Display
print(master_figure)

# Save
ggsave("Figure_Master_Heterogeneity_Aligned.png", master_figure, 
       width = 12, height = 10, dpi = 600, bg = "white")

ggsave("Figure_Master_Heterogeneity_Aligned.pdf", master_figure, 
       width = 12, height = 10, dpi = 600, bg = "white")


