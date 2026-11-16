leaf_soil_data <- read.xlsx("leaf_soil_data.xlsx", sheet = 1)
head(leaf_soil_data)
#delete wrong samples of Lotus corniculatus
leaf_soil_data <- leaf_soil_data %>%
  filter(!sample_ID %in% c("H3_LC", "H40_LC"))
#leaf_soil_data <- leaf_soil_data %>%
  #rename(
    #leaf_Si_mg_g = Si_mg_gr,
    #leaf_Ca_mg_g = ca_mg_per_g
#  )
unique(leaf_soil_data$group)#there is NA, and it belongs to Dactylis

# Replace NA values in the 'group' column with 'Grass'
leaf_soil_data <- leaf_soil_data %>%
  mutate(group = ifelse(is.na(group), "Grass", group))
unique(leaf_soil_data$group)#perfect! 
# 1. Overall dimensions and structure
str(leaf_soil_data)
dim(leaf_soil_data)

# 2. How many observations per species
table(leaf_soil_data$species, leaf_soil_data$group)

# 3. How many plots per region
table(leaf_soil_data$exploratorium, leaf_soil_data$plot_ID) %>% rowSums()
length(unique(leaf_soil_data$plot_ID))

# 4. Missing data per column
sapply(leaf_soil_data, function(x) sum(is.na(x)))

# 5. Summary stats for the key variables
summary(leaf_soil_data[, c("leaf_Si_mg_gr", "leaf_ca_mg_per_g", "SLA", 
                           "Soil_Si_mg_g", "Soil_Ca_mg_g")])

# 6. How many observations per plot
table(leaf_soil_data$plot_ID) %>% summary()
# =============================================================================
# Manuscript: Leaf Si and Ca accumulation across plant functional groups
# Analysis pipeline — clean version
# Author: Nahid R. Paeenroudposhti
# =============================================================================

# ---- 1.1 Load packages ----
# Data handling
library(openxlsx)
library(dplyr)
library(tidyr)

# Modeling
library(lme4)        # linear mixed-effects models
library(lmerTest)    # adds p-values to lmer output (Satterthwaite df)
library(car)         # Anova() function for Type II/III tests
library(MuMIn)       # r.squaredGLMM() for marginal/conditional R²
library(emmeans)     # post-hoc comparisons and estimated marginal means

# Visualization (we'll use these later, loading now for completeness)
library(ggplot2)
library(patchwork)
library(ggeffects)


# ---- 1.2 Load raw data ----
raw_data <- read.xlsx("leaf_soil_data.xlsx", sheet = 1)

cat("Raw data loaded:", nrow(raw_data), "observations,", 
    ncol(raw_data), "columns\n")


# ---- 1.3 Document and apply exclusions ----
# Two samples were excluded because they were mislabeled in the field 
# (the species recorded was not actually Lotus corniculatus).
# This was determined during quality control after sampling.
#(2) H13_TR: the leaf chemistry values (leaf Si = 77 mg/g, leaf Ca = 4.2 
#     mg/g) are inconsistent with Trifolium repens and all other legume 
#     samples (leaf Si: 1-7 mg/g, leaf Ca: 13-22 mg/g), indicating a 
#     likely sample swap with another species at the same plot. As the 
#     correct species identity could not be unambiguously determined, 
#     the sample was excluded from all analyses.

mislabeled_samples <- c("H3_LC", "H40_LC", "H13_TR")

cat("Excluding", length(mislabeled_samples), 
    "mislabeled samples:", paste(mislabeled_samples, collapse = ", "), "\n")

data <- raw_data %>%
  filter(!sample_ID %in% mislabeled_samples)


# ---- 1.4 Fix the missing functional group label ----
# One Dactylis glomerata sample has NA in the 'group' column.
# Dactylis glomerata is a grass, so we assign 'Grass'.

cat("Fixing functional group label for Dactylis glomerata...\n")

data <- data %>%
  mutate(group = ifelse(is.na(group) & species == "Dactylis glomerata", 
                        "Grass", 
                        group))


# ---- 1.5 Select and rename relevant columns only ----
# Remove pre-computed log columns (we will recompute) and Treatment column 
# (all values are "RP-SN" — uninformative).

data <- data %>%
  dplyr::select(sample_ID, 
         exploratorium, 
         plot_ID, 
         species, 
         group,
         leaf_Si  = leaf_Si_mg_gr,        # rename for clarity
         leaf_Ca  = leaf_ca_mg_per_g,
         SLA,
         soil_Si  = Soil_Si_mg_g,
         soil_Ca  = Soil_Ca_mg_g)


# ---- 1.6 Set factor levels explicitly ----
# This makes results reproducible across machines and ensures correct 
# reference levels in models (Forb as reference for group).

data <- data %>%
  mutate(
    group         = factor(group, levels = c("Forb", "Grass", "Legume")),
    exploratorium = factor(exploratorium, levels = c("AEG", "HEG", "SEG")),
    species       = factor(species),
    plot_ID       = factor(plot_ID)
  )


# ---- 1.7 Verify the cleaned dataset ----
cat("\n=== CLEANED DATA SUMMARY ===\n")
cat("Total observations:", nrow(data), "\n")
cat("Number of plots:", length(unique(data$plot_ID)), "\n")
cat("Number of species:", length(unique(data$species)), "\n")
cat("Number of regions:", length(unique(data$exploratorium)), "\n\n")

cat("Observations per functional group:\n")
print(table(data$group))

cat("\nObservations per region:\n")
print(table(data$exploratorium))

cat("\nMissing values per variable:\n")
print(sapply(data[, c("leaf_Si", "leaf_Ca", "SLA", "soil_Si", "soil_Ca")],
             function(x) sum(is.na(x))))

# Save cleaned data for later use
# write.xlsx(data, "cleaned_data.xlsx")  # uncomment if you want to save
# =============================================================================
# Section 2: Inspect distributions and decide on transformations
# =============================================================================

# ---- 2.1 Histograms and skewness of key variables ----

# Set up a panel of histograms for visual inspection
par(mfrow = c(2, 3))

hist(data$leaf_Si, main = "Leaf Si (mg/g)", xlab = "leaf Si", col = "lightgray")
hist(data$leaf_Ca, main = "Leaf Ca (mg/g)", xlab = "leaf Ca", col = "lightgray")
hist(data$SLA, main = "SLA", xlab = "SLA", col = "lightgray")
hist(data$soil_Si, main = "Soil Si (mg/g)", xlab = "soil Si", col = "lightgray")
hist(data$soil_Ca, main = "Soil Ca (mg/g)", xlab = "soil Ca", col = "lightgray")

par(mfrow = c(1, 1))  # reset


# ---- 2.2 Numerical skewness check ----
# As a rule of thumb:
#   |skewness| < 0.5     → approximately symmetric
#   0.5 < |skewness| < 1 → moderately skewed
#   |skewness| > 1       → highly skewed (transformation usually recommended)

library(e1071)  # for skewness function

skew_check <- data.frame(
  variable = c("leaf_Si", "leaf_Ca", "SLA", "soil_Si", "soil_Ca"),
  skewness = sapply(data[, c("leaf_Si", "leaf_Ca", "SLA", "soil_Si", "soil_Ca")], 
                    function(x) round(skewness(x, na.rm = TRUE), 2))
)

cat("Skewness of variables (raw scale):\n")
print(skew_check)


# ---- 2.3 Test log transformation for leaf_Si ----
# Leaf Si values range from 0 to 282 mg/g with strong right skew.
# We use log(x + 1) because some values are 0 (legumes can have near-zero Si).

cat("\nSkewness of log(leaf_Si + 1):", 
    round(skewness(log(data$leaf_Si + 1), na.rm = TRUE), 2), "\n")

par(mfrow = c(1, 2))
hist(data$leaf_Si, main = "Leaf Si raw", xlab = "leaf Si", col = "lightgray")
hist(log(data$leaf_Si + 1), main = "Leaf Si log(x+1)", 
     xlab = "log(leaf Si + 1)", col = "lightblue")
par(mfrow = c(1, 1))


# ---- 2.4 Test log transformation for leaf_Ca ----
# Leaf Ca is moderately skewed — check if log helps

cat("Skewness of log(leaf_Ca):", 
    round(skewness(log(data$leaf_Ca), na.rm = TRUE), 2), "\n")

par(mfrow = c(1, 2))
hist(data$leaf_Ca, main = "Leaf Ca raw", xlab = "leaf Ca", col = "lightgray")
hist(log(data$leaf_Ca), main = "Leaf Ca log", 
     xlab = "log(leaf Ca)", col = "lightblue")
par(mfrow = c(1, 1))


# ---- 2.5 Apply transformations ----
# Decision rules:
# - leaf_Si: log(x + 1) because it's strongly right-skewed and contains zeros
# - leaf_Ca: check skewness after log — apply only if it improves substantially
# - SLA, soil_Si, soil_Ca: keep raw unless they show strong skewness

data <- data %>%
  mutate(
    leaf_Si_log = log(leaf_Si + 1)
    # leaf_Ca_log = log(leaf_Ca)   # only add if skewness improves substantially
  )

cat("\nTransformations applied:\n")
cat("- leaf_Si → log(leaf_Si + 1) stored as 'leaf_Si_log'\n")
# ---- 2.6 Investigate soil_Ca extreme skewness ----

# View the distribution of unique soil Ca values per plot
soil_ca_check <- data %>%
  dplyr::select(plot_ID, exploratorium, soil_Ca) %>%
  distinct() %>%
  arrange(desc(soil_Ca))

print(head(soil_ca_check, 10))   # top 10 highest soil Ca values
print(summary(soil_ca_check$soil_Ca))

# Histogram of plot-level soil_Ca (one value per plot, not repeated)
hist(soil_ca_check$soil_Ca, 
     main = "Soil Ca per plot (n = 45)", 
     xlab = "Soil Ca (mg/g)",
     col = "lightgray",
     breaks = 20)

# Also test log transformation
cat("\nSkewness of log(soil_Ca) at plot level:", 
    round(skewness(log(soil_ca_check$soil_Ca), na.rm = TRUE), 2), "\n")
# Check the two suspicious plots in detail
suspect_plots <- data %>%
  filter(plot_ID %in% c("SEG14", "SEG3")) %>%
  dplyr::select(plot_ID, species, leaf_Ca, soil_Ca, soil_Si)

print(suspect_plots)

# Compare to other SEG plots
seg_plots <- data %>%
  filter(exploratorium == "SEG") %>%
  group_by(plot_ID) %>%
  summarise(soil_Ca = first(soil_Ca),
            soil_Si = first(soil_Si)) %>%
  arrange(desc(soil_Ca))

print(seg_plots)
# ---- 2.6 Quality control: implausibly high soil Ca in two SEG plots ----
# SEG14 and SEG3 show soil Ca values (54.1 and 41.4 mg/g) that are
# 6-7x higher than all other plots in the Schorfheide region, which is
# characterized by sandy, acidic, low-Ca soils. The corresponding leaf
# Ca values do not reflect calcareous soil conditions, and the soil_Si
# values for these plots are consistent with the regional geology.
# These soil_Ca measurements are therefore considered probable analytical
# errors and set to NA for the soil-leaf Ca analyses.

implausible_soil_Ca_plots <- c("SEG14", "SEG3")

data <- data %>%
  mutate(
    soil_Ca = ifelse(plot_ID %in% implausible_soil_Ca_plots, NA, soil_Ca)
  )

# Verify
cat("Soil Ca values set to NA for plots:", 
    paste(implausible_soil_Ca_plots, collapse = ", "), "\n")

cat("\nUpdated soil Ca summary:\n")
print(summary(data$soil_Ca))

cat("\nUpdated soil Ca skewness:", 
    round(skewness(data$soil_Ca, na.rm = TRUE), 2), "\n")
# =============================================================================
# Section 3: Fit all linear mixed-effects models
# =============================================================================
# Model structure for all models:
#   response ~ predictor (* group) + exploratorium + (1|plot_ID) + (1|species)
#
# This structure accounts for:
#   - Fixed effect of region (exploratorium, 3 levels)
#   - Plot-level non-independence (1|plot_ID)  [addresses pseudoreplication]
#   - Species-level non-independence (1|species)
# =============================================================================


# ---- 3.1 Functional group effects (M1, M2, M3) ----
# These test whether functional groups differ in leaf chemistry / SLA,
# controlling for  plot and species variation.

M1 <- lmer(leaf_Si_log ~ group  + 
             (1|plot_ID) + (1|species), 
           data = data)

M2 <- lmer(leaf_Ca ~ group  + 
             (1|plot_ID) + (1|species), 
           data = data)

M3 <- lmer(SLA ~ group  + 
             (1|plot_ID) + (1|species), 
           data = data)


# ---- 3.2 Si–Ca relationship within functional groups (M4) ----

M4 <- lmer(leaf_Si_log ~ leaf_Ca * group  + 
             (1|plot_ID)+(1|species), 
           data = data)


# ---- 3.3 SLA–chemistry relationships (M5, M6) ----

M5 <- lmer(SLA ~ leaf_Ca * group  + 
             (1|plot_ID) + (1|species), 
           data = data)

M6 <- lmer(SLA ~ leaf_Si_log * group  + 
             (1|plot_ID) + (1|species), 
           data = data)


# ---- 3.4 Soil–leaf relationships (M7, M8) ----

M7 <- lmer(leaf_Si_log ~ soil_Si * group  + 
             (1|plot_ID) + (1|species), 
           data = data)

M8 <- lmer(leaf_Ca ~ soil_Ca * group  + 
             (1|plot_ID) + (1|species), 
           data = data)


# ---- 3.5 Print summaries and ANOVAs for all models ----
# Use Type III ANOVA for models with interactions (M4-M8)
# Use Type II ANOVA for models without interactions (M1-M3)

cat("\n================ M1: leaf_Si ~ group ================\n")
print(summary(M1))
cat("\n--- ANOVA Type II ---\n")
print(Anova(M1, type = "II"))
cat("\n--- R² (marginal/conditional) ---\n")
print(r.squaredGLMM(M1))

cat("\n================ M2: leaf_Ca ~ group ================\n")
print(summary(M2))
cat("\n--- ANOVA Type II ---\n")
print(Anova(M2, type = "II"))
cat("\n--- R² (marginal/conditional) ---\n")
print(r.squaredGLMM(M2))

cat("\n================ M3: SLA ~ group ================\n")
print(summary(M3))
cat("\n--- ANOVA Type II ---\n")
print(Anova(M3, type = "II"))
cat("\n--- R² (marginal/conditional) ---\n")
print(r.squaredGLMM(M3))

cat("\n================ M4: leaf_Si ~ leaf_Ca * group ================\n")
print(summary(M4))
cat("\n--- ANOVA Type III ---\n")
print(Anova(M4, type = "III"))
cat("\n--- R² (marginal/conditional) ---\n")
print(r.squaredGLMM(M4))

cat("\n================ M5: SLA ~ leaf_Ca * group ================\n")
print(summary(M5))
cat("\n--- ANOVA Type III ---\n")
print(Anova(M5, type = "III"))
cat("\n--- R² (marginal/conditional) ---\n")
print(r.squaredGLMM(M5))

cat("\n================ M6: SLA ~ leaf_Si * group ================\n")
print(summary(M6))
cat("\n--- ANOVA Type III ---\n")
print(Anova(M6, type = "III"))
cat("\n--- R² (marginal/conditional) ---\n")
print(r.squaredGLMM(M6))

cat("\n================ M7: leaf_Si ~ soil_Si * group ================\n")
print(summary(M7))
cat("\n--- ANOVA Type III ---\n")
print(Anova(M7, type = "III"))
cat("\n--- R² (marginal/conditional) ---\n")
print(r.squaredGLMM(M7))

cat("\n================ M8: leaf_Ca ~ soil_Ca * group ================\n")
print(summary(M8))
cat("\n--- ANOVA Type III ---\n")
print(Anova(M8, type = "III"))
cat("\n--- R² (marginal/conditional) ---\n")
print(r.squaredGLMM(M8))
# =============================================================================
# M4 reduced: leaf_Si_log ~ leaf_Ca + group + (1|plot_ID) + (1|species)
# =============================================================================
# Full M4 had a non-significant Ca × group interaction (P = 0.271),
# so we drop the interaction and refit with additive main effects only.
# =============================================================================

# ---- Fit reduced model ----
M4_reduced <- lmer(leaf_Si_log ~ leaf_Ca + group +
                     (1|plot_ID) + (1|species),
                   data = data)

# ---- Print summary and Type II ANOVA (since no interaction, Type II is appropriate) ----
cat("\n========== M4 reduced: model summary ==========\n")
print(summary(M4_reduced))

cat("\n========== M4 reduced: Type II ANOVA ==========\n")
print(Anova(M4_reduced, type = "II"))

# ---- R² ----
cat("\n========== M4 reduced: R² ==========\n")
print(r.squaredGLMM(M4_reduced))

# ---- Pairwise group contrasts (for consistency with full-model reporting) ----
cat("\n========== M4 reduced: Tukey pairwise group comparisons ==========\n")
emm_M4red <- emmeans(M4_reduced, ~ group)
print(pairs(emm_M4red, adjust = "tukey"))
print(confint(pairs(emm_M4red, adjust = "tukey")))

# ---- Common Ca slope (single slope across groups) ----
cat("\n========== M4 reduced: common leaf_Ca slope ==========\n")
slope_M4red <- emtrends(M4_reduced, ~ 1, var = "leaf_Ca")
print(summary(slope_M4red, infer = c(TRUE, TRUE)))

# ---- Predictions for figure (three parallel lines) ----
library(ggeffects)
pred_M4red <- ggpredict(M4_reduced, terms = c("leaf_Ca", "group")) %>%
  as.data.frame()

cat("\n========== M4 reduced: prediction ranges per group ==========\n")
pred_M4red %>%
  group_by(group) %>%
  summarise(
    n_pred    = n(),
    min_pred  = min(predicted),
    max_pred  = max(predicted),
    .groups   = "drop"
  ) %>%
  print()
# =============================================================================
# M6 reduced: SLA ~ leaf_Si_log + group + (1|plot_ID) + (1|species)
# =============================================================================
# Full M6 had a non-significant Si × group interaction (P = 0.211),
# so we drop the interaction and refit with additive main effects only.
# =============================================================================

# ---- Fit reduced model ----
M6_reduced <- lmer(SLA ~ leaf_Si_log + group +
                     (1|plot_ID) + (1|species),
                   data = data)

# ---- Model summary and Type II ANOVA ----
cat("\n========== M6 reduced: model summary ==========\n")
print(summary(M6_reduced))

cat("\n========== M6 reduced: Type II ANOVA ==========\n")
print(Anova(M6_reduced, type = "II"))

cat("\n========== M6 reduced: R² ==========\n")
print(r.squaredGLMM(M6_reduced))

# ---- Common leaf Si slope ----
cat("\n========== M6 reduced: common leaf_Si slope ==========\n")
slope_M6red <- emtrends(M6_reduced, ~ 1, var = "leaf_Si_log")
print(summary(slope_M6red, infer = c(TRUE, TRUE)))

# ---- Predictions for figure ----
pred_M6red <- ggpredict(M6_reduced, terms = c("leaf_Si_log", "group")) %>%
  as.data.frame()

cat("\n========== M6 reduced: prediction summary ==========\n")
print(head(pred_M6red))
# =============================================================================
# M7 and M8 reduced: drop the non-significant interactions
# =============================================================================

# ---- Fit reduced models ----
M7_reduced <- lmer(leaf_Si_log ~ soil_Si + group +
                     (1|plot_ID) + (1|species),
                   data = data)

M8_reduced <- lmer(leaf_Ca ~ soil_Ca + group +
                     (1|plot_ID) + (1|species),
                   data = data)

# ---- M7 reduced ----
cat("\n========== M7 reduced: leaf_Si_log ~ soil_Si + group ==========\n")
print(summary(M7_reduced))

cat("\n========== M7 reduced: Type II ANOVA ==========\n")
print(Anova(M7_reduced, type = "II"))

cat("\n========== M7 reduced: R² ==========\n")
print(r.squaredGLMM(M7_reduced))

# ---- Common soil Si slope for M7 (for the figure/results text if needed) ----
cat("\n========== M7 reduced: common soil_Si slope ==========\n")
slope_M7red <- emtrends(M7_reduced, ~ 1, var = "soil_Si")
print(summary(slope_M7red, infer = c(TRUE, TRUE)))

# ---- M8 reduced ----
cat("\n========== M8 reduced: leaf_Ca ~ soil_Ca + group ==========\n")
print(summary(M8_reduced))

cat("\n========== M8 reduced: Type II ANOVA ==========\n")
print(Anova(M8_reduced, type = "II"))

cat("\n========== M8 reduced: R² ==========\n")
print(r.squaredGLMM(M8_reduced))

# ---- Common soil Ca slope for M8 ----
cat("\n========== M8 reduced: common soil_Ca slope ==========\n")
slope_M8red <- emtrends(M8_reduced, ~ 1, var = "soil_Ca")
print(summary(slope_M8red, infer = c(TRUE, TRUE)))

# Check if the soil_Ca effect was real or driven by region
# Run M8 WITHOUT region to see what the original model was capturing
M8_no_region <- lmer(leaf_Ca ~ soil_Ca * group + 
                       (1|plot_ID) + (1|species), 
                     data = data)
Anova(M8_no_region, type = "III")
# Within-region soil-leaf Ca relationships
data_aeg <- data %>% filter(exploratorium == "AEG")
data_heg <- data %>% filter(exploratorium == "HEG")
data_seg <- data %>% filter(exploratorium == "SEG")

M8_aeg <- lmer(leaf_Ca ~ soil_Ca * group + (1|plot_ID) + (1|species), data = data_aeg)
M8_heg <- lmer(leaf_Ca ~ soil_Ca * group + (1|plot_ID) + (1|species), data = data_heg)
M8_seg <- lmer(leaf_Ca ~ soil_Ca * group + (1|plot_ID) + (1|species), data = data_seg)

Anova(M8_aeg, type = "III")
Anova(M8_heg, type = "III")
Anova(M8_seg, type = "III")
# =============================================================================
# Section 4: Post-hoc pairwise comparisons (Tukey-adjusted)
# =============================================================================
# For models with significant group effects (M1, M2, M3), test which 
# functional groups differ from each other.

# ---- 4.1 Leaf Si by functional group (M1) ----
cat("\n========== M1: Pairwise comparisons of leaf Si ==========\n")
emm_M1 <- emmeans(M1, pairwise ~ group, adjust = "tukey")
print(emm_M1)

# Compact letter display (a, b, c notation for plots)
cat("\nCompact letter display:\n")
print(multcomp::cld(emm_M1$emmeans, Letters = letters, adjust = "tukey"))


# ---- 4.2 Leaf Ca by functional group (M2) ----
cat("\n========== M2: Pairwise comparisons of leaf Ca ==========\n")
emm_M2 <- emmeans(M2, pairwise ~ group, adjust = "tukey")
print(emm_M2)

cat("\nCompact letter display:\n")
print(multcomp::cld(emm_M2$emmeans, Letters = letters, adjust = "tukey"))


# ---- 4.3 SLA by functional group (M3) ----
cat("\n========== M3: Pairwise comparisons of SLA ==========\n")
emm_M3 <- emmeans(M3, pairwise ~ group, adjust = "tukey")
print(emm_M3)

cat("\nCompact letter display:\n")
print(multcomp::cld(emm_M3$emmeans, Letters = letters, adjust = "tukey"))
#install.packages("multcomp")
library(multcomp)
#install.packages("multcompView")

library(multcompView)

# Then rerun:
cat("\n========== Compact letter displays ==========\n")
cat("\nM1 (leaf Si):\n")
print(multcomp::cld(emm_M1$emmeans, Letters = letters, adjust = "tukey"))

cat("\nM2 (leaf Ca):\n")
print(multcomp::cld(emm_M2$emmeans, Letters = letters, adjust = "tukey"))

cat("\nM3 (SLA):\n")
print(multcomp::cld(emm_M3$emmeans, Letters = letters, adjust = "tukey"))
# =============================================================================
# Section 5: Pearson correlations — overall and within functional groups
# =============================================================================
# Reviewer #2 asked us to distinguish between:
# (a) the overall Si–Ca correlation across all observations
# (b) interspecific correlations (between-species)
# (c) within-functional-group correlations
# (d) within-species correlations

# ---- 5.1 Overall correlation (all observations) ----
cat("\n========== Overall correlation: log(leaf_Si+1) vs leaf_Ca ==========\n")
cor_overall <- cor.test(data$leaf_Si_log, data$leaf_Ca, method = "pearson")
print(cor_overall)


# ---- 5.2 Within-functional-group correlations ----
cat("\n========== Within-functional-group correlations ==========\n")
for (g in levels(data$group)) {
  d_sub <- data %>% filter(group == g)
  cat("\n--- Functional group:", g, "---\n")
  print(cor.test(d_sub$leaf_Si_log, d_sub$leaf_Ca, method = "pearson"))
}


# ---- 5.3 Between-species correlation (species means) ----
cat("\n========== Between-species correlation (species means) ==========\n")
species_means <- data %>%
  group_by(species, group) %>%
  summarise(
    mean_leaf_Si_log = mean(leaf_Si_log, na.rm = TRUE),
    mean_leaf_Ca     = mean(leaf_Ca, na.rm = TRUE),
    .groups = "drop"
  )

print(species_means)

cor_species <- cor.test(species_means$mean_leaf_Si_log, 
                        species_means$mean_leaf_Ca, 
                        method = "pearson")
print(cor_species)


# ---- 5.4 Within-species correlations ----
# For species with enough observations, test Si-Ca within each species
cat("\n========== Within-species correlations ==========\n")
for (sp in levels(data$species)) {
  d_sub <- data %>% 
    filter(species == sp) %>%
    filter(!is.na(leaf_Si_log) & !is.na(leaf_Ca))
  
  if (nrow(d_sub) >= 10) {  # Only test species with enough data
    ct <- cor.test(d_sub$leaf_Si_log, d_sub$leaf_Ca, method = "pearson")
    cat(sprintf("%-25s n=%2d  r=%6.3f  p=%.3f\n", 
                sp, nrow(d_sub), ct$estimate, ct$p.value))
  }
}
# =============================================================================
# Section 6: Group-specific slopes for M5 (SLA ~ leaf_Ca * group)
# =============================================================================
# Reviewer #2 asked for slopes, SE, 95% CI, and p-values for each functional 
# group's Ca-SLA relationship.

cat("\n========== M5: Group-specific Ca-SLA slopes ==========\n")
slopes_M5 <- emtrends(M5, ~ group, var = "leaf_Ca")
print(summary(slopes_M5, infer = c(TRUE, TRUE)))


# Also for M4 (Si-Ca relationship) — though interaction was NS
cat("\n========== M4: Group-specific Ca-Si slopes ==========\n")
slopes_M4 <- emtrends(M4, ~ group, var = "leaf_Ca")
print(summary(slopes_M4, infer = c(TRUE, TRUE)))


# Also for M6 (Si-SLA relationship) — though interaction was NS
cat("\n========== M6: Group-specific Si-SLA slopes ==========\n")
slopes_M6 <- emtrends(M6, ~ group, var = "leaf_Si_log")
print(summary(slopes_M6, infer = c(TRUE, TRUE)))


# Also for M7 and M8 (soil-leaf relationships)
cat("\n========== M7: Group-specific Soil Si-Leaf Si slopes ==========\n")
slopes_M7 <- emtrends(M7, ~ group, var = "soil_Si")
print(summary(slopes_M7, infer = c(TRUE, TRUE)))

cat("\n========== M8: Group-specific Soil Ca-Leaf Ca slopes ==========\n")
slopes_M8 <- emtrends(M8, ~ group, var = "soil_Ca")
print(summary(slopes_M8, infer = c(TRUE, TRUE)))
# =============================================================================
# Section 7: Sensitivity analysis — leave-one-species-out (Achillea millefolium)
# =============================================================================
# Reviewer #2 raised the concern that the high Si concentrations in forbs are 
# largely driven by Achillea millefolium. To test the robustness of our 
# findings, we refit all 8 models excluding A. millefolium and compare the 
# key statistics with the full-dataset models.
#
# IMPORTANT: For M4, M6, M7, and M8, we refit the REDUCED (additive) models 
# rather than the original full-interaction versions, because those are our 
# primary analysis specifications (interactions were not significant in the 
# full dataset). M5 keeps its full interaction structure (interaction was 
# significant). M1, M2, M3 have no interaction and are unchanged.
# =============================================================================

# ---- 7.1 Create dataset without Achillea millefolium ----
data_noAchillea <- data %>%
  filter(species != "Achillea millefolium") %>%
  droplevels()

cat("\nFull dataset:    ", nrow(data), "observations,", 
    length(unique(data$species)), "species\n")
cat("Without Achillea:", nrow(data_noAchillea), "observations,", 
    length(unique(data_noAchillea$species)), "species\n")


# ---- 7.2 Refit models without Achillea, matching FINAL specifications ----

# M1, M2, M3 — no interaction, unchanged
M1_noA <- lmer(leaf_Si_log ~ group + (1|plot_ID) + (1|species), 
               data = data_noAchillea)
M2_noA <- lmer(leaf_Ca ~ group + (1|plot_ID) + (1|species), 
               data = data_noAchillea)
M3_noA <- lmer(SLA ~ group + (1|plot_ID) + (1|species), 
               data = data_noAchillea)

# M4, M6, M7, M8 — REDUCED (additive) models
M4_noA <- lmer(leaf_Si_log ~ leaf_Ca + group + (1|plot_ID) + (1|species), 
               data = data_noAchillea)
M6_noA <- lmer(SLA ~ leaf_Si_log + group + (1|plot_ID) + (1|species), 
               data = data_noAchillea)
M7_noA <- lmer(leaf_Si_log ~ soil_Si + group + (1|plot_ID) + (1|species), 
               data = data_noAchillea)
M8_noA <- lmer(leaf_Ca ~ soil_Ca + group + (1|plot_ID) + (1|species), 
               data = data_noAchillea)

# M5 — FULL interaction model retained (interaction was significant)
M5_noA <- lmer(SLA ~ leaf_Ca * group + (1|plot_ID) + (1|species), 
               data = data_noAchillea)


# ---- 7.3 Print ANOVA results for each model ----
# M1-M4, M6-M8: Type II (no interaction)
# M5: Type III (has interaction)

cat("\n========== M1 (leaf_Si ~ group) — WITHOUT Achillea ==========\n")
print(Anova(M1_noA, type = "II"))
cat("R²:\n"); print(r.squaredGLMM(M1_noA))

cat("\n========== M2 (leaf_Ca ~ group) — WITHOUT Achillea ==========\n")
print(Anova(M2_noA, type = "II"))
cat("R²:\n"); print(r.squaredGLMM(M2_noA))

cat("\n========== M3 (SLA ~ group) — WITHOUT Achillea ==========\n")
print(Anova(M3_noA, type = "II"))
cat("R²:\n"); print(r.squaredGLMM(M3_noA))

cat("\n========== M4 REDUCED (leaf_Si ~ leaf_Ca + group) — WITHOUT Achillea ==========\n")
print(Anova(M4_noA, type = "II"))
cat("R²:\n"); print(r.squaredGLMM(M4_noA))

cat("\n========== M5 (SLA ~ leaf_Ca * group) — WITHOUT Achillea ==========\n")
print(Anova(M5_noA, type = "III"))
cat("R²:\n"); print(r.squaredGLMM(M5_noA))

cat("\n========== M6 REDUCED (SLA ~ leaf_Si_log + group) — WITHOUT Achillea ==========\n")
print(Anova(M6_noA, type = "II"))
cat("R²:\n"); print(r.squaredGLMM(M6_noA))

cat("\n========== M7 REDUCED (leaf_Si ~ soil_Si + group) — WITHOUT Achillea ==========\n")
print(Anova(M7_noA, type = "II"))
cat("R²:\n"); print(r.squaredGLMM(M7_noA))

cat("\n========== M8 REDUCED (leaf_Ca ~ soil_Ca + group) — WITHOUT Achillea ==========\n")
print(Anova(M8_noA, type = "II"))
cat("R²:\n"); print(r.squaredGLMM(M8_noA))


# ---- 7.4 Side-by-side comparison table (full vs. without Achillea) ----

cat("\n========== COMPARISON: full dataset vs without Achillea ==========\n\n")

# Helper function to extract p-values from Anova output
get_p <- function(model, type) {
  a <- car::Anova(model, type = type)
  setNames(a$`Pr(>Chisq)`, rownames(a))
}

# M1: leaf_Si ~ group
p1_full <- get_p(M1, "II")
p1_noA  <- get_p(M1_noA, "II")
cat("M1: leaf_Si ~ group\n")
cat(sprintf("  group:  full p = %.4f | without Achillea p = %.4f\n\n", 
            p1_full["group"], p1_noA["group"]))

# M2: leaf_Ca ~ group
p2_full <- get_p(M2, "II")
p2_noA  <- get_p(M2_noA, "II")
cat("M2: leaf_Ca ~ group\n")
cat(sprintf("  group:  full p = %.4f | without Achillea p = %.4f\n\n", 
            p2_full["group"], p2_noA["group"]))

# M3: SLA ~ group
p3_full <- get_p(M3, "II")
p3_noA  <- get_p(M3_noA, "II")
cat("M3: SLA ~ group\n")
cat(sprintf("  group:  full p = %.4f | without Achillea p = %.4f\n\n", 
            p3_full["group"], p3_noA["group"]))

# M4 reduced: leaf_Si ~ leaf_Ca + group
# Note: compare to M4_reduced (full dataset) not M4 (full interaction)
p4_full <- get_p(M4_reduced, "II")   # M4_reduced from primary analysis
p4_noA  <- get_p(M4_noA,     "II")
cat("M4 REDUCED: leaf_Si ~ leaf_Ca + group\n")
cat(sprintf("  leaf_Ca: full p = %.4f | without Achillea p = %.4f\n", 
            p4_full["leaf_Ca"], p4_noA["leaf_Ca"]))
cat(sprintf("  group:   full p = %.4f | without Achillea p = %.4f\n\n", 
            p4_full["group"], p4_noA["group"]))

# M5: SLA ~ leaf_Ca * group (full interaction, unchanged)
p5_full <- get_p(M5, "III")
p5_noA  <- get_p(M5_noA, "III")
cat("M5: SLA ~ leaf_Ca * group  (full interaction)\n")
cat(sprintf("  leaf_Ca:        full p = %.4f | without Achillea p = %.4f\n", 
            p5_full["leaf_Ca"], p5_noA["leaf_Ca"]))
cat(sprintf("  group:          full p = %.4f | without Achillea p = %.4f\n", 
            p5_full["group"], p5_noA["group"]))
cat(sprintf("  leaf_Ca:group:  full p = %.4f | without Achillea p = %.4f\n\n", 
            p5_full["leaf_Ca:group"], p5_noA["leaf_Ca:group"]))

# M6 reduced: SLA ~ leaf_Si + group
p6_full <- get_p(M6_reduced, "II")   # M6_reduced from primary analysis
p6_noA  <- get_p(M6_noA,     "II")
cat("M6 REDUCED: SLA ~ leaf_Si_log + group\n")
cat(sprintf("  leaf_Si_log: full p = %.4f | without Achillea p = %.4f\n", 
            p6_full["leaf_Si_log"], p6_noA["leaf_Si_log"]))
cat(sprintf("  group:       full p = %.4f | without Achillea p = %.4f\n\n", 
            p6_full["group"], p6_noA["group"]))

# M7 reduced: leaf_Si ~ soil_Si + group
p7_full <- get_p(M7_reduced, "II")   # M7_reduced from primary analysis
p7_noA  <- get_p(M7_noA,     "II")
cat("M7 REDUCED: leaf_Si ~ soil_Si + group\n")
cat(sprintf("  soil_Si: full p = %.4f | without Achillea p = %.4f\n", 
            p7_full["soil_Si"], p7_noA["soil_Si"]))
cat(sprintf("  group:   full p = %.4f | without Achillea p = %.4f\n\n", 
            p7_full["group"], p7_noA["group"]))

# M8 reduced: leaf_Ca ~ soil_Ca + group
p8_full <- get_p(M8_reduced, "II")   # M8_reduced from primary analysis
p8_noA  <- get_p(M8_noA,     "II")
cat("M8 REDUCED: leaf_Ca ~ soil_Ca + group\n")
cat(sprintf("  soil_Ca: full p = %.4f | without Achillea p = %.4f\n", 
            p8_full["soil_Ca"], p8_noA["soil_Ca"]))
cat(sprintf("  group:   full p = %.4f | without Achillea p = %.4f\n\n", 
            p8_full["group"], p8_noA["group"]))
# =============================================================================
# Section 8: Model diagnostics
# =============================================================================
# Reviewer #2 asked for: "rationale for log-transforming leaf Si, residual 
# normality and homoscedasticity, exact sample sizes for each model, random-
# effect variance components, possible singular fits"
# =============================================================================


# ---- 8.1 Diagnostic plots for all 8 models ----
# For each model, we check:
#   (1) Residuals vs Fitted (look for homoscedasticity — random scatter)
#   (2) Q-Q plot (look for normality — points along the diagonal)

# Set up panel for 4 rows x 4 cols (2 plots per model x 8 models)
par(mfrow = c(4, 4), mar = c(3.5, 3.5, 2, 1), mgp = c(2, 0.7, 0))
library(e1071)
models_list <- list(
  M1 = M1, M2 = M2, M3 = M3, M4 = M4,
  M5 = M5, M6 = M6, M7 = M7, M8 = M8
)

for (name in names(models_list)) {
  m <- models_list[[name]]
  # Residuals vs Fitted
  plot(fitted(m), resid(m),
       xlab = "Fitted values", ylab = "Residuals",
       main = paste(name, "Residuals vs Fitted"),
       pch = 16, col = rgb(0, 0, 0, 0.3))
  abline(h = 0, col = "red", lty = 2)
  
  # Q-Q plot
  qqnorm(resid(m), main = paste(name, "Q-Q plot"),
         pch = 16, col = rgb(0, 0, 0, 0.3))
  qqline(resid(m), col = "red")
}

n <- length(models_list)

tiff("model_diagnostics.tiff", 
     width = 8, height = 3 * n, units = "in", res = 150, compression = "lzw")
par(mfrow = c(n, 2), mar = c(4, 4, 2, 1))

for (name in names(models_list)) {
  m <- models_list[[name]]
  
  plot(fitted(m), resid(m),
       xlab = "Fitted values", ylab = "Residuals",
       main = paste(name, "Residuals vs Fitted"),
       pch = 16, col = rgb(0, 0, 0, 0.3))
  abline(h = 0, col = "red", lty = 2)
  
  qqnorm(resid(m), main = paste(name, "Q-Q plot"),
         pch = 16, col = rgb(0, 0, 0, 0.3))
  qqline(resid(m), col = "red")
}

dev.off()
par(mfrow = c(1, 1))  # reset


# ---- 8.2 Numerical summary of residuals ----
# Skewness and kurtosis tell us about normality
# |skewness| < 1 is acceptable; |kurtosis - 3| < 2 is acceptable

cat("\n========== Residual diagnostics summary ==========\n")
cat(sprintf("%-5s %-12s %-12s %-12s\n", 
            "Model", "Skewness", "Kurtosis", "Singular?"))
cat(strrep("-", 50), "\n")

for (name in names(models_list)) {
  m <- models_list[[name]]
  r <- resid(m)
  sk <- round(skewness(r), 2)
  kr <- round(kurtosis(r), 2)
  sing <- ifelse(isSingular(m), "YES", "no")
  cat(sprintf("%-5s %-12s %-12s %-12s\n", name, sk, kr, sing))
}


# ---- 8.3 Sample sizes for each model ----
cat("\n========== Sample sizes per model ==========\n")
cat(sprintf("%-5s %-15s %-15s %-15s\n", 
            "Model", "Observations", "Plots", "Species"))
cat(strrep("-", 55), "\n")

for (name in names(models_list)) {
  m <- models_list[[name]]
  n_obs <- nobs(m)
  n_plots <- length(unique(m@frame$plot_ID))
  n_sp    <- length(unique(m@frame$species))
  cat(sprintf("%-5s %-15s %-15s %-15s\n", name, n_obs, n_plots, n_sp))
}


# ---- 8.4 Random effect variance components ----
cat("\n========== Random effect variance components ==========\n")
cat(sprintf("%-5s %-15s %-15s %-15s %-15s\n",
            "Model", "Var(plot_ID)", "Var(species)", "Var(Residual)", 
            "Species/Total"))
cat(strrep("-", 70), "\n")

for (name in names(models_list)) {
  m <- models_list[[name]]
  vc <- as.data.frame(VarCorr(m))
  v_plot <- round(vc$vcov[vc$grp == "plot_ID"], 3)
  v_sp   <- round(vc$vcov[vc$grp == "species"], 3)
  v_res  <- round(vc$vcov[vc$grp == "Residual"], 3)
  total  <- v_plot + v_sp + v_res
  prop_sp <- round(v_sp / total, 2)
  cat(sprintf("%-5s %-15s %-15s %-15s %-15s\n",
              name, v_plot, v_sp, v_res, prop_sp))
}
# =============================================================================
# Section 9: Species summary table (Table S1)
# =============================================================================
# Mean ± SD per species for leaf Si, leaf Ca, and SLA
# Reviewer #1 requested raw Si values for reference.

species_summary <- data %>%
  group_by(species, group) %>%
  summarise(
    n = n(),
    
    # Leaf Si — raw mg/g
    leaf_Si_mean = mean(leaf_Si, na.rm = TRUE),
    leaf_Si_sd   = sd(leaf_Si, na.rm = TRUE),
    
    # Leaf Ca — raw mg/g
    leaf_Ca_mean = mean(leaf_Ca, na.rm = TRUE),
    leaf_Ca_sd   = sd(leaf_Ca, na.rm = TRUE),
    
    # SLA
    SLA_mean = mean(SLA, na.rm = TRUE),
    SLA_sd   = sd(SLA, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(group, species)

cat("\n========== Species summary table (raw values) ==========\n")
print(species_summary, n = Inf)


# ---- 9.1 Publication-formatted version (mean ± SD) ----
species_table_pub <- data %>%
  group_by(species, group) %>%
  summarise(
    n = n(),
    leaf_Si = sprintf("%.1f ± %.1f", 
                      mean(leaf_Si, na.rm = TRUE), 
                      sd(leaf_Si, na.rm = TRUE)),
    leaf_Ca = sprintf("%.2f ± %.2f", 
                      mean(leaf_Ca, na.rm = TRUE), 
                      sd(leaf_Ca, na.rm = TRUE)),
    SLA     = sprintf("%.2f ± %.2f", 
                      mean(SLA, na.rm = TRUE), 
                      sd(SLA, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(group, species)

cat("\n========== Publication-formatted species table ==========\n")
print(species_table_pub, n = Inf)


# ---- 9.2 Save to Excel file for the supplementary material ----
write.xlsx(species_table_pub, "Table_S1_species_summary.xlsx")
cat("\nSaved: Table_S1_species_summary.xlsx\n")


# ---- 9.3 Group-level summary (means across species within each group) ----
group_summary <- data %>%
  group_by(group) %>%
  summarise(
    n = n(),
    leaf_Si_mean = round(mean(leaf_Si, na.rm = TRUE), 1),
    leaf_Si_sd   = round(sd(leaf_Si, na.rm = TRUE), 1),
    leaf_Ca_mean = round(mean(leaf_Ca, na.rm = TRUE), 2),
    leaf_Ca_sd   = round(sd(leaf_Ca, na.rm = TRUE), 2),
    SLA_mean     = round(mean(SLA, na.rm = TRUE), 2),
    SLA_sd       = round(sd(SLA, na.rm = TRUE), 2),
    .groups = "drop"
  )

cat("\n========== Group-level summary ==========\n")
print(group_summary)
data %>%
  filter(species == "Trifolium repens") %>%
  dplyr::select(plot_ID, leaf_Si, leaf_Ca, SLA) %>%
  arrange(desc(leaf_Si))
# Check what the data looks like for HEG13 (other species in that plot)
data %>%
  filter(plot_ID == "HEG13") %>%
  dplyr::select(sample_ID, species, group, leaf_Si, leaf_Ca, SLA)
# =============================================================================
# Section 10: Figures
# =============================================================================

# ---- 10.1 Common setup ----

# Okabe-Ito colorblind-friendly palette
group_colors <- c("Forb"   = "#E69F00",   # orange
                  "Grass"  = "#009E73",   # green
                  "Legume" = "#56B4E9")   # blue

# Custom minimal theme for all figures (consistent style)
my_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    plot.title       = element_blank(),
    axis.title       = element_text(family = "sans", size = 10, face = "plain"),
    axis.text        = element_text(family = "sans", size = 9),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_text(family = "sans", size = 9, face = "plain"),
    legend.text      = element_text(family = "sans", size = 9),
    legend.key.size  = unit(0.5, "cm")
  )


# =============================================================================
# Figure 1a — Leaf Si by functional group
# =============================================================================

# Letters from Tukey post-hoc (based on M1 results we ran earlier)
letters_si <- data.frame(
  group = factor(c("Forb", "Grass", "Legume"),
                 levels = c("Forb", "Grass", "Legume")),
  label = c("a", "b", "a"),
  y     = max(data$leaf_Si_log, na.rm = TRUE) + 0.3
)


# Define breaks at reasonable raw values, converted to log-scale positions
raw_breaks_Si <- c(0, 1, 5, 20, 100, 250)
log_breaks_Si <- log(raw_breaks_Si + 1)

Fig1a <- ggplot(data %>% filter(!is.na(leaf_Si_log)),
                aes(x = group, y = leaf_Si_log, fill = group)) +
  
  geom_violin(alpha = 0.8, color = NA, trim = FALSE) +
  
  geom_jitter(aes(color = group),
              width = 0.15, size = 1.8, alpha = 0.6, show.legend = FALSE) +
  
  # Median bar
  stat_summary(fun = median, geom = "point",
               shape = 95, size = 6, color = "black") +
  
  # Tukey letters
  geom_text(data = letters_si,
            aes(x = group, y = y, label = label),
            inherit.aes = FALSE,
            size = 4, fontface = "bold") +
  
  scale_y_continuous(breaks = log_breaks_Si,
                     labels = raw_breaks_Si) +
  
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  
  labs(x = "Functional group",
       y = expression("Leaf Si (mg g"^{-1}*")")) +
  
  my_theme +
  theme(legend.position = "none")
Fig1a

# =============================================================================
# Figure 1b — Leaf Ca by functional group
# =============================================================================

# Letters from Tukey post-hoc (based on M2 results)
letters_ca <- data.frame(
  group = factor(c("Forb", "Grass", "Legume"),
                 levels = c("Forb", "Grass", "Legume")),
  label = c("b", "a", "b"),
  y     = max(data$leaf_Ca, na.rm = TRUE) + 2
)

Fig1b <- ggplot(data %>% filter(!is.na(leaf_Ca)),
                aes(x = group, y = leaf_Ca, fill = group)) +
  
  geom_violin(alpha = 0.8, color = NA, trim = FALSE) +
  
  geom_jitter(aes(color = group),
              width = 0.15, size = 1.8, alpha = 0.6, show.legend = FALSE) +
  
  stat_summary(fun = median, geom = "point",
               shape = 95, size = 6, color = "black") +
  
  geom_text(data = letters_ca,
            aes(x = group, y = y, label = label),
            inherit.aes = FALSE,
            size = 4, fontface = "bold") +
  
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  
  labs(x = "Functional group",
       y = expression("Leaf Ca (mg g"^{-1}*")")) +
  
  my_theme +
  theme(legend.position = "none")


# =============================================================================
# Combine Fig 1a + Fig 1b
# =============================================================================

Figure1 <- Fig1a + Fig1b +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = "sans", size = 12, face = "bold"),
        plot.tag.position = c(0.02, 0.98))

# Display
Figure1

# Save as TIFF (publication quality)
ggsave("Figure1.tiff",
       plot = Figure1,
       width = 168,   # 84 + 84 mm (two single-column panels)
       height = 84,
       units = "mm",
       dpi = 300,
       compression = "lzw",
       bg = "white")
# ============================================================================

# =============================================================================
# Figure 2 — Soil Si → Leaf Si (a) and Soil Ca → Leaf Ca (b) 
# Comparison: intrinsic species/group trait vs. environmentally-tracked trait
# =============================================================================

# ----- Panel A: Soil Si → Leaf Si (based on M7) -----------------------------

# Group-specific slopes from M7 emtrends
slope_M7 <- emtrends(M7, ~ group, var = "soil_Si") %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as.data.frame()

slope_M7$linetype <- ifelse(slope_M7$p.value < 0.05, "solid", "dashed")
slope_M7$group    <- as.character(slope_M7$group)

print(slope_M7)  # Sanity check

# Predictions
pred_M7 <- ggpredict(M7, terms = c("soil_Si", "group")) %>%
  as.data.frame()

pred_M7$group <- as.character(pred_M7$group)
pred_M7 <- pred_M7 %>%
  left_join(slope_M7[, c("group", "linetype")], by = "group")
pred_M7$group <- factor(pred_M7$group, levels = c("Forb", "Grass", "Legume"))

# Define breaks at reasonable raw values, converted to log-scale positions
raw_breaks <- c(0, 1, 5, 20, 100, 250)
log_breaks <- log(raw_breaks + 1)   # this is where the ticks go

Fig2a <- ggplot(data %>% filter(!is.na(leaf_Si_log) & !is.na(soil_Si)),
                aes(x = soil_Si, y = leaf_Si_log)) +
  geom_point(aes(color = group),
             alpha = 0.5, size = 2, na.rm = TRUE) +
  geom_ribbon(data = pred_M7,
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.15,
              inherit.aes = FALSE) +
  geom_line(data = pred_M7,
            aes(x = x, y = predicted, color = group, linetype = linetype),
            linewidth = 1.1,
            inherit.aes = FALSE) +
  scale_y_continuous(
    breaks = log_breaks,                          # tick positions on log scale
    labels = raw_breaks                           # labels showing raw values
  ) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_linetype_identity() +
  labs(x = expression("Soil Si (mg g"^{-1}*")"),
       y = expression("Leaf Si (mg g"^{-1}*")"),   # raw units!
       color = "Functional group",
       fill  = "Functional group") +
  my_theme
Fig2a
# ----- Panel B: Soil Ca → Leaf Ca (based on M8) -----------------------------

slope_M8 <- emtrends(M8, ~ group, var = "soil_Ca") %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as.data.frame()

slope_M8$linetype <- ifelse(slope_M8$p.value < 0.05, "solid", "dashed")
slope_M8$group    <- as.character(slope_M8$group)

print(slope_M8)  # Sanity check

pred_M8 <- ggpredict(M8, terms = c("soil_Ca", "group")) %>%
  as.data.frame()

pred_M8$group <- as.character(pred_M8$group)
pred_M8 <- pred_M8 %>%
  left_join(slope_M8[, c("group", "linetype")], by = "group")
pred_M8$group <- factor(pred_M8$group, levels = c("Forb", "Grass", "Legume"))

Fig2b <- ggplot(data %>% filter(!is.na(leaf_Ca) & !is.na(soil_Ca)),
                aes(x = soil_Ca, y = leaf_Ca)) +
  geom_point(aes(color = group),
             alpha = 0.5, size = 2, na.rm = TRUE) +
  geom_ribbon(data = pred_M8,
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.15,
              inherit.aes = FALSE) +
  geom_line(data = pred_M8,
            aes(x = x, y = predicted, color = group, linetype = linetype),
            linewidth = 1.1,
            inherit.aes = FALSE) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_linetype_identity() +
  labs(x = expression("Soil Ca (mg g"^{-1}*")"),   # ← confirm unit
       y = expression("Leaf Ca (mg g"^{-1}*")"),
       color = "Functional group",
       fill  = "Functional group") +
  my_theme
Fig2b
# ----- Combine panels ------------------------------------------------------

Figure2 <- Fig2a + Fig2b +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = "sans", size = 12, face = "bold"),
        plot.tag.position = c(0.02, 0.98),
        legend.position = "bottom")

Figure2

# ----- Save ----------------------------------------------------------------

ggsave("Figure2.tiff",
       plot = Figure2,
       width = 168, height = 90, units = "mm",
       dpi = 300, compression = "lzw", bg = "white")

# =============================================================================
# =============================================================================
# Figure S1 — Species-level leaf Si and Ca by functional group
# =============================================================================
# Reviewer #1 requested species-level breakdown of leaf traits.
# We're keeping this in supplementary as a species-level companion to Fig. 1.
# =============================================================================

# Reorder species by functional group, then alphabetically within group
data_for_fig <- data %>%
  mutate(group = factor(group, levels = c("Forb", "Grass", "Legume"))) %>%
  arrange(group, species) %>%
  mutate(species = factor(species, levels = unique(species)))

# ---- Define log-axis breaks for Si (raw values shown on log-scaled axis) ----
raw_breaks_Si <- c(0, 1, 5, 20, 100, 250)
log_breaks_Si <- log(raw_breaks_Si + 1)

# ---- Panel a: Leaf Si by species ----
FigS1a <- ggplot(data_for_fig %>% filter(!is.na(leaf_Si_log)),
                 aes(x = species, y = leaf_Si_log, fill = group)) +
  
  geom_violin(alpha = 0.8, color = NA, trim = TRUE) +
  
  geom_jitter(aes(color = group),
              width = 0.15, size = 1.5, alpha = 0.5,
              show.legend = FALSE) +
  
  stat_summary(fun = median, geom = "crossbar",
               width = 0.5, fatten = 2, color = "black") +
  
  scale_y_continuous(breaks = log_breaks_Si,
                     labels = raw_breaks_Si) +
  
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  
  labs(x = NULL,
       y = expression("Leaf Si (mg g"^{-1}*")"),
       fill = "Functional group") +
  
  my_theme +
  theme(
    axis.text.x = element_text(family = "sans", size = 9,
                               angle = 45, hjust = 1, face = "italic"),
    legend.position = "none"
  )

# ---- Panel b: Leaf Ca by species ----
FigS1b <- ggplot(data_for_fig %>% filter(!is.na(leaf_Ca)),
                 aes(x = species, y = leaf_Ca, fill = group)) +
  
  geom_violin(alpha = 0.8, color = NA, trim = TRUE) +
  
  geom_jitter(aes(color = group),
              width = 0.15, size = 1.5, alpha = 0.5,
              show.legend = FALSE) +
  
  stat_summary(fun = median, geom = "crossbar",
               width = 0.5, fatten = 2, color = "black") +
  
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  
  labs(x = "Species",
       y = expression("Leaf Ca (mg g"^{-1}*")"),
       fill = "Functional group") +
  
  my_theme +
  theme(
    axis.text.x = element_text(family = "sans", size = 9,
                               angle = 45, hjust = 1, face = "italic")
  )

# =============================================================================
# Combine vertically so x-axis labels align
# =============================================================================
FigureS1 <- FigS1a / FigS1b +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = "sans", size = 12, face = "bold"),
        plot.tag.position = c(0.02, 0.98))

FigureS1

# Save
ggsave("FigureS1.tiff",
       plot = FigureS1,
       width = 174,
       height = 200,
       units = "mm",
       dpi = 300,
       compression = "lzw",
       bg = "white")
######################################################################################
# =============================================================================
# Figure S2 — SLA by plant functional group (based on M3)
# =============================================================================

letters_sla <- data.frame(
  group = factor(c("Forb", "Grass", "Legume"),
                 levels = c("Forb", "Grass", "Legume")),
  label = c("a", "a", "a"),
  y     = max(data$SLA, na.rm = TRUE) + 2
)

FigS2 <- ggplot(data %>% filter(!is.na(SLA)),
                aes(x = group, y = SLA, fill = group)) +
  
  geom_violin(alpha = 0.8, color = NA, trim = TRUE) +
  
  geom_jitter(aes(color = group),
              width = 0.15, size = 1.8, alpha = 0.6, show.legend = FALSE) +
  
  stat_summary(fun = median, geom = "crossbar",
               width = 0.4, fatten = 2, color = "black") +
  
  geom_text(data = letters_sla,
            aes(x = group, y = y, label = label),
            inherit.aes = FALSE,
            size = 4, fontface = "bold") +
  
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  
  labs(x = "Functional group",
       y = expression("SLA (mm"^2*" mg"^{-1}*")")) +   # ← updated
  
  my_theme +
  theme(legend.position = "none")

FigS2

ggsave("FigureS2.tiff",
       plot = FigS2,
       width = 84,
       height = 90,
       units = "mm",
       dpi = 300,
       compression = "lzw",
       bg = "white")

# =============================================================================
# Fig. S3 — Leaf Ca vs. Leaf Si (based on M4 reduced)
# =============================================================================
# M4 reduced: leaf_Si_log ~ leaf_Ca + group (additive, no interaction)
# The full-model interaction (Ca × group) was not significant (P = 0.271),
# so the model was simplified. Three lines share a common Ca slope
# (-0.011, P = 0.14) but differ in intercept by group.
# All three lines are dashed because the common Ca slope is non-significant.
# =============================================================================

# Define log-axis breaks for Si (raw values on log-scaled axis)
raw_breaks_Si <- c(0, 1, 5, 20, 100, 250)
log_breaks_Si <- log(raw_breaks_Si + 1)

# Predictions from reduced model
pred_M4red <- ggpredict(M4_reduced, terms = c("leaf_Ca", "group")) %>%
  as.data.frame()

pred_M4red$group <- factor(pred_M4red$group,
                           levels = c("Forb", "Grass", "Legume"))

# Common Ca slope is non-significant (P = 0.14), so all three lines are dashed
common_line_type <- "dashed"

FigS3 <- ggplot(data %>% filter(!is.na(leaf_Si_log) & !is.na(leaf_Ca)),
                aes(x = leaf_Ca, y = leaf_Si_log)) +
  geom_point(aes(color = group),
             alpha = 0.5, size = 2, na.rm = TRUE) +
  geom_ribbon(data = pred_M4red,
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.15,
              inherit.aes = FALSE) +
  geom_line(data = pred_M4red,
            aes(x = x, y = predicted, color = group),
            linetype = common_line_type,
            linewidth = 1.1,
            inherit.aes = FALSE) +
  scale_y_continuous(breaks = log_breaks_Si,
                     labels = raw_breaks_Si) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(x = expression("Leaf Ca (mg g"^{-1}*")"),
       y = expression("Leaf Si (mg g"^{-1}*")"),
       color = "Functional group",
       fill  = "Functional group") +
  my_theme +
  theme(legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6))

FigS3

ggsave("FigureS3.tiff",
       plot = FigS3,
       width = 84, height = 90, units = "mm",
       dpi = 300, compression = "lzw", bg = "white")
# =============================================================================
# Fig. S4 — SLA vs. leaf Ca by functional group (based on M5)
# =============================================================================
# M5: SLA ~ leaf_Ca * group + (1|plot_ID) + (1|species)
# Interaction P = 0.012 (significant) — retained.
# Group-specific slopes (from emtrends):
#   Forb:   +0.332, P = 0.065 (marginal)
#   Grass:  +0.275, P = 0.321 (NS)
#   Legume: -0.230, P = 0.053 (marginal)
# =============================================================================

# Get group-specific slopes and their P-values
slope_M5 <- emtrends(M5, ~ group, var = "leaf_Ca") %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as.data.frame()

# Line type: solid if group slope P < 0.05, dashed otherwise
# (All three end up dashed here; see caption for explanation)
slope_M5$linetype <- ifelse(slope_M5$p.value < 0.05, "solid", "dashed")
slope_M5$group    <- as.character(slope_M5$group)

print(slope_M5)  # sanity check

# Predictions
pred_M5 <- ggpredict(M5, terms = c("leaf_Ca", "group")) %>%
  as.data.frame()

pred_M5$group <- as.character(pred_M5$group)

pred_M5 <- pred_M5 %>%
  left_join(slope_M5[, c("group", "linetype")], by = "group")

pred_M5$group <- factor(pred_M5$group, levels = c("Forb", "Grass", "Legume"))

FigS4 <- ggplot(data %>% filter(!is.na(SLA) & !is.na(leaf_Ca)),
                aes(x = leaf_Ca, y = SLA)) +
  geom_point(aes(color = group),
             alpha = 0.5, size = 2, na.rm = TRUE) +
  geom_ribbon(data = pred_M5,
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.15,
              inherit.aes = FALSE) +
  geom_line(data = pred_M5,
            aes(x = x, y = predicted, color = group, linetype = linetype),
            linewidth = 1.1,
            inherit.aes = FALSE) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_linetype_identity() +
  labs(x = expression("Leaf Ca (mg g"^{-1}*")"),
       y = expression("SLA (mm"^2*" mg"^{-1}*")"),
       color = "Functional group",
       fill  = "Functional group") +
  my_theme +
  theme(legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6))

FigS4

ggsave("FigureS4.tiff",
       plot = FigS4,
       width = 84, height = 90, units = "mm",
       dpi = 300, compression = "lzw", bg = "white")
# =============================================================================
# Fig. S5 — SLA vs. leaf Si (based on M6 reduced)
# =============================================================================
# M6 reduced: SLA ~ leaf_Si_log + group (additive, no interaction)
# Full-model interaction (Si × group) was not significant (P = 0.211),
# so we simplified to the additive model. Common Si slope is marginal
# (β = +1.13, P = 0.098), so all three lines are dashed.
# =============================================================================

# Define log-axis breaks for Si (raw values shown on log-scaled x-axis)
raw_breaks_Si <- c(0, 1, 5, 20, 100, 250)
log_breaks_Si <- log(raw_breaks_Si + 1)

# Predictions from reduced model
pred_M6red <- ggpredict(M6_reduced, terms = c("leaf_Si_log", "group")) %>%
  as.data.frame()

pred_M6red$group <- factor(pred_M6red$group,
                           levels = c("Forb", "Grass", "Legume"))

# Common Si slope is P = 0.098 — not below 0.05, so all lines dashed
common_line_type <- "dashed"

FigS5 <- ggplot(data %>% filter(!is.na(SLA) & !is.na(leaf_Si_log)),
                aes(x = leaf_Si_log, y = SLA)) +
  geom_point(aes(color = group),
             alpha = 0.5, size = 2, na.rm = TRUE) +
  geom_ribbon(data = pred_M6red,
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.15,
              inherit.aes = FALSE) +
  geom_line(data = pred_M6red,
            aes(x = x, y = predicted, color = group),
            linetype = common_line_type,
            linewidth = 1.1,
            inherit.aes = FALSE) +
  scale_x_continuous(breaks = log_breaks_Si,
                     labels = raw_breaks_Si) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(x = expression("Leaf Si (mg g"^{-1}*")"),
       y = expression("SLA (mm"^2*" mg"^{-1}*")"),
       color = "Functional group",
       fill  = "Functional group") +
  my_theme +
  theme(legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6))

FigS5

ggsave("FigureS5.tiff",
       plot = FigS5,
       width = 84, height = 90, units = "mm",
       dpi = 300, compression = "lzw", bg = "white")





# Check VIFs for the soil-leaf models
cat("\n========== VIF for M7 (Leaf Si ~ Soil Si × group + exploratorium) ==========\n")
print(vif(M7))

cat("\n========== VIF for M8 (Leaf Ca ~ Soil Ca × group + exploratorium) ==========\n")
print(vif(M8))
# =============================================================================
# Test 3-way interactions with exploratorium
# =============================================================================
# Following  suggestion: test if the shape of the relationship differs 
# among regions. Start with the soil-leaf models where this makes most 
# ecological sense.

# ---- M7 with 3-way interaction ----
M7_3way <- lmer(leaf_Si_log ~ soil_Si * group * exploratorium + 
                  (1|plot_ID) + (1|species),
                data = data)

cat("\n========== M7 with 3-way interaction ==========\n")
print(Anova(M7_3way, type = "III"))
cat("R²:\n"); print(r.squaredGLMM(M7_3way))


# ---- M8 with 3-way interaction ----
M8_3way <- lmer(leaf_Ca ~ soil_Ca * group * exploratorium + 
                  (1|plot_ID) + (1|species),
                data = data)

cat("\n========== M8 with 3-way interaction ==========\n")
print(Anova(M8_3way, type = "III"))
cat("R²:\n"); print(r.squaredGLMM(M8_3way))


# ---- M4 with 3-way (leaf Si vs leaf Ca) ----
M4_3way <- lmer(leaf_Si_log ~ leaf_Ca * group * exploratorium + 
                  (1|plot_ID) + (1|species),
                data = data)

cat("\n========== M4 with 3-way interaction ==========\n")
print(Anova(M4_3way, type = "III"))


# ---- M5 with 3-way (SLA vs leaf Ca) ----
M5_3way <- lmer(SLA ~ leaf_Ca * group * exploratorium + 
                  (1|plot_ID) + (1|species),
                data = data)

cat("\n========== M5 with 3-way interaction ==========\n")
print(Anova(M5_3way, type = "III"))
# =============================================================================
# Refit models keeping only significant interactions with region
# =============================================================================
# Rule (following Joana): drop 3-way interactions that are not significant

# ---- M7: keep soil_Si × exploratorium (2-way sig), drop 3-way ----
M7_final <- lmer(leaf_Si_log ~ soil_Si * group + soil_Si * exploratorium + 
                   (1|plot_ID) + (1|species),
                 data = data)

cat("\n========== M7 FINAL: with soil_Si × region ==========\n")
print(Anova(M7_final, type = "III"))
cat("R²:\n"); print(r.squaredGLMM(M7_final))


# ---- M8: keep group × exploratorium, drop 3-way and soil_Ca:exploratorium ----
M8_final <- lmer(leaf_Ca ~ soil_Ca * group + group * exploratorium + 
                   (1|plot_ID) + (1|species),
                 data = data)

cat("\n========== M8 FINAL: with group × region ==========\n")
print(Anova(M8_final, type = "III"))
cat("R²:\n"); print(r.squaredGLMM(M8_final))


# ---- M4: keep the full 3-way (all significant) ----
M4_final <- lmer(leaf_Si_log ~ leaf_Ca * group * exploratorium + 
                   (1|plot_ID) + (1|species),
                 data = data)

cat("\n========== M4 FINAL: full 3-way interaction ==========\n")
print(Anova(M4_final, type = "III"))
cat("R²:\n"); print(r.squaredGLMM(M4_final))


# ---- M5: keep group × exploratorium and leaf_Ca × group, drop 3-way ----
M5_final <- lmer(SLA ~ leaf_Ca * group + group * exploratorium + 
                   (1|plot_ID) + (1|species),
                 data = data)

cat("\n========== M5 FINAL: with group × region ==========\n")
print(Anova(M5_final, type = "III"))
cat("R²:\n"); print(r.squaredGLMM(M5_final))
# =============================================================================
# Updated Figure 2 — Si-Ca (M4) and SLA-Ca (M5) with region interactions
# =============================================================================
# Now that M4 has a significant 3-way interaction (leaf_Ca × group × region),
# we need to show these relationships FACETED by region.
# For M5, we keep single-panel (no 3-way, but keep group × region additive).

# =========================================================
# Figure 2a — Leaf Si vs Leaf Ca × group × exploratorium (M4)
# =========================================================

# Get group-specific slopes per region
slope_M4_final <- emtrends(M4_final, ~ group | exploratorium, var = "leaf_Ca") %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as.data.frame()

slope_M4_final$linetype <- ifelse(slope_M4_final$p.value < 0.05, "solid", "dashed")
slope_M4_final$group    <- as.character(slope_M4_final$group)

# View the region-specific slopes
print(slope_M4_final)

# Get predictions from M4_final
pred_M4_final <- ggpredict(M4_final, 
                           terms = c("leaf_Ca", "group", "exploratorium")) %>%
  as.data.frame()

pred_M4_final$group         <- as.character(pred_M4_final$group)
pred_M4_final$exploratorium <- as.character(pred_M4_final$facet)

# Match linetype based on group AND region
pred_M4_final <- pred_M4_final %>%
  left_join(slope_M4_final[, c("group", "exploratorium", "linetype")],
            by = c("group", "exploratorium"))

pred_M4_final$group <- factor(pred_M4_final$group, 
                              levels = c("Forb", "Grass", "Legume"))
pred_M4_final$exploratorium <- factor(pred_M4_final$exploratorium, 
                                      levels = c("AEG", "HEG", "SEG"))


Fig2a <- ggplot(data %>% filter(!is.na(leaf_Si_log) & !is.na(leaf_Ca)),
                aes(x = leaf_Ca, y = leaf_Si_log)) +
  
  geom_point(aes(color = group), alpha = 0.5, size = 1.5, na.rm = TRUE) +
  
  geom_ribbon(data = pred_M4_final,
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.20,
              inherit.aes = FALSE) +
  
  geom_line(data = pred_M4_final,
            aes(x = x, y = predicted, color = group, linetype = linetype),
            linewidth = 1,
            inherit.aes = FALSE) +
  
  facet_wrap(~ exploratorium, nrow = 1) +   # 3 panels: AEG | HEG | SEG
  
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_linetype_identity() +
  
  labs(x = "Leaf calcium (mg/g)",
       y = "Leaf silicon (log-transformed, mg/g)",
       color = "Functional group",
       fill  = "Functional group") +
  
  my_theme +
  theme(strip.text = element_text(face = "bold", size = 10))


# =========================================================
# Figure 2b — SLA vs Leaf Ca × group (M5) — SINGLE PANEL
# =========================================================
# The 3-way was NOT significant, so we keep a single-panel figure showing
# the Ca × group interaction (which is significant).

slope_M5_final <- emtrends(M5_final, ~ group, var = "leaf_Ca") %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as.data.frame()

slope_M5_final$linetype <- ifelse(slope_M5_final$p.value < 0.05, "solid", "dashed")
slope_M5_final$group    <- as.character(slope_M5_final$group)

print(slope_M5_final)

pred_M5_final <- ggpredict(M5_final, terms = c("leaf_Ca", "group")) %>%
  as.data.frame()

pred_M5_final$group <- as.character(pred_M5_final$group)

pred_M5_final <- pred_M5_final %>%
  left_join(slope_M5_final[, c("group", "linetype")], by = "group")

pred_M5_final$group <- factor(pred_M5_final$group, 
                              levels = c("Forb", "Grass", "Legume"))


Fig2b <- ggplot(data %>% filter(!is.na(SLA) & !is.na(leaf_Ca)),
                aes(x = leaf_Ca, y = SLA)) +
  
  geom_point(aes(color = group), alpha = 0.5, size = 2, na.rm = TRUE) +
  
  geom_ribbon(data = pred_M5_final,
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.20,
              inherit.aes = FALSE) +
  
  geom_line(data = pred_M5_final,
            aes(x = x, y = predicted, color = group, linetype = linetype),
            linewidth = 1.2,
            inherit.aes = FALSE) +
  
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_linetype_identity() +
  
  labs(x = "Leaf calcium (mg/g)",
       y = expression("Specific leaf area (cm"^2*"/g)"),
       color = "Functional group",
       fill  = "Functional group") +
  
  my_theme


# =========================================================
# Combine Fig 2a (3 panels) on top, Fig 2b (1 panel) below
# =========================================================

Figure2 <- Fig2a / Fig2b +
  plot_layout(guides = "collect", heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = "sans", size = 12, face = "bold"),
        plot.tag.position = c(0.02, 0.98),
        legend.position = "bottom")

Figure2

# Save
ggsave("Figure2.tiff",
       plot = Figure2,
       width = 174,
       height = 160,
       units = "mm",
       dpi = 300,
       compression = "lzw",
       bg = "white")
# =============================================================================
# Updated Figure 3 — Soil-leaf Si (M7) and Soil-leaf Ca (M8)
# =============================================================================
# M7 now has soil_Si × region interaction → 3-panel faceted plot
# M8 has group × region interaction → we can either facet or show single plot
#     Let's do 3-panel for consistency

# =========================================================
# Figure 3a — Leaf Si vs Soil Si × group (M7)
# =========================================================
# soil_Si × group NS, but soil_Si × region SIGNIFICANT
# So we facet by region and show ONE line per region (no group-specific)
# OR we can show per group per region for symmetry with Figure 2

# For consistency with Fig 2, let's show groups within regions
slope_M7_final <- emtrends(M7_final, ~ group | exploratorium, var = "soil_Si") %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as.data.frame()

slope_M7_final$linetype <- ifelse(slope_M7_final$p.value < 0.05, "solid", "dashed")
slope_M7_final$group    <- as.character(slope_M7_final$group)

print(slope_M7_final)

pred_M7_final <- ggpredict(M7_final, 
                           terms = c("soil_Si", "group", "exploratorium")) %>%
  as.data.frame()

pred_M7_final$group         <- as.character(pred_M7_final$group)
pred_M7_final$exploratorium <- as.character(pred_M7_final$facet)

pred_M7_final <- pred_M7_final %>%
  left_join(slope_M7_final[, c("group", "exploratorium", "linetype")],
            by = c("group", "exploratorium"))

pred_M7_final$group <- factor(pred_M7_final$group, 
                              levels = c("Forb", "Grass", "Legume"))
pred_M7_final$exploratorium <- factor(pred_M7_final$exploratorium, 
                                      levels = c("AEG", "HEG", "SEG"))


Fig3a <- ggplot(data %>% filter(!is.na(leaf_Si_log) & !is.na(soil_Si)),
                aes(x = soil_Si, y = leaf_Si_log)) +
  
  geom_point(aes(color = group), alpha = 0.5, size = 1.5, na.rm = TRUE) +
  
  geom_ribbon(data = pred_M7_final,
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.20,
              inherit.aes = FALSE) +
  
  geom_line(data = pred_M7_final,
            aes(x = x, y = predicted, color = group, linetype = linetype),
            linewidth = 1,
            inherit.aes = FALSE) +
  
  facet_wrap(~ exploratorium, nrow = 1) +
  
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_linetype_identity() +
  
  labs(x = "Soil silicon (mg/g)",
       y = "Leaf silicon (log-transformed, mg/g)",
       color = "Functional group",
       fill  = "Functional group") +
  
  my_theme +
  theme(strip.text = element_text(face = "bold", size = 10))


# =========================================================
# Figure 3b — Leaf Ca vs Soil Ca × group (M8)
# =========================================================
# soil_Ca × group NS, group × region SIG
# Show groups within regions for consistency

slope_M8_final <- emtrends(M8_final, ~ group | exploratorium, var = "soil_Ca") %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as.data.frame()

slope_M8_final$linetype <- ifelse(slope_M8_final$p.value < 0.05, "solid", "dashed")
slope_M8_final$group    <- as.character(slope_M8_final$group)

print(slope_M8_final)

pred_M8_final <- ggpredict(M8_final, 
                           terms = c("soil_Ca", "group", "exploratorium")) %>%
  as.data.frame()

pred_M8_final$group         <- as.character(pred_M8_final$group)
pred_M8_final$exploratorium <- as.character(pred_M8_final$facet)

pred_M8_final <- pred_M8_final %>%
  left_join(slope_M8_final[, c("group", "exploratorium", "linetype")],
            by = c("group", "exploratorium"))

pred_M8_final$group <- factor(pred_M8_final$group, 
                              levels = c("Forb", "Grass", "Legume"))
pred_M8_final$exploratorium <- factor(pred_M8_final$exploratorium, 
                                      levels = c("AEG", "HEG", "SEG"))


Fig3b <- ggplot(data %>% filter(!is.na(leaf_Ca) & !is.na(soil_Ca)),
                aes(x = soil_Ca, y = leaf_Ca)) +
  
  geom_point(aes(color = group), alpha = 0.5, size = 1.5, na.rm = TRUE) +
  
  geom_ribbon(data = pred_M8_final,
              aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.20,
              inherit.aes = FALSE) +
  
  geom_line(data = pred_M8_final,
            aes(x = x, y = predicted, color = group, linetype = linetype),
            linewidth = 1,
            inherit.aes = FALSE) +
  
  facet_wrap(~ exploratorium, nrow = 1, scales = "free_x") +   # free_x since soil Ca ranges differ hugely
  
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_linetype_identity() +
  
  labs(x = "Soil calcium (mg/g)",
       y = "Leaf calcium (mg/g)",
       color = "Functional group",
       fill  = "Functional group") +
  
  my_theme +
  theme(strip.text = element_text(face = "bold", size = 10))


# =========================================================
# Combine Fig 3a on top, Fig 3b below
# =========================================================

Figure3 <- Fig3a / Fig3b +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = "sans", size = 12, face = "bold"),
        plot.tag.position = c(0.02, 0.98),
        legend.position = "bottom")

Figure3

# Save
ggsave("Figure3.tiff",
       plot = Figure3,
       width = 174,
       height = 160,
       units = "mm",
       dpi = 300,
       compression = "lzw",
       bg = "white")
# =============================================================================
# Section 11: Soil chemistry patterns across regions
# =============================================================================
# Explore how soil Si and soil Ca vary among AEG, HEG, and SEG regions.
# This helps interpret why region absorbs soil signals in the models.


# ---- 11.1 Descriptive statistics per region ----
# Get plot-level values (one per plot, not repeated across species)

soil_by_region <- data %>%
  dplyr::select(plot_ID, exploratorium, soil_Si, soil_Ca) %>%
  distinct() %>%
  group_by(exploratorium) %>%
  summarise(
    n_plots      = n(),
    # Soil Si
    soil_Si_mean = round(mean(soil_Si, na.rm = TRUE), 3),
    soil_Si_sd   = round(sd(soil_Si, na.rm = TRUE), 3),
    soil_Si_min  = round(min(soil_Si, na.rm = TRUE), 3),
    soil_Si_max  = round(max(soil_Si, na.rm = TRUE), 3),
    # Soil Ca
    soil_Ca_mean = round(mean(soil_Ca, na.rm = TRUE), 2),
    soil_Ca_sd   = round(sd(soil_Ca, na.rm = TRUE), 2),
    soil_Ca_min  = round(min(soil_Ca, na.rm = TRUE), 2),
    soil_Ca_max  = round(max(soil_Ca, na.rm = TRUE), 2)
  )

cat("\n========== Soil chemistry summary per region ==========\n")
print(soil_by_region)


# ---- 11.2 Statistical test: does soil chemistry differ by region? ----

# Get one soil value per plot (not repeated)
soil_plot_level <- data %>%
  dplyr::select(plot_ID, exploratorium, soil_Si, soil_Ca) %>%
  distinct()

# Test soil Si by region
cat("\n========== ANOVA: Soil Si by region ==========\n")
soil_Si_model <- lm(soil_Si ~ exploratorium, data = soil_plot_level)
print(Anova(soil_Si_model, type = "II"))

# Tukey pairwise for soil Si
cat("\nTukey pairwise for soil Si:\n")
emm_soilSi <- emmeans(soil_Si_model, ~ exploratorium)
print(pairs(emm_soilSi, adjust = "tukey"))
print(multcomp::cld(emm_soilSi, Letters = letters, adjust = "tukey"))


# Test soil Ca by region
cat("\n========== ANOVA: Soil Ca by region ==========\n")
soil_Ca_model <- lm(soil_Ca ~ exploratorium, data = soil_plot_level)
print(Anova(soil_Ca_model, type = "II"))

# Tukey pairwise for soil Ca
cat("\nTukey pairwise for soil Ca:\n")
emm_soilCa <- emmeans(soil_Ca_model, ~ exploratorium)
print(pairs(emm_soilCa, adjust = "tukey"))
print(multcomp::cld(emm_soilCa, Letters = letters, adjust = "tukey"))


# ---- 11.3 Correlation between soil Si and soil Ca ----
cat("\n========== Correlation between soil Si and soil Ca ==========\n")
cor_test <- cor.test(soil_plot_level$soil_Si, soil_plot_level$soil_Ca,
                     method = "pearson", use = "complete.obs")
print(cor_test)

# Within each region
cat("\nWithin-region correlations:\n")
for (r in levels(soil_plot_level$exploratorium)) {
  sub <- soil_plot_level %>% filter(exploratorium == r)
  ct <- cor.test(sub$soil_Si, sub$soil_Ca, method = "pearson", 
                 use = "complete.obs")
  cat(sprintf("%s: n=%d  r=%6.3f  p=%.4f\n",
              r, nrow(sub), ct$estimate, ct$p.value))
}


# ---- 11.4 Visualize soil Si and soil Ca by region ----

# Panel a: soil Si by region (violin plot)
# Updated soil-region colors
# AEG (calcareous, Swabian Alb) → dark brown
# HEG (limestone, Hainich)      → sienna
# SEG (sandy acidic, Schorfheide) → gold

region_colors <- c("AEG" = "#3E2C1C",   # very dark brown
                   "HEG" = "#A0522D",   # sienna (unchanged)
                   "SEG" = "#DAA520")   # gold (unchanged)


# Panel a: soil Si by region
FigSoil_a <- ggplot(soil_plot_level %>% filter(!is.na(soil_Si)),
                    aes(x = exploratorium, y = soil_Si, fill = exploratorium)) +
  geom_violin(alpha = 0.7, color = NA, trim = FALSE) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.7,
              aes(color = exploratorium), show.legend = FALSE) +
  stat_summary(fun = median, geom = "point",
               shape = 95, size = 8, color = "black") +
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = region_colors) +
  labs(x = "Region", y = "Soil Si (mg/g)") +
  my_theme +
  theme(legend.position = "none")


# Panel b: soil Ca by region
FigSoil_b <- ggplot(soil_plot_level %>% filter(!is.na(soil_Ca)),
                    aes(x = exploratorium, y = soil_Ca, fill = exploratorium)) +
  geom_violin(alpha = 0.7, color = NA, trim = FALSE) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.7,
              aes(color = exploratorium), show.legend = FALSE) +
  stat_summary(fun = median, geom = "point",
               shape = 95, size = 8, color = "black") +
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = region_colors) +
  labs(x = "Region", y = "Soil Ca (mg/g)") +
  my_theme +
  theme(legend.position = "none")


# Panel c: soil Si vs soil Ca (scatter)
FigSoil_c <- ggplot(soil_plot_level %>% filter(!is.na(soil_Si) & !is.na(soil_Ca)),
                    aes(x = soil_Si, y = soil_Ca, color = exploratorium)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = region_colors) +
  labs(x = "Soil Si (mg/g)", y = "Soil Ca (mg/g)",
       color = "Region") +
  my_theme


# Combine
FigureSoil <- (FigSoil_a + FigSoil_b) / FigSoil_c +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = "sans", size = 12, face = "bold"),
        plot.tag.position = c(0.02, 0.98))

FigureSoil

# Save
ggsave("Figure_soil_by_region.tiff",
       plot = FigureSoil,
       width = 174,
       height = 180,
       units = "mm",
       dpi = 300,
       compression = "lzw",
       bg = "white")
# =============================================================================
# Section 12: Leaf Si and Ca variation by region (violin plots)
# =============================================================================
# Check how leaf chemistry varies across the three regions,
# ignoring functional group for now.

# Region colors (matching the soil-by-region figures)
region_colors <- c("AEG" = "#3E2C1C",   # very dark brown
                   "HEG" = "#A0522D",   # sienna
                   "SEG" = "#DAA520")   # gold


# ---- 12.1 Statistical test: does leaf chemistry differ by region? ----

# Simple mixed-effects models: response ~ region + (1|species) + (1|plot_ID)
# We keep the random effects to correctly account for sampling structure

# Leaf Si by region
leafSi_by_region <- lmer(leaf_Si_log ~ exploratorium + (1|plot_ID) + (1|species),
                         data = data)

cat("\n========== Leaf Si by region (ANOVA) ==========\n")
print(Anova(leafSi_by_region, type = "II"))

cat("\nTukey pairwise for leaf Si:\n")
emm_leafSi <- emmeans(leafSi_by_region, ~ exploratorium)
print(pairs(emm_leafSi, adjust = "tukey"))
print(multcomp::cld(emm_leafSi, Letters = letters, adjust = "tukey"))

# Leaf Ca by region
leafCa_by_region <- lmer(leaf_Ca ~ exploratorium + (1|plot_ID) + (1|species),
                         data = data)

cat("\n========== Leaf Ca by region (ANOVA) ==========\n")
print(Anova(leafCa_by_region, type = "II"))

cat("\nTukey pairwise for leaf Ca:\n")
emm_leafCa <- emmeans(leafCa_by_region, ~ exploratorium)
print(pairs(emm_leafCa, adjust = "tukey"))
print(multcomp::cld(emm_leafCa, Letters = letters, adjust = "tukey"))


# ---- 12.2 Descriptive statistics per region ----

leaf_by_region <- data %>%
  group_by(exploratorium) %>%
  summarise(
    n            = n(),
    leaf_Si_mean = round(mean(leaf_Si, na.rm = TRUE), 1),
    leaf_Si_sd   = round(sd(leaf_Si, na.rm = TRUE), 1),
    leaf_Si_min  = round(min(leaf_Si, na.rm = TRUE), 1),
    leaf_Si_max  = round(max(leaf_Si, na.rm = TRUE), 1),
    leaf_Ca_mean = round(mean(leaf_Ca, na.rm = TRUE), 2),
    leaf_Ca_sd   = round(sd(leaf_Ca, na.rm = TRUE), 2),
    leaf_Ca_min  = round(min(leaf_Ca, na.rm = TRUE), 2),
    leaf_Ca_max  = round(max(leaf_Ca, na.rm = TRUE), 2)
  )

cat("\n========== Leaf chemistry summary per region ==========\n")
print(leaf_by_region)


# ---- 12.3 Violin plots ----

# We need Tukey letters — extract them from the CLD output
# After running the code above, manually update these based on your output:
# (I'll use placeholders — replace after seeing the letters)

# Letters for leaf Si — you'll need to check output and update
letters_leafSi <- data.frame(
  exploratorium = factor(c("AEG", "HEG", "SEG"),
                         levels = c("AEG", "HEG", "SEG")),
  label = c("a", "b", "a"),  # UPDATE with actual Tukey letters
  y     = max(data$leaf_Si_log, na.rm = TRUE) + 0.3
)

# Letters for leaf Ca — you'll need to check output and update
letters_leafCa <- data.frame(
  exploratorium = factor(c("AEG", "HEG", "SEG"),
                         levels = c("AEG", "HEG", "SEG")),
  label = c("b", "a", "a"),  # UPDATE with actual Tukey letters
  y     = max(data$leaf_Ca, na.rm = TRUE) + 2
)


# Panel a: Leaf Si by region
FigLeaf_a <- ggplot(data %>% filter(!is.na(leaf_Si_log)),
                    aes(x = exploratorium, y = leaf_Si_log, fill = exploratorium)) +
  
  geom_violin(alpha = 0.7, color = NA, trim = FALSE) +
  
  geom_jitter(aes(color = exploratorium),
              width = 0.15, size = 1.5, alpha = 0.5,
              show.legend = FALSE) +
  
  stat_summary(fun = median, geom = "point",
               shape = 95, size = 6, color = "black") +
  
  geom_text(data = letters_leafSi,
            aes(x = exploratorium, y = y, label = label),
            inherit.aes = FALSE,
            size = 4, fontface = "bold") +
  
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = region_colors) +
  
  labs(x = "Region",
       y = "Leaf silicon (log-transformed, mg/g)") +
  
  my_theme +
  theme(legend.position = "none")


# Panel b: Leaf Ca by region
FigLeaf_b <- ggplot(data %>% filter(!is.na(leaf_Ca)),
                    aes(x = exploratorium, y = leaf_Ca, fill = exploratorium)) +
  
  geom_violin(alpha = 0.7, color = NA, trim = FALSE) +
  
  geom_jitter(aes(color = exploratorium),
              width = 0.15, size = 1.5, alpha = 0.5,
              show.legend = FALSE) +
  
  stat_summary(fun = median, geom = "point",
               shape = 95, size = 6, color = "black") +
  
  geom_text(data = letters_leafCa,
            aes(x = exploratorium, y = y, label = label),
            inherit.aes = FALSE,
            size = 4, fontface = "bold") +
  
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = region_colors) +
  
  labs(x = "Region",
       y = "Leaf calcium (mg/g)") +
  
  my_theme +
  theme(legend.position = "none")


# Combine
FigureLeaf <- FigLeaf_a + FigLeaf_b +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = "sans", size = 12, face = "bold"),
        plot.tag.position = c(0.02, 0.98))

FigureLeaf

# Save
ggsave("Figure_leaf_by_region.tiff",
       plot = FigureLeaf,
       width = 168,
       height = 84,
       units = "mm",
       dpi = 300,
       compression = "lzw",
       bg = "white")
# =============================================================================
# Simple ANOVA (raw data, no random effects) for regional letters
# =============================================================================

# Leaf Si — simple model
leafSi_simple <- lm(leaf_Si_log ~ exploratorium, data = data)

cat("\n========== Leaf Si by region — simple ANOVA ==========\n")
print(Anova(leafSi_simple, type = "II"))

cat("\nTukey pairwise (simple model):\n")
emm_leafSi_simple <- emmeans(leafSi_simple, ~ exploratorium)
print(pairs(emm_leafSi_simple, adjust = "tukey"))
print(multcomp::cld(emm_leafSi_simple, Letters = letters, adjust = "tukey"))


# Leaf Ca — simple model
leafCa_simple <- lm(leaf_Ca ~ exploratorium, data = data)

cat("\n========== Leaf Ca by region — simple ANOVA ==========\n")
print(Anova(leafCa_simple, type = "II"))

cat("\nTukey pairwise (simple model):\n")
emm_leafCa_simple <- emmeans(leafCa_simple, ~ exploratorium)
print(pairs(emm_leafCa_simple, adjust = "tukey"))
print(multcomp::cld(emm_leafCa_simple, Letters = letters, adjust = "tukey"))
# Within-species Ca-SLA correlations
cat("\n========== Within-species Ca-SLA relationships ==========\n")
for (sp in levels(data$species)) {
  d_sub <- data %>% 
    filter(species == sp & !is.na(SLA) & !is.na(leaf_Ca))
  
  if (nrow(d_sub) >= 10) {  # only species with enough data
    ct <- cor.test(d_sub$leaf_Ca, d_sub$SLA, method = "pearson")
    cat(sprintf("%-25s (%s) n=%2d  r=%+.3f  p=%.3f\n",
                sp, unique(as.character(d_sub$group)), 
                nrow(d_sub), ct$estimate, ct$p.value))
  }
}
# ---------------------------------------------------------
# Raw (untransformed) group means ± SE for leaf_Si and leaf_Ca
# by functional group (Grass / Legume / Forb)
# ---------------------------------------------------------

library(dplyr)

# Standard error function
se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

group_summary <- data %>%
  group_by(group) %>%
  summarise(
    n_Si       = sum(!is.na(leaf_Si)),
    mean_Si    = mean(leaf_Si, na.rm = TRUE),
    se_Si      = se(leaf_Si),
    n_Ca       = sum(!is.na(leaf_Ca)),
    mean_Ca    = mean(leaf_Ca, na.rm = TRUE),
    se_Ca      = se(leaf_Ca),
    .groups = "drop"
  )

print(group_summary)

# Optional: also get overall n per group (rows), useful for the Methods/Table
n_per_group <- data %>%
  group_by(group) %>%
  summarise(n_rows = n(), .groups = "drop")

print(n_per_group)

# ---------------------------------------------------------
# If you'd rather have model-based estimated marginal means
# (accounting for the random effects structure) instead of
# raw arithmetic means, use emmeans on your fitted models:
# ---------------------------------------------------------

# library(emmeans)
#
# For Ca (untransformed model, M2) - gives means directly on raw scale:
# emmeans(M2, ~ group)
#
# For Si (M1 is fit on leaf_Si_log) - back-transform to get geometric means:
# emmeans(M1, ~ group, type = "response")
#   Note: type = "response" back-transforms automatically since M1's
#   response is leaf_Si_log; check that emmeans detects the log link
#   correctly (it works automatically if you fit it as log(leaf_Si) in
#   the formula, but if leaf_Si_log is a pre-computed column, you may
#   need to specify the transformation manually, e.g.:
#
# emm_Si <- emmeans(M1, ~ group)
# summary(emm_Si, type = "response")  # may not back-transform automatically
#
# Safer manual back-transform if needed:
# emm_Si_df <- as.data.frame(emmeans(M1, ~ group))
# emm_Si_df$response_mean <- exp(emm_Si_df$emmean)
# emm_Si_df$lower <- exp(emm_Si_df$emmean - emm_Si_df$SE)
# emm_Si_df$upper <- exp(emm_Si_df$emmean + emm_Si_df$SE)
# print(emm_Si_df)
# ---------------------------------------------------------
# Leave-one-species-out (LOSO) sensitivity check for M5
# Refits SLA ~ leaf_Ca * group + (1|plot_ID) + (1|species)
# 12 times, each time excluding one species.
# For each refit:
#   - Records the Ca × group interaction P-value (Type III)
#   - Records the group-specific slopes for Forb, Grass, Legume
#     with their SEs and P-values
# ---------------------------------------------------------

library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(dplyr)

species_list <- levels(data$species)
loso_results <- data.frame()

for (sp in species_list) {
  d_loso <- data %>% filter(species != sp) %>% droplevels()
  
  # Refit M5 without species sp
  m5_loso <- tryCatch(
    lmer(SLA ~ leaf_Ca * group + (1|plot_ID) + (1|species), 
         data = d_loso),
    error = function(e) NULL
  )
  
  if (is.null(m5_loso)) next
  
  # Type III ANOVA — get interaction P
  anova_out <- Anova(m5_loso, type = "III")
  interaction_p <- anova_out["leaf_Ca:group", "Pr(>Chisq)"]
  main_group_p  <- anova_out["group",         "Pr(>Chisq)"]
  main_Ca_p     <- anova_out["leaf_Ca",       "Pr(>Chisq)"]
  
  # Group-specific slopes with SE, CI, P
  slopes <- summary(emtrends(m5_loso, ~ group, var = "leaf_Ca"),
                    infer = c(TRUE, TRUE))
  
  # Extract slope info per group (handle if a group is missing)
  get_group <- function(grp, col) {
    if (grp %in% slopes$group) slopes[slopes$group == grp, col] else NA
  }
  
  loso_results <- rbind(loso_results, data.frame(
    excluded_species = sp,
    n_obs            = nobs(m5_loso),
    n_species        = length(unique(d_loso$species)),
    interaction_p    = interaction_p,
    group_p          = main_group_p,
    Ca_p             = main_Ca_p,
    Forb_slope       = get_group("Forb",   "leaf_Ca.trend"),
    Forb_SE          = get_group("Forb",   "SE"),
    Forb_p           = get_group("Forb",   "p.value"),
    Grass_slope      = get_group("Grass",  "leaf_Ca.trend"),
    Grass_SE         = get_group("Grass",  "SE"),
    Grass_p          = get_group("Grass",  "p.value"),
    Legume_slope     = get_group("Legume", "leaf_Ca.trend"),
    Legume_SE        = get_group("Legume", "SE"),
    Legume_p         = get_group("Legume", "p.value")
  ))
}

slopes_full <- summary(emtrends(M5, ~ group, var = "leaf_Ca"), infer = c(TRUE, TRUE))
anova_full  <- Anova(M5, type = "III")

full_row <- data.frame(
  excluded_species = "None (full model)",
  n_obs            = nobs(M5),
  n_species        = length(unique(data$species)),
  interaction_p    = anova_full["leaf_Ca:group", "Pr(>Chisq)"],
  group_p          = anova_full["group", "Pr(>Chisq)"],
  Ca_p             = anova_full["leaf_Ca", "Pr(>Chisq)"],
  Forb_slope       = slopes_full[slopes_full$group == "Forb", "leaf_Ca.trend"],
  Forb_SE          = slopes_full[slopes_full$group == "Forb", "SE"],
  Forb_p           = slopes_full[slopes_full$group == "Forb", "p.value"],
  Grass_slope      = slopes_full[slopes_full$group == "Grass", "leaf_Ca.trend"],
  Grass_SE         = slopes_full[slopes_full$group == "Grass", "SE"],
  Grass_p          = slopes_full[slopes_full$group == "Grass", "p.value"],
  Legume_slope     = slopes_full[slopes_full$group == "Legume", "leaf_Ca.trend"],
  Legume_SE        = slopes_full[slopes_full$group == "Legume", "SE"],
  Legume_p         = slopes_full[slopes_full$group == "Legume", "p.value"]
)
loso_results <- rbind(loso_results, full_row)

# Print full LOSO table
cat("\n========== Leave-one-species-out results for M5 ==========\n")
print(loso_results, row.names = FALSE, digits = 3)

# Optional: save for supplementary
# write.csv(table_S1, "Table_S1_species_means.csv", row.names = FALSE)
## M4 simplified: drop the leaf_Ca × group interaction
M4_reduced <- lmer(leaf_Si_log ~ leaf_Ca + group + 
                     (1|plot_ID) + (1|species), 
                   data = data)

# Model summary
summary(M4_reduced)

# Type II ANOVA (Type II is appropriate for main-effects-only models,
# since there's no interaction — Type II is more powerful than Type III here)
Anova(M4_reduced, type = "II")

# R²
r.squaredGLMM(M4_reduced)

# Optional: pairwise group contrasts and Ca slope
emmeans(M4_reduced, pairwise ~ group, adjust = "tukey")
emtrends(M4_reduced, ~ 1, var = "leaf_Ca")  # single slope, common across groups


# ---- Compute raw Si:Ca ratio for every sample ----
data <- data %>%
  mutate(
    Si_Ca_raw = leaf_Si / leaf_Ca   # simple raw ratio
  )

# ---- Print all raw ratios per sample (to check for zeros, extremes) ----
cat("\n== Raw Si:Ca ratio distribution ==\n")
summary(data$Si_Ca_raw)

# Full print of ratios (with group and species for context)
options(na.print = "NA")

data %>%
  dplyr::select(species, group, leaf_Si, leaf_Ca, Si_Ca_raw) %>%
  arrange(group, species) %>%
  print(n = Inf)

# ---- Distribution by group ----
cat("\n== Raw ratio summary by group ==\n")
data %>%
  group_by(group) %>%
  summarise(
    n       = sum(!is.na(Si_Ca_raw)),
    min     = min(Si_Ca_raw, na.rm = TRUE),
    q25     = quantile(Si_Ca_raw, 0.25, na.rm = TRUE),
    median  = median(Si_Ca_raw, na.rm = TRUE),
    mean    = mean(Si_Ca_raw, na.rm = TRUE),
    q75     = quantile(Si_Ca_raw, 0.75, na.rm = TRUE),
    max     = max(Si_Ca_raw, na.rm = TRUE),
    sd      = sd(Si_Ca_raw, na.rm = TRUE)
  )

# ---- Check for zeros or infinities ----
cat("\n== Special values check ==\n")
cat("Number of NA in leaf_Si:", sum(is.na(data$leaf_Si)), "\n")
cat("Number of NA in leaf_Ca:", sum(is.na(data$leaf_Ca)), "\n")
cat("Number of leaf_Si = 0:", sum(data$leaf_Si == 0, na.rm = TRUE), "\n")
cat("Number of leaf_Ca = 0:", sum(data$leaf_Ca == 0, na.rm = TRUE), "\n")
cat("Number of Si_Ca_raw = 0 (Si=0):", sum(data$Si_Ca_raw == 0, na.rm = TRUE), "\n")
cat("Number of Si_Ca_raw = Inf (Ca=0):", sum(is.infinite(data$Si_Ca_raw), na.rm = TRUE), "\n")
cat("Number of Si_Ca_raw = NaN:", sum(is.nan(data$Si_Ca_raw)), "\n")

# ---- Histogram of raw ratio ----
p1 <- ggplot(data %>% filter(!is.na(Si_Ca_raw)),
             aes(x = Si_Ca_raw, fill = group)) +
  geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = group_colors) +
  labs(x = "Raw Si:Ca ratio", y = "Count",
       title = "Raw Si:Ca distribution (untransformed)") +
  theme_minimal()

# ---- Histogram of log-transformed ratio for comparison ----
p2 <- ggplot(data %>% filter(!is.na(Si_Ca_raw) & Si_Ca_raw > 0),
             aes(x = log(Si_Ca_raw), fill = group)) +
  geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = group_colors) +
  labs(x = "log(Si:Ca ratio)", y = "Count",
       title = "Log-transformed Si:Ca distribution") +
  theme_minimal()

# Show both side-by-side
library(patchwork)
p1 / p2   # stacked vertically for easy comparison

# ---- Shapiro-Wilk normality test (for raw vs log) ----
# Only informative for smaller samples; approximate here
cat("\n== Shapiro-Wilk normality tests (approximate) ==\n")
cat("Raw ratio:\n")
print(shapiro.test(sample(data$Si_Ca_raw[!is.na(data$Si_Ca_raw)],
                          min(5000, sum(!is.na(data$Si_Ca_raw))))))
cat("\nLog ratio (excluding zero and negative Si):\n")
valid_log <- data$Si_Ca_raw[!is.na(data$Si_Ca_raw) & data$Si_Ca_raw > 0]
print(shapiro.test(sample(log(valid_log), min(5000, length(valid_log)))))
# =============================================================================
# =============================================================================
# Fig. S3 (or Fig. 3) — Leaf Ca vs. Leaf Si with pooled correlation line
# =============================================================================
# Shows all observations color-coded by functional group, with a single
# pooled OLS regression line (log Si ~ leaf Ca) to visualize the overall
# negative correlation (r = −0.70).
# =============================================================================

# ---- Prepare y-axis breaks for log-scaled Si (raw labels) ----
raw_breaks_Si <- c(0, 1, 5, 20, 100, 250)
log_breaks_Si <- log(raw_breaks_Si + 1)

# ---- Filter data to complete cases (matching what M4 uses) ----
plot_data <- data %>%
  filter(!is.na(leaf_Si_log) & !is.na(leaf_Ca))

# ---- Build the plot ----
FigS10 <- ggplot(plot_data,
                aes(x = leaf_Ca, y = leaf_Si_log)) +
  
  # Color-coded points by functional group
  geom_point(aes(color = group),
             alpha = 0.6, size = 2, na.rm = TRUE) +
  
  # Pooled OLS regression line (single line across all data)
  # Note: this is NOT the mixed-model line — it visualizes the pooled
  # Pearson correlation only. See caption note.
  geom_smooth(method = "lm", 
              se = TRUE,               # 95% CI band
              color = "black",         # neutral color to distinguish from group lines
              linetype = "solid",      # pooled correlation is significant (r = -0.70, P < 0.001)
              linewidth = 0.9,
              fill = "gray60",         # gray CI ribbon
              alpha = 0.25) +
  
  # Log-scaled y-axis with raw-value labels
  scale_y_continuous(breaks = log_breaks_Si,
                     labels = raw_breaks_Si) +
  
  scale_color_manual(values = group_colors) +
  
  labs(x = expression("Leaf Ca (mg g"^{-1}*")"),
       y = expression("Leaf Si (mg g"^{-1}*")"),
       color = "Functional group") +
  
  my_theme

FigS10
########################
M4 <- lmer(leaf_Si_log ~ leaf_Ca+ 
             (1|plot_ID)+(1|species), 
           data = data)
cat("\n================ M4: leaf_Si ~ leaf_Ca * group ================\n")
print(summary(M4))
cat("\n--- ANOVA Type III ---\n")
print(Anova(M4, type = "II"))
cat("\n--- R² (marginal/conditional) ---\n")
print(r.squaredGLMM(M4))
###############################
fig_3 <- ggplot(data %>% filter(!is.na(leaf_Si_log) & !is.na(leaf_Ca)),
                aes(x = leaf_Ca, y = leaf_Si, colour = group)) +
  
  geom_point(alpha = 0.5, size = 2, na.rm = TRUE) +
  
  scale_colour_manual(values = group_colors, name = "Functional group") +
  scale_y_continuous(
    trans = "log1p",
    breaks = c(0, 1, 5, 20, 60, 120),
    labels = c("0", "1", "5", "20", "60", "120")
  ) +
  
  labs(
    x = expression("Leaf Ca (mg g"^{-1}*")"),
    y = expression("Leaf Si (mg g"^{-1}*")")
  ) +
  
  my_theme

fig_3
fig_3 <- fig_3 +
  theme(
    legend.text  = element_text(family = "sans", size = 7),
    legend.title = element_text(family = "sans", size = 7),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.x = unit(0.1, "cm")
  )
ggsave("Figure3.tiff",
       plot = fig_3,
       width = 84,
       height = 90,
       units = "mm",
       dpi = 300,
       compression = "lzw",
       bg = "white")
# Within-species correlations: leaf Si vs SLA
library(dplyr)

si_sla_by_species <- data %>%
  filter(!is.na(leaf_Si_log) & !is.na(SLA)) %>%
  group_by(species, group) %>%
  filter(n() >= 5) %>%   # minimum sample size
  summarise(
    n     = n(),
    r     = cor.test(leaf_Si_log, SLA)$estimate,
    ci_lo = cor.test(leaf_Si_log, SLA)$conf.int[1],
    ci_hi = cor.test(leaf_Si_log, SLA)$conf.int[2],
    p     = cor.test(leaf_Si_log, SLA)$p.value,
    .groups = "drop"
  ) %>%
  mutate(sig = ifelse(p < 0.05, "*", "")) %>%
  arrange(group, species)

print(si_sla_by_species, n = 20)

# Also within-group correlations
si_sla_by_group <- data %>%
  filter(!is.na(leaf_Si_log) & !is.na(SLA)) %>%
  group_by(group) %>%
  summarise(
    n     = n(),
    r     = cor.test(leaf_Si_log, SLA)$estimate,
    ci_lo = cor.test(leaf_Si_log, SLA)$conf.int[1],
    ci_hi = cor.test(leaf_Si_log, SLA)$conf.int[2],
    p     = cor.test(leaf_Si_log, SLA)$p.value,
    .groups = "drop"
  )

print(si_sla_by_group)

# Overall correlation
cor_si_sla_overall <- cor.test(data$leaf_Si_log, data$SLA, use = "complete.obs")
print(cor_si_sla_overall)
# Within-species correlations: leaf Ca vs SLA
ca_sla_by_species <- data %>%
  filter(!is.na(leaf_Ca) & !is.na(SLA)) %>%
  group_by(species, group) %>%
  filter(n() >= 5) %>%
  summarise(
    n     = n(),
    r     = cor.test(leaf_Ca, SLA)$estimate,
    ci_lo = cor.test(leaf_Ca, SLA)$conf.int[1],
    ci_hi = cor.test(leaf_Ca, SLA)$conf.int[2],
    p     = cor.test(leaf_Ca, SLA)$p.value,
    .groups = "drop"
  ) %>%
  mutate(sig = ifelse(p < 0.05, "*", "")) %>%
  arrange(group, species)

print(ca_sla_by_species, n = 20)

# Within-group correlations
ca_sla_by_group <- data %>%
  filter(!is.na(leaf_Ca) & !is.na(SLA)) %>%
  group_by(group) %>%
  summarise(
    n     = n(),
    r     = cor.test(leaf_Ca, SLA)$estimate,
    ci_lo = cor.test(leaf_Ca, SLA)$conf.int[1],
    ci_hi = cor.test(leaf_Ca, SLA)$conf.int[2],
    p     = cor.test(leaf_Ca, SLA)$p.value,
    .groups = "drop"
  )

print(ca_sla_by_group)

# Overall correlation
cor_ca_sla_overall <- cor.test(data$leaf_Ca, data$SLA, use = "complete.obs")
print(cor_ca_sla_overall)
#######################
getOption("contrasts")
contrasts(data$group) <- contr.sum(3)
M5 <- lmer(SLA ~ leaf_Ca * group + (1|plot_ID) + (1|species), data = data)
Anova(M5, type = "III")
contrasts(data$group) <- contr.sum(3)
M5_sum <- lmer(SLA ~ leaf_Ca * group + (1|plot_ID) + (1|species), data = data)
Anova(M5_sum, type = "III")
contrasts(data$group) <- contr.treatment(3)
