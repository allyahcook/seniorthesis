### Senior Thesis Analysis R-code ###
# Ally Ah Cook
# Completed March 9, 2026

# Set WD & load necessary packages
setwd("/Users/ally/Downloads/Senior Thesis/R Files") # Change to your own directory
library(vegan)
library(tidyr)
library(dplyr)
library(ggplot2)
library(nlme)
library(readxl)

# set.seed() for PERMANOVAs
set.seed(123)

# Load full data file, replace with your own pathname
full_data <- read_excel('/Users/ally/Downloads/Senior Thesis/Processed Data/20260129_Hawaii_Thesis_FullData.xlsx')

# For summing isomers
pfas_order <- names(full_data)[14:77] %>% 
  str_remove("^(Br-|L-)") %>% 
  unique()


### 3.1 Detection Frequencies ####################################################################################################
# For detection frequencies only, change <LOQ to 0 in new DF
full_data_detect <- full_data %>% 
  mutate(across(14:77, ~ as.numeric(replace(as.character(.), . == "<LOQ", 0)))) %>%  # Columns 14-77 are all PFAS 
  pivot_longer(14:77, names_to = "compound", values_to = "conc") %>%
  mutate(compound = str_remove(compound, "^(Br-|L-)")) %>%
  group_by(across(1:13), compound) %>%   # group by metadata columns only
  summarise(conc = sum(conc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = compound, values_from = conc) %>%
  relocate(all_of(pfas_order), .after = 13)

# Make object with list of all PFAS 
pfas_cols <- colnames(full_data_detect)[14:70] # Columns 14-70 are PFAS after summing isomers

### Water detection frequencies by site
detect_freq_water <- full_data_detect %>%
  filter(type == "water") %>% 
  group_by(site) %>% 
  summarise(across(all_of(pfas_cols),
                   ~ mean(. > 0, na.rm = TRUE) * 100), # Marks compound as "TRUE" if [PFAS] > 0, 
            .groups = "drop") %>%                      # averages # of "TRUE" values by site
  pivot_longer(
    cols = -site,
    names_to = "Compound",
    values_to = "Detection_Frequency_Water"
  )

# Extra DF for plotting
detect_freq_water_plot <- detect_freq_water %>%
  group_by(Compound) %>%
  filter(any(Detection_Frequency_Water > 0)) %>%
  ungroup() # Only keep compounds with non-zero detection frequencies
detect_freq_water_plot$Compound <- factor(
  detect_freq_water_plot$Compound,
  levels = original_order
) # Keep compounds in original order
detect_freq_water_plot$site <- factor(
  detect_freq_water_plot$site,
  levels = c("PF", "WR", "SB", "MB")
) # Keep sites in order

### Biota detection frequencies by site (All biota types)
detect_freq_biota <- full_data_detect %>%
  filter(type == "biota") %>% 
  group_by(site) %>% 
  summarise(across(all_of(pfas_cols),
                   ~ mean(. > 0, na.rm = TRUE) * 100),
            .groups = "drop") %>%
  pivot_longer(
    cols = -site,
    names_to = "Compound",
    values_to = "Detection_Frequency_Biota"
  )

# Extra DF for plotting
detect_freq_biota_plot <- detect_freq_biota %>%
  group_by(Compound) %>%
  filter(any(Detection_Frequency_Biota > 0)) %>%
  ungroup() # Only keep compounds with non-zero detection frequencies
detect_freq_biota_plot$Compound <- factor(
  detect_freq_biota_plot$Compound,
  levels = original_order
) # Keep compounds in original order
detect_freq_biota_plot$site <- factor(
  detect_freq_biota_plot$site,
  levels = c("PF", "SB", "OO", "MB")
)

### Biota detection frequencies by site, separated by biota type
detect_freq_biota_type <- full_data_detect %>%
  filter(type == "biota") %>% 
  group_by(biota_type, site) %>% 
  summarise(across(all_of(pfas_cols),
                   ~ mean(. > 0, na.rm = TRUE) * 100),
            .groups = "drop") %>%
  pivot_longer(
    cols = -(biota_type:site),
    names_to = "Compound",
    values_to = "Detection_Frequency_Biota"
  )



### Compounds/DFs for Analysis ###################################################################################################
# Filter out compounds with less than 70% (67% for n = 3) detection at any site - Water
finalcmpds_w <- detect_freq_water_plot %>%
  group_by(Compound) %>% 
  filter(any(Detection_Frequency_Water > 66)) %>%
  ungroup() # Only keep compounds with detection frequencies > 67
finalcmpds_w$Compound <- factor(
  finalcmpds_w$Compound,
  levels = pfas_cols # Keeps compounds in original order
)

# Filter out compounds with less than 70% detection at any site - Biota
finalcmpds_b <- detect_freq_biota_plot %>%
  group_by(Compound) %>% 
  filter(any(Detection_Frequency_Biota > 70)) %>%
  ungroup() # Only keep compounds with detection frequencies > 70
finalcmpds_b$Compound <- factor(
  finalcmpds_b$Compound,
  levels = pfas_cols
)

# Lists of final compounds for analysis
cmpds_w <- unique(finalcmpds_w$Compound) # List of final 8 PFAS in water
cmpds_b <- unique(finalcmpds_b$Compound) # List of final 11 PFAS in biota
all_cmpds <- union(cmpds_w, cmpds_b) # List of final 14 PFAS in all samples
# Object with list of final 14 PFAS in order
pfas_cols_12 <- pfas_cols[pfas_cols %in% all_cmpds]  



### 3.2/3.3 Concentrations #######################################################################################################
### Impute values for samples <MDL for compounds with ≥70% detection frequency using MDL/√2
# Read in MDL sheet, replace with your own pathname
mdl_df <- read_excel('/Users/ally/Downloads/Senior Thesis/Processed Data/20260129_Hawaii_Thesis_MDL.xlsx') 
mdl_w <- mdl_df %>% filter(Media == "Water") # MDL values for water only
mdl_b <- mdl_df %>% filter(Media == "Biota") # MDL values for biota only
mdl_cols <- colnames(mdl_df)[2:65] # Names of PFAS columns with Br- and L-
# Replace <LOQ with MDL/sqrt(2)
full_data_imputed <- full_data %>%
  mutate(across(all_of(mdl_cols), ~{
    col_name <- cur_column()
    mdl_w_val <- as.numeric(mdl_w[[col_name]])
    mdl_b_val <- as.numeric(mdl_b[[col_name]])
    as.numeric(
      ifelse(. == "<LOQ",
             ifelse(type == "water",
                    mdl_w_val / sqrt(2),
                    mdl_b_val / sqrt(2)),
             .)
    )
  }))
# Sum Br- and L- isomers
full_data_isomers <- full_data_imputed %>%
  pivot_longer(all_of(mdl_cols), names_to = "compound", values_to = "conc") %>%
  mutate(compound = str_remove(compound, "^(Br-|L-)")) %>%
  group_by(across(-c(compound, conc)), compound) %>%
  summarise(conc = sum(conc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = compound, values_from = conc) %>%
  relocate(all_of(pfas_order))
# Keep only final 12 compounds
full_data_mdl <- full_data_isomers %>%
  select(-all_of(setdiff(pfas_order, all_cmpds))) %>%
  mutate(site = factor(site,
                       levels = c("PF", "WR", "SB", "OO", "MB")))

### Object to add PFAS class to DFs
pfclass_lookup <- data.frame(
  compound = c(
  # PFCAs (9)
    "PFEA", "PFBA", "PFPeA", "PFHpA", "PFOA", "PFNA", "PFDA", "PFUnDA", "PFDoDA",
  # PFSAs (2)
    "PFHxS", "PFOS",
  # Fluorotelomer (1)
    "8:2 FTUCA"),
  class = c(
    rep("PFCA", 9),
    rep("PFSA", 2),
    "Fluorotelomer"),
  stringsAsFactors = FALSE)
pfclass_lookup$compound <- factor(
  pfclass_lookup$compound,
  levels = pfas_cols_12)

### Combine PFAS matrix and metadata
full_data_long <- full_data_mdl %>% # Long format for graphing
  pivot_longer(cols = all_of(all_cmpds), 
               names_to = "compound", 
               values_to = "Concentration") %>%
  left_join(pfclass_lookup, by = "compound")  # Add PFAS class to DF
# DF without TFA
pfas_long <- full_data_mdl %>% # Long format for graphing
  pivot_longer(cols = all_of(all_cmpds), 
               names_to = "compound", 
               values_to = "concentration") %>%
  filter(compound != "PFEA") %>% # Removed TFA from analysis
  left_join(pfclass_lookup, by = "compound")  # Add PFAS class to DF

### Calculate concentrations for each compound, averaged by site & type (water/biota)
pfas_summary <- pfas_long %>%
  group_by(site, type, compound, class) %>%
  summarise(avg_conc = mean(concentration, na.rm = TRUE), .groups = "drop")
pfas_summary$compound <- factor(
  pfas_summary$compound,
  levels = pfas_cols_12)

### Water Concentration Summary Stats
# TFA summary stats
full_data_long %>% 
  filter(type == "water", compound == "PFEA") %>% # Replace "water" with "biota" for biota data
  summarise(
    n    = sum(!is.na(Concentration)),
    mean = ifelse(n > 0, mean(Concentration, na.rm = TRUE), NA),
    median = ifelse(n > 0, median(Concentration, na.rm = TRUE), NA),
    min  = ifelse(n > 0, min(Concentration, na.rm = TRUE), NA),
    max  = ifelse(n > 0, max(Concentration, na.rm = TRUE), NA)
  ) %>% 
  print()
# Summary stats by type without TFA
full_data_long %>% 
  filter(type == "water", compound != "PFEA") %>% 
  summarise(
    n    = sum(!is.na(Concentration)),
    mean = ifelse(n > 0, mean(Concentration, na.rm = TRUE), NA),
    median = ifelse(n > 0, median(Concentration, na.rm = TRUE), NA),
    min  = ifelse(n > 0, min(Concentration, na.rm = TRUE), NA),
    max  = ifelse(n > 0, max(Concentration, na.rm = TRUE), NA)
  ) %>% 
  print()
# Summary stats for individual sites
pfas_long %>% 
  filter(type == "water", site == "SB") %>% # Replace with site of interest
  summarise(avg_sum = sum(concentration)/3) %>% # Replace with appropriate number of water samples
  print()

### Biota Concentration Summary Stats
# Summary stats by type without TFA
full_data_long %>% 
  filter(type == "biota", compound != "PFEA") %>% 
  summarise(
    n    = sum(!is.na(Concentration)),
    mean = ifelse(n > 0, mean(Concentration, na.rm = TRUE), NA),
    median = ifelse(n > 0, median(Concentration, na.rm = TRUE), NA),
    min  = ifelse(n > 0, min(Concentration, na.rm = TRUE), NA),
    max  = ifelse(n > 0, max(Concentration, na.rm = TRUE), NA)
  ) %>% 
  print()
# Summary stats for individual sites
pfas_long %>% 
  filter(type == "biota", site == "MB") %>% # Replace with site of interest
  summarise(avg_sum = sum(concentration)/7) %>% # Replace with appropriate number of biota samples
  print()



### (pair-wise) PERMANOVA & PCA DFs & functions ##################################################################################
# Water matrix for PERMANOVA/PCA (compounds & concentrations only)
water_matrix <- full_data_mdl %>%
  filter(type == "water") %>% 
  select(all_of(all_cmpds)) %>%
  as.matrix()
# Water metadata (sample ID and site)
water_meta <- full_data_mdl %>%
  filter(type == "water") %>% 
  select(sample, site)

# Biota matrix for PERMANOVA/PCA (compounds & concentrations only)
biota_matrix <- full_data_mdl %>%
  filter(type == "biota") %>% 
  select(all_of(all_cmpds)) %>%
  as.matrix()
# Biota metadata (sample ID and site)
biota_meta <- full_data_mdl %>%
  filter(type == "biota") %>% 
  select(sample, site)
# Biota matrix for PERMANOVA/PCA excluding sample SB-B10
biota_matrix10 <- full_data_mdl %>%
  filter(type == "biota", sample != "SB-B10") %>% 
  select(all_of(all_cmpds)) %>%
  as.matrix()
# Biota metadata excluding SB-B10
biota_meta10 <- full_data_mdl %>%
  filter(type == "biota", sample != "SB-B10") %>% 
  select(sample, site)

# Water & biota matrix for PERMANOVA/PCA (compounds & concentrations only)
waterbiota_matrix <- full_data_mdl %>%
  select(all_of(all_cmpds)) %>%
  as.matrix()
# Water & biota metadata (sample ID and site)
waterbiota_meta <- full_data_mdl %>%
  select(sample, site, type, biota_type)

# pairwise PERMANOVA function, includes BH p-value adjustment for multiple testing
pairwise.adonis2 <- function(X, factors, permutations=999, p.adjust.method = "BH"){
  library(vegan)
  levels <- unique(factors)
  results <- list()
  for(i in 1:(length(levels)-1)){
    for(j in (i+1):length(levels)){
      subset_idx <- factors %in% c(levels[i], levels[j])
      ad <- adonis2(X[subset_idx,] ~ factors[subset_idx], permutations=permutations, method="bray")
      results[[paste(levels[i], "vs", levels[j])]] <- ad
    }
  }
  return(results)
}



### 3.4 Comparing water composition between sites ################################################################################
### PERMANOVA - not included in report b/c sample size is too small
adonis2(water_matrix ~ site,
        data = water_meta,
        method = "bray")
# Check dispersion
dist_matrix_w <- vegdist(water_matrix, method = "bray")
disp <- betadisper(dist_matrix_w, water_meta$site)
anova(disp)

### PCA
water_matrix_filtered <- water_matrix[, apply(water_matrix, 2, sd) > 0]
pca_w <- prcomp(water_matrix_filtered,
                center = TRUE,
                scale. = TRUE)
scores_w <- as.data.frame(pca_w$x)
scores_w$site <- water_meta$site
loadings_w <- as.data.frame(pca_w$rotation)
loadings_w$Compound <- rownames(loadings_w)
var_explained_w <- summary(pca_w)$importance[2, ] * 100



### 3.5 Comparing biota composition between sites ################################################################################
### PERMANOVA
adonis2(biota_matrix ~ site,
        data = biota_meta,
        method = "bray")
# Excluding SB-B10
adonis2(biota_matrix10 ~ site,
        data = biota_meta10,
        method = "bray")
# Check dispersion
dist_matrix_b <- vegdist(biota_matrix, method = "bray")
disp_b <- betadisper(dist_matrix_b, biota_meta$site)
anova(disp_b)
# Check dispersion excluding SB-B10
dist_matrix_b10 <- vegdist(biota_matrix10, method = "bray")
disp_b10 <- betadisper(dist_matrix_b10, biota_meta10$site)
anova(disp_b10)

### PCA
biota_matrix_filtered <- biota_matrix[, apply(biota_matrix, 2, sd) > 0]
pca_b<- prcomp(biota_matrix_filtered,
               center = TRUE,
               scale. = TRUE)
scores_b <- as.data.frame(pca_b$x)
scores_b$site <- biota_meta$site
loadings_b <- as.data.frame(pca_b$rotation)
loadings_b$Compound <- rownames(loadings_b)
var_explained_b <- summary(pca_b)$importance[2, ] * 100
# PCA without SB-B10
biota_matrix_filtered10 <- biota_matrix10[, apply(biota_matrix10, 2, sd) > 0]
pca_b10<- prcomp(biota_matrix_filtered10,
                 center = TRUE,
                 scale. = TRUE)
scores_b10 <- as.data.frame(pca_b10$x)
scores_b10$site <- biota_meta10$site
loadings_b10 <- as.data.frame(pca_b10$rotation)
loadings_b10$Compound <- rownames(loadings_b10)
var_explained_b10 <- summary(pca_b10)$importance[2, ] * 100

### Pairwise PERMANOVA
pairwise_results_b <- pairwise.adonis2(biota_matrix, biota_meta$site)
# Pairwise PERMANOVA without SB-B10
pairwise_results_b10 <- pairwise.adonis2(biota_matrix_filtered10, biota_meta10$site)



### 3.6 Comparing composition between sites (water & biota) ######################################################################
### Visualize concentration distributions to choose best test for difference
# Water histogram
ggplot(water_matrix, aes(x = PFBA)) +
  geom_histogram(bins = 30) +
  theme_minimal() # Not normally distributed -> Kruskal-Wallis
# Biota histogram
ggplot(biota_matrix, aes(x = PFHpA)) +
  geom_histogram(bins = 20) +
  theme_minimal() # Not normally distributed -> Kruskal-Wallis

### Water individual compounds Kruskal-Wallis & BH p-value adjustment
# Filter to water samples
water_data <- full_data_mdl %>%
  filter(type == "water", site != "WR") # Exclude WR (n = 1)
# Run Kruskal-Wallis for each compound
kw_list_w <- lapply(all_cmpds, function(compound) {
  formula <- as.formula(paste0("`", compound, "` ~ site"))
  test <- kruskal.test(formula, data = water_data)
  data.frame(
    Compound = compound,
    statistic = as.numeric(test$statistic),
    df = as.numeric(test$parameter),
    p_value = test$p.value
  )
})
kw_results_w <- do.call(rbind, kw_list_w)
# Adjust p-values (BH method)
kw_results_w$p_adj_BH   <- p.adjust(kw_results_w$p_value, method = "BH")
kw_results_w

### Biota individual compounds Kruskal-Wallis & BH p-value adjustment
# Filter to biota samples
biota_data <- full_data_mdl %>%
  filter(type == "biota", site != "OO")
# Run Kruskal–Wallis for each compound
kw_list_b <- lapply(all_cmpds, function(compound) {
  formula <- as.formula(paste0("`", compound, "` ~ site"))
  test <- kruskal.test(formula, data = biota_data)
  data.frame(
    Compound = compound,
    statistic = as.numeric(test$statistic),
    df = as.numeric(test$parameter),
    p_value = test$p.value
  )
})
kw_results_b <- do.call(rbind, kw_list_b)
# Adjust p-values (BH method)
kw_results_b$p_adj_BH   <- p.adjust(kw_results_b$p_value, method = "BH")
kw_results_b

### Biota Dunn-Bonferroni Post Hoc
library(FSA)
# List compounds with significant differences from Kruskal-Wallis testing
sig_cmpds_b <- c("PFEA", "PFOA", 
                 "PFOS", "PFNA", 
                 "PFDA", "PFDoDA")
# Test run (run first)
dunnTest(PFOS ~ site,
         data = full_data_mdl %>% filter(type == "biota"),
         method = "bonferroni") 
# Dunn-Bonferroni post hoc
dunn_results <- lapply(sig_cmpds_b, function(compound) {
  formula <- as.formula(paste0("`", compound, "` ~ site"))
  test <- dunnTest(formula,
                   data = biota_data,
                   method = "bonferroni")
  out <- test$res
  out$Compound <- compound
  out
})
dunn_results <- do.call(rbind, dunn_results)
dunn_results

### Water PFAS class Kruskal-Wallis & BH p-value adjustment
library(rstatix)
# Summarize into PFAS classes
class_totals_w <- pfas_long %>%
  filter(type == "water") %>%
  group_by(sample, site, class) %>%
  summarise(class_conc_w = sum(concentration, na.rm = TRUE), .groups = "drop")
# Calculate fraction of composition for each PFAS class
class_frac_w <- class_totals_w %>%
  group_by(sample) %>%
  mutate(fraction = class_conc_w / sum(class_conc_w)) %>%
  ungroup()
# Kruskal-Wallis with BH p-value adjustment
kw_results_wclass <- class_frac_w %>%
  group_by(class) %>%
  summarise(
    p_value = kruskal.test(fraction ~ site)$p.value
  ) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
kw_results_wclass

### Biota PFAS class Kruskal-Wallis & BH p-value adjustment
# Summarize into PFAS classes
class_totals_b <- pfas_long %>%
  filter(type == "biota") %>%
  group_by(sample, site, class) %>%
  summarise(class_conc_b = sum(concentration, na.rm = TRUE), .groups = "drop")
# Calculate fraction of composition for each PFAS class
class_frac_b <- class_totals_b %>%
  group_by(sample) %>%
  mutate(fraction = class_conc_b / sum(class_conc_b)) %>%
  ungroup()
# Kruskal-Wallis with BH p-value adjustment
kw_results_bclass <- class_frac_b %>%
  group_by(class) %>%
  summarise(
    p_value = kruskal.test(fraction ~ site)$p.value
  ) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
kw_results_bclass

# Dunn-Bonferroni post hoc
dunn_results_b <- class_frac_b %>%
  group_by(class) %>%
  dunn_test(fraction ~ site, p.adjust.method = "bonferroni")
dunn_results_b



### 3.7 BAFs #####################################################################################################################
### New DF with all info
pfas_data <- full_data_detect[, 14:70]
# Add metadata columns
pfas_df <- cbind(
  pfas_data,
  site = full_data_detect$site,
  type = full_data_detect$type,
  biota_type = full_data_detect$biota_type,
  sample = full_data_detect$sample
)
# Convert to data.frame if needed
pfas_df <- as.data.frame(pfas_df)
# Keep only PFAS compounds in pfas_cols_12 + metadata
pfas_df <- full_data_detect[, c(pfas_cols_12, "site", "type", "biota_type", "sample")]
# Replace 0 with NA in the PFAS columns
pfas_df[, pfas_cols_12][pfas_df[, pfas_cols_12] == 0] <- NA

# Select only water samples
water_baf <- pfas_df %>%
  filter(type == "water") %>%
  group_by(site) %>%
  summarise(across(all_of(pfas_cols_12), mean, na.rm = TRUE)) %>%
  ungroup()

# Select only biota samples
biota_baf <- pfas_df %>% filter(type == "biota") %>% 
  mutate(across(all_of(pfas_cols_12), ~ . * 1000))

# Join average water concentrations by site
baf_df <- biota_baf %>%
  left_join(water_baf, by = "site", suffix = c("", "_water")) %>%
  # Divide each PFAS concentration by corresponding water average
  rowwise() %>%
  mutate(across(all_of(all_cmpds), 
                ~ .x / get(paste0(cur_column(), "_water")), 
                .names = "{.col}_BAF")) %>%
  ungroup()

# Convert to long for plotting
baf_long <- baf_df %>%
  select(sample, site, biota_type, ends_with("_BAF")) %>%
  filter(site != "OO") %>% 
  pivot_longer(cols = ends_with("_BAF"),
               names_to = "compound",
               values_to = "BAF") %>%
  mutate(compound = sub("_BAF$", "", compound),
         compound = factor(compound, levels = pfas_cols_12)) %>%
  group_by(compound) %>%
  filter(!all(is.na(BAF))) %>%
  ungroup()

### One-way ANOVA: Difference in BAFs between biota types?
baf_long_log <- baf_long %>%
  filter(BAF > 0) %>% 
  mutate(log10_BAF = log10(BAF))
kruskal.test(log10_BAF ~ biota_type, data = baf_long_log)



### 4.2 Health Guidelines ################################################################################################
# PFOA and PFOS summaries for WR
full_data_mdl %>% 
  filter(site == "WR", # Fill in site of interest
         type == "water") %>% 
  select(PFOS, PFOA)

### Fish consumption summary
# Calculated local weekly exposure
full_data_mdl %>% 
  filter(site == "PF", # Fill in site
         biota_type == "fish") %>% 
  summarise(avg_pfos = mean(PFOS), 
            avg_pfoa = mean(PFOA),
            avg_pfba = mean(PFBA),
            avg_pfna = mean(PFNA),
            avg_pfda = mean(PFDA),
            avg_pfhxs = mean(PFHxS),
            avg_ftuca = mean(`8:2 FTUCA`)) %>% 
  summarise(guide_pfos = avg_pfos*73.34/70, # Replace 73.34 with 248.83 for OO calculations
            guide_pfoa = avg_pfoa*73.34/70,
            guide_pfba = avg_pfba*73.34/70,
            guide_pfna = avg_pfna*73.34/70,
            guide_pfda = avg_pfda*73.34/70,
            guide_pfhxs = avg_pfhxs*73.34/70,
            guide_ftuca = avg_ftuca*73.34/70) 
# Calculated 1 meal exposure
full_data_mdl %>% 
  filter(site == "PF", # Fill in site
         biota_type == "fish") %>% 
  summarise(avg_pfos = mean(PFOS), 
            avg_pfoa = mean(PFOA),
            avg_pfba = mean(PFBA),
            avg_pfna = mean(PFNA),
            avg_pfda = mean(PFDA),
            avg_pfhxs = mean(PFHxS),
            avg_ftuca = mean(`8:2 FTUCA`)) %>% 
  summarise(.guide_pfos = avg_pfos*277/70,
            .guide_pfoa = avg_pfoa*277/70,
            .guide_pfba = avg_pfba*277/70,
            .guide_pfna = avg_pfna*277/70,
            .guide_pfda = avg_pfda*277/70,
            .guide_pfhxs = avg_pfhxs*277/70,
            .guide_ftuca = avg_ftuca*277/70)



