
fullargs <- commandArgs(trailingOnly = FALSE)
args <- commandArgs(trailingOnly = TRUE)

script_name <- normalizePath(
    sub("--file=", "", fullargs[grep("--file=", fullargs)])
)
script_dir <- dirname(script_name)

repo_dir <- dirname(dirname(script_dir))
analysis_name <- basename(script_dir)
print(paste("analysis_name:", analysis_name))

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

print(paste("start source", repo_dir, "paths.sh"))
source(file.path(repo_dir, 'paths.sh'))
print(paste("start source", repo_dir, "/helpers/functions.R"))
source(file.path(repo_dir, 'helpers', '/functions.R'))
print("finish source")

####################################################################
# output files
####################################################################

print(paste("script_dir:", script_dir))
pc_plot_name <- file.path(script_dir, "pc_plot.png")
#sample_wts_f <- file.path(
#    script_dir, "sample_wts_train_v2.rds")

#print(paste("sample_wts_f", sample_wts_f))

# all 1s
#dummy_wts_f <- file.path(
#        script_dir, "dummy_wts_train_v2.rds")

####################################################################
# input files
####################################################################

phe_f <- file.path(
    ukb21942_user, "pheno",
    "brain_mri.tsv.gz"
)

phe_name <- dirname(script_dir)
print(paste("phe_name:", phe_name))
phe_long_f   <- file.path(phe_name, "brain_mri.long.tsv.gz")

####################################################################
# main
####################################################################

# Read SQC file and clean data
sqc_f %>%
  fread(select = c("#FID", "IID", "population", "PC1", "PC2", "split")) %>%
  rename_with(~str_replace(., '#', ''), starts_with("#")) %>%
  mutate(
    across(all_of(c("FID", "IID")), as.character),
    population = factor(population, levels = c("WB", "NBW", "Afr", "SA", "others"))
  ) -> sqc_df

sqc_df %>%
  filter(split == "train") -> train_sqc_df

phe_long_f %>%
  fread() %>%
  rename(FID = `#FID`) %>%
  mutate(
    across(all_of(c("FID", "IID")), as.character)) %>%
  inner_join(train_sqc_df, by = c("FID", "IID")) %>%
  select(FID, IID, population, PC1, PC2, GBE_ID, value) -> train_with_phenotypes_long

train_with_phenotypes_long %>%
  filter(population == "SA") -> train_with_phenotypes_long

print(head(train_with_phenotypes_long))
print(tail(train_with_phenotypes_long))

print("train with phenotypes summary")

train_with_phenotypes_long %>%
  group_by(GBE_ID) %>%
  summarise(
    PC1_raw_mean = mean(PC1, na.rm = TRUE),
    PC1_raw_sd = sd(PC1, na.rm = TRUE),
    PC2_raw_mean = mean(PC2, na.rm = TRUE),
    PC2_raw_sd = sd(PC2, na.rm = TRUE)
  ) %>%
  ungroup()

# distribution of phenotype data
summary(train_with_phenotypes_long$value)

## PCA

# Center & Scale
train_with_phenotypes_long %>%
  mutate(log_value = ifelse(value > 0, log10(value + 1), NA)) %>%
  group_by(GBE_ID) %>%
  mutate(log_value_scaled = scale(log_value, center = TRUE, scale = TRUE)[, 1],
        PC1_scale=scale(PC1, center = TRUE, scale = TRUE)[, 1],
	PC2_scale=scale(PC2, center = TRUE, scale = TRUE)[, 1]
  ) %>%
  ungroup() -> scaled_train_with_phenotypes_long

print(head(scaled_train_with_phenotypes_long))
print(tail(scaled_train_with_phenotypes_long))

scaled_train_with_phenotypes_long %>%
  group_by(GBE_ID) %>%
  summarise(
    PC1_mean = mean(PC1_scale, na.rm = TRUE),
    PC1_sd = sd(PC1_scale, na.rm = TRUE),
    PC2_mean = mean(PC2_scale, na.rm = TRUE),
    PC2_sd = sd(PC2_scale, na.rm = TRUE)
  ) %>%
  ungroup()

print("variance of pca components")
summary(scaled_train_with_phenotypes_long$PC1_scale)
summary(scaled_train_with_phenotypes_long$PC2_scale)

scaled_train_with_phenotypes_long <- scaled_train_with_phenotypes_long %>%
  mutate(GBE_ID = as.factor(GBE_ID))

print(unique(scaled_train_with_phenotypes_long$GBE_ID))

scaled_train_with_phenotypes_long %>%
  count(GBE_ID)

scaled_train_with_phenotypes_long %>%
  group_by(GBE_ID) %>%
  summarise(unique_PC1 = n_distinct(PC1_scale), unique_PC2 = n_distinct(PC2_scale)) 

pc_plot <- scaled_train_with_phenotypes_long %>%
  mutate(GBE_ID = as.factor(GBE_ID)) %>%
  ggplot(aes(x=PC1_scale, y=PC2_scale, color=GBE_ID))+
  geom_point(alpha = 0.7, size = 2) +
  geom_jitter(alpha = 0.7, size = 2, width = 0.05, height = 0.05) +
  theme_minimal() +
  labs(title="PC1 v. PC2 Phenotypes", xlab="PC1", ylab="PC2", color="Phenotype")+
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

pc_plot_1 <- scaled_train_with_phenotypes_long %>%
  filter(GBE_ID %in% c("INI25014", "INI25781")) %>%
  ggplot(aes(x = PC1_scale, y = PC2_scale, color = GBE_ID)) +
  geom_jitter(width = 0.05, height = 0.05, alpha = 0.7, size = 2) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(
    title = "Subset: PC1 vs PC2 Phenotypes", 
    x = "PC1", 
    y = "PC2", 
    color = "Phenotype"
  )

print(paste("pc_plot_name", pc_plot_name))
ggsave(pc_plot_name, pc_plot, width = 8, height = 6)

pc1 <- ggplot(scaled_train_with_phenotypes_long, aes(x = PC1_scale)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of PC1", x = "PC1", y = "Count")

pc2 <- ggplot(scaled_train_with_phenotypes_long, aes(x = PC2_scale)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of PC2", x = "PC2", y = "Count")
ggsave("pc1_distb.png", pc1, width = 8, height = 6)
ggsave("pc2_distb.png", pc2, width = 8, height = 6)
