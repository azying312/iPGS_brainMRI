fullargs <- commandArgs(trailingOnly = FALSE)
args <- commandArgs(trailingOnly = TRUE)
print(paste("fullargs:", fullargs))
print(paste("args:", args))
script_name <- normalizePath(
  sub("--file=", "", fullargs[grep("--file=", fullargs)])
)
script_dir <- dirname(script_name)
suppressWarnings(suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
}))
print(paste("script name:", script_name))
print(paste("dir name:", script_dir))
####################################################################
source(file.path(dirname(dirname(script_dir)), 'paths.sh'))
####################################################################

# from 2_define_phenotypes.R: save image - pop_counts.RData
#load(file.path(script_dir, "pop_counts.RData"))

####################################################################
# output files
####################################################################
svd_results_dir <- file.path(script_dir, "svd_results")
scree_plot_name <- file.path(script_dir, "scree_plot.png")
svd_plot_name12_pop <- file.path(script_dir, "svd1_2_pop.png")
svd_plot_name12 <- file.path(script_dir, "svd1_2.png")
svd_plot_name23 <- file.path(script_dir, "svd2_3.png")
svd_plot_name34 <- file.path(script_dir, "svd3_4.png")
svd_plot_name45 <- file.path(script_dir, "svd4_5.png")
svd_plot_name56 <- file.path(script_dir, "svd5_6.png")
svd_plot_name67 <- file.path(script_dir, "svd6_7.png")
svd_plot_name78 <- file.path(script_dir, "svd7_8.png")

####################################################################
# input files
####################################################################
keep_file <- "/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942/sqc/population/train.keep"
sqc_f <- file.path(ukb21942_d, "sqc", "sqc.20220316.tsv.gz")
fields_f <- file.path(script_dir, "misc", "mri_traits_449_names.tsv") #"fields.tsv")
#SQL_dump_f <- file.path("${ukb21942_user}/pheno/tmp/${analysis_name}/tab_dump.tsv.gz")
SQL_dump_f <- file.path(ukb21942_user, "pheno", "tmp", basename(script_dir), "tab_dump.tsv.gz")

phe_f <- file.path(
    ukb21942_user, "pheno",
    "brain_mri.tsv.gz"
)

fields_f <- file.path(
    script_dir, "misc",
    "mri_traits_449_names.tsv"
)

phe_name <- dirname(script_dir)
phe_long_f   <- file.path(phe_name, "brain_mri.long.tsv.gz")
phe_counts_f <- file.path(phe_name, "brain_mri.counts.tsv.gz")

####################################################################
# main
####################################################################

sqc_f %>%
  fread(
    select = c("#FID", "IID", "population", "population_MM", "pass_QC_filter"),
    colClasses = c('#FID'='character', 'IID'='character')
  ) %>%
  rename_with(
    function(x){str_replace(x, '#', '')}, starts_with("#")
  ) -> sqc_df

# set missing population to NA
sqc_df <- sqc_df %>%
   mutate(population=ifelse(is.na(population), "NA", population))

print("sqc_df")
print(table(sqc_df$pass_QC_filter))
print(table(sqc_df$population, sqc_df$pass_QC_filter))

fields_f %>%
  fread() %>%
  rename_with(~str_replace(., "#", ""), starts_with("#")) -> fields_df

fields_df %>%
   rename("GBE_ID" = "Field.ID", "GBE_NAME" = "Description") %>%
   mutate(
      GBE_ID = paste0("INI", GBE_ID),
      GBE_short_name = GBE_NAME,
      GBE_short_name = str_replace(GBE_short_name, "Volume", "Vol."),
      GBE_short_name = str_replace(GBE_short_name, "\\(left\\)", "L"),
      GBE_short_name = str_replace(GBE_short_name, "\\(right\\)", "R"),
      GBE_short_name = str_replace(GBE_short_name, "percentage", "%"),
      GBE_short_name = str_replace(GBE_short_name, "Impedance", "Impd."),) -> field_names_df

SQL_dump_f %>%
  fread(colClasses = c("#eid" = "character", "value" = "character")) %>%
  rename("IID" = "#eid") -> tab_dump_df

####################################################################

# filter out QC terms from field names
# phenotypes
field_names_df <- field_names_df %>%
   filter(GBE_ID != "INI25000") %>%
   filter(GBE_ID != "INI25734")

tab_dump_df <- tab_dump_df %>%
   mutate(GBE_ID=paste0("INI", field))

# sqc_df: IID  and population; 488377 - total number of ppl from our data
# tab dump: IID, phenotype,  GBE ID, value; 22225418 - each person and trait

print(sum(is.na(tab_dump_df$value)))
print(sum(tab_dump_df$value == 0))
# nothing has no value, 3268 have 0 value

print("merge sqc and tab dumb")
merged.df <- sqc_df %>%
   left_join(tab_dump_df, by = "IID")
#merged.df <- merged.df %>%
#   mutate(population=ifelse(is.na(population), "NA", population))

table(merged.df$pass_QC_filter)
table(merged.df$population)

# 488377 individual and 450 phenotypes

population_counts <- merged.df %>%
  #filter(pass_QC_filter) %>% # 21553397 original
  #filter(!is.na(value)) %>%
  #filter(!is.na(population)) %>%
  filter(GBE_ID != "INI25000") %>% # head size
  filter(GBE_ID != "INI25734") %>% # Inverted_signal-to-noise_ratio_in_T1
  mutate(
    phenotype_group = ifelse(GBE_ID %in% field_names_df$GBE_ID, "has_brain_measure", "other_measure")
  )

ukbb_summary <- population_counts %>%
   #filter(pass_QC_filter == TRUE) %>%
   #filter(!is.na(value)) %>%
   #filter(!is.na(population)) %>%
   group_by(IID, population, pass_QC_filter) %>%
   summarise(hasBM = any(phenotype_group == "has_brain_measure")) %>%
   ungroup()

print("ukbb")
print(tail(ukbb_summary))
print(table(ukbb_summary$pass_QC_filter))
print("qc filter, hasBM")
print(table(ukbb_summary$pass_QC_filter, ukbb_summary$hasBM))
print("population v. brain measure table")
print(table(ukbb_summary$population, ukbb_summary$hasBM))

# Keep ID
keep_data <- read.table(keep_file, header = FALSE, sep = "", stringsAsFactors = FALSE)
keepID <- keep_data[,1]

phe_f %>%
   fread() -> pheno_df
print("pheno")
dim(pheno_df)

# pop counts
phe_counts_f %>%
   fread() -> phe_counts_df

phe_counts_df <- phe_counts_df %>%
   select("#GBE_ID", "total") %>%
   rename(GBE_ID = "#GBE_ID") %>%
   filter(GBE_ID != "INI25000") %>%
   filter(GBE_ID != "INI25734")

# Filter for BM traits
brain_measure_ids <- fields_df %>%
  mutate(GBE_ID = paste0("INI", Field.ID)) %>%
  filter(GBE_ID != "INI25000", GBE_ID != "INI25734") %>%
  pull(GBE_ID)

bm_traits_in_pheno <- intersect(brain_measure_ids, names(pheno_df))
print("bm_traits_in_pheno")
length(bm_traits_in_pheno)
length(names(pheno_df))

print("bm traits")
print(dim(pheno_df))

pheno_df <- pheno_df %>%
  filter(if_any(all_of(bm_traits_in_pheno), ~ !is.na(.)))
print(dim(pheno_df))

summary(rowSums(is.na(pheno_df[, ..bm_traits_in_pheno])))
summary(colSums(is.na(pheno_df[, ..bm_traits_in_pheno])))

train_keepID <- intersect(pheno_df$IID, keepID)
print(length(train_keepID))
print(length(keepID))

########################################################################

## Normalize volumetric traits
print("traits_to_norm")
normalized_traits <- c("INI25001", "INI25003", "INI25005", "INI25007", "INI25009")

print(head(fields_df))

traits_to_norm <- field_names_df %>%
   filter(GBE_ID != "INI25000", GBE_ID != "INI25734")
   #filter(grepl("Vol\\.|volume", Description, ignore.case = TRUE)) %>%
   #filter(GBE_ID != "INI25000") #%>%
   #filter(!(GBE_ID %in% normalized_traits))
traits_to_norm <- traits_to_norm %>%
   filter(!(GBE_ID %in% normalized_traits))

norm_traits <- c("INI25000", normalized_traits, unique(traits_to_norm$GBE_ID))
norm_pheno_df <- pheno_df %>%
   select(all_of(c("#FID", "IID", norm_traits)))

filter_traits <- c("INI25000", unique(traits_to_norm$GBE_ID))

norm_pheno_df <- norm_pheno_df %>% 
   select(all_of(c("#FID", "IID", filter_traits)))

norm_pheno_df <- norm_pheno_df %>%
   mutate(across(-c("#FID", "IID", "INI25000"), ~ ifelse(!is.na(.), . * INI25000, NA)))

print("norm df")
print(dim(norm_pheno_df))

# remake pheno_df
pheno_df <- pheno_df %>%
   select(-all_of(filter_traits))
dim(filter_traits)

pheno_df <- pheno_df %>%
   left_join(norm_pheno_df, by=c("#FID", "IID")) %>%
   select(-INI25000)

print(dim(pheno_df))

ini_cols <- grep("^INI", names(pheno_df), value = TRUE)
length(ini_cols)

########################################################################

# Test set
test_set <- pheno_df %>%
   filter(!(IID %in% keepID))
print("test set")
print(dim(test_set))
write.csv(test_set, file.path(script_dir, "test_set.csv"))

# Filter for training set
pheno_df <- pheno_df %>%
   filter(IID %in% keepID)
print("train")
print(dim(pheno_df))

norm_cols <- c("INI25001", "INI25003", "INI25005", "INI25007", "INI25009","INI25002", "INI25004", "INI25006", "INI25008", "INI25010")

pop_counts <- pheno_df %>%
   summarise(across(all_of(ini_cols),  ~ sum(!is.na(.))))

# Center and scale the selected columns
scaling_params <- pheno_df %>%
  summarise(across(all_of(ini_cols), ~ list(mean = mean(., na.rm=TRUE), sd = sd(., na.rm=TRUE)))) %>%
  unnest(cols=everything()) %>%
  as.data.frame()
rownames(scaling_params) <- c("mean", "sd")

#print(scaling_params[, 1:5])
dim(scaling_params)

# save scaling params
write.csv(scaling_params,"scaling_params.csv", row.names = TRUE)

scaling_params <- pheno_df %>%
  summarise(across(all_of(ini_cols), ~ list(mean = mean(., na.rm=TRUE), sd = sd(., na.rm=TRUE))))

phe_scaled <- pheno_df %>%
  mutate(across(all_of(ini_cols), ~ (.-scaling_params[[cur_column()]][["mean"]]) / scaling_params[[cur_column()]][["sd"]]))

# Dealing with NA
na_count <- sum(is.na(phe_scaled %>% select(all_of(ini_cols))))
nan_count <- sum(is.nan(as.matrix(phe_scaled %>% select(all_of(ini_cols)))))
inf_count <- sum(is.infinite(as.matrix(phe_scaled %>% select(all_of(ini_cols)))))
print(paste("Count of NA values before imputation:", na_count, "\n"))
print(paste("Count of NaN values before imputation:", nan_count, "\n"))
print(paste("Count of Inf values before imputation:", inf_count, "\n"))

# total count
total_count <- prod(dim(phe_scaled %>% select(all_of(ini_cols))))
print(paste("Count of total values before imputation:", total_count, "\n"))

na_vals <- phe_scaled %>%
  summarise(across(all_of(ini_cols), ~ sum(is.na(.))))
#print(unique(na_vals))

# Impute with col means
phe_imputed <- phe_scaled %>%
   mutate(across(all_of(ini_cols), ~ replace(., is.na(.), mean(., na.rm=TRUE))))

# Remove columns with all NA
all_na <- phe_imputed %>%
  select(where(~ all(is.na(.))))
#print(paste("All NA col number:", dim(all_na)))

phe_imputed <- phe_imputed %>%
  select(where(~ !all(is.na(.))))

ini_cols <- grep("^INI", names(phe_imputed), value = TRUE)

fields_df <- field_names_df %>%
   filter(GBE_ID %in% ini_cols) %>%
   rename(Phenotype=GBE_ID)

pop_counts <- phe_imputed %>%
   summarise(across(all_of(ini_cols),  ~ sum(!is.na(.))))

print("pop_counts")
print(head(pop_counts[, 1:5]))

scaling_metadata <- list(
  original_summary = pheno_df %>%
    summarise(across(all_of(ini_cols), ~ list(mean = mean(., na.rm=TRUE), sd = sd(., na.rm=TRUE), min = min(., na.rm=TRUE), max = max(., na.rm=TRUE)))),
  scaling_params = scaling_params,
  missing_counts = list(
    na_count = na_count,
    nan_count = nan_count,
    inf_count = inf_count
  ),
  removed_columns = names(all_na),
  imputation_method = "Mean imputation"
)

saveRDS(scaling_metadata, "scaling_metadata.rds")

save.image("cleaned_svd.RData")

####################################################################
# SVD
####################################################################

# SVD
phe_matrix <- phe_imputed %>%
  select(all_of(ini_cols)) %>%
  as.matrix()

# SVD calculation
svd_result <- svd(phe_matrix, nu = 8, nv = 8)
v_matrix <- svd_result$v[, 1:8]
d_vec <- svd_result$d
svd_v <- d_vec^2/sum(d_vec^2)
scree_df <- data.frame(
   Component=c(1:length(svd_v)),
   PVE=svd_v,
   PVE_sum=cumsum(svd_v)
)
write.csv(svd_result$v, "training_svd_v_mat.csv")
write.csv(scree_df, "training_scree_components.csv")
if (!dir.exists(svd_results_dir)) {
  dir.create(svd_results_dir)}

training_projection <- phe_matrix %*% v_matrix

print("training_projection")
print(dim(training_projection))
print(head(training_projection[, 1:4]))
print("phe_matrix")
print(dim(phe_matrix))
print(head(phe_matrix[, 1:4]))
print("v_matrix")
print(dim(v_matrix))
print(head(v_matrix[, 1:4]))

### Plotting
# Scree plot
thresholds=c(0.5, 0.6, 0.8, 0.9, 0.95)
thres_idx <- sapply(thresholds, function(threshold) {which(scree_df$PVE_sum >= threshold)[1] })
scree_plt <- ggplot(scree_df, aes(x=Component, y=PVE)) +
   geom_bar(stat="identity", fill="orchid") +
   geom_line(aes(group=1), color="red", size=1) +
   geom_point(color="black", size=0.5) +
   geom_vline(xintercept=thres_idx,
        linetype="solid", color="blue") +
   labs(
    title = "Scree Plot",
    x = "Singular Value",
    y = "Proportion of Variance Explained"
  ) +
  theme_minimal()
ggsave(scree_plot_name, scree_plt, width = 8, height = 6, dpi=300)

#############################################################################
# Create v_df
create_vdf <- function(comp1, comp2) {
   print("create vdf")
   v_df <- data.frame(Phenotype = ini_cols, SVD1 = v_matrix[, comp1], SVD2 = v_matrix[, comp2])
   v_df <- v_df %>%
      left_join(fields_df, by="Phenotype")
   v_df <- v_df %>%
      filter(Phenotype != "INI25000")
   return(v_df)
}
print("12")
v_df <- create_vdf(1,2)

#print(head(v_df))
print(v_df[which(is.na(v_df$Category)),]) # INI25734 is the inverted signal to noise, measurement factor
print(which(is.na(v_df$Category)))

pve_svd <- round(scree_df$PVE * 100, 2)
pve_svd1 <- round(scree_df$PVE[1] * 100, 2)
pve_svd2 <- round(scree_df$PVE[2] * 100, 2)
pve_svd3 <- round(scree_df$PVE[3] * 100, 2)

#############################################################################

#############################################################################
# Components 1 and 2

pve1<-pve_svd[1]
pve2<-pve_svd[2]
svd_plot <- ggplot(v_df, aes(x = SVD1, y = SVD2, color = Category)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid") +
   geom_vline(xintercept=0, color="black", linetype="solid") +
   geom_text(aes(label = Phenotype), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD1 (", pve1, "%)"),
        y = paste0("SVD2 (", pve2, "%)"),
        colors="MRI Type") +
   theme_minimal()
#   guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))
ggsave(svd_plot_name12, svd_plot, width = 8, height = 6, dpi=300)

#############################################################################
# population
#############################################################################

# Components 1 and 2
svd1 <- training_projection[, 1]
svd2 <- training_projection[, 2]
svd3 <- training_projection[, 3]
svd4 <- training_projection[, 4]
svd5 <- training_projection[, 5]
svd6 <- training_projection[, 6]
svd7 <- training_projection[, 7]
svd8 <- training_projection[, 8]
plot_df <- data.frame(
  IID = pheno_df$IID,
  SVD1 = svd1,
  SVD2 = svd2,
  SVD3 = svd3,
  SVD4 = svd4,
  SVD5 = svd5,
  SVD6 = svd6,
  SVD7 = svd7,
  SVD8 = svd8
)
sqc_df$IID <- as.numeric(sqc_df$IID)
print("pop counts")
print(dim(sqc_df))
print(dim(plot_df))

plot_df <- plot_df %>%
   left_join(sqc_df, by="IID")
pve1<-pve_svd[1]
pve2<-pve_svd[2]

dim(plot_df)

print(head(plot_df[,2:9]))

print(table(duplicated(plot_df[,2:9])))

svd_plot <- ggplot(plot_df, aes(x = SVD1, y = SVD2, color = population)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid") +
   geom_vline(xintercept=0, color="black", linetype="solid") +
   #geom_text(aes(label = population), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD1 (", pve1, "%)"),
        y = paste0("SVD2 (", pve2, "%)") ,
        colors="Population") +
   theme_minimal()
ggsave(svd_plot_name12_pop, svd_plot, width = 8, height = 6, dpi=300)
##########################
print("done")

