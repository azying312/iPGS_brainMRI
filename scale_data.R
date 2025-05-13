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

####################################################################
# input files
####################################################################
keep_file <- "/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942/sqc/population/train.keep"

phe_f <- file.path(
    ukb21942_user, "pheno",
    "brain_mri.tsv.gz"
)


fields_f <- file.path(
    script_dir, "misc",
    "mri_traits_449_names.tsv"
)

print(paste("fields_f:", fields_f))

phe_name <- dirname(script_dir)
print(paste("phe_name:", phe_name))
phe_long_f   <- file.path(phe_name, "brain_mri.long.tsv.gz")
phe_counts_f <- file.path(phe_name, "brain_mri.counts.tsv.gz")

####################################################################
# main
####################################################################
#genotype.pfile = "/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942/geno/ukb_genoHM3/ukb_genoHM3"
#ids <- list()
#ids[['psam']] <- readIDsFromPsam(sprintf('%s.psam', genotype.pfile))

keep_data <- read.table(keep_file, header = FALSE, sep = "", stringsAsFactors = FALSE)
keepID <- keep_data[,1]

print(length(keepID))

phe_f %>%
   fread() -> pheno_df

# pop counts
phe_counts_f %>%
   fread() -> phe_counts_df

phe_counts_df <- phe_counts_df %>%
   select("#GBE_ID", "total") %>%
   rename(GBE_ID = "#GBE_ID") #%>%
   #filter(GBE_ID != "INI25000") %>%
   #filter(GBE_ID!="INI25734")

fields_df <- read.table(fields_f, sep='\t', header=TRUE)

# Filter out QC traits
fields_df <- fields_df %>%
   mutate(GBE_ID=paste0("INI", Field.ID))# %>%
   #filter(GBE_ID!="INI25000") %>%
   #filter(GBE_ID!="INI25734")

# Filter for BM traits
brain_measure_ids <- fields_df %>%
  mutate(GBE_ID = paste0("INI", Field.ID)) %>%
  filter(GBE_ID != "INI25000", GBE_ID != "INI25734") %>%
  pull(GBE_ID)

bm_traits_in_pheno <- intersect(brain_measure_ids, names(pheno_df))
print("bm traits")
print(dim(pheno_df))

pheno_df <- pheno_df %>%
  filter(if_any(all_of(bm_traits_in_pheno), ~ !is.na(.)))
print(dim(pheno_df))

print("keep ID")
print(length(keepID))
train_keepID <- intersect(keepID, pheno_df$IID)
print(length(keepID))

## Normalize volumetric traits
print("traits_to_norm")

normalized_traits <- c("INI25001", "INI25003", "INI25005", "INI25007", "INI25009")

#print(names(fields_df))

traits_to_norm <- fields_df %>%
   filter(grepl("Vol\\.|volume", Description, ignore.case = TRUE)) %>%
   filter(GBE_ID != "INI25000") #%>%
   #filter(!(GBE_ID %in% normalized_traits))
print("filter2")
traits_to_norm <- traits_to_norm %>%
   filter(!(GBE_ID %in% normalized_traits))


norm_traits <- c("INI25000", normalized_traits, unique(traits_to_norm$GBE_ID))
#print(norm_traits)

print("filtering")
norm_pheno_df <- pheno_df %>%
   select(all_of(c("#FID", "IID", norm_traits)))

filter_traits <- c("INI25000", unique(traits_to_norm$GBE_ID))

print("norm_pheno_df")
norm_pheno_df <- norm_pheno_df %>% 
   select(all_of(c("#FID", "IID", filter_traits)))
   #filter(GBE_ID %in% filter_traits)

norm_pheno_df <- norm_pheno_df %>%
   mutate(across(-c("#FID", "IID", "INI25000"), ~ ifelse(!is.na(.), . * INI25000, NA)))

# remake pheno_df
pheno_df <- pheno_df %>%
   select(-all_of(norm_traits))

pheno_df <- pheno_df %>%
   left_join(norm_pheno_df, by=c("#FID", "IID")) %>%
   select(-INI25000)

#print("pheno_df")
#dim(pheno_df)

ini_cols <- grep("^INI", names(pheno_df), value = TRUE)

# Test set
test_set <- pheno_df %>%
   filter(!(IID %in% keepID))
print("test set")
print(dim(test_set))
write.csv(test_set, file.path(script_dir, "test_set.csv"))

#print("test")
#print(dim(test_set))
#print(length(unique(test_set$IID)))
#print(tail(test_set))

# Filter for training set
pheno_df <- pheno_df %>%
   filter(IID %in% keepID)
print("train")
print(dim(pheno_df))
#print("check training set")
#print(dim(pheno_df))
#print(names(pheno_df))
#print(table(pheno_df$population))

norm_cols <- c("INI25001", "INI25003", "INI25005", "INI25007", "INI25009","INI25002", "INI25004", "INI25006", "INI25008", "INI25010")

pop_counts <- pheno_df %>%
   summarise(across(all_of(ini_cols),  ~ sum(!is.na(.))))

# Center and scale the selected columns
scaling_params <- pheno_df %>%
  summarise(across(all_of(ini_cols), ~ list(mean = mean(., na.rm=TRUE), sd = sd(., na.rm=TRUE)))) %>%
  unnest(cols=everything()) %>%
  as.data.frame()
rownames(scaling_params) <- c("mean", "sd")

#str(scaling_params)

print("save scaling params")
# save scaling params
write.csv(scaling_params,"scaling_params.csv", row.names = TRUE)

scaling_params <- pheno_df %>%
  summarise(across(all_of(ini_cols), ~ list(mean = mean(., na.rm=TRUE), sd = sd(., na.rm=TRUE))))

phe_scaled <- pheno_df %>%
  mutate(across(all_of(ini_cols), ~ (.-scaling_params[[cur_column()]][["mean"]]) / scaling_params[[cur_column()]][["sd"]]))
#phe_scaled <- pheno_df %>%
#  mutate(across(all_of(ini_cols), 
#                ~ (. - scaling_params["mean", cur_column()]) / scaling_params["sd", cur_column()]))

# Dealing with NA
na_count <- sum(is.na(phe_scaled %>% select(all_of(ini_cols))))
nan_count <- sum(is.nan(as.matrix(phe_scaled %>% select(all_of(ini_cols)))))
inf_count <- sum(is.infinite(as.matrix(phe_scaled %>% select(all_of(ini_cols)))))
print(paste("Count of NA values before imputation:", na_count, "\n"))
print(paste("Count of NaN values before imputation:", nan_count, "\n"))
print(paste("Count of Inf values before imputation:", inf_count, "\n"))

print("NA counts")
na_vals <- phe_scaled %>%
  summarise(across(all_of(ini_cols), ~ sum(is.na(.))))
print(unique(na_vals))

# Impute with col means
phe_imputed <- phe_scaled %>%
   mutate(across(all_of(ini_cols), ~ replace(., is.na(.), mean(., na.rm=TRUE))))

# Remove columns with all NA
all_na <- phe_imputed %>%
  select(where(~ all(is.na(.))))
print(paste("All NA col number:", dim(all_na)))

phe_imputed <- phe_imputed %>%
  select(where(~ !all(is.na(.))))

ini_cols <- grep("^INI", names(phe_imputed), value = TRUE)

fields_df <- fields_df %>%
   filter(GBE_ID %in% ini_cols) %>%
   rename(Phenotype=GBE_ID)

pop_counts <- phe_imputed %>%
   summarise(across(all_of(ini_cols),  ~ sum(!is.na(.))))

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
#print(scaling_metadata)

save.image("cleaned_svd.RData")
