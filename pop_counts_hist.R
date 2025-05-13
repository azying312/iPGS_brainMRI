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

source(file.path(repo_dir, 'paths.sh'))
source(file.path(repo_dir, 'helpers', '/functions.R'))

####################################################################
# output files
####################################################################
pop_count_plt <- file.path(script_dir, "pop_count_hist.png")
norm_plot_name12 <- file.path(script_dir, "norm_traits1_2.png")
norm_plot_name23 <- file.path(script_dir, "norm_traits2_3.png")
norm_plot_name34 <- file.path(script_dir, "norm_traits3_4.png")
scree_plot_name <- file.path(script_dir, "scree_plot.png")
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

phe_f %>%
   fread() -> pheno_df
print(head(pheno_df))
phe_counts_f %>%
   fread() -> phe_counts_df
phe_counts_df <- phe_counts_df %>%
   select("#GBE_ID", "total") %>%
   rename(GBE_ID = "#GBE_ID") %>%
   filter(GBE_ID != "INI25000")
fields_df <- read.table(fields_f, sep='\t', header=TRUE)
fields_df <- fields_df %>%
   mutate(GBE_ID=paste0("INI", Field.ID)) %>%
   filter(GBE_ID!="INI25000")
ini_cols <- grep("^INI", names(pheno_df), value = TRUE)

pop_counts <- pheno_df %>%
   summarise(across(all_of(ini_cols),  ~ sum(!is.na(.))))

# Histogram of pop counts
pop_counts_df <- pivot_longer(pop_counts, cols = everything(), names_to = "Phenotype", values_to = "pop_counts")
pop_counts_df <- pop_counts_df %>%
   filter(!(Phenotype %in% c("#FID", "IID", "INI25000")))

hist_plt <- ggplot(pop_counts_df, aes(x = pop_counts)) +
  geom_histogram(fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Population Counts", x = "Population Count", y = "Frequency") +
  theme_minimal()

ggsave(pop_count_plt, hist_plt, width = 8, height = 6, dpi=300)


