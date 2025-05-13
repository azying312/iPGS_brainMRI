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

####################################################################
# output files
####################################################################

phe_f        <- file.path(ukb21942_user, "pheno", "brain_mri.tsv_keep.gz")
phe_long_f   <- file.path(ukb21942_user, "pheno", "brain_mri.long_keep.tsv.gz")
phe_counts_f <- file.path(ukb21942_user, "pheno", "brain_mri.counts_keep.tsv.gz")
pop_counts_f <- file.path(script_dir, "pop_counts.tsv.gz")
print("output files")
print(paste("phe_f:", phe_f))
print(paste("phe_long_f:", phe_long_f))
print(paste("phe_counts_f:", phe_counts_f))

####################################################################
# input files
####################################################################
keep_file <- "/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942/sqc/population/train.keep"
sqc_f <- file.path(ukb21942_d, "sqc", "sqc.20220316.tsv.gz")
fields_f <- file.path(script_dir, "misc", "mri_traits_449_names.tsv") #"fields.tsv")
#SQL_dump_f <- file.path("${ukb21942_user}/pheno/tmp/${analysis_name}/tab_dump.tsv.gz")
SQL_dump_f <- file.path(ukb21942_user, "pheno", "tmp", basename(script_dir), "tab_dump.tsv.gz")

print(sqc_f)
print(SQL_dump_f)

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

print("sqc_df")
print(table(sqc_df$pass_QC_filter))
print(table(sqc_df$population, sqc_df$pass_QC_filter))
print(dim(sqc_df))
print(head(sqc_df))
print(length(unique(sqc_df$IID)))

fields_f %>%
  fread() %>%
  rename_with(~str_replace(., "#", ""), starts_with("#")) -> fields_df

print("fields_df")
print(tail(fields_df))
print(length(unique(fields_df$`Field ID`)))

fields_df %>%
rename("GBE_ID" = "Field.ID", "GBE_NAME" = "Description") %>%
mutate(
    GBE_ID = paste0("INI", GBE_ID),
    GBE_short_name = GBE_NAME,
    GBE_short_name = str_replace(GBE_short_name, "Volume", "Vol."),
    GBE_short_name = str_replace(GBE_short_name, "\\(left\\)", "L"),
    GBE_short_name = str_replace(GBE_short_name, "\\(right\\)", "R"),
    GBE_short_name = str_replace(GBE_short_name, "percentage", "%"),
    GBE_short_name = str_replace(GBE_short_name, "Impedance", "Impd."),
) -> field_names_df

print("brain GBE IDs")
print(length(unique(field_names_df$GBE_ID)))

SQL_dump_f %>%
  fread(colClasses = c("#eid" = "character", "value" = "character")) %>%
  rename("IID" = "#eid") -> tab_dump_df
print("tab dump df")
print(dim(tab_dump_df))
print(head(tab_dump_df))
print(length(unique(tab_dump_df$IID)))
tab_dump_df %>% count(value < 0) %>% print()
# we confirm that all the values are non-negative

tab_dump_df %>% count(array) %>% print()
# we confirm that array == 0 for all lines

# aggregate multiple measurements at multiple instances
# we take non-NA median

save.image(file.path(script_dir, "pop_counts.RData"))

print("save image")

tab_dump_df %>%
  select(-array) %>%
  mutate(value = as.numeric(value)) %>%
  group_by(IID, field) %>%
  summarise(value = median(value), .groups = "drop") %>%
  mutate(GBE_ID = paste0("INI", field)) %>%
  select(IID, GBE_ID, value) -> tab_aggregated_df

# Filter for keep file
keep_data <- read.table(keep_file, header = FALSE, sep = "", stringsAsFactors = FALSE)
keepID <- keep_data[,1]
tab_aggregated_df <- tab_aggregated_df %>%
   filter(IID %in% keepID)
print(paste("tab agg dim:", dim(tab_aggregated_df)))

# save phe file

tab_aggregated_df %>%
  select(GBE_ID) %>%
  unique %>%
  mutate(
    GBE_num = as.integer(str_replace_all(GBE_ID, "[a-zA-Z]", "")),
    GBE_cat = str_replace_all(GBE_ID, "[0-9]", "")
  ) %>%
  arrange(GBE_cat, GBE_num) %>%
  pull(GBE_ID) -> GBE_ID_sorted

tab_aggregated_df %>%
  mutate(FID = IID) %>%
  select(FID, IID, GBE_ID, value) %>%
  unique() %>%
  mutate(
    GBE_ID = factor(GBE_ID, levels = GBE_ID_sorted)
  ) %>%
  arrange(GBE_ID, FID, IID) -> phe_long_df

phe_long_df %>%
  rename("#FID" = "FID") %>%
  fwrite(phe_long_f, sep = "\t", na = "NA", quote = F)

print("phe_long_df")
print(head(phe_long_df))

sqc_df %>%
  select(FID, IID) %>%
  left_join(
    phe_long_df %>%
      pivot_wider(
        names_from = GBE_ID,
        values_from = value
      ),
    by = c("FID", "IID")
  ) -> phe_wide_df

print("phe_wide_df")
print(dim(phe_wide_df))

phe_wide_df %>%
  rename("#FID" = "FID") %>%
  fwrite(phe_f, sep = "\t", na = "NA", quote = F)

pops <- c("WB", "NBW", "Afr", "SA", "others")

bind_rows(
  phe_long_df %>%
    inner_join(
      sqc_df %>%
        select(FID, IID, population) %>%
        drop_na(population),
      by = c("FID", "IID")
    ) %>% count(population, GBE_ID),
  
  phe_long_df %>%
    inner_join(
      sqc_df %>%
        select(FID, IID, population_MM) %>%
        drop_na(population_MM),
      by = c("FID", "IID")
    ) %>% count(population_MM, GBE_ID) %>%
    rename("population" = "population_MM")
) %>%
  drop_na(GBE_ID) %>%
  pivot_wider(
    names_from = population, values_from = n, values_fill = 0
  ) %>%
  mutate(
    total = WB + NBW + Afr + SA + others,
    total_MM = WB_MM + NBW_MM + Afr_MM + SA_MM + others_MM
  ) %>%
  left_join(
    field_names_df %>%
      select(GBE_ID, GBE_NAME, GBE_short_name),
    by = "GBE_ID"
  ) %>%
  select(all_of(c(
    "GBE_ID", "GBE_NAME", "GBE_short_name",
    pops, paste0(pops, "_MM"), "total", "total_MM"
  ))) %>%
  mutate(
    GBE_ID = factor(GBE_ID, levels = GBE_ID_sorted)
  ) %>%
  arrange(GBE_ID) -> phe_counts_df

print("phe_counts_df")
print(head(phe_counts_df))

phe_counts_df %>%
  rename("#GBE_ID" = "GBE_ID") %>%
  fwrite(phe_counts_f, sep = "\t", na = "NA", quote = F)

print("WB")
print(unique(phe_counts_df$WB))
print("NBW")
print(unique(phe_counts_df$NBW))
print("Afr")
print(unique(phe_counts_df$Afr))
print("SA")
print(unique(phe_counts_df$SA))
print("others")
print(unique(phe_counts_df$others))

print("done")
