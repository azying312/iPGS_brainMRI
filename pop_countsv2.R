
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

load(file.path(script_dir, "pop_counts.RData"))

print(tail(field_names_df))
print(dim(field_names_df))

# filter out QC terms from field names
field_names_df <- field_names_df %>%
   filter(GBE_ID != "INI25000") %>%
   filter(GBE_ID != "INI25734")

print(dim(field_names_df))

tab_dump_df <- tab_dump_df %>%
   mutate(GBE_ID=paste0("INI", field))

print("sqc_df")
print(dim(sqc_df))
print(tail(sqc_df))
print(length(unique(sqc_df$IID)))

print("tab dump df")
print(dim(tab_dump_df))
print(tail(tab_dump_df))
print(length(unique(tab_dump_df$IID)))

print("merge sqc and tab dumb")
merged.df <- sqc_df %>%
   left_join(tab_dump_df, by = "IID")
merged.df <- merged.df %>%
   mutate(population=ifelse(is.na(population), "NA", population))

print(dim(merged.df))
print(tail(merged.df))
print(length(unique(merged.df$IID)))
print(length(unique(merged.df$field)))


print("pass qc filter")
table(merged.df$pass_QC_filter)

print(table(merged.df$pass_QC_filter, merged.df$population))
print(length(unique(merged.df$IID)))

population_counts <- merged.df %>%
  #filter(pass_QC_filter) %>%
  #filter(!is.na(value)) %>%
  #filter(!is.na(population)) %>%
  filter(GBE_ID != "INI25000") %>%
  filter(GBE_ID != "INI25734") %>%
  mutate(
    phenotype_group = ifelse(GBE_ID %in% field_names_df$GBE_ID, "has_brain_measure", "other_measure")
  )
print(tail(population_counts))

print("qc and phenotype group")
print(table(population_counts$pass_QC_filter, population_counts$phenotype_group))

ukbb_summary <- population_counts %>%
   #filter(pass_QC_filter == TRUE) %>%
   #filter(!is.na(value)) %>%
   #filter(!is.na(population)) %>%
   group_by(IID, population, pass_QC_filter) %>%
   summarise(hasBM = any(phenotype_group == "has_brain_measure")) %>%
   ungroup()
#print("phenotypes")
#print(unique(length(ukbb_summary$field)))

print("ukbb")
print(tail(ukbb_summary))
print(table(ukbb_summary$pass_QC_filter))

print("qc filter, hasBM")
print(table(ukbb_summary$pass_QC_filter, ukbb_summary$hasBM))

print("population v. brain measure table")
print(table(ukbb_summary$population, ukbb_summary$hasBM))

max_population_phenotype <- population_counts %>%
  # filter out the non brain measure traits
  filter(!is.na(field)) %>%
  filter(field != 25000) %>%
  filter(field != 25734) %>%
  distinct(IID, field, .keep_all = TRUE) %>%
  group_by(field) %>%
  summarise(total_population = n_distinct(IID)) %>%
  arrange(desc(total_population)) %>%
  slice(1)

print(max_population_phenotype)
print("INI25001")
df_25001 <- merged.df %>%
  #distinct(IID, .keep_all = TRUE) %>%
  filter(field == 25001)

dim(df_25001)
print("population numbers")
print(length(unique(df_25001$IID)))
sum((df_25001$population == "NA"))
table(df_25001$population)
table(df_25001$pass_QC_filter)

brain_ukbb <- gsub("INI", "", field_names_df$GBE_ID)
print(length(brain_ukbb))
print(brain_ukbb[400:436])
print(class(brain_ukbb[1]))

print("25883")

print("25883" %in% brain_ukbb)
print(25883 %in% brain_ukbb)

max_population_phenotype <- population_counts %>%
  #filter out the non brain measure traits
  filter(!(field %in% brain_ukbb))# %>%
  #filter(field != 25000) %>%
  #filter(field != 25734) %>%
  #distinct(IID, field, .keep_all = TRUE) %>%
  #group_by(field) %>%
  #summarise(total_population = n_distinct(IID)) %>%
  #arrange(desc(total_population)) %>%
  #slice(1)
dim(population_counts)
dim(max_population_phenotype)
#max_population_phenotype
#table(max_population_phenotype$population)

tail(max_population_phenotype)

save.image("pop_counts.RData")
#population_counts %>% fwrite(pop_counts_f, sep = "\t", na = "NA", quote = F)
