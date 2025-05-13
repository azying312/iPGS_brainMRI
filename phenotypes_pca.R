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

hist1_name <- file.path(script_dir, "initial_pop_count_df.png")
hist2_name <- file.path(script_dir, "filtered_pop_count_df.png")
pop_counts_name <- file.path(script_dir, "pop_counts_scatterplt.png")

norm_plot_name12 <- file.path(script_dir, "norm_traits1_2.png")
norm_plot_name23 <- file.path(script_dir, "norm_traits2_3.png")
norm_plot_name34 <- file.path(script_dir, "norm_traits3_4.png")

svd_plot_name12 <- file.path(script_dir, "svd1_2.png")
svd_plot_name23 <- file.path(script_dir, "svd2_3.png")
svd_plot_name34 <- file.path(script_dir, "svd3_4.png")
svd_plot_name45 <- file.path(script_dir, "svd4_5.png")
svd_plot_name56 <- file.path(script_dir, "svd5_6.png")
svd_plot_name67 <- file.path(script_dir, "svd6_7.png")
svd_plot_name78 <- file.path(script_dir, "svd7_8.png")

scree_plot_name <- file.path(script_dir, "scree_plot.png") 

#bi_plot_name <- file.path(script_dir, "biplot.png")
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

print(paste("phe long f:", phe_long_f))

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
   rename(GBE_ID = "#GBE_ID") %>%
   filter(GBE_ID != "INI25000")
#print(head(phe_counts_df))

fields_df <- read.table(fields_f, sep='\t', header=TRUE)

fields_df <- fields_df %>%
   mutate(GBE_ID=paste0("INI", Field.ID)) %>%
   filter(GBE_ID!="INI25000")

ini_cols <- grep("^INI", names(pheno_df), value = TRUE)

# Initial pop counts
pop_counts <- pheno_df %>%
   summarise(across(all_of(ini_cols),  ~ sum(!is.na(.))))
print(head(pop_counts))

# hist of initial pop df
pheno_initial <- pivot_longer(pop_counts, cols = everything(), names_to = "Phenotype", values_to = "pop_counts_initial")

print("hist1")

pheno_initial_filtered <- pheno_initial[!is.na(pheno_initial$pop_counts_initial), ]
summary(pheno_initial_filtered$pop_counts_initial)
table(is.na(pheno_initial_filtered$pop_counts_initial))
range(pheno_initial_filtered$pop_counts_initial)
print("pheno initial")
head(pheno_initial_filtered)
print(tail(pheno_initial_filtered))

pheno_initial_filtered <- pheno_initial_filtered %>%
   filter(!(Phenotype %in% c("#FID", "IID", "INI25000")))

hist1 <- ggplot(pheno_initial_filtered, aes(x = pop_counts_initial)) +
  geom_histogram(fill = "blue", color = "black", alpha = 0.7) +
  #scale_x_log10() +
  labs(title = "Histogram of Initial Population Counts", x = "Population Count", y = "Frequency") +
  theme_minimal()
ggsave(hist1_name, hist1, width = 8, height = 6, dpi=300)

# Filter for training set
print("pheno")
print(tail(pheno_df))
dim(pheno_df)
pheno_df <- pheno_df %>%
   filter(IID %in% keepID)
dim(pheno_df)

print("hist2")
# hist of after filtered pop df
pop_counts <- pheno_df %>%
   summarise(across(all_of(ini_cols),  ~ sum(!is.na(.))))

pheno_filtered <- pivot_longer(pop_counts, cols = everything(), names_to = "Phenotype", values_to = "pop_counts_filtered")
pheno_filtered <- pheno_filtered %>%
   filter(!(Phenotype %in% c("#FID", "IID", "INI25000")))

print("filtered")
head(pheno_filtered)
tail(pheno_filtered)

hist2 <- ggplot(pheno_filtered, aes(x = pop_counts_filtered)) +
  geom_histogram(fill = "blue", color = "black", alpha = 0.7) +
  #scale_x_log10() +
  labs(title = "Histogram of Filtered Population Counts", x = "Population Count", y = "Frequency") +
  theme_minimal()

#ggsave(hist1_name, hist1, width = 8, height = 6, dpi=300)
ggsave(hist2_name, hist2, width = 8, height = 6, dpi=300)

pheno_initial_filtered <- pheno_initial %>%
  semi_join(pheno_filtered, by = "Phenotype")
#pheno_joined <- pheno_initial_filtered %>%
#  left_join(pheno_filtered, by = "Phenotype")

dim(pheno_initial)
dim(pheno_filtered)

library(data.table)

setDT(pheno_initial_filtered)
setDT(pheno_filtered)

class(pheno_initial_filtered)
class(pheno_filtered)

pheno_initial_filtered[, .N, by = Phenotype][N > 1]
pheno_filtered[, .N, by = Phenotype][N > 1]
pheno_initial_filtered <- pheno_initial_filtered[, .SD[1], by = Phenotype]
pheno_filtered <- pheno_filtered[, .SD[1], by = Phenotype]
gc()
pheno_joined <- pheno_initial_filtered[pheno_filtered, on = "Phenotype"]

print("summary")
summary(pheno_joined$pop_counts_initial)
summary(pheno_joined$pop_counts_filtered)

print(head(pheno_joined))

print("pop counts plt")
pop_counts_plt <- ggplot(pheno_joined, aes(x=pop_counts_initial, y=pop_counts_filtered, color=Phenotype))+
   geom_point(alpha=0.7, size=3)+
   geom_abline(slope=1, intercept=0, linetype="dashed", color="black")+
   labs(title = "Scatterplot of Initial vs. Filtered Population Counts",
        x = "Initial Population Count",
        y = "Filtered Population Count") +
   theme_minimal() +
   theme(legend.position = "none")

#ggsave(hist1_name, hist1, width = 8, height = 6, dpi=300)
#ggsave(hist2_name, hist2, width = 8, height = 6, dpi=300)
ggsave(pop_counts_name, pop_counts_plt, width = 8, height = 6, dpi=300)

####

norm_cols <- c("INI25001", "INI25003", "INI25005", "INI25007", "INI25009","INI25002", "INI25004", "INI25006", "INI25008", "INI25010")
print("matches")
print(match(norm_cols, ini_cols))

pop_counts <- pheno_df %>%
   summarise(across(all_of(ini_cols),  ~ sum(!is.na(.))))

print("num brain MRI people")
#print(max(pop_counts))
#print(unique(pop_counts))
unique_values <- unique(unlist(pop_counts))
print(unique_values)

# Center and scale the selected columns
scaling_params <- pheno_df %>%
  summarise(across(all_of(ini_cols), ~ list(mean = mean(., na.rm=TRUE), sd = sd(., na.rm=TRUE)))) %>%
  unnest(cols=everything()) %>%
  as.data.frame()
rownames(scaling_params) <- c("mean", "sd")

# save scaling params
write.csv(scaling_params,"scaling_params.csv", row.names = TRUE)
#print("write scaling params")
#print(head(scaling_params))

scaling_params <- pheno_df %>%
  summarise(across(all_of(ini_cols), ~ list(mean = mean(., na.rm=TRUE), sd = sd(., na.rm=TRUE))))

phe_scaled <- pheno_df %>%
  mutate(across(all_of(ini_cols), ~ (.-scaling_params[[cur_column()]][["mean"]]) / scaling_params[[cur_column()]][["sd"]]))
#phe_scaled <- pheno_df %>%
#  mutate(across(all_of(ini_cols), scale)) 
#print(head(phe_scaled))

# Dealing with NA
na_count <- sum(is.na(phe_scaled %>% select(all_of(ini_cols))))
nan_count <- sum(is.nan(as.matrix(phe_scaled %>% select(all_of(ini_cols)))))
inf_count <- sum(is.infinite(as.matrix(phe_scaled %>% select(all_of(ini_cols)))))
print(paste("Count of NA values before imputation:", na_count, "\n"))
print(paste("Count of NaN values before imputation:", nan_count, "\n"))
print(paste("Count of Inf values before imputation:", inf_count, "\n"))

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

na_count <- sum(is.na(phe_imputed %>% select(all_of(ini_cols))))
nan_count <- sum(is.nan(as.matrix(phe_imputed %>% select(all_of(ini_cols)))))
inf_count <- sum(is.infinite(as.matrix(phe_imputed %>% select(all_of(ini_cols)))))
print(paste("Count of NA values after imputation:", na_count, "\n"))
print(paste("Count of NaN values after imputation:", nan_count, "\n"))
print(paste("Count of Inf values after imputation:", inf_count, "\n"))

cols_with_na <- phe_imputed %>%
  select(all_of(ini_cols)) %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "na_count") %>%
  filter(na_count > 0)

ini_cols <- grep("^INI", names(phe_imputed), value = TRUE)

fields_df <- fields_df %>%
   filter(GBE_ID %in% ini_cols) %>%
   rename(Phenotype=GBE_ID)
print("fields df dim:")
print(dim(fields_df))
print(names(fields_df))
pop_counts <- phe_imputed %>%
   summarise(across(all_of(ini_cols),  ~ sum(!is.na(.))))
print("new pop counts")
print(unique(as.numeric(pop_counts)))

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

#v_df <- v_df %>%
#        filter(Phenotype != "INI25000")

# Create v_df
create_vdf <- function(comp1, comp2) {
   #print("create vdf")
   v_df <- data.frame(Phenotype = ini_cols, SVD1 = v_matrix[, comp1], SVD2 = v_matrix[, comp2])
   #print(head(v_df))
   v_df <- v_df %>%
      left_join(fields_df, by="Phenotype")
   # row 407 didn't read in properly
   #print(v_df[407,])
   v_df$Field.ID[407] <- 25871
   v_df$Description[407] <- "Vol._of_grey_matter_in_Heschl's_Gyrus_(includes_H1_and_H2)_(R)"
   v_df$Category[407] <- "T1_structural_brain_MRI_Regional_grey_matter_volumes_(FAST)"
   v_df$Category_Name[407] <- "T1_structural_brain_MRI"

   v_df <- v_df %>% 
      filter(Phenotype != "INI25000")
      
   # normalized phenotypes
   #ini_norm <- c("INI25001", "INI25003", "INI25005", "INI25007", "INI25009")
   #print("ini norm")
   #print(match(ini_norm, ini_cols))
    
   # not normalized
   #ini_reg <- c("INI25002", "INI25004", "INI25006", "INI25008", "INI25010")
   #print("ini reg")
   #print(match(ini_reg, ini_cols))
   # LR Experiment
   #LR_pheno <- c("INI25062", "INI25063", "INI25080", "INI25081","INI25571","INI25572","INI25026","INI25027","INI25011","INI25012","INI25870","INI25871")
   
   #v_df <- v_df %>%
   # mutate(norm_v_not = case_when(
   #     Phenotype %in% ini_norm ~ "normalised",
   #     Phenotype %in% ini_reg ~ "not_normalised",
   #     TRUE ~ NA_character_
   # )) #%>%
	#mutate(Category_Name=ifelse(Phenotype %in% LR_pheno, "LR Experiment", Category_Name),
	#Phenotype=ifelse(Category_Name=="LR Experiment", Phenotype, ""))
   return(v_df)
}

v_df <- create_vdf(1,2)
print(unique(v_df$Category_Name))
print(unique(v_df$Phenotype))
print(sum(is.na(v_df$Phenotype)))

#norm_df <- v_df %>%
#   filter(!is.na(norm_v_not))

#norm_plot1 <- ggplot(norm_df, aes(x=SVD1, y=SVD2, color=as.factor(Field.ID))) +
#   geom_point()+
#   geom_hline(yintercept=0, color="black", linetype="solid")+
#   geom_vline(xintercept=0, color="black", linetype="solid")+
#   geom_text(aes(label = as.factor(Field.ID)), hjust = 0, vjust = 0, size = 3, max.overlaps = Inf) + #, check_overlap = TRUE) +
#   labs(x = "SVD1",
#        y = "SVD2") +
#   theme_minimal()

#ggsave(norm_plot_name12, norm_plot1, width = 8, height = 6, dpi=300)

print(table(v_df$Category))
print("no category")
no_category <- v_df %>%
   filter(is.na(Category))
print(no_category)

print("svd plots")

pve_svd1 <- round(scree_df$PVE[1] * 100, 2)
pve_svd2 <- round(scree_df$PVE[2] * 100, 2)
pve_svd <- round(scree_df$PVE * 100, 2)

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

print("svd 2, 3")

# Component 2 and 3
pve2<-pve_svd[2]
pve3<-pve_svd[3]
print("pve")
print(pve2)
print(pve3)

print("V_DF for 2 and 3")
v_df <- create_vdf(2,3)
print(dim(v_df))
#norm_df <- v_df %>%
#   filter(!is.na(norm_v_not))

#norm_plot1 <- ggplot(norm_df, aes(x=SVD1, y=SVD2, color=as.factor(Field.ID))) +
#   geom_point()+
#   geom_hline(yintercept=0, color="black", linetype="solid")+
#   geom_vline(xintercept=0, color="black", linetype="solid")+
#   geom_text(aes(label = as.factor(Field.ID)), hjust = 0, vjust = 0, size = 3, max.overlaps = Inf) + #, check_overlap = TRUE) +
#   labs(x = "SVD2",
#        y = "SVD3") +
#   theme_minimal()
#
#ggsave(norm_plot_name23, norm_plot1, width = 8, height = 6, dpi=300)

#no_category <- v_df %>%
#   filter(is.na(Category_Name))

print(svd_plot_name23)

svd_plot <- ggplot(v_df, aes(x = SVD1, y = SVD2, color = Category)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid") +
   geom_vline(xintercept=0, color="black", linetype="solid") +
   geom_text(aes(label = Phenotype), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD2 (", pve2, "%)"),
        y = paste0("SVD3 (", pve3, "%)"),
        colors="MRI Type") +
   theme_minimal()

print("plt name")
print(svd_plot)

ggsave(svd_plot_name23, svd_plot, width = 8, height = 6, dpi=300)

# Components 3 and 4
pve4 <- pve_svd[4]
#v_df <- data.frame(Phenotype = ini_cols, SVD1 = v_matrix[, 1], SVD2 = v_matrix[, 4])
v_df <- create_vdf(3, 4)

#norm_df <- v_df %>%
#   filter(!is.na(norm_v_not))

#norm_plot1 <- ggplot(norm_df, aes(x=SVD1, y=SVD2, color=as.factor(Field.ID))) +
#   geom_point()+
#   geom_hline(yintercept=0, color="black", linetype="solid")+
#   geom_vline(xintercept=0, color="black", linetype="solid")+
#   geom_text(aes(label = as.factor(Field.ID)), hjust = 0, vjust = 0, size = 3, max.overlaps = Inf) + #, check_overlap = TRUE) +
#   labs(x = "SVD3",
#        y = "SVD4") +
#   theme_minimal()

#ggsave(norm_plot_name34, norm_plot1, width = 8, height = 6, dpi=300)

svd_plot <- ggplot(v_df, aes(x = SVD1, y = SVD2, color = Category)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid")+
   geom_vline(xintercept=0, color="black", linetype="solid")+
   geom_text(aes(label = Phenotype), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD3 (", pve_svd[3], "%)"),
        y = paste0("SVD4 (", pve_svd[4], "%)",
	colors="MRI Types") +
   theme_minimal()
#svd_plot <- svd_plot +
#   guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))

ggsave(svd_plot_name34, svd_plot, width = 8, height = 6, dpi=300)

print("svd 4 5")
# Components 4 and 5
pve5 <- pve_svd[5]
#v_df <- data.frame(Phenotype = ini_cols, SVD1 = v_matrix[, 2], SVD2 = v_matrix[, 4])
v_df <- create_vdf(4, 5)
svd_plot <- ggplot(v_df, aes(x = SVD1, y = SVD2, color = Category)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid")+
   geom_vline(xintercept=0, color="black", linetype="solid")+
   #geom_text(aes(label = Phenotype), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD4 (", pve_svd[4], "%)"),
        y = paste0("SVD5 (", pve_svd[5], "%)",
	colors="MRI Types") +
   theme_minimal()

#svd_plot <- svd_plot +
#   guides(color = guide_legend(ncol = 5, override.aes = list(size = 2)))

ggsave(svd_plot_name45, svd_plot, width = 8, height = 6, dpi=300)

print("svd 5 6")
pve6 <- pve_svd[6]
# Components 5 and 6
#v_df <- data.frame(Phenotype = ini_cols, SVD1 = v_matrix[, 3], SVD2 = v_matrix[, 4])
v_df <- create_vdf(5,6)
svd_plot <- ggplot(v_df, aes(x = SVD1, y = SVD2, color = Category)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid")+
   geom_vline(xintercept=0, color="black", linetype="solid")+
   #geom_text(aes(label = Phenotype), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD5 (", pve_svd[5], "%)"),
        y = paste0("SVD6 (", pve_svd[6], "%)"
	colors="MRI Types") +
   theme_minimal()
#svd_plot <- svd_plot +
#   guides(color = guide_legend(ncol = 5, override.aes = list(size = 2)))

ggsave(svd_plot_name56, svd_plot, width = 8, height = 6, dpi=300)

print("svd 6 7")
pve7 <- pve_svd[7]
# Components 6 and 7
#v_df <- data.frame(Phenotype = ini_cols, SVD1 = v_matrix[, 2], SVD2 = v_matrix[, 3])
v_df <- create_vdf(6, 7)
svd_plot <- ggplot(v_df, aes(x = SVD1, y = SVD2, color = Category)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid")+
   geom_vline(xintercept=0, color="black", linetype="solid")+
   #geom_text(aes(label = Phenotype), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD6 (", pve_svd[6], "%)"),
        y = paste0("SVD7 (", pve_svd[7], "%)",
	colors="MRI Types") +
   theme_minimal()
#svd_plot <- svd_plot +
#   guides(color = guide_legend(ncol = 5, override.aes = list(size = 2)))

ggsave(svd_plot_name67, svd_plot, width = 8, height = 6, dpi=300)

print("svd 7 8")
# Comp 7 and 8
pve8 <- pve_svd[8]
v_df <- create_vdf(7,8)
svd_plot <- ggplot(v_df, aes(x = SVD1, y = SVD2, color = Category)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid")+
   geom_vline(xintercept=0, color="black", linetype="solid")+
   #geom_text(aes(label = Phenotype), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD7 (", pve_svd[7], "%)"),
        y = paste0("SVD8 (", pve_svd[8], "%)",
	colors="MRI Types") +
   theme_minimal()
ggsave(svd_plot_name78, svd_plot, width = 8, height = 6, dpi=300)


# Biplot
#sample_points <- data.frame(
#  SVD1 = svd_result$u[,comp1] * svd_result$d[comp1],
#  SVD2 = svd_result$u[,comp2] * svd_result$d[comp2],
#  Type = "Sample"
#)

#phenotype_points <- data.frame(
#  SVD1 = svd_result$v[,comp1] * svd_result$d[comp1],
#  SVD2 = svd_result$v[,comp2] * svd_result$d[comp2],
#  Type = "Phenotype"
#)

#biplot_plt <- ggplot() +
#  # scatter plot points (individuals)
#  geom_point(data = sample_points, aes(x = SVD1, y = SVD2, color = Type)) +
#  # feature vectors (phenotypes)
#  geom_segment(data = v_df, aes(x = 0, y = 0, xend = SVD1, yend = SVD2, color = Phenotype),
#               arrow = arrow(length = unit(0.2, "cm"))) +
#  geom_text(check_overlap = TRUE, size = 3, vjust = 1.5) +
#  geom_hline(yintercept = 0, linetype = "solid") +
#  geom_vline(xintercept = 0, linetype = "solid") +
#  labs(x = x_name, y = y_name, title = "SVD Biplot: Samples and Phenotypes") +
#  scale_x_continuous(
#    breaks = seq(-10, 10, by = 2),  
#    labels = seq(-10, 10, by = 2)   
#  ) +
#  scale_y_continuous(
#    breaks = seq(-10, 10, by = 2), 
#    labels = seq(-10, 10, by = 2)  
#  ) +
#  theme_minimal()

#biplot_plt <- biplot_plt +
#   guides(color = guide_legend(ncol = 5, override.aes = list(size = 2)))

#ggsave(bi_plot_name, biplot_plt, width = 8, height = 6, dpi=300)

print("done")
