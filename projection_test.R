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
svd_results_dir <- file.path(script_dir, "svd_results")
svd_plot_name12 <- file.path(script_dir, "svd1_2_test.png")
svd_plot_name23 <- file.path(script_dir, "svd2_3_test.png")
svd_plot_name34 <- file.path(script_dir, "svd3_4_test.png")
svd_plot_name45 <- file.path(script_dir, "svd4_5_test.png")

####################################################################
# input files
####################################################################

training_svd_v <-  file.path(script_dir, "training_svd_v_mat.csv")
scaling_params <- file.path(script_dir, "scaling_params.csv")
test_set <- file.path(script_dir, "test_set.csv")
fields_f <- file.path(
    script_dir, "misc",
    "mri_traits_449_names.tsv"
)
scree_f <- file.path(script_dir, "training_scree_components.csv")

####################################################################
# main
####################################################################

fields_df <- read.table(fields_f, sep='\t', header=TRUE)
scree_df <- read.csv(scree_f)

print(head(scree_df))

# Filter out QC traits
fields_df <- fields_df %>%
   mutate(Phenotype=paste0("INI", Field.ID))   
#mutate(GBE_ID=paste0("INI", Field.ID))
#fields_df <- fields_df %>%
# rename(Phenotype = Field.ID)

print(training_svd_v)
training_svd_v_df <- read.csv(training_svd_v)
training_svd_v_df <- training_svd_v_df %>%
   select(-X)

scaling_params_df <- read.csv(scaling_params)
scaling_params_df <- scaling_params_df %>%
   select(-X)

test_set_df <- read.csv(test_set)
test_set_df <- test_set_df %>%
   select(-X) %>%
   rename(FID = X.FID)

rownames(scaling_params_df) <- c("mean", "sd")

ini_cols_scaling <- grep("^INI", names(scaling_params_df), value = TRUE)
ini_cols <- grep("^INI", names(test_set_df), value = TRUE)
scaling_params_t <- as.data.frame(t(scaling_params_df))

scaling_params_list <- lapply(seq_len(nrow(scaling_params_t)), function(i) {
  list(mean = scaling_params_t$mean[i], sd = scaling_params_t$sd[i])
})
names(scaling_params_list) <- rownames(scaling_params_t)

# Center and Scale
phe_scaled <- test_set_df
phe_scaled[ini_cols] <- (test_set_df[ini_cols] - 
  as.numeric(scaling_params_df["mean", ini_cols_scaling])) /
  as.numeric(scaling_params_df["sd", ini_cols_scaling])

# Impute with col means
#phe_imputed <- phe_scaled %>%
#   mutate(across(all_of(ini_cols), ~ replace(., is.na(.), scaling_params_list[[cur_column()]][["mean"]])))
phe_imputed <- phe_scaled %>%
   mutate(across(all_of(ini_cols), ~ replace(., is.na(.), mean(., na.rm = TRUE))))

#print(head(phe_imputed))
#phe_imputed <- phe_imputed %>%
#   select(-any_of(c("INI25000", "INI25734")))

# Remove columns with all NA
all_na <- phe_imputed %>%
  select(where(~ all(is.na(.))))

phe_imputed <- phe_imputed %>%
  select(where(~ !all(is.na(.))))

ini_cols <- grep("^INI", names(phe_imputed), value = TRUE)

phe_matrix_test <- phe_imputed %>%
  select(all_of(ini_cols)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()
training_svd_v_df <- as.matrix(sapply(training_svd_v_df, as.numeric))


load(file.path(script_dir, "pop_counts.RData"))

# test pop projection
test_projection <- phe_matrix_test %*% training_svd_v_df
print("test proj")
print(head(test_projection))

# SVD
phe_matrix <- phe_imputed %>%
  select(all_of(ini_cols)) %>%
  as.matrix()

# SVD calculation
#svd_result <- svd(phe_matrix, nu = 8, nv = 8)
#u_matrix <- svd_result$u[, 1:8]
#d_vec <- svd_result$d
#svd_u <- d_vec^2/sum(d_vec^2)
#scree_df <- data.frame(
#   Component=c(1:length(svd_u)),
#   PVE=svd_u,
#   PVE_sum=cumsum(svd_u)
#)

#test_proj <- u_matrix %*% training_svd_v[, 1:8]

if (!dir.exists(svd_results_dir)) {
  dir.create(svd_results_dir)}

#############
print("issues in test set")
summary(apply(phe_matrix_test, 1, function(x) sd(x)))
#tail(phe_matrix_test)
str(training_svd_v_df)
summary(training_svd_v_df)
dim(phe_matrix_test)
dim(training_svd_v_df)
print("unique rows")
unique_rows <- !duplicated(phe_matrix_test)
sum(!unique_rows)
#print("non 0 var")
#apply(phe_matrix_test, 2, function(col) sum(is.na(col)))
#apply(phe_matrix_test, 2, sd)

sum(apply(phe_matrix_test, 1, function(row) length(unique(row)) == 1))

sum(duplicated(test_set_df$IID))

print("projected")
projected_phe <- as.matrix(phe_matrix_test) %*% as.matrix(training_svd_v_df)
sum(duplicated(round(projected_phe, digits = 20)))
sum(duplicated(round(projected_phe, digits = 4)))


############

# Create v_df
create_vdf <- function(comp1, comp2) {
   print("create vdf")
   #v_df <- data.frame(IID = test_set_df$IID , SVD1 = u_matrix[, comp1], SVD2 = u_matrix[, comp2])
   v_df <- data.frame(IID = test_set_df$IID , SVD1 = test_projection[comp1], SVD2 = test_projection[, comp2])
   return(v_df)
}

print("12")
v_df <- create_vdf(1,2)

pve_svd <- round(scree_df$PVE * 100, 2)

# Components 1 and 2
svd1 <- test_projection[, 1] 
svd2 <- test_projection[, 2]
svd3 <- test_projection[, 3]
svd4 <- test_projection[, 4]
svd5 <- test_projection[, 5]
svd6 <- test_projection[, 6]
svd7 <- test_projection[, 7]
svd8 <- test_projection[, 8]

plot_df <- data.frame(
  IID = test_set_df$IID,
  SVD1 = svd1,
  SVD2 = svd2,
  SVD3 = svd3,
  SVD4 = svd4,
  SVD5 = svd5,
  SVD6 = svd6,
  SVD7 = svd7,
  SVD8 = svd8
)

print("sqc")
print(head(sqc_df))
sqc_df$IID <- as.numeric(sqc_df$IID)

print("plot df dim")
print(dim(plot_df))
print(head(plot_df))
plot_df <- plot_df %>%
   left_join(sqc_df, by="IID") #%>%
   #filter(pass_QC_filter==TRUE) %>%
   #filter(population != "WB") %>%
   #filter(population != "others")

print("plot df dims")
print(dim(plot_df))

print(table(plot_df$population))
print(table(plot_df$pass_QC_filter))

pve1<-pve_svd[1]
pve2<-pve_svd[2]

svd_plot <- ggplot(plot_df, aes(x = SVD1, y = SVD2, color = population)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid") +
   geom_vline(xintercept=0, color="black", linetype="solid") +
   #geom_text(aes(label = population), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD1 (", pve1, "%)"),
        y = paste0("SVD2 (", pve2, "%)") ,
        colors="Population") +
   theme_minimal()
ggsave(svd_plot_name12, svd_plot, width = 8, height = 6, dpi=300)

# Components 2 and 3
pve3<-pve_svd[3]
svd_plot <- ggplot(plot_df, aes(x = SVD2, y = SVD3, color = population)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid") +
   geom_vline(xintercept=0, color="black", linetype="solid") +
   #geom_text(aes(label = population), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD2 (", pve2, "%)"),
        y = paste0("SVD3 (", pve3, "%)"),
        colors="Population") +
   theme_minimal()

ggsave(svd_plot_name23, svd_plot, width = 8, height = 6, dpi=300)

# Components 3 and 4
pve4<-pve_svd[4]

svd_plot <- ggplot(plot_df, aes(x = SVD3, y = SVD4, color = population)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid") +
   geom_vline(xintercept=0, color="black", linetype="solid") +
   #geom_text(aes(label = population), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD3 (", pve3, "%)"),
        y = paste0("SVD4 (", pve4, "%)"),
        colors="Population") +
   theme_minimal()

ggsave(svd_plot_name34, svd_plot, width = 8, height = 6, dpi=300)


# Components 4 and 5
pve5<-pve_svd[5]
#v_df <- create_vdf(4,5)
svd_plot <- ggplot(plot_df, aes(x = SVD4, y = SVD5, color = population)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid") +
   geom_vline(xintercept=0, color="black", linetype="solid") +
   #geom_text(aes(label = population), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD4 (", pve4, "%)"),
        y = paste0("SVD5 (", pve5, "%)"),
        colors="Population") +
   theme_minimal()

ggsave(svd_plot_name45, svd_plot, width = 8, height = 6, dpi=300)



