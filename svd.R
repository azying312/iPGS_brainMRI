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
load("cleaned_svd.RData")

####################################################################
# main
####################################################################
# SVD
phe_matrix <- phe_imputed %>%
  select(all_of(ini_cols)) %>%
  as.matrix()

dim(phe_matrix)

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

print("check training")
#dim(training_projection)
apply(phe_matrix, 2, sd)
sum(duplicated(phe_matrix))
dim(phe_matrix)
sum(duplicated(phe_matrix))

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

###########

# Create v_df
create_vdf <- function(comp1, comp2) {
   print("create vdf")
   v_df <- data.frame(Phenotype = ini_cols, SVD1 = v_matrix[, comp1], SVD2 = v_matrix[, comp2])
   v_df <- v_df %>%
      left_join(fields_df, by="Phenotype")
   # row 407 didn't read in properly
   #print(v_df[407,])
   v_df$Field.ID[407] <- 25871
   v_df$Description[407] <- "Vol._of_grey_matter_in_Heschl's_Gyrus_(includes_H1_and_H2)_(R)"
   v_df$Category[407] <- "T1_structural_brain_MRI_Regional_grey_matter_volumes_(FAST)"
   v_df$Category_Name[407] <- "T1_structural_brain_MRI"
   # row 314 didn't read in properly
   v_df$Field.ID[314] <- 25871
   v_df$Description[314] <- "Vol._of_grey_matter_in_Heschl's_Gyrus_(includes_H1_and_H2)_(R)"
   v_df$Category[314] <- "T1_structural_brain_MRI_Regional_grey_matter_volumes_(FAST)"
   v_df$Category_Name[314] <- "T1_structural_brain_MRI"
   v_df <- v_df %>%
      filter(Phenotype != "INI25000")
   return(v_df)
}

print("12")
v_df <- create_vdf(1,2)
print(head(v_df))
print(v_df[which(is.na(v_df$Category)),])
print(which(is.na(v_df$Category)))
pve_svd <- round(scree_df$PVE * 100, 2)
pve_svd1 <- round(scree_df$PVE[1] * 100, 2)
pve_svd2 <- round(scree_df$PVE[2] * 100, 2)
pve_svd3 <- round(scree_df$PVE[3] * 100, 2)

##########################

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

print("23")

# Components 2 and 3
pve1<-pve_svd[2]
pve2<-pve_svd[3]
v_df <- create_vdf(2,3)
svd_plot <- ggplot(v_df, aes(x = SVD1, y = SVD2, color = Category)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid") +
   geom_vline(xintercept=0, color="black", linetype="solid") +
   geom_text(aes(label = Phenotype), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD2 (", pve1, "%)"),
        y = paste0("SVD3 (", pve2, "%)"),
        colors="MRI Type") +
   theme_minimal()
#   guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))

ggsave(svd_plot_name23, svd_plot, width = 8, height = 6, dpi=300)

# Components 3 and 4
pve1<-pve_svd[3]
pve2<-pve_svd[4]
v_df <- create_vdf(3,4)
svd_plot <- ggplot(v_df, aes(x = SVD1, y = SVD2, color = Category)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid") +
   geom_vline(xintercept=0, color="black", linetype="solid") +
   geom_text(aes(label = Phenotype), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD3 (", pve1, "%)"),
        y = paste0("SVD4 (", pve2, "%)"),
        colors="MRI Type") +
   theme_minimal()
#   guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))

ggsave(svd_plot_name34, svd_plot, width = 8, height = 6, dpi=300)

# Components 4 and 5
pve1<-pve_svd[4]
pve2<-pve_svd[5]
v_df <- create_vdf(4,5)
svd_plot <- ggplot(v_df, aes(x = SVD1, y = SVD2, color = Category)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid") +
   geom_vline(xintercept=0, color="black", linetype="solid") +
   geom_text(aes(label = Phenotype), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD4 (", pve1, "%)"),
        y = paste0("SVD5 (", pve2, "%)"),
        colors="MRI Type") +
   theme_minimal()
#   guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))

ggsave(svd_plot_name45, svd_plot, width = 8, height = 6, dpi=300)

# Components 5 and 6
pve1<-pve_svd[5]
pve2<-pve_svd[6]
v_df <- create_vdf(5,6)
svd_plot <- ggplot(v_df, aes(x = SVD1, y = SVD2, color = Category)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid") +
   geom_vline(xintercept=0, color="black", linetype="solid") +
   geom_text(aes(label = Phenotype), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD5 (", pve1, "%)"),
        y = paste0("SVD6 (", pve2, "%)"),
        colors="MRI Type") +
   theme_minimal()
#   guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))

ggsave(svd_plot_name56, svd_plot, width = 8, height = 6, dpi=300)

# Components 6 and 7
pve1<-pve_svd[6]
pve2<-pve_svd[7]
v_df <- create_vdf(6,7)
svd_plot <- ggplot(v_df, aes(x = SVD1, y = SVD2, color = Category)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid") +
   geom_vline(xintercept=0, color="black", linetype="solid") +
   geom_text(aes(label = Phenotype), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD6 (", pve1, "%)"),
        y = paste0("SVD7 (", pve2, "%)"),
        colors="MRI Type") +
   theme_minimal()
#   guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))

ggsave(svd_plot_name67, svd_plot, width = 8, height = 6, dpi=300)

# Components 7 and 8
pve1<-pve_svd[7]
pve2<-pve_svd[8]
v_df <- create_vdf(7,8)
svd_plot <- ggplot(v_df, aes(x = SVD1, y = SVD2, color = Category)) +
   geom_point() +
   geom_hline(yintercept=0, color="black", linetype="solid") +
   geom_vline(xintercept=0, color="black", linetype="solid") +
   geom_text(aes(label = Phenotype), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) + # Add labels
   labs(x = paste0("SVD7 (", pve1, "%)"),
        y = paste0("SVD8 (", pve2, "%)"),
        colors="MRI Type") +
   theme_minimal()
#   guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))

ggsave(svd_plot_name78, svd_plot, width = 8, height = 6, dpi=300)

##########################
load(file.path(script_dir, "pop_counts.RData"))
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

print(head(sqc_df))
print(head(plot_df))

plot_df <- plot_df %>%
   left_join(sqc_df, by="IID")
pve1<-pve_svd[1]
pve2<-pve_svd[2]

dim(plot_df)

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
