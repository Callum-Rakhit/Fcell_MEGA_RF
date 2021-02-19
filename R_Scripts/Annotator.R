# Function to install packages
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

# Install required packages
GetPackages(c("tidyverse", "randomForest", "randomForestExplainer", "dplyr", "igraph", "caret", 
              "reshape", "rsample", "reshape2", "devtools", "PerformanceAnalytics", "ggplot2",
              "lubridate","bootstrap", "corrplot", "ggraph", "doParallel", "ranger", "data.table", 
              "h2o", "sparsio", "qqman"))

options(scipen = 999)

# Get the positional information for all predictors in the over 10% f-cell group
sickle_model_over10fcell_all_predictors <- read.csv("~/Desktop/sickle_fcell_rf/Sickle_All_over10fcell_NonLog_1000.csv")
sickle_model_over10fcell_all_predictors$variable <- 
  substr(sickle_model_over10fcell_all_predictors$variable, 1, 
         nchar(sickle_model_over10fcell_all_predictors$variable)-2) # Remove last 2 characters
positional_info <- read.delim("~/Desktop/sickle_fcell_rf/SickleMEGA_Pruned_NoSexMismatch.bim", sep = "\t", header = F)
colnames(positional_info) <- c("Chr", "variable", "Chr_Position", "position", "ref", "alt")
merged <- merge(sickle_model_over10fcell_all_predictors, positional_info, by = "variable")
merged$Chr_Position <- paste(merged$Chr, merged$position, sep = ":")

# Same for the healthy group
healthy_model_top_1000_predictors <- read.csv("~/Desktop/sickle_fcell_rf/Healthy_Top_100.csv.csv")
healthy_model_top_1000_predictors$variable <- 
  substr(healthy_model_top_1000_predictors$variable, 1, 
         nchar(healthy_model_top_1000_predictors$variable)-2) # Remove last 2 characters

files <- list.files(path = "~/Twins/Data", pattern = "*.bim", full.names = T)
healthy_positional_info <-  read.delim(files[1], sep = "\t", header = F)
healthy_positional_info2 <-  read.delim(files[2], sep = "\t", header = F)
healthy_positional_info3 <-  read.delim(files[3], sep = "\t", header = F)
healthy_positional_info <- rbind(healthy_positional_info, healthy_positional_info2)
healthy_positional_info <- rbind(healthy_positional_info, healthy_positional_info3)

colnames(healthy_positional_info) <- c("Chr", "variable", "Chr_Position", "position", "ref", "alt")
merged_healthy <- merge(healthy_model_top_1000_predictors, healthy_positional_info, by = "variable")
# merged_healthy$Chr_Position <- paste(merged_healthy$Chr, merged_healthy$position, sep = ":")

# Clear old dataframes
rm(healthy_model_top_1000_predictors, healthy_positional_info, healthy_positional_info2, healthy_positional_info3)

# # Load in John's rs table
# dbSNP <- read.csv("~/Desktop/sickle_fcell_rf/dbSNP_IDs/All_RSid_matches.txt", header = F, sep = " ")
# colnames(dbSNP) <- c("Chr_Position", "rsID")
# merged <- merge(merged, dbSNP, by = "Chr_Position", all.x = T)

# Generate and save over 10% f-cell positions to file if needed
postions <- as.data.frame(merged$Chr)
postions$start <- merged$position
postions$end <- merged$position+1
postions$ref <- merged$ref
postions$alt <- merged$alt
postions$rsID <- merged$variable
postions$relative_importance <- merged$relative_importance
colnames(postions) <- c("Chr", "Start", "End", "Ref", "Alt", "rsID", "relative_importance")
rownames(postions) <- NULL
write.table(postions, file = "~/Desktop/sickle_fcell_rf/dbSNP_IDs/Over10fcell_NonLog_Sickle_rsID.txt", sep = "\t", row.names = F, col.names = F)              

# Generate and save healthy positions to file if needed
postions_healthy <- as.data.frame(merged_healthy$Chr)
postions_healthy$start <- merged_healthy$position
postions_healthy$end <- merged_healthy$position+1
postions_healthy$ref <- merged_healthy$ref
postions_healthy$alt <- merged_healthy$alt
postions_healthy$rsID <- merged_healthy$variable
postions_healthy$relative_importance <- merged_healthy$relative_importance
colnames(postions_healthy) <- c("Chr", "Start", "End", "Ref", "Alt", "rsID", "relative_importance")
rownames(postions_healthy) <- NULL
write.table(postions_healthy, file = "~/Desktop/sickle_fcell_rf/dbSNP_IDs/Healthy_Predictor_Positions.txt", sep = "\t", row.names = F, col.names = F)              

# Manhatten plot for the over 10% f-cell model
pdf("~/Desktop/sickle_fcell_rf/Figures/over10fcell_all_relativeimportance_mahatten.pdf", width = 13.33, height = 7.5)
manhattan(postions, chr = "Chr", bp = "Start", snp = "rsID",
          p = "relative_importance", logp = F, ylab = "Relative Importance",
          col = c("black", "red"), genomewideline = F, suggestiveline = F, 
          main = "Manhattan plot of the relative importance of the SNPs in the Model")
graphics.off()

# Manhatten plot for the healthy model
pdf("~/Desktop/sickle_fcell_rf/Figures/healthy_relativeimportance_mahatten.pdf", width = 13.33, height = 7.5)
manhattan(postions_healthy, chr = "Chr", bp = "Start", snp = "rsID",
          p = "relative_importance", logp = F, ylab = "Relative Importance",
          col = c("black", "red"), genomewideline = F, suggestiveline = F,
          main = "Manhattan plot of the relative importance of the SNPs in the Model")
graphics.off()

#### Annotate with dbSNP IDs ####
if (!requireNamespace("BiocManager", quietly = T))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

results <- c() # Initialise storage vector
snp_mart <- useMart(biomart = "ENSEMBL_MART_SNP", host = "grch37.ensembl.org", path = "/biomart/martservice", dataset = "hsapiens_snp")
for (snp in 1:dim(postions)[1]) { 
  # chrom <- gsub(pattern = 'chr', replacement = '', x = df[snp, 1]) # Remove 'chr' for searching
  temp <- getBM(attributes = c('refsnp_id', 'allele', 'chrom_start', 'chrom_strand'), 
                filters = c('chr_name', 'start', 'end'), 
                values = list(chrom, df[snp, 2], df[snp, 3]), mart )
  temp[, 5] <- snp  # Store SNP file index row next to each answer
  results <- rbind(results, temp)
}

# Using h2o to load in random forest models
h2o.init(max_mem_size = "4G", nthreads = -1)

# Plot 100 most important variables for the over 10% f-cell model
over10_model <- h2o.loadModel(path = "/home/callum/Desktop/sickle_fcell_rf/Models/sickle_model_over10fcell_nonlog/DRF_model_R_1608035413926_1")

pdf("~/Desktop/sickle_fcell_rf/Figures/over10_top100_relativeimportance.pdf", width = 13.33, height = 7.5)
as.data.frame(over10_model@model[["variable_importances"]]) %>%
  dplyr::arrange(desc(relative_importance)) %>%
  dplyr::top_n(n = 100, wt = relative_importance) %>%
  ggplot(aes(reorder(variable, relative_importance), relative_importance)) +
  geom_col() +
  coord_flip() +
  xlab("SNPs") +
  ylab("Relative Importance") +
  ggtitle("The 100 variables which most reduce the OOB RMSE") +
  theme(
    # Lengends to the top
    plot.title = element_text(hjust = 0.5),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 45, hjust = 1))
graphics.off()

#### Compare to the healthy model ####
healthy_model <- readRDS(file = "~/Desktop/sickle_fcell_rf/Models/model_1.rds")

# Plot 100 most important variables for the healthy model
pdf("~/Desktop/sickle_fcell_rf/Figures/heathy_top100_relativeimportance.pdf", width = 13.33, height = 7.5)
as.data.frame(healthy_model@model[["variable_importances"]]) %>%
  dplyr::arrange(desc(relative_importance)) %>%
  dplyr::top_n(n = 100, wt = relative_importance) %>%
  ggplot(aes(reorder(variable, relative_importance), relative_importance)) +
  geom_col() +
  coord_flip() +
  xlab("SNPs") +
  ylab("Relative Importance") +
  ggtitle("The 100 variables which most reduce the OOB RMSE") +
  theme(
    # Lengends to the top
    plot.title = element_text(hjust = 0.5),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 45, hjust = 1))
graphics.off()

# Pick out the top 100 rsIDs (ranked by relative importance)
as.data.frame(healthy_model@model[["variable_importances"]]) %>%
  dplyr::arrange(desc(relative_importance)) %>%
  dplyr::top_n(n = 1000, wt = relative_importance) %>% 
  {. ->> healthy_model_top_1000_predictors }

# Write results to output
write.csv(healthy_model_top_1000_predictors, file = "~/Desktop/sickle_fcell_rf/Healthy_Top_1000.csv")
healthy_model_top_1000_predictors <- read.csv(file = "~/Desktop/sickle_fcell_rf/Healthy_Top_1000.csv")
colnames(healthy_model_top_1000_predictors) <- c("rownames", "rsID", "relative_importance", "scaled_importance", "percentage")
healthy_model_top_1000_predictors$rsID <- substr(healthy_model_top_1000_predictors$rsID, 1, nchar(healthy_model_top_1000_predictors$rsID)-2) # Remove last 2 characters
combined <- merge(healthy_model_top_1000_predictors, merged, by = "rsID", all.x = T)

healthy_model_top_1000_predictors$rsID
sum(merged$percentage)

# Shutdown h2o instance
h2o.shutdown(prompt = F)

