#!/usr/bin/env Rscript

##### Install required libraries

GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

GetPackages(c("tidyverse", "randomForest", "dplyr", "igraph", "caret", "reshape", "rsample",
              "reshape2", "devtools", "PerformanceAnalytics", "ggplot2", "lubridate",
              "bootstrap", "corrplot", "ggraph", "doParallel", "ranger", "data.table", "h2o",
              "sparsio"))

devtools::install_github("krlmlr/ulimit")

##### Load in the genotype matrices #####

# The input arguments should supply the PLINK .raw file location 
args <- commandArgs(trailingOnly = T)
# file_loaction <- args[1]

file_loaction <- "/media/callum/New Volume/Fcell_SNPs_RF/Data/Twins/1-4-T23-gwa610K.output.raw"

SNP_list <- fread(input = file_loaction)  # Load the plink data
SNP_list <- SNP_list_1_4[, grep("HET", colnames(SNP_list_1_4)):=NULL]  # Remove the HET columns (cut dataframe size in half)
SNP_Fcell_test <- fread(input = "~/Twins/FID_ID_FCellFACS.csv")  # Load fcell level information
SNP_merged_test <- merge(x = SNP_list_1_4, y = SNP_Fcell_test, by = "IID")  # Add the fcell information to the SNP data

rm(SNP_list_1_4, SNP_Fcell_test)
gc()

checkpoint_1 <- paste(dirname(file_loaction), "/temp_noHET_incFACS.rds", sep = "")

saveRDS(SNP_merged_test, file = checkpoint_1)
SNP_merged_test <- readRDS(file = checkpoint_1)

SNP_merged_test <- SNP_merged_test[, grep("AT", colnames(SNP_merged_test)):=NULL]
SNP_merged_test <- SNP_merged_test[, grep("SEX", colnames(SNP_merged_test)):=NULL]
SNP_merged_test <- SNP_merged_test[, grep("PHENOTYPE", colnames(SNP_merged_test)):=NULL]

# SNP_merged_test[,1:5][1,] # Checking columns names
# SNP_merged_test[,dim(SNP_merged_test_1)[2]-5:dim(SNP_merged_test_1)][1,] # Checking columns names

IDs <- SNP_merged_test[,1:2]
FACS_info <- SNP_merged_test[,dim(SNP_merged_test)[2]:dim(SNP_merged_test)[2]]
SNP_merged_test_1 <- SNP_merged_test[,3:round(dim(SNP_merged_test)[2]/3)]
SNP_merged_test_2 <- SNP_merged_test[,(round(dim(SNP_merged_test)[2]/3)+1):(2*round(dim(SNP_merged_test)[2]/3))]
SNP_merged_test_3 <- SNP_merged_test[,((2*round(dim(SNP_merged_test)[2]/3))+1):((dim(SNP_merged_test)[2])-1)]

# Part 1
SNP_merged_test_1 <- cbind(IDs, SNP_merged_test_1, FACS_info)

# Part 2
SNP_merged_test_1 <- cbind(IDs, SNP_merged_test_2, FACS_info)
rm(SNP_merged_test_2)

# Part 3
SNP_merged_test_1 <- cbind(IDs, SNP_merged_test_3, FACS_info)

SNP_merged_test_Twin1 <- SNP_merged_test_1[!duplicated(SNP_merged_test_1[, 2])]  # Separate the twins out, one in group 1
SNP_merged_test_Twin2 <- SNP_merged_test_1[duplicated(SNP_merged_test_1[, 2])]  # The other in group 2

SNP_merged_test_Twin1 <- SNP_merged_test_Twin1[, grep("ID", colnames(SNP_merged_test_Twin1)):=NULL]
SNP_merged_test_Twin2 <- SNP_merged_test_Twin2[, grep("ID", colnames(SNP_merged_test_Twin2)):=NULL]

rm(IDs)
rm(FACS_info)
rm(SNP_merged_test)
rm(SNP_merged_test_1)
gc()

NAs.to.zero <- function(DT) {
  # Call columns by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT, which(is.na(DT[[j]])), j, 0)
}

NAs.to.zero(SNP_merged_test_Twin1)
NAs.to.zero(SNP_merged_test_Twin2)

# Save/Load checkpoint 2
saveRDS(SNP_merged_test_Twin1, file = "~/Documents/twins_ML_project/plink/ML_project/1-4_SNP_merged_test_Twin1_part3.rds")  # Save locally as compact rds file
saveRDS(SNP_merged_test_Twin2, file = "~/Documents/twins_ML_project/plink/ML_project/1-4_SNP_merged_test_Twin2_part3.rds")  # Save locally as compact rds file

SNP_merged_test_Twin1 <- readRDS(file = "~/Documents/twins_ML_project/plink/ML_project/1-4_SNP_merged_test_Twin1_part1.rds")
SNP_merged_test_Twin2 <- readRDS(file = "~/Documents/twins_ML_project/plink/ML_project/1-4_SNP_merged_test_Twin2_part1.rds")

m2 <- randomForest::tuneRF(
  x = SNP_merged_test_Twin1,
  y = SNP_merged_test_Twin1$geom_mean_FCFACS,
  ntreeTry = 500,
  mtryStart = 100,
  stepFactor = 1.5,
  improve = 0.01,
  trace = T  # show real-time progress
)

# Split the glycomics data into training and test sets
# Twins1_SNPs_fcell_split <- rsample::initial_split(SNP_merged_test_Twin1, prop = 0.8)
# Twins1_SNPs_fcell_train <- rsample::training(Twins1_SNPs_fcell_split)
# Twins1_SNPs_fcell_test  <- rsample::testing(Twins1_SNPs_fcell_split)

Twins1_SNPs_fcell_train <- SNP_merged_test_Twin1
Twins1_SNPs_fcell_test <- SNP_merged_test_Twin2

rm(SNP_merged_test_Twin1)
rm(SNP_merged_test_Twin2)
rm(Twins1_SNPs_fcell_split)
gc()

# Using h2o for grid optimisation and random forest generation, setup a cluster with 16Gb of RAM and all cores
h2o.init(max_mem_size = "16G", nthreads = -1)

# Use all cores
# registerDoParallel(parallel::detectCores())

# Turn training set into h2o object (slow)
# train.h2o <- as.h2o(Twins1_SNPs_fcell_train)
# test.h2o <- as.h2o(Twins1_SNPs_fcell_test)

# as.h2o very slo, use h2o.importFile instead (better to write to disk then load into h2o)
# Is the data being clipped?
data.table::fwrite(x = Twins1_SNPs_fcell_train, file = "/home/callumrakhit/Documents/twins_ML_project/plink/ML_project/Twins1_SNPs_fcell_train_part3.csv")
data.table::fwrite(x = Twins1_SNPs_fcell_test, file = "/home/callumrakhit/Documents/twins_ML_project/plink/ML_project/Twins2_SNPs_fcell_test_part3.csv")

train.h2o <- h2o.importFile(path = "/home/callumrakhit/Documents/twins_ML_project/plink/ML_project/Twins1_SNPs_fcell_train_part3.csv")
test.h2o <- h2o.importFile(path = "/home/callumrakhit/Documents/twins_ML_project/plink/ML_project/Twins2_SNPs_fcell_test_part3.csv")

rm(Twins1_SNPs_fcell_test)
rm(Twins1_SNPs_fcell_train)
gc()

# head(names(test.h2o))
# tail(names(test.h2o))  # Check if geom_mean_FCFACS is the last column

# Hyperparameter grid
hyper_grid.h2o <- list(
  ntrees = seq(501, 801, by = 100),
  mtries = seq((ncol(Twins1_SNPs_fcell_test)/5), (ncol(Twins1_SNPs_fcell_test)/5) + 50000, by = 10000),
  sample_rate = c(.55, .6, .65, .7)
)

# Crappy hyperparameter grid
hyper_grid.h2o <- list(
  ntrees = 601,
  mtries = 50000,
  sample_rate = c(0.6)
)

# Random grid search criteria
search_criteria <- list(
  strategy = "RandomDiscrete",
  stopping_metric = "mse",
  stopping_tolerance = 0.005,
  stopping_rounds = 10,
  max_runtime_secs = 30*60
)

# Set values for x and y
y <- "geom_mean_FCFACS"
x <- setdiff(names(train.h2o), y)

# Need to enlarge the stack size in bash (ulimit -s 16384) - enlarge stack limit to 16Mb
ulimit::memory_limit(size = 16384)

# Identify stack limit
Cstack_info()["size"]

# Build the grid search 
grid_3 <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_grid_3",
  x = x, 
  y = y,  
  training_frame = train.h2o,
  hyper_params = hyper_grid.h2o,
  search_criteria = search_criteria
)

# Save the model locally
saveRDS(object = grid_2, file = "~/Documents/twins_ML_project/plink/ML_project/rf_3.rds")

# Load model
grid <- readRDS("~/Documents/twins_ML_project/plink/ML_project/rf_3")

# Collect the results and sort by our model performance metric of choice
best_grid <- h2o.getGrid(
  grid_id = "rf_grid_3", 
  sort_by = "mse",
  decreasing = F
)

# Print results
print(best_grid)

# Try to get some plots of the optimisation metrics
plot(best_grid@summary_table[["mtries"]], best_grid@summary_table[["mse"]])
plot(best_grid@summary_table[["ntrees"]], best_grid@summary_table[["mse"]])
plot(best_grid@summary_table[["sample_rate"]], best_grid@summary_table[["mse"]])

# Get the model_id for the top model, chosen by validation error
best_model_id <- best_grid@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)

# Evaluate the model performance on a test set
best_model_performance <- h2o.performance(model = best_model, newdata = test.h2o)

# RMSE of best model
h2o.mse(best_model_performance) %>% sqrt()
best_model_performance

# Make a prediction from the h2o random forest on the test set
pred_h2o <- predict(best_model, test.h2o)
head(pred_h2o)

best_model@model[["variable_importances"]]

# Plot 100 most important variables
as.data.frame(best_model@model[["variable_importances"]]) %>%
  dplyr::arrange(desc(relative_importance)) %>%
  dplyr::top_n(n = 100, wt = relative_importance) %>%
  ggplot(aes(reorder(variable, relative_importance), relative_importance)) +
  geom_col() +
  coord_flip() +
  ggtitle("The 100 variables which most reduce the OOB RMSE") +
  theme_minimal()

# Pick out the top 10,000 rsIDs (ranked by relative importance)
as.data.frame(best_model@model[["variable_importances"]]) %>%
  dplyr::arrange(desc(relative_importance)) %>%
  dplyr::top_n(n = 100, wt = relative_importance) %>% 
  {. ->> rf_1_top_100_predictors }

best_model@model[["variable_importances"]]

rf_1_top_100_predictors <- rf_1_top_100_predictors[,1]

View(rf_1_top_100_predictors)

write_csv(x = as.data.frame(rf_1_top_100_predictors), path = "~/Documents/twins_ML_project/plink/ML_project/rf_3_top100.csv")

h2o.shutdown(prompt = F)  # Shuts down the Java cluster

