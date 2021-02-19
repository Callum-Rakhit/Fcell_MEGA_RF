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
              "h2o", "sparsio"))

devtools::install_github("krlmlr/ulimit")

# Get lastest package versions
update.packages()

# Troubleshooting "Error in readRDS(pfile)". List the library paths, the issue is likely to be in the first directory
paths = .libPaths()

## Try and detect bad files
list.files(paths, pattern = "^00LOCK*|*\\.rds$|*\\.RDS$", full.names = T)

## List files of size 0
l = list.files(paths, full.names = T)
l[sapply(l, file.size) == 0]

# Load in the genotypes and FACS data
SNP_list <- fread(input = "~/Desktop/sickle_fcell_rf/Dataset_2.raw")  # Load the plink data
SNP_list <- SNP_list[, grep("HET", colnames(SNP_list)):=NULL]  # Remove the HET columns (cut dataframe size in half)
# SNP_list$IID <- as.numeric(gsub("ID-([0-9]+)", "\\1", SNP_list$IID)) # Convert ID-001 to 001 for matching with FACS
SNP_Fcell <- fread(input = "~/Desktop/sickle_fcell_rf/HbF.txt")  # Load fcell level information
# SNP_Fcell$IID <- as.numeric(gsub("ID-([0-9]+)", "\\1", SNP_list$IID)) # Convert ID-001 to 001 format

# Invert the natural log for the LnHbF
SNP_Fcell <- SNP_Fcell[(SNP_Fcell$LnHbF != -9)]
SNP_Fcell$LnHbF <- exp(SNP_Fcell$LnHbF)
hist(SNP_Fcell$LnHbF)

# Merge and remove unwanted columns
SNP_merged <- merge(x = SNP_list, y = SNP_Fcell, by = "IID")  # Add the fcell information to the SNP data
SNP_merged <- SNP_merged[, grep("AT", colnames(SNP_merged)):=NULL]
SNP_merged <- SNP_merged[, grep("SEX", colnames(SNP_merged)):=NULL]
SNP_merged <- SNP_merged[, grep("PHENOTYPE", colnames(SNP_merged)):=NULL]

# Remove IID and FID
SNP_merged$FID.x <- NULL  # Remove FID.x
SNP_merged <- subset(SNP_merged, select = -c(IID, FID.y))

# Change NAs to 0 (missing genotypes) as not tollerated by random forest
NAs.to.zero <- function(DT) {
  # Call columns by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT, which(is.na(DT[[j]])), j, 0)
}

NAs.to.zero(SNP_merged)

# Save data for easy loading
saveRDS(SNP_merged, file = "~/Desktop/sickle_fcell_rf/Dataset_2_nonlog.rds")
SNP_merged <- readRDS(file = "~/Desktop/sickle_fcell_rf/Dataset_2_nonlog.rds")

# Remove missing LnHbF values
SNP_merged <- SNP_merged[(SNP_merged$LnHbF != -9)]

# Create train and test datasets
Sickle_data_train <- SNP_merged[1:760,]  # All samples
Sickle_data_test <- sample_n(SNP_merged, 304) # Random 40%

# Clear memory
gc()

# Basic random forest
m2 <- randomForest::tuneRF(
  x = Sickle_data_train,
  y = Sickle_data_train$LnHbF,
  ntreeTry = 501,
  mtryStart = 100,
  stepFactor = 1.5,
  improve = 0.01,
  trace = T  # show real-time progress
)

# An mtry = 757 gave an OOB error of 9.972194 

# install.packages("randomForestExplainer")
library(randomForestExplainer)
library(randomForest)

# Subset the data into smaller packages
SNP_merged_test_1 <- SNP_merged_test[,3:round(dim(SNP_merged_test)[2]/3)] # First 1/3
Sickle_data_train_1 <- Sickle_data_train[,3:round(dim(Sickle_data_train)[2]/3)]
Sickle_data_train_1$LnHbF <- Sickle_data_train$LnHbF

set.seed(2020)
forest <- randomForest(LnHbF ~ ., mtry = 1000, ntree = 501, data = Sickle_data_train, do.trace = T)
save(forest, file = "~/Stephan_SNP_HbSS_Dataset/Sickle_Model.RData")
forest_load_test = get(load("~/Stephan_SNP_HbSS_Dataset/Sickle_Model.RData"))
explain_forest(forest, interactions = T, data = Sickle_data_train)

# Using h2o for grid optimisation and random forest generation, setup a cluster with 16Gb of RAM and all cores
h2o.init(max_mem_size = "4G", nthreads = -1)

# as.h2o very slo, use h2o.importFile instead (better to write to disk then load into h2o). Data clipped - maximum of ~265,000 columns?
data.table::fwrite(x = Sickle_data_train, file = "/home/callum/Desktop/sickle_fcell_rf/Sickle_SNPs_train_nonlog.csv")
data.table::fwrite(x = Sickle_data_test, file = "/home/callum/Desktop/sickle_fcell_rf/Sickle_SNPs_test_nonlog.csv")

train.h2o <- h2o.importFile(path = "/home/callum/Desktop/sickle_fcell_rf/Sickle_SNPs_train_nonlog.csv")
test.h2o <- h2o.importFile(path = "/home/callum/Desktop/sickle_fcell_rf/Sickle_SNPs_test_nonlog.csv")

y <- "LnHbF"
x <- setdiff(names(train.h2o), y)

# Need to enlarge the stack size in bash (ulimit -s 16384) - enlarge stack limit to 16Mb
ulimit::memory_limit(size = 4000)

# Create the model
sickle_model_nonlog <- h2o.randomForest(
  x = x, 
  y = y,
  ntrees = 501,
  mtries = 1000,
  # sample_rate = c(0.6),
  training_frame = train.h2o,
  # calibrate_model = T,
  # calibration_frame = test.h2o
)

# Save the model locally
h2o.saveModel(object = sickle_model_nonlog, path = "/home/callum/Desktop/sickle_fcell_rf/Models/sickle_model_nonlog")
sickle_model_nonlog <- h2o.loadModel(path = "/home/callum/Desktop/sickle_fcell_rf/Models/sickle_model_nonlog/DRF_model_R_1605876137250_1")

# Predict accuracy
predicted_fcell_twins <- h2o.predict(sickle_model_nonlog, newdata = test.h2o)
test_values <- read.csv(file = "/home/callum/Desktop/sickle_fcell_rf/Sickle_SNPs_test_nonlog.csv")
pred_vs_actual_test <- as.data.frame(h2o.cbind(test.h2o$LnHbF, predicted_fcell_twins))

# Save predicted vs actual
write.csv(pred_vs_actual, file = "~/Desktop/sickle_fcell_rf/Sickle_Top_NonLog_Test_PredvsAct.csv")

hist(test_values$LnHbF)
hist(pred_vs_actual$LnHbF)
hist(pred_vs_actual$predict)
plot(pred_vs_actual$LnHbF, pred_vs_actual$predict)
plot(pred_vs_actual_test$LnHbF, pred_vs_actual_test$predict)
pred_vs_actual_test$rescaled <- rescale(x = pred_vs_actual_test$predict, to = c(0, 25))
plot(pred_vs_actual_test$LnHbF, rescale(x = pred_vs_actual_test$predict, to = c(0, 25)))

library(ggplot2)
ggplot(pred_vs_actual_test, aes(x = LnHbF, y = rescaled)) +
  xlab("HbF Actual") +
  ylab("HbF Predicted") +
  geom_point() + 
  geom_smooth() +
  theme_bw() +
  theme(legend.position = "top")

over_ten <- pred_vs_actual_test[pred_vs_actual_test$rescaled >= 10,]
wrong <- over_ten[over_ten$LnHbF < 10,]
(1-(21/74))*100


install.packages("scales")
library(scales)
rescale(x = pred_vs_actual_test$predict, to = c(0, 25))

pred_vs_actual_test$predict

performance_fcell_twins <- h2o.performance(sickle_model_nonlog, test.h2o)
performance_fcell_twins
test.h2o
mean(predicted_fcell_twins)

hist(SNP_merged$LnHbF)

# Get the best model
best_grid <- h2o.getGrid(grid_id = "sickle_model_nonlog", sort_by = "mse", decreasing = F)
best_model_id <- sickle_model_nonlog@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)

# Plot 100 most important variables
as.data.frame(sickle_model_nonlog@model[["variable_importances"]]) %>%
  dplyr::arrange(desc(relative_importance)) %>%
  dplyr::top_n(n = 100, wt = relative_importance) %>%
  ggplot(aes(reorder(variable, relative_importance), relative_importance)) +
  geom_col() +
  coord_flip() +
  ggtitle("The 100 variables which most reduce the OOB RMSE") +
  theme_minimal()

# Pick out the top 100 rsIDs (ranked by relative importance)
as.data.frame(sickle_model_nonlog@model[["variable_importances"]]) %>%
  dplyr::arrange(desc(relative_importance)) %>%
  dplyr::top_n(n = 1000, wt = relative_importance) %>% 
  {. ->> sickle_model_top_1000_predictors }

# Write results to output
write.csv(sickle_model_top_1000_predictors, file = "~/Desktop/sickle_fcell_rf/Sickle_Top_NonLog_1000.csv")

# Shutdown h2o instance
h2o.shutdown()

