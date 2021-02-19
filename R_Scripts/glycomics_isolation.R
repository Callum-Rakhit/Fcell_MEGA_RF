# TODO(Callum)
#  Produce the explain random forest plots manually
#  Import a SNP subset/any part of the SNP data
#  Sync rf "template" across all scripts

##### Load/install the required packages #####

# for tidyverse needed libxml2-dev, libssl-dev, libcurl4-openssl-dev in the OS
GetPackages <- function(required.packages) {
  packages.not.installed <- 
    required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){
    install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

options(pkgType = "source")
install.packages("tidyverse", dependencies = T) 
install.packages("caret", dependencies = T)
install.packages("ggpubr", dependencies = T)

library(tidyverse)

GetPackages(c("tidyverse", "randomForest", "dplyr", "igraph", "caret", "reshape", "reshape2",
              "devtools", "PerformanceAnalytics", "ggplot2", "car", "ggpubr", "lubridate",
              "bootstrap", "corrplot", "ggraph", "doParallel", "ranger", "data.table"))

# install_github is a devtools package, so load this first
devtools::install_github(c("kassambara/easyGgplot2", "MI2DataLab/randomForestExplainer"))
lapply(c("randomForestExplainer", "easyGgplot2"), require, character.only = T)

##### Import the datasets #####
twins_data_glycomics_raw <- readr::read_delim(
  "~/Documents/twins_ML_project/Glycomics/glycans.igg.global.combat.scale.processed.txt", delim = "\t")

# Need for DOB
twins_data_metabolomics_raw <- readr::read_delim(
  "~/Documents/twins_ML_project/metabolon_2015_scaleRunDayMedian_imputeRunDayMin_normInverse/metabolon_2015_scaleRunDayMedian_imputeRunDayMin_normInverse.txt",
  delim = "\t")

# Split the data into twins
SNP_Fcell_test <- fread(input = "~/Twins/FID_ID_FCellFACS.csv")  # Load fcell level information
SNP_merged_test <- merge(x = twins_data_glycomics_raw, y = SNP_Fcell_test, by = "IID")  # Add the fcell information to the SNP data

rm(SNP_list_1_4)
rm(SNP_Fcell_test)
gc()

saveRDS(SNP_merged_test, file = "~/Documents/twins_ML_project/plink/ML_project/1-4_dataset_noHET_FACSinc.rds")
SNP_merged_test <- readRDS(file = "~/Documents/twins_ML_project/plink/ML_project/1-4_dataset_noHET_FACSinc.rds")

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


##### Adding additional age information to the datasets #####

# Extract year of birth information to new data frame
twins_data_demographics <- base::data.frame(twins_data_metabolomics_raw$PublicID, twins_data_metabolomics_raw$month_year_birth)
colnames(twins_data_demographics) <- c("PublicID", "BirthYear")
twins_data_demographics$BirthYear <- lubridate::dmy(twins_data_demographics$BirthYear)
twins_data_demographics$BirthYear <- format(as.Date(twins_data_demographics[["BirthYear"]]), "%Y")
twins_data_demographics$BirthYear <- twins_data_demographics$BirthYear - 2012
twins_data_demographics$YearOfTest <- 2012
twins_data_demographics$Age <- twins_data_demographics$YearOfTest - as.numeric(twins_data_demographics$BirthYear)

ggplot(data=twins_data_demographics, aes(twins_data_demographics$Age)) + 
  geom_histogram(breaks=seq(45, 90, by = 1), 
                 col="black", 
                 fill="red", 
                 alpha = .2) + 
  labs(title="Histogram for Age") +
  labs(x="Age", y="Count") +
  theme(# Lengends to the top
    legend.position = "none",
    # Remove the y-axis
    axis.title.y = element_blank(),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major.x = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size = .25, color = "black"),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 0, hjust = 0))

ks.test(twins_data_demographics$Age, pnorm)

# Remove duplicates from the list
twins_data_demographics <- twins_data_demographics[!base::duplicated(twins_data_demographics$PublicID), ]

# Merge information into glycomics, calculate individuals age at time of testing
twins_data_glycomics <- base::merge(twins_data_glycomics_raw, twins_data_demographics, by = "PublicID")
twins_data_glycomics$date <- base::as.Date(twins_data_glycomics$date, format = "%d/%m/%Y")
twins_data_glycomics$age_at_test <- base::as.numeric(format(as.Date(twins_data_glycomics[["date"]]), "%Y")) - base::as.numeric(twins_data_glycomics$BirthYear)
drops <- c("PublicID", "batch", "plate", "date", "Year.Reported", "Month.Reported", "month_year_birth", "BirthYear", "ID_visityear")
twins_data_glycomics <- twins_data_glycomics[,!(names(twins_data_glycomics) %in% drops)]

##### Visually testing for normailty with several plots #####

# Get variable names
glycans_only <- melt(twins_data_glycomics[1:76])
age_only <- melt(twins_data_glycomics[77])

# Plot (quartile values)
layout(matrix(c(1:1), 1, 1))
qqPlot(glycans_only$value)
qqPlot(age_only$value)

# Glycan values histogram
ggplot(glycans_only, aes(x = value, fill = variable)) + theme_minimal() + 
  geom_density(alpha = 0.1) +
  guides(fill = F) +
  scale_color_brewer(palette = "Dark2") +
  xlim(-5, 5)

# Histogram of ages
ggplot(age_only, aes(x = value, fill = variable)) + theme_minimal() + 
  geom_density(alpha = 0.1) +
  guides(fill = F) +
  scale_color_brewer(palette = "Dark2")

# Probablity of data occuring from a randomly sampled normal distribution
# Works for distrubtions with >5000 variables unlike some other algorithms
stats::ks.test(glycans_only$value, y = pnorm, alternative = "two.sided")
stats::ks.test(age_only$value, y = pnorm, alternative = "two.sided")

# Testing "normality" of variables being used as predictors in a regression model before the fit is unwarranted. 
# It does make sense to test the normality of the residuals since that is what is assumed in the modeling theory.

##### Random Forest Glycomics #####

# Convert to numeric (via character to avoid coverting factors to numeric storage values), remove NAs
twins_data_glycomics <- dplyr::mutate_all(twins_data_glycomics, function(x) as.numeric(as.character(x)))
twins_data_glycomics <- na.omit(twins_data_glycomics)

# Age as a numerical factor, standard random forest, default mtry value is max(floor(ncol(x)/3), nodesize is 1 or 5 if factor
rf_gly_1 <- randomForest::randomForest(age_at_test ~ ., data = twins_data_glycomics, importance = T)

# Make the html reports and pull out other useful factors (broken)
# randomForestExplainer::explain_forest(rf_met_1, interactions = T, data = twins_data_glycomics)
# randomForest::getTree(rf_met_1, labelVar = T)
# top_25_met_predictors_rf <- randomForestExplainer::important_variables(rf_met_1, k = 25)
# plot(rf_met_1)  # Plot cumulative mse vs number of trees
# which.min(rf_met_1$mse)  # Find tree with the lowest mse
# sqrt(rf_met_1$mse[which.min(rf_met_1$mse)])  # Root mean square error (SD of residuals) of this optimal random forest

# Split the glycomics data into training and test sets
glycomics_split <- rsample::initial_split(twins_data_glycomics, prop = .7)

# Training data
glycomics_train <- rsample::training(twins_data_glycomics)

# Test data
glycomics_test  <- rsample::testing(twins_data_glycomics)

# Set age_at_test as the variable to be predicted
x_test <- glycomics_train[setdiff(names(glycomics_train), "age_at_test")]
y_test <- glycomics_train$age_at_test

rf_oob_comp <- randomForest(
  formula = age_at_test ~ .,
  data    = glycomics_train,
  xtest   = x_test,
  ytest   = y_test
)

# Extract OOB & validation errors
oob <- sqrt(rf_oob_comp$mse)
validation <- sqrt(rf_oob_comp$test$mse)

# Compare error rates
tibble::tibble(
  `Out of Bag Error` = oob,
  `Test error` = validation,
  ntrees = 1:rf_oob_comp$ntree) %>%
  gather(Metric, RMSE, -ntrees) %>%
  ggplot(aes(ntrees, RMSE, color = Metric)) +
  geom_line() +
  scale_y_continuous() +
  xlab("Number of trees")

# Tuning features using random forest
features <- setdiff(names(twins_data_glycomics), "age_at_test")

m2 <- randomForest::tuneRF(
  x = twins_data_glycomics[features],
  y = twins_data_glycomics$age_at_test,
  ntreeTry = 500,
  mtryStart = 5,
  stepFactor = 1.5,
  improve = 0.01,
  trace = T  # show real-time progress 
)

# Hyperparameter grid search (use input from tuneRF)
hyper_grid <- expand.grid(
  mtry       = seq(4, 7, by = 1),
  node_size  = seq(2, 9, by = 1),
  sampe_size = c(.55, .60, .65, .70, .75, .80) #,
  # OOB_RMSE  = 0
)

# Multiple cores
registerDoParallel(parallel::detectCores())

# Perform another grid search using ranger random forest
for(i in 1:nrow(hyper_grid)) {
  
  # train model
  test_gly_model <- ranger(
    formula         = age_at_test ~ ., 
    data            = twins_data_glycomics, 
    num.trees       = 500,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$node_size[i],
    sample.fraction = hyper_grid$sampe_size[i]
  )
  
  # add OOB error to grid
  hyper_grid$OOB_RMSE[i] <- sqrt(test_gly_model$prediction.error)
}

hyper_grid %>% 
  dplyr::arrange(OOB_RMSE) %>%
  head(10)

# Plot histogram of the error rates on a tuned ranger rf
hist(hyper_grid$OOB_RMSE, breaks = 20)

# Now have the best perameters, use these to construct a ranger random forest

OOB_RMSE <- vector(mode = "numeric", length = 1)

for(i in seq_along(OOB_RMSE)) {
  
  optimal_gly_model <- ranger(
    formula         = age_at_test ~ ., 
    data            = twins_data_glycomics, 
    num.trees       = 500,
    mtry            = 4,
    min.node.size   = 3,
    sample.fraction = .55,
    importance      = 'impurity'
  )
  
  OOB_RMSE[i] <- sqrt(optimal_gly_model$prediction.error)
}

hist(OOB_RMSE, breaks = 20)

print(optimal_gly_model)

# ROC curve

require(pROC)
rf.roc <- roc(twins_data_glycomics$age_at_test, optimal_gly_model$votes[,2])
plot(rf.roc)
auc(rf.roc)

# Then plot the top 25 most important variables across the random forests
optimal_gly_model$variable.importance %>% 
  tidy()%>%
  dplyr::arrange(desc(x)) %>%
  dplyr::top_n(25) %>%
  ggplot(aes(reorder(names, x), x)) +
  geom_col() +
  coord_flip() +
  ggtitle("The 25 variables which most reduce the OOB RMSE") +
  theme_minimal()

##### Using h2o for grid optimisation and random forest generation

h2o.init(max_mem_size = "32G")

# Turn training set into h2o object
train.h2o <- as.h2o(metabolomics_train)

# Hyperparameter grid
hyper_grid.h2o <- list(
  ntrees      = seq(100, 500, by = 100),
  mtries      = seq(10, 110, by = 10),
  sample_rate = c(.55, .625, .70, .775)
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
y <- "age_at_test"
x <- setdiff(names(twins_data_metabolomics), y)

# Build the grid search 
grid <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_grid",
  x = x, 
  y = y, 
  training_frame = train.h2o,
  hyper_params = hyper_grid.h2o,
  search_criteria = search_criteria
)

# Collect the results and sort by our model performance metric of choice
grid_perf <- h2o.getGrid(
  grid_id = "rf_grid_2", 
  sort_by = "mse", 
  decreasing = F
)

# Print results
print(grid_perf)

# Try to get some plots of the optimisation metrics
plot(grid_perf@summary_table[["mtries"]], grid_perf@summary_table[["mse"]])
plot(grid_perf@summary_table[["ntrees"]], grid_perf@summary_table[["mse"]])
plot(grid_perf@summary_table[["sample_rate"]], grid_perf@summary_table[["mse"]])

# Grab the model_id for the top model, chosen by validation error
best_model_id <- grid_perf@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)

# Now let’s evaluate the model performance on a test set
metabolomics_train.h2o <- as.h2o(metabolomics_train)
metabolomics_test.h2o <- as.h2o(metabolomics_test)
best_model_perf <- h2o.performance(model = best_model, newdata = metabolomics_test.h2o)

# RMSE of best model
h2o.mse(best_model_perf) %>% sqrt()

# Make a prediction from the h2o random forest on the test set
pred_h2o <- predict(best_model, metabolomics_test.h2o)
head(pred_h2o)

# Plot 25 most important variables
as.data.frame(best_model@model[["variable_importances"]]) %>%
  dplyr::arrange(desc(relative_importance)) %>%
  dplyr::top_n(n = 25, wt = relative_importance) %>%
  ggplot(aes(reorder(variable, relative_importance), relative_importance)) +
  geom_col() +
  coord_flip() +
  ggtitle("The 25 variables which most reduce the OOB RMSE") +
  theme_minimal()

h2o.shutdown()  # Shuts down the Java cluster

##### Glycomics Regression #####

# Multiple linear regression 
glycans <- paste(colnames(twins_data_glycomics)[-length(twins_data_glycomics)], collapse = " + ")
glycomic_regression <- stats::lm(paste("age_at_test ~ ", glycans, sep = ""),
                                 data = twins_data_glycomics, na.action = na.exclude)
base::summary(glycomic_regression) # show results

# Diagnostic plots
layout(matrix(c(1:4), 2, 2)) # optional 4 graphs/page
plot(glycomic_regression)

# Split lm models into positive and negative, pick out top predictors
glycans <- paste(row.names(as.data.frame(sort(optimal_gly_model$variable.importance, decreasing = T)))[1:10], collapse = " + ")
glycomic_regression <- stats::lm(paste("age_at_test ~ ", glycans, sep = ""),
                                 data = twins_data_glycomics, na.action = na.exclude)
base::summary(glycomic_regression) # show results

# Diagnostic plots
layout(matrix(c(1:4), 2, 2)) # optional 4 graphs/page
plot(glycomic_regression)

# Pick out positive and negitive predictors
a <- summary(glycomic_regression)
a <- as.data.frame(a$coefficients)
positives <- a[a$`t value` > 0, ]
negitives <- a[a$`t value` <= 0, ]

# Assessing R2 shrinkage using 10-Fold cross-validation

# Create matrices of predictors
x <- as.matrix(twins_data_glycomics[colnames(twins_data_glycomics)[-length(twins_data_glycomics)]])
y <- as.matrix(twins_data_glycomics[c("age_at_test")]) 

# Define theta fit/predict functions 
ThetaFit <- function(x, y){lsfit(x, y)}
ThetaPredict <- function(fit, x){cbind(1, x) %*% fit$coef}
crossval_results <- bootstrap::crossval(x, y, ThetaFit, ThetaPredict, ngroup = 10)
stats::cor(y, glycomic_regression$fitted.values)**2 # raw R^2 
stats::cor(y, crossval_results$cv.fit)**2 # cross-validated R2

