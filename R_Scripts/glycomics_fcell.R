##### Load/install the required packages #####

# for tidyverse needed libxml2-dev, libssl-dev, libcurl4-openssl-dev in the OS
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T)) }

GetPackages(c("tidyverse", "randomForest", "dplyr", "igraph", "caret", "reshape", "reshape2",
              "devtools", "PerformanceAnalytics", "ggplot2", "car", "ggpubr", "lubridate",
              "bootstrap", "corrplot", "ggraph", "doParallel", "ranger", "data.table"))

# install_github is a devtools package, so load this first
devtools::install_github(c("kassambara/easyGgplot2", "MI2DataLab/randomForestExplainer"))
lapply(c("randomForestExplainer", "easyGgplot2"), require, character.only = T)

##### Import the datasets #####

Fcell_levels <- readr::read_csv("~/Twins/twin_IDs_Fcell_FACs.csv")
twins_fcell_glycomics_raw <- utils::read.delim(
  "~/Documents/twins_ML_project/Glycomics/glycans.igg.global.combat.scale.processed.txt")

##### Adding additional age information to the datasets #####

# Extract year of birth information to new data frame
colnames(Fcell_levels) <- c("PublicID", "FAM", "TwinType", "Sex", "YEAR_BIRTH", "ethnicity", "geom_mean_FCFACS",
                            "Sanger_Genome_Seq", "Glycomics", "metabolomic", "overlapping")
twins_fcell_glycomics <- base::merge(twins_fcell_glycomics_raw, Fcell_levels, by = "PublicID")
twins_fcell_glycomics <- twins_fcell_glycomics[!base::duplicated(twins_fcell_glycomics$PublicID), ]

# Remove duplicates from the list
drops <- c("PublicID", "batch", "plate", "date", "Year.Reported", "Month.Reported", "month_year_birth", 
           "BirthYear", "ID_visityear", "FAM", "TwinType", "Sex", "YEAR_BIRTH", "ethnicity", 
           "Sanger_Genome_Seq", "Glycomics", "metabolomic", "overlapping")
twins_fcell_glycomics <- twins_fcell_glycomics[,!(names(twins_fcell_glycomics) %in% drops)]

##### Visually testing for normailty with several plots #####

# Get variable names
glycans_only <- melt(twins_fcell_glycomics[1:6])
fcell_only <- melt(twins_fcell_glycomics[77])

# Plot (quartile values)
layout(matrix(c(1:1), 1, 1))
qqPlot(glycans_only$value)
qqPlot(fcell_only$value)

ggplot(glycans_only, aes(x = value, fill = variable)) + theme_minimal() + 
  geom_density(alpha = 0.1) +
  guides(fill = F) +
  scale_color_brewer(palette = "Dark2") +
  xlim(-5, 5)

ggplot(fcell_only, aes(x = value, fill = variable)) + theme_minimal() + 
  geom_density(alpha = 0.1) +
  guides(fill = F) +
  scale_color_brewer(palette = "Dark2")

# Probablity of your data occuring in a randomly sampled normal distribution
# Works for distrubtions with >5000 variables unlike some other algorithms
stats::ks.test(glycans_only$value, y = pnorm, alternative = "two.sided")
stats::ks.test(fcell_only$value, y = pnorm, alternative = "two.sided")

# Testing "normality" of variables being used as predictors in a regression model before the fit is unwarranted. 
# It does make sense to test the normality of the residuals since that is what is assumed in the modeling theory.

##### Random Forest Glycomics #####

# Convert to numeric (via character to avoid coverting factors to numeric storage values), remove NAs
twins_fcell_glycomics <- dplyr::mutate_all(twins_fcell_glycomics, function(x) as.numeric(as.character(x)))
twins_fcell_glycomics <- na.omit(twins_fcell_glycomics)

# Setup training and test sets

# Split the glycomics data into training and test sets
glycomics_fcell_split <- rsample::initial_split(twins_fcell_glycomics, prop = 0.8)
glycomics_fcell_train <- rsample::training(glycomics_fcell_split)
glycomics_fcell_test  <- rsample::testing(glycomics_fcell_split)

# Using h2o for grid optimisation and random forest generation, setup a cluster with 32Gb of RAM
h2o.init(max_mem_size = "32G")

# Turn training set into h2o object
train.h2o <- as.h2o(glycomics_fcell_train)

# Hyperparameter grid
hyper_grid.h2o <- list(
  ntrees = seq(200, 600, by = 100),
  mtries = seq(2, 5, by = 1),
  sample_rate = c(.55, .6, .65, .7)
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
x <- setdiff(names(twins_fcell_glycomics), y)

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
best_grid <- h2o.getGrid(
  grid_id = "rf_grid", 
  sort_by = "mse",
  decreasing = F
)

# Print results
print(best_grid)

# Try to get some plots of the optimisation metrics
plot(best_grid@summary_table[["mtries"]], best_grid@summary_table[["mse"]])
plot(best_grid@summary_table[["ntrees"]], best_grid@summary_table[["mse"]])
plot(best_grid@summary_table[["sample_rate"]], best_grid@summary_table[["mse"]])

# Grab the model_id for the top model, chosen by validation error
best_model_id <- best_grid@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)

# Now let’s evaluate the model performance on a test set
test.h2o <- as.h2o(glycomics_fcell_test)
best_model_performance <- h2o.performance(model = best_model, newdata = test.h2o)

# RMSE of best model
h2o.mse(best_model_performance) %>% sqrt()
best_model_performance

# Make a prediction from the h2o random forest on the test set
pred_h2o <- predict(best_model, test.h2o)
head(pred_h2o)

# Plot 25 most important variables
as.data.frame(best_model@model[["variable_importances"]]) %>%
  dplyr::arrange(desc(relative_importance)) %>%
  dplyr::top_n(n = 76, wt = relative_importance) %>%
  ggplot(aes(reorder(variable, relative_importance), relative_importance)) +
  geom_col() +
  coord_flip() +
  ggtitle("The 25 variables which most reduce the OOB RMSE") +
  theme_minimal()

h2o.shutdown(prompt = F)  # Shuts down the Java cluster

##### Glycomics Regression #####

# Multiple linear regression 
glycans <- paste(colnames(twins_fcell_glycomics)[-length(twins_fcell_glycomics)], collapse = " + ")
glycomic_regression <- stats::lm(paste("geom_mean_FCFACS ~ ", glycans, sep = ""),
                                 data = twins_fcell_glycomics, na.action = na.exclude)

# Diagnostic plots
layout(matrix(c(1:4), 2, 2)) # optional 4 graphs/page
plot(glycomic_regression)

# Split lm models into positive and negative, pick out top predictors
glycans <- paste(row.names(as.data.frame(sort(optimal_gly_model$variable.importance, decreasing = T)))[1:10], collapse = " + ")
glycomic_regression <- stats::lm(paste("geom_mean_FCFACS ~ ", glycans, sep = ""),
                                 data = twins_fcell_glycomics, na.action = na.exclude)
base::summary(glycomic_regression) # show results

# Diagnostic plots
layout(matrix(c(1:4), 2, 2)) # optional 4 graphs/page
plot(glycomic_regression)

# Pick out positive and negative predictors
a <- summary(glycomic_regression)
a <- as.data.frame(a$coefficients)
positives <- a[a$`t value` > 0, ]
negatives <- a[a$`t value` <= 0, ]

# Assessing R2 shrinkage using 10-Fold cross-validation

# Create matrices of predictors
x <- as.matrix(twins_fcell_glycomics[colnames(twins_fcell_glycomics)[-length(twins_fcell_glycomics)]])
y <- as.matrix(twins_fcell_glycomics[c("geom_mean_FCFACS")]) 

# Define theta fit/predict functions 
ThetaFit <- function(x, y){lsfit(x, y)}
ThetaPredict <- function(fit, x){cbind(1, x) %*% fit$coef}
crossval_results <- bootstrap::crossval(x, y, ThetaFit, ThetaPredict, ngroup = 10)
stats::cor(y, glycomic_regression$fitted.values)**2 # raw R^2 
stats::cor(y, crossval_results$cv.fit)**2 # cross-validated R2
