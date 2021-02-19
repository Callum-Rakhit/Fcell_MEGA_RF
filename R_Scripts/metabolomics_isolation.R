##### Load/install the required packages #####

# for tidyverse needed libxml2-dev, libssl-dev, libcurl4-openssl-dev in the OS
GetPackages <- function(required.packages) {
  packages.not.installed <- 
    required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){
    install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}


GetPackages(c("tidyverse", "randomForest", "dplyr", "igraph", "caret", "reshape", "reshape2", "devtools", 
              "rsample", "PerformanceAnalytics", "ggplot2", "car", "ggpubr", "lubridate", "broom", 
              "bootstrap", "h2o", "corrplot", "ggraph", "doParallel", "ranger", "data.table"))

# install_github is a devtools package, so load this first
devtools::install_github(c("kassambara/easyGgplot2", "MI2DataLab/randomForestExplainer", "cardiomoon/ggiraphExtra"))
base::lapply(c("easyGgplot2", "randomForestExplainer", "ggiraphExtra"), base::require, character.only = T)

##### Import the datasets #####

twins_data_metabolomics <- read.delim(
  "~/Documents/twins_ML_project/metabolon_2015_scaleRunDayMedian_imputeRunDayMin_normInverse/metabolon_2015_scaleRunDayMedian_imputeRunDayMin_normInverse.txt")

##### Adding additional age information to the datasets #####

# Extract year of birth information to new data frame
twins_data_demographics <- base::data.frame(twins_data_metabolomics$PublicID, twins_data_metabolomics$month_year_birth)
colnames(twins_data_demographics) <- c("PublicID", "BirthYear")
twins_data_demographics$BirthYear <- dmy(twins_data_demographics$BirthYear)
twins_data_demographics$BirthYear <- format(as.Date(twins_data_demographics[["BirthYear"]]), "%Y")

# Remove duplicates from the list
twins_data_demographics <- twins_data_demographics[!duplicated(twins_data_demographics$PublicID), ]

# Merge information into metabolomics
twins_data_metabolomics <- merge(twins_data_metabolomics, twins_data_demographics, by = "PublicID")
twins_data_metabolomics$age_at_test <- as.numeric(twins_data_metabolomics$year_hli) - as.numeric(twins_data_metabolomics$BirthYear)
drops <- c("PublicID", "month_year_birth.x", "age_hli", "sex", "zygosity", "bmi_hli", "date_hli", "year_hli", 
           "visit_hli", "batch_hli", "Year.Reported", "Month.Reported", "month_year_birth", "BirthYear")
twins_data_metabolomics <- twins_data_metabolomics[ , !(names(twins_data_metabolomics) %in% drops)]

##### Visually testing for normailty with several plots #####

# Visual
metabolites_only <- melt(twins_data_metabolomics[1:756])
age_only <- melt(twins_data_metabolomics[757])
layout(matrix(c(1:1), 1, 1))
qqPlot(metabolites_only$value)
qqPlot(age_only$value)

ggplot(metabolites_only, aes(x = value, fill = variable)) + theme_minimal() + 
  geom_density(alpha = 0.1) +
  guides(fill = F) +
  scale_color_brewer(palette = "Dark2") +
  xlim(-4, 4)

ggplot(age_only, aes(x = value, fill = variable)) + theme_minimal() + 
  geom_density(alpha = 0.1) +
  guides(fill = F) +
  scale_color_brewer(palette = "Dark2")

# Statistical (probablity of your data occuring from a randomly sampled normal distribution)
stats::ks.test(metabolites_only$value, y = pnorm, alternative = "two.sided") # Metabolics is normal
stats::ks.test(age_only$value, y = pnorm, alternative = "two.sided")

# Testing "normality" of variables being used as predictors in a regression model before the fit is unwarranted. 
# It does make sense to test the normality of the residuals since that is what is assumed in the modeling theory.

##### Random Forest Metabolomics #####

# Convert to numeric (via character to avoid coverting factors to numeric storage values), remove NAs
twins_data_metabolomics_1 <- dplyr::mutate_all(twins_data_metabolomics, function(x) as.numeric(as.character(x)))
twins_data_metabolomics <- na.omit(twins_data_metabolomics)

# Age as a discrete factor, standard random forest
rf_met_1 <- randomForest::randomForest(as.factor(age_at_test) ~ ., data = twins_data_metabolomics, importance = TRUE) 

# Make the html reports and pull out other useful factors (broken)
# randomForestExplainer::explain_forest(rf_met_1, interactions = T, data = twins_data_metabolomics)
# randomForest::getTree(rf_met_1, labelVar = T)
# top_25_met_predictors_rf <- randomForestExplainer::important_variables(rf_met_1, k = 25)
# plot(rf_met_1)  # Plot cumulative mse vs number of trees
# which.min(rf_met_1$mse)  # Find tree with the lowest mse
# sqrt(rf_met_1$mse[which.min(rf_met_1$mse)])  # Root mean square error (SD of residuals) of this optimal random forest

# Split the glycomics data into training and test sets
metabolomics_split <- rsample::initial_split(twins_data_metabolomics, prop = .7)

# Training data
metabolomics_train <- rsample::training(metabolomics_split)

# Test data
metabolomics_test  <- rsample::testing(metabolomics_split)

# Set age_at_test as the variable to be predicted
x_test <- metabolomics_train[setdiff(names(metabolomics_train), "age_at_test")]
y_test <- metabolomics_train$age_at_test

rf_oob_comp <- randomForest(
  formula = age_at_test ~ .,
  data    = metabolomics_train_v1,
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
  ntrees = 1:rf_oob_comp$ntree
  ) %>%
  gather(Metric, RMSE, -ntrees) %>%
  ggplot(aes(ntrees, RMSE, color = Metric)) +
  geom_line() +
  scale_y_continuous() +
  xlab("Number of trees")

# Tuning features using random forest
features <- setdiff(names(twins_data_metabolomics), "age_at_test")

m2 <- randomForest::tuneRF(
  x = twins_data_metabolomics[features],
  y = twins_data_metabolomics$age_at_test,
  ntreeTry = 500,
  mtryStart = 5,
  stepFactor = 1.5,
  improve = 0.01,
  trace = T  # show real-time progress 
)

# Hyperparameter grid search (use input from tuneRF)
hyper_grid_met <- expand.grid(
  mtry       = seq(150, 180, by = 5),
  node_size  = seq(4, 12, by = 2),
  sampe_size = c(.55, .60, .65, .70, .75, .80),
  OOB_RMSE  = 0  # Set this as empty to fill during the subsequent for loop
)

# Multiple cores
registerDoParallel(parallel::detectCores())

# Perform another grid search using ranger random forest
for(i in 1:nrow(hyper_grid_met)) {
  
  # train model
  test_met_model <- ranger(
    formula         = age_at_test ~ ., 
    data            = twins_data_metabolomics, 
    num.trees       = 500,
    mtry            = hyper_grid_met$mtry[i],
    min.node.size   = hyper_grid_met$node_size[i],
    sample.fraction = hyper_grid_met$sampe_size[i]
  )
  
  # add OOB error to the hyperparameter grid
  hyper_grid_met$OOB_RMSE[i] <- sqrt(test_met_model$prediction.error)
}

hyper_grid_met %>% 
  dplyr::arrange(OOB_RMSE) %>%
  head(10)

# Plot histogram of the error rates on a tuned ranger rf
hist(hyper_grid_met$OOB_RMSE, breaks = 20)

# Now have the best perameters, use these to construct a ranger random forest
OOB_RMSE <- vector(mode = "numeric", length = 1)

for(i in seq_along(OOB_RMSE)) {
  
  optimal_met_model <- ranger(
    formula         = age_at_test ~ ., 
    data            = twins_data_metabolomics, 
    num.trees       = 500,
    mtry            = 180,
    min.node.size   = 10,
    sample.fraction = .75,
    importance      = 'impurity'
  )
  
  OOB_RMSE[i] <- sqrt(optimal_met_model$prediction.error)
}

hist(OOB_RMSE, breaks = 20)

# Then plot the top 25 most important variables across the random forests
optimal_met_model$variable.importance %>% 
  tidy()%>%
  dplyr::arrange(desc(x)) %>%
  dplyr::top_n(25) %>%
  ggplot(aes(reorder(names, x), x)) +
  geom_col() +
  coord_flip() +
  ggtitle("The 25 variables which most reduce the OOB RMSE") +
  theme_minimal()

# one-hot encode our categorical variables
one_hot <- dummyVars(~ ., ames_train, fullRank = FALSE)
ames_train_hot <- predict(one_hot, ames_train) %>% as.data.frame()

# make ranger compatible names
names(ames_train_hot) <- make.names(names(ames_train_hot), allow_ = FALSE)

# hyperparameter grid search --> same as above but with increased mtry values
hyper_grid_2 <- expand.grid(
  mtry       = seq(50, 200, by = 25),
  node_size  = seq(3, 9, by = 2),
  sampe_size = c(.55, .632, .70, .80),
  OOB_RMSE  = 0
)

# perform grid search
for(i in 1:nrow(hyper_grid_2)) {
  
  # train model
  model <- ranger(
    formula         = Sale.Price ~ ., 
    data            = ames_train_hot, 
    num.trees       = 500,
    mtry            = hyper_grid_2$mtry[i],
    min.node.size   = hyper_grid_2$node_size[i],
    sample.fraction = hyper_grid_2$sampe_size[i],
    seed            = 123
  )
  
  # add OOB error to grid
  hyper_grid_2$OOB_RMSE[i] <- sqrt(model$prediction.error)
}

hyper_grid_2 %>% 
  dplyr::arrange(OOB_RMSE) %>%
  head(10)

##### Using h2o for grid optimisation and random forest generation

h2o.init(max_mem_size = "32G")

y <- "age_at_test"
x <- setdiff(names(twins_data_metabolomics), y)

# Turn training set into h2o object
train.h2o <- as.h2o(metabolomics_train)

# Hyperparameter grid
hyper_grid.h2o <- list(
  ntrees      = seq(100, 500, by = 100),
  mtries      = seq(10, 110, by = 10),
  sample_rate = c(.55, .625, .70, .775)
)

# random grid search criteria
search_criteria <- list(
  strategy = "RandomDiscrete",
  stopping_metric = "mse",
  stopping_tolerance = 0.005,
  stopping_rounds = 10,
  max_runtime_secs = 30*60
)

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








##### Metabolomics Regression #####

# Multiple Linear Regression Example
metabolites <- paste(colnames(twins_data_metabolomics)[-length(twins_data_metabolomics)], collapse = " + ")
fit_met <- stats::lm(paste("age_at_test ~ ", metabolites, sep = ""),
                     data = twins_data_metabolomics, na.action = na.exclude)
base::summary(fit_met) # show results

##### Split lm models into positive and negative, pick out top predictors #####

as.data.frame(best_model@model[["variable_importances"]]) %>%
  dplyr::arrange(desc(relative_importance)) %>%
  dplyr::top_n(n = 25, wt = relative_importance) %>%
  ggplot(aes(reorder(variable, relative_importance), relative_importance)) +
  geom_col() +
  coord_flip() +
  ggtitle("The 25 variables which most reduce the OOB RMSE") +
  theme_minimal()


as.data.frame(a[["coefficients"]])
a <- summary(fit_met)
a <- as.data.frame(a$coefficients)
positives <- a[a$`t value` > 0, ]
negitives <- a[a$`t value` <= 0, ]
positives <- data.table::setDT(positives[order(positives$`Pr(>|t|)`),][1:4,], keep.rownames = TRUE)[]
positive_predictors <- paste(positives$rn[2:4], collapse = " + ")
positives <- lm(paste("age_at_test ~ ", positive_predictors, sep = ""), 
                data = twins_data_metabolomics,
                na.action = na.exclude)
negitives <- data.table::setDT(negitives[order(negitives$`Pr(>|t|)`),][1:4,], keep.rownames = TRUE)[]
negitive_predictors <- paste(negitives$rn[2:4], collapse = " + ")
negitives <- lm(paste("age_at_test ~ ", negitive_predictors, sep = ""), 
                data = twins_data_metabolomics,
                na.action = na.exclude)
layout(matrix(c(1:4), 2, 2)) # optional 4 graphs/page
plot(positives)
plot(negitives)
ggPredict(positives, se = T,  interactive = T)
ggPredict(negitives, se = T,  interactive = T)
d <- twins_data_metabolomics

# Obtain predicted and residual values
d$predicted <- predict(positives)
d$residuals <- residuals(positives)

# Create plot
d %>% 
  gather(key = "iv", value = "x", -age_at_test, -predicted, -residuals) %>%
  ggplot(aes(x = x, y = age_at_test)) +
  geom_segment(aes(xend = x, yend = predicted), alpha = .2) +
  geom_point(aes(color = residuals)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  guides(color = FALSE) +
  geom_point(aes(y = predicted), shape = 1) +
  facet_grid(~ iv, scales = "free_x") +
  theme_minimal()

##### Assessing R2 shrinkage using 10-Fold Cross-Validation  #####

# matrix of predictors
x <- as.matrix(twins_data_met_dob[c(names(twins_data_met_dob)[2:757])])

# vector of predicted values
y <- as.matrix(twins_data_met_dob[c("age_at_test")]) 

# define functions 
theta.fit.met <- function(x, y){lsfit(x, y)}
theta.predict.met <- function(fit_met, x){cbind(1, x) %*% fit_met$coef} 

# Doesn't work?
results_met <- crossval(x, y, theta.fit.met, theta.predict.met, ngroup = 10)
results_met <- crossval(x, y, theta.fit.met, theta.predict.met, ngroup = 700)

cor(y, fit_met$fitted.values)**2 # raw R^2 
cor(y, results_met$cv.fit)**2 # cross-validated R2
