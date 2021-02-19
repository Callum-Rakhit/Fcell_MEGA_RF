# Function to install packages
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

# update.packages()
# 
# # List the library paths
# # The issue is likely to be in the first directory
# paths = .libPaths()
# 
# ## Try and detect bad files
# list.files(paths,
#            pattern = "^00LOCK*|*\\.rds$|*\\.RDS$",
#            full.names = T)
# 
# ## List files of size 0
# l = list.files(paths, full.names = T)
# l[sapply(l, file.size) == 0]


# Install required packages
required.packages <- c("tidyverse", "randomForest", "dplyr", "igraph", "caret", "reshape", "rsample", "reshape2",
                       "devtools", "PerformanceAnalytics", "ggplot2", "lubridate", "bootstrap", "corrplot",
                       "ggraph", "doParallel", "ranger", "data.table", "h2o", "sparsio", "lattice", "rsnps")

GetPackages(required.packages)
 
install.packages("tidyverse")
install.packages("rsnps")
install.packages("tidyverse", dependencies = T, INSTALL_opts = c('--no-lock'))
library(tidyverse)
library(rsnps)

# Get filenames for first twins and individuals without a twin
twin1_1 <- as.data.frame(Sys.glob(
  paths = "~/Documents/twins_ML_project/plink/ML_project/*1-4_SNP_merged_train_Twin1*"))
twin1_1$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*1-4_SNP_merged_train_Twin1*")` <- as.character(
  twin1_1$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*1-4_SNP_merged_train_Twin1*")`)
colnames(twin1_1) <- c("Twin1")

twin1_2 <- as.data.frame(Sys.glob(
  paths = "~/Documents/twins_ML_project/plink/ML_project/*5-10_SNP_merged_train_Twin1*"))
twin1_2$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*5-10_SNP_merged_train_Twin1*")` <- as.character(
  twin1_2$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*5-10_SNP_merged_train_Twin1*")`)
colnames(twin1_2) <- c("Twin1")

twin1_3 <- as.data.frame(Sys.glob(
  paths = "~/Documents/twins_ML_project/plink/ML_project/*11-22_SNP_merged_train_Twin1*"))
twin1_3$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*11-22_SNP_merged_train_Twin1*")` <- as.character(
  twin1_3$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*11-22_SNP_merged_train_Twin1*")`)
colnames(twin1_3) <- c("Twin1")

twin_1 <- as.data.frame(rbind(twin1_1, twin1_2, twin1_3))

# Get filenames for second twins and individuals without a twin
twin2_1 <- as.data.frame(Sys.glob(
  paths = "~/Documents/twins_ML_project/plink/ML_project/*1-4_SNP_merged_test_Twin2*"))
twin2_1$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*1-4_SNP_merged_test_Twin2*")` <- as.character(
  twin2_1$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*1-4_SNP_merged_test_Twin2*")`)
colnames(twin2_1) <- c("Twin2")

twin2_2 <- as.data.frame(Sys.glob(
  paths = "~/Documents/twins_ML_project/plink/ML_project/*5-10_SNP_merged_test_Twin2*"))
twin2_2$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*5-10_SNP_merged_test_Twin2*")` <- as.character(
  twin2_2$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*5-10_SNP_merged_test_Twin2*")`)
colnames(twin2_2) <- c("Twin2")

twin2_3 <- as.data.frame(Sys.glob(
  paths = "~/Documents/twins_ML_project/plink/ML_project/*11-22_SNP_merged_test_Twin2*"))
twin2_3$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*11-22_SNP_merged_test_Twin2*")` <- as.character(
  twin2_3$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*11-22_SNP_merged_test_Twin2*")`)
colnames(twin2_3) <- c("Twin2")

twin_2 <- as.data.frame(rbind(twin2_1, twin2_2, twin2_3))

# Load in the SNP files
files <- as.data.frame(Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*top100.*"))
files$Twin1 <- twin_1
files$Twin2 <- twin_2
colnames(files) <- c("top100", "Twin1", "Twin2")
files$top100 <- as.character(files$top100)

# Remove unwanted data
rm(list = ls(pattern = "^twin")) 

# Select the top SNPs
SNPselector <- function(files){
  top_100 <- read.csv(file = files[,1][1], header = T, stringsAsFactors = F)
  top_100[nrow(top_100) + 1,] <- c("geom_mean_FCFACS")
  SNPs_and_FACS_Twins1 <- readRDS(file = files[,2][1,])
  SNPs_and_FACS_Twins2 <- readRDS(file = files[,3][1,])
  SNPs_and_FACS_train_Twin1_filtered <<- subset(SNPs_and_FACS_Twins1, select = top_100[,1])
  SNPs_and_FACS_test_Twin2_filtered <<- subset(SNPs_and_FACS_Twins2, select = top_100[,1])
  for(i in 2:9){
    top_100 <- read.csv(file = files[,1][i], header = T, stringsAsFactors = F)
    SNPs_and_FACS_Twins1 <- readRDS(file = files[,2][i,])
    SNPs_and_FACS_Twins2 <- readRDS(file = files[,3][i,])
    SNPs_and_FACS_train_Twin1_filtered_i <- subset(SNPs_and_FACS_Twins1, select = top_100[,1])
    SNPs_and_FACS_train_Twin1_filtered <<- cbind(SNPs_and_FACS_train_Twin1_filtered, SNPs_and_FACS_train_Twin1_filtered_i)
    SNPs_and_FACS_test_Twin2_filtered_i <- subset(SNPs_and_FACS_Twins2, select = top_100[,1])
    SNPs_and_FACS_test_Twin2_filtered <<- cbind(SNPs_and_FACS_test_Twin2_filtered, SNPs_and_FACS_test_Twin2_filtered_i)
  }
}

SNPselector(files)

FACS_data <- rbind(SNPs_and_FACS_train_Twin1_filtered, SNPs_and_FACS_test_Twin2_filtered)

output <- ggplot(FACS_data, aes(x = geom_mean_FCFACS)) +
  geom_histogram(color = "darkblue", fill = "lightblue", binwidth = 1) +
  theme(
    # Centre title
    plot.title = element_text(hjust = 0.5),
    # Lengends to the top
    legend.position = "top",
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major.x = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size = .25, color = "black"),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank())

# output  # If you want a graph of the F-cell distributions
# ggsave("~/Dropbox/STP/MSc Bioinformatics/MSc Project/FACS_Fcell_histogram.pdf", output, width = 16*0.75, height = 9*0.75)

rm(list = c("SNPs_and_FACS_Twins1", "SNPs_and_FACS_Twins2", "top_100", "files"))
gc()  # Helps to manage PC resources

saveRDS(SNPs_and_FACS_train_Twin1_filtered, file = "~/Documents/twins_ML_project/plink/ML_project/SNPs_and_FACS_train_Twin1_filtered.rds")  # Save locally as compact rds file
saveRDS(SNPs_and_FACS_test_Twin2_filtered, file = "~/Documents/twins_ML_project/plink/ML_project/SNPs_and_FACS_test_Twin2_filtered.rds")  # Save locally as compact rds file

SNPs_and_FACS_train_Twin1_filtered.rds <- readRDS(file = "~/Documents/twins_ML_project/plink/ML_project/SNPs_and_FACS_train_Twin1_filtered.rds")
SNPs_and_FACS_test_Twin2_filtered.rds <- readRDS(file = "~/Documents/twins_ML_project/plink/ML_project/SNPs_and_FACS_test_Twin2_filtered.rds")

# Determine the best mtry value to use in the later analysis
m2 <- randomForest::tuneRF(
  x = SNPs_and_FACS_train_Twin1_filtered,
  y = SNPs_and_FACS_train_Twin1_filtered$geom_mean_FCFACS,
  ntreeTry = 500,
  mtryStart = 100,
  stepFactor = 1.25,
  improve = 0.01,
  trace = T  # show real-time progress
)

saveRDS(m2, file = "~/Documents/twins_ML_project/plink/ML_project/mtry_testing_900SNPs.rds")

m3 <- readRDS(file = "~/Documents/twins_ML_project/plink/ML_project/mtry_testing_900SNPs.rds")
m3 <- as.data.frame(m3)
plot(m3$mtry, m3$OOBError, type = "l")
m3[2,]

plot(m3)

# Split the glycomics data into training and test sets
# Twins1_SNPs_fcell_split <- rsample::initial_split(SNP_merged_test_Twin1, prop = 0.8)
# Twins1_SNPs_fcell_train <- rsample::training(Twins1_SNPs_fcell_split)
# Twins1_SNPs_fcell_test  <- rsample::testing(Twins1_SNPs_fcell_split)

Twins1_SNPs_fcell_train <- SNPs_and_FACS_train_Twin1_filtered
Twins1_SNPs_fcell_test <- SNPs_and_FACS_test_Twin2_filtered

rm(list = c("SNPs_and_FACS_train_Twin1_filtered", "SNPs_and_FACS_test_Twin2_filtered", "Twins1_SNPs_fcell_split", "m2"))
gc()

# Using h2o for grid optimisation and random forest generation, setup a cluster with 16Gb of RAM and all cores
h2o.init(max_mem_size = "16G", nthreads = -1)

# Use all cores
# registerDoParallel(parallel::detectCores())

# Turn training set into h2o object (slow)
# train.h2o <- as.h2o(Twins1_SNPs_fcell_train)
# test.h2o <- as.h2o(Twins1_SNPs_fcell_test)

# as.h2o very slo, use h2o.importFile instead (better to write to disk then load into h2o)
# data clipped? Maximum of ~265,000 columns?
data.table::fwrite(x = Twins1_SNPs_fcell_train, file = "/home/callumrakhit/Documents/twins_ML_project/plink/ML_project/Twins1_900_SNPs_fcell_train.csv")
data.table::fwrite(x = Twins1_SNPs_fcell_test, file = "/home/callumrakhit/Documents/twins_ML_project/plink/ML_project/Twins2_900_SNPs_fcell_test.csv")

train.h2o <- h2o.importFile(path = "/home/callumrakhit/Documents/twins_ML_project/plink/ML_project/Twins1_900_SNPs_fcell_train.csv")
test.h2o <- h2o.importFile(path = "/home/callumrakhit/Documents/twins_ML_project/plink/ML_project/Twins2_900_SNPs_fcell_test.csv")

rm(Twins1_SNPs_fcell_test)
rm(Twins1_SNPs_fcell_train)
gc()

# head(names(test.h2o))
# tail(names(test.h2o))  # Check if geom_mean_FCFACS is the last column

# Hyperparameter grid
hyper_grid.h2o <- list(
  ntrees = seq(101, 401, by = 50),
  mtries = seq(251, 501, by = 50),
  sample_rate = c(.55, .6, .65, .7)
)

# # Crappy hyperparameter grid based on the previous best values
# hyper_grid.h2o <- list(
#   ntrees = 251,
#   mtries = 251,
#   sample_rate = c(0.6)
# )

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
grid_900 <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_grid_900",
  x = x, 
  y = y,  
  training_frame = train.h2o,
  hyper_params = hyper_grid.h2o,
  search_criteria = search_criteria
)

# Collect the results and sort by our model performance metric of choice
best_grid <- h2o.getGrid(
  grid_id = "rf_grid_900", 
  sort_by = "mse",
  decreasing = F
)

model_1 <- h2o.randomForest(
  x = x, 
  y = y,
  ntrees = 251,
  mtries = 251,
  sample_rate = c(0.6),
  training_frame = train.h2o,
  # calibrate_model = T,
  # calibration_frame = test.h2o
)

saveRDS(object = model_1, file = "~/Twins/Models/model_1.rds")

# Save the model locally
h2o.saveModel(object = model_1, path = "~/Twins/Models/Model_1")
performance_fcell_twins <- h2o.performance(model_1)

test_88 <- as.data.frame(test.h2o)
test_88 <- test_88[1:88,]

data.table::fwrite(x = test_88, file = "~/Twins/Data/test_88.csv")
test_88.h2o <- h2o.importFile(path = "/home/callumrakhit/Twins/Data/test_88.csv")

predicted_fcell_twins <- h2o.predict(model_1, test.h2o)
predicted_fcell_twins_88 <- h2o.predict(model_1, test_88.h2o)

saveRDS(object = model_1, file = "~/Documents/twins_ML_project/plink/ML_project/model_1.rds")
model_1 <- h2o.loadModel
model_1
which(colnames(SNPs_and_FACS_test_Twin2_filtered.rds)=="geom_mean_FCFACS")

plot(as.data.frame(predicted_fcell_twins)[,1], as.data.frame(SNPs_and_FACS_test_Twin2_filtered.rds[,101])[,1])

plot(as.data.frame(predicted_fcell_twins_88)[,1], as.data.frame(SNPs_and_FACS_test_Twin2_filtered.rds[,101])[,1][1:88])

for_plot <- as.data.frame(cbind(as.data.frame(predicted_fcell_twins_88)[,1], as.data.frame(SNPs_and_FACS_test_Twin2_filtered.rds[,101])[,1][1:88]))

ggplot(data = for_plot, aes(x = V1, y = V2)) +
  geom_point() +
  geom_smooth(method = "loess") +
  theme_minimal() +
  xlab("Predicted") +
  ylab("Observed")



as.data.frame(predicted_fcell_twins_88)[,1]

as.data.frame(predicted_fcell_twins)[,1]
as.data.frame(SNPs_and_FACS_test_Twin2_filtered.rds[,101])[,1]
# Load model
# grid <- readRDS("~/Documents/twins_ML_project/plink/ML_project/rf_900.rds")

## Load in new disease data2
library(data.table)

Dataset_2 <- readRDS(file = "~/Twins/Data/Dataset_2.rds")
names(Dataset_2) <- gsub("JHU_1.", "rs", names(Dataset_2))
names(Dataset_2) <- gsub("exm", "rs", names(Dataset_2))
names(Dataset_2) <- gsub("1:", "rs", names(Dataset_2))
names(Dataset_2) <- gsub("-T", "", names(Dataset_2))
names(Dataset_2) <- gsub("-C", "", names(Dataset_2))
names(Dataset_2) <- gsub("-G", "", names(Dataset_2))
names(Dataset_2) <- gsub("-A", "", names(Dataset_2))
test <- subset(Dataset_2, select = colnames(SNPs_and_FACS_train_Twin1_filtered.rds))
hist(Dataset_2[,98269]$geom_mean_FCFACS, breaks = 100)

names(Dataset_2)

Dataset_2 %>% filter(names(Dataset_2) %in% names(SNPs_and_FACS_train_Twin1_filtered.rds))

names(SNPs_and_FACS_train_Twin1_filtered.rds) %in% names(Dataset_2)

names(SNPs_and_FACS_train_Twin1_filtered.rds)
names(Dataset_2)

class(names(SNPs_and_FACS_train_Twin1_filtered.rds))
class(names(Dataset_2))

colnames(Dataset_2)
Dataset_2$rs67765457_A
colnames(SNPs_and_FACS_train_Twin1_filtered.rds)
setDT(Dataset_2)[]










































# Print results
print(best_grid)

# Try to get some plots of the optimisation metrics
plot(best_grid@summary_table[["mtries"]], best_grid@summary_table[["mse"]])
plot(best_grid@summary_table[["ntrees"]], best_grid@summary_table[["mse"]])
plot(best_grid@summary_table[["sample_rate"]], best_grid@summary_table[["mse"]])

# Get the model_id for the top model, chosen by validation error
# best_model_id <- grid@model_ids[[1]]
best_model_id <- best_grid@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)

saveRDS(object = best_grid, file = "~/Documents/twins_ML_project/plink/ML_project/rf_900_bestmodel.rds")

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
  dplyr::top_n(n = 50, wt = relative_importance) %>%
  ggplot(aes(reorder(variable, relative_importance), relative_importance)) +
  geom_col() +
  coord_flip() +
  ggtitle("The 50 variables which most reduce the OOB RMSE") +
  theme_minimal()

# Pick out the top 10,000 rsIDs (ranked by relative importance)
as.data.frame(best_model@model[["variable_importances"]]) %>%
  dplyr::arrange(desc(relative_importance)) %>%
  dplyr::top_n(n = 100, wt = relative_importance) %>% 
  {. ->> rf_1_top_100_predictors }

best_model@model[["variable_importances"]]

rf_1_top_100_predictors[,1]

rf_1_top_100_predictors <- rf_1_top_100_predictors[,1]

View(rf_1_top_100_predictors)

rf_1_top_100_predictors == "rs2071348"

write_csv(x = as.data.frame(rf_1_top_100_predictors), path = "~/Documents/twins_ML_project/plink/ML_project/rf_900_top100_alt.csv")

h2o.shutdown(prompt = F)  # Shuts down the Java cluster

plot_multi_way_importance(best_model@model[["variable_importances"]][["variable"]], 
                          x_measure = best_model@model[["variable_importances"]][["scaled_importance"]], 
                          y_measure = best_model@model[["variable_importances"]][["relative_importance"]], 
                          size_measure = best_model@model[["variable_importances"]][["percentage"]], 
                          no_of_labels = 5)

# Get a list of the SNPs, relative importance and their locations
SNP_list <- as.data.frame(gsub('.{2}$', '', noquote(best_model@model[["variable_importances"]][["variable"]])))
colnames(SNP_list) <- c("rsID")
SNP_list$relative_importance <- gsub('.{2}$', '', noquote(best_model@model[["variable_importances"]][["relative_importance"]]))
SNP_list$location <- "1"
SNP_list$top_relevant_publication <- "1"

lapply(1:length(SNP_list$location), function(i){
  SNP_location <- annotations(snp = SNP_list$rsID[i], output = 'metadata')
  ifelse(dim(SNP_location)[1] == 0,
         SNP_list$location[i] <<- paste("Not Found"),
         SNP_list$location[i] <<- paste(SNP_location[2,][2], ":", SNP_location[3,][2], sep = "")
         )
  }
)

SNP_publication <- annotations(snp = SNP_list$rsID[850], output = 'all')
dim(SNP_publication)[1] == 0
SNP_publication$title[1]

lapply(1:length(SNP_list$location), function(i){
  SNP_publication <- annotations(snp = SNP_list$rsID[i], output = 'all')
  ifelse(dim(SNP_publication)[1] == 0,
         SNP_list$top_relevant_publication[i] <<- paste("Not Found"),
         SNP_list$top_relevant_publication[i] <<- SNP_publication$title[1]
         )
  }
)

write_csv(x = SNP_list, path = "~/Desktop/SNP_list.csv", col_names = T)

# This one may work better for getting locations
BiocManager::install()
library(BiocManager)
install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

snp_mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

snp_ids <- SNPs$`gsub(".{2}$", "", noquote(best_model@model[["variable_importances"]][["variable"]]))`

snp_attributes <- c("refsnp_id", "chr_name", "chrom_start")

snp_locations <- getBM(attributes = snp_attributes, filters = "snp_filter", values = SNP_list$rsID, mart = snp_mart)

manhatten <- as.data.frame(snp_locations)

SNPs <- as.data.frame(gsub('.{2}$', '', noquote(best_model@model[["variable_importances"]][["variable"]])))
SNPs$variable.importance <- best_model@model[["variable_importances"]][["percentage"]]
SNPs$`gsub(".{2}$", "", noquote(best_model@model[["variable_importances"]][["variable"]]))` <- as.character(SNPs$`gsub(".{2}$", "", noquote(best_model@model[["variable_importances"]][["variable"]]))`)
colnames(SNPs) <- c("refsnp_id", "variable.importance")

SNPs <- merge(manhatten, y = SNPs, by = "refsnp_id")
SNPs <- SNPs[!is.na(as.numeric(as.character(SNPs$chr_name))),]
SNPs$chr_name <- as.numeric(SNPs$chr_name)
SNPs$chrom_start <- as.numeric(SNPs$chrom_start)

options(scipen = 999)

View(reorder(x = SNPs$chrom_start, X = SNPs$chr_name))

annotations(snp = 'rs766432', output = 'all')

SNP_info <- annotations(snp = SNP_list$rsID[1], output = 'all')
SNP_info$title[1]

SNP_info <- annotations(snp = SNP_list$SNPs[1], output = 'metadata')

SNP_info <- rbind(SNP_info, annotations(snp = SNP_list$SNPs[2], output = 'all'))

manhattan.plot <- function(chr, pos, pvalue, 
                         sig.level = NA, 
                         annotate = NULL,
                         ann.default = list(),
                         should.thin = T, 
                         thin.pos.places = 2, 
                         thin.logp.places = 2, 
                         xlab = "Chromosome", 
                         # ylab = expression(-log[10](p-value)),
                         ylab = expression("Relative Importance in Random Forest Model"),
                         col = c("gray","darkgray"), 
                         # col = c("blue", "lightblue"),
                         panel.extra = NULL, pch = 20, cex = 0.8,...) {
  
  if (length(chr) == 0) stop("chromosome vector is empty")
  if (length(pos) == 0) stop("position vector is empty")
  if (length(pvalue) == 0) stop("pvalue vector is empty")
  
  # make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }
  
  # make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;
  
  # calculate absolute genomic position from relative chromosomal positions
  posmin <- tapply(pos, chr, min);
  posmax <- tapply(pos, chr, max);
  posshift <- head(c(0, cumsum(posmax)), - 1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos <- function(cchr, cpos) {
    p <- posshift[as.character(cchr)] + cpos
    return(p)
  }
  
  # parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default <- list(x = "peak", y = "peak", adj = NULL, pos = 3, offset = 0.5, 
                      col = NULL, fontface = NULL, fontsize = NULL, show = F)
  parse.label <- function(rawval, groupname) {
    r <- list(text = groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval) >= 1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times = length(pvalue)))
  }
  
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  
  if (length(ann.settings)>1) { 
    lcols <- trellis.par.get("superpose.symbol")$col 
    lfills <- trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch = pch, 
                              col = lcols[(i-2) %% length(lcols) + 1 ], 
                              fill = lfills[(i-2) %% length(lfills) + 1 ], 
                              cex = cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings) <- levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i > 1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate) > 1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1]) != "")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols - 1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  # reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      # logp = round(-log10(pvalue), thin.logp.places),
      logp = pvalue,
      pos = round(genpos, thin.pos.places), 
      chr = chr,
      grp = grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    # logp <- -log10(pvalue)
    logp <- pvalue
  }
  rm(pos, pvalue)
  gc()
  
  # custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side == "bottom") {
      panel.axis(side = side, outside = T,
                 at = ((posmax+posmin)/2+posshift),
                 labels = levels(chr), 
                 ticks = F, rot = 0,
                 check.overlap = F
      )
    } else if (side == "top" || side == "right") {
      panel.axis(side = side, draw.labels = F, ticks = F);
    }
    else {
      axis.default(side = side,...);
    }
  }
  
  # make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) { 
    A <- list();
    # maxy <- ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0))) + 0.5;
    maxy <- ceiling(max(y, ifelse(!is.na(sig.level), sig.level, 0))) - 0.85;
    A$ylim = c(0, maxy);
    A;
  }
  
  xyplot(logp~genpos, chr = chr, groups = grp,
         axis = axis.chr, ann.settings=ann.settings, 
         prepanel = prepanel.chr, scales=list(axs="i"),
         panel = function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             # add significance line (if requested)
             panel.abline(h = -log10(sig.level), lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           # allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt) == "text")] <- "labels"
             gt$show <- NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x == "peak") {gt$x<-x[peak]}
                 if(gt$x == "center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y == "peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x <- A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab = xlab, ylab = ylab, 
         panel.extra = panel.extra, getgenpos=getGenPos, ...
  );
  }

ann <- rep(1, length(SNPs$variable.importance))
ann[with(SNPs, chr_name == 2)] <- 2
ann[with(SNPs, chr_name == 6)] <- 3
ann[with(SNPs, chr_name == 11)] <- 4
ann <- factor(ann, levels = 1:4, labels=c("", "Chromsome 2", "Chromsome 6", "Chromsome 11"))

output <- manhattan.plot(factor(SNPs$chr_name, levels = c(1:22)), SNPs$chrom_start, SNPs$variable.importance, annotate = ann)

pdf("~/Desktop/Manhatten.pdf", width = 16*1.25, height = 9*1.25)
print(output)
dev.off()

View(ann)

