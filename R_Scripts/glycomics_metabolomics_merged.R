##### Assessing the overlap in the datasets #####

count(
  merge(
    as.data.frame(sort(unique(twins_data_glycomics$PublicID))),
    as.data.frame(sort(unique(twins_data_metabolomics$PublicID))),
    by.x = 'sort(unique(twins_data_glycomics$PublicID))', 
    by.y = 'sort(unique(twins_data_metabolomics$PublicID))'
  )
)

##### Attempting to merge glycomics and metabolomics #####

# Merge information from 2 into 1
twins_data_bmi_alldemographics <- merge(twins_data_bmi, twins_data_demographics, by = "PublicID")
twins_data_bmi_alldemographics <- merge(twins_data_bmi_alldemographics, twins_dob, by = "PublicID")

# Split visit date into month and year of visit for demographics
twins_data_bmi_alldemographics$Visit_Date <- as.Date(twins_data_bmi_alldemographics$Visit_Date, "%d/%m/%Y")
twins_data_bmi_alldemographics$Year.Reported <- format(as.Date(twins_data_bmi_alldemographics[["Visit_Date"]]), "%Y")
twins_data_bmi_alldemographics$Month.Reported <- format(as.Date(twins_data_bmi_alldemographics[["Visit_Date"]]), "%m")

# Drop unwanted columns
drops <- c("STUDY_TYPE", "month_year_birth", "Visit_Date")
twins_data_bmi_alldemographics <- twins_data_bmi_alldemographics[ , !(names(twins_data_bmi_alldemographics) %in% drops)]

# Split visit date into month and year of visit for glycomics
twins_data_glycomics$date <- as.Date(twins_data_glycomics$date, "%d/%m/%Y")
twins_data_glycomics$Year.Reported <- format(as.Date(twins_data_glycomics[["date"]]), "%Y")
twins_data_glycomics$Month.Reported <- format(as.Date(twins_data_glycomics[["date"]]), "%m")

# Split visit date into month and year of visit for metabolobics
twins_data_metabolomics$date_hli <- as.Date(twins_data_metabolomics$date_hli, "%d/%m/%Y")
twins_data_metabolomics$Year.Reported <- format(as.Date(twins_data_metabolomics[["date_hli"]]), "%Y")
twins_data_metabolomics$Month.Reported <- format(as.Date(twins_data_metabolomics[["date_hli"]]), "%m")

# Create a unique PublicID/year column, can use to merge datasets with
twins_data_bmi_alldemographics$ID_visityear <- paste(twins_data_bmi_alldemographics$PublicID, twins_data_bmi_alldemographics$Year.Reported, sep = "_")
twins_data_glycomics$ID_visityear <- paste(twins_data_glycomics$PublicID, twins_data_glycomics$Year.Reported, sep = "_")
twins_data_metabolomics$ID_visityear <- paste(twins_data_metabolomics$PublicID, twins_data_metabolomics$Year.Reported, sep = "_")

# NAs not tollerated in random forest methodology, need to identify and remove them
sapply(twins_data_met_dob, function(x) sum(is.na(x)))
which(is.na(twins_data_met_dob$m52689))
twins_data_met_dob <- twins_data_met_dob[-c(which(is.na(twins_data_met_dob$m52689))),]
twins_data_met_dob <- twins_data_met_dob[-c(which(is.na(twins_data_met_dob$m6146))),]

# Convert to numeric (via character to avoid coverting factors to numeric storage values) 
twins_data_gly_dob <- mutate_all(twins_data_gly_dob, function(x) as.numeric(as.character(x)))
twins_data_met_dob[is.na(twins_data_met_dob)] <- 0
twins_data_met_dob <- mutate_all(twins_data_met_dob, function(x) as.numeric(as.character(x)))

# Generate the random forest
rf_gly_1 <- randomForest(age_at_test ~ ., data = twins_data_gly_dob, importance = T)  
rf_met_1 <- randomForest(age_at_test ~ ., data = twins_data_met_dob, importance = T)
tune <- tuneRF(twins_data_gly_dob, y = as.factor(twins_data_gly_dob$age_at_test), doBest = T)

rf_gly_2 <- randomForest(as.factor(age_at_test) ~ ., data = twins_data_gly_dob, importance = T)
rf_gly_smol <- randomForest(as.factor(age_at_test) ~ ., data = twins_data_gly_dob, importance = T, nodesize = 1500)

getTree(rf, k = 1, labelVar = T)
getTree(rf_gly_2, k = 1, labelVar = T)
getTree(rf_gly_smol, k = 1, labelVar = T)

ID1 <- twins_data_glycomics$PublicID
ID2 <- twins_data_metabolomics$PublicID
overlap <- as.data.frame(intersect(ID1, ID2))

write.csv(overlap, "~/Desktop/overlapping.csv")