myPackages <- c("lme4", "rsq", "MuMIn", "optimx", "MASS", "dplyr", "ggplot2", 
                "HLMdiag", "lmtest", "lmerTest", "psy", "PowerUpR","nlme", "MVN", 
                "skimr", "tidyverse", "pacman", "lavaan", "semTools", "gsl", "graphics", "psych", "GPArotation")

# then check if all the needed packages are installed by checking against the list 
# of already installed packages in R:
installed <- sapply(myPackages, function(p) p %in% rownames(installed.packages()))

# and install any missing packages
if (any(!installed)) {
  install.packages(myPackages[!installed])
}

# then initiate the packages management library:
library(pacman)

# to load the rest of the needed libraries in this lab at once:
p_load(lme4, rsq, MuMIn, optimx, MASS, dplyr, ggplot2, HLMdiag, lmtest, lmerTest, 
       psy, PowerUpR, nlme, MVN, skimr, tidyverse, lavaan, semTools, gsl, graphics, psych, GPArotation)

#------------------------------------------------------------------------------
Cancer.Data <- Cancer.export[Cancer.export$YEAR %in% c("2006-2010", "2010-2014", 
                                                       "2012-2016", "2014-2018", 
                                                       "2016-2020"),]

Cancer.Data <- Cancer.Data %>%
  mutate(YEAR = case_when(
    YEAR == "2006-2010" ~ "2010",
    YEAR == "2010-2014" ~ "2014",
    YEAR == "2012-2016" ~ "2016",
    YEAR == "2014-2018" ~ "2018", 
    YEAR == "2016-2020" ~ "2020",
    # Add more conditions as needed
    TRUE ~ YEAR  # Default case to keep the original value
  ))
#====================The problem of having negative values======================
# Columns to process
columns <- c("COUNTS", "AADJ_RATE", "LOWR_LIM", "UPPR_LIM", "STD_ERR")

# Function to calculate mean excluding -5
calculate_mean_excluding_neg5 <- function(column_name) {
  values <- Cancer.export[[column_name]][Cancer.export$COUNTY == "Menominee"]
  filtered_values <- values[values != -5]
  mean(filtered_values, na.rm = TRUE)
}

# Apply the function to each column and collect results
means <- table(sapply(columns, calculate_mean_excluding_neg5))

#======================Solve the problem by replacing the means=================

# Columns to process
columns <- c("COUNTS", "AADJ_RATE", "LOWR_LIM", "UPPR_LIM", "STD_ERR")

# Iterate over each column and replace the values for "Menominee" county, excluding -5, with the mean
for (col in columns) {
  # Filter out -5 values and calculate the mean for "Menominee" county
  mean_value <- mean(Cancer.Data[Cancer.Data$COUNTY == "Menominee" & Cancer.Data[[col]] != -5, col], na.rm = TRUE)
  
  # Replace values in the column for "Menominee" county, excluding -5, with the mean value
  indices <- which(Cancer.Data$COUNTY == "Menominee" & Cancer.Data[[col]] == -5)
  Cancer.Data[indices, col] <- mean_value
}

#==========================Furhter Cleaning===================================

Cancer.Data <- subset(Cancer.Data, select = -c(OBJECTID, NAME))

#=============================================================================

Household.composition.disability.percentile.rank <- 
  Household.composition.disability.percentile.rank[, c("year", "COUNT", "OBJECTID", "FIPS", "COUNTY")]

Household.composition.disability.percentile.rank <- Household.composition.disability.percentile.rank %>% 
  rename(HCDPR = COUNT)

Household.composition.disability.percentile.rank$YEAR <- Household.composition.disability.percentile.rank$year
Household.composition.disability.percentile.rank$year <- NULL

Housing.transportation.percentile.rank <- 
  Housing.transportation.percentile.rank[, c("year", "COUNT", "OBJECTID", "FIPS", "COUNTY")]

Housing.transportation.percentile.rank <- Housing.transportation.percentile.rank %>% 
  rename(HTPR = COUNT)

Housing.transportation.percentile.rank$YEAR <- Housing.transportation.percentile.rank$year
Housing.transportation.percentile.rank$year <- NULL

Minority.status.language.percentile.rank <-
  Minority.status.language.percentile.rank[, c("year", "COUNT", "OBJECTID", "FIPS", "COUNTY")]

Minority.status.language.percentile.rank <- Minority.status.language.percentile.rank %>% 
  rename(MSLPR = COUNT)

Minority.status.language.percentile.rank$YEAR <- Minority.status.language.percentile.rank$year
Minority.status.language.percentile.rank$year<- NULL

Overall.percentile.vulnerability.rank <-
  Overall.percentile.vulnerability.rank[, c("year", "COUNT", "OBJECTID", "FIPS", "COUNTY")]

Overall.percentile.vulnerability.rank <- Overall.percentile.vulnerability.rank %>% 
  rename(OPVR = COUNT)

Overall.percentile.vulnerability.rank$YEAR <- Overall.percentile.vulnerability.rank$year
Overall.percentile.vulnerability.rank$year <- NULL

Socioeconomic.percentile.vulnerability.rank <-
  Socioeconomic.percentile.vulnerability.rank[, c("year", "COUNT", "OBJECTID", "FIPS", "COUNTY")]

Socioeconomic.percentile.vulnerability.rank <- Socioeconomic.percentile.vulnerability.rank %>% 
  rename(SPVR = COUNT)

Socioeconomic.percentile.vulnerability.rank$YEAR <- Socioeconomic.percentile.vulnerability.rank$year
Socioeconomic.percentile.vulnerability.rank$year <- NULL

#=================================Standardization===============================
# Function to standardize a variable
standardize_variable <- function(variable) {
  (variable - mean(variable, na.rm = TRUE)) / sd(variable, na.rm = TRUE)
}

# Standardize and add as new columns
Household.composition.disability.percentile.rank$HCDPR_standardized <- standardize_variable(Household.composition.disability.percentile.rank$HCDPR)
Housing.transportation.percentile.rank$HTPR_standardized <- standardize_variable(Housing.transportation.percentile.rank$HTPR)
Minority.status.language.percentile.rank$MSLPR_standardized <- standardize_variable(Minority.status.language.percentile.rank$MSLPR)
Overall.percentile.vulnerability.rank$OPVR_standardized <- standardize_variable(Overall.percentile.vulnerability.rank$OPVR)
Socioeconomic.percentile.vulnerability.rank$SPVR_standardized <- standardize_variable(Socioeconomic.percentile.vulnerability.rank$SPVR)

#===============================================================================
# Calculate quantiles for each standardized variable
Household.composition.disability.percentile.rank$HCDPR_standardized_quantiles <- cut(Household.composition.disability.percentile.rank$HCDPR_standardized, breaks = quantile(Household.composition.disability.percentile.rank$HCDPR_standardized, probs = 0:4/4), include.lowest = TRUE, labels = FALSE)
Housing.transportation.percentile.rank$HTPR_standardized_quantiles <- cut(Housing.transportation.percentile.rank$HTPR_standardized, breaks = quantile(Housing.transportation.percentile.rank$HTPR_standardized, probs = 0:4/4), include.lowest = TRUE, labels = FALSE)
Minority.status.language.percentile.rank$MSLPR_standardized_quantiles <- cut(Minority.status.language.percentile.rank$MSLPR_standardized, breaks = quantile(Minority.status.language.percentile.rank$MSLPR_standardized, probs = 0:4/4), include.lowest = TRUE, labels = FALSE)
Overall.percentile.vulnerability.rank$OPVR_standardized_quantiles <- cut(Overall.percentile.vulnerability.rank$OPVR_standardized, breaks = quantile(Overall.percentile.vulnerability.rank$OPVR_standardized, probs = 0:4/4), include.lowest = TRUE, labels = FALSE)
Socioeconomic.percentile.vulnerability.rank$SPVR_standardized_quantiles <- cut(Socioeconomic.percentile.vulnerability.rank$SPVR_standardized, breaks = quantile(Socioeconomic.percentile.vulnerability.rank$SPVR_standardized, probs = 0:4/4), include.lowest = TRUE, labels = FALSE)

#===========================Merge, Diagnose, & Clean Outlier Data==========================================
# Step 1: Ensure key columns are consistent
# Convert to the same case if necessary (e.g., both to uppercase)
#Household.composition.disability.percentile.rank$COUNTY <- toupper(Household.composition.disability.percentile.rank$COUNTY)
#Cancer.Data$COUNTY <- toupper(Cancer.Data$COUNTY)
# Ensure YEAR columns are of the same format and represent the same period
# For simplicity, let's assume both datasets already align on this

Household.composition.disability.percentile.rank <- as.data.frame(Household.composition.disability.percentile.rank)

Cancer.Data <- as.data.frame(Cancer.Data)
# Step 2: Merge the datasets
merged_data <- merge(Household.composition.disability.percentile.rank, Cancer.Data, by = c("FIPS", "COUNTY", "YEAR"), all = TRUE)
#------------------step in between-----------
skim(merged_data)
plot(density(merged_data$COUNTS))

# Calculate the 5th and 95th percentiles of the COUNTS variable
percentile_5 <- quantile(merged_data$COUNTS, probs = 0.05, na.rm = TRUE)
percentile_95 <- quantile(merged_data$COUNTS, probs = 0.95, na.rm = TRUE)

# Filter the dataset to only include values within these percentiles
Quantiled_data <- merged_data[merged_data$COUNTS > percentile_5 & merged_data$COUNTS < percentile_95, ]
skim(Quantiled_data)

#Transformation
Log_Counts <- log(Quantiled_data$COUNTS)
plot(density(Log_Counts))
#------------------------------------------
# Step 3: Analysis and Comparison
# For example, to compare the COUNTS between counties based on their HCDPR_standardized_quantiles
comparison <- merged_data %>%
  group_by(HCDPR_standardized_quantiles, COUNTY) %>%
  summarize(Mean_COUNTS = mean(COUNTS, na.rm = TRUE),
            Median_COUNTS = median(COUNTS, na.rm = TRUE),
            Quantile_COUNTS = quantile(COUNTS, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)) %>%
  arrange(HCDPR_standardized_quantiles, Mean_COUNTS)

# This code snippet assumes that 'COUNTS' is numeric and suitable for these operations.
# Adjust the 'summarize' function to fit your specific analytical goals.

#====================2010===========================================================
# Filter for the year 2010
data_2010 <- merged_data %>% 
  filter(YEAR == "2010" | YEAR == 2010) # Adjust based on whether YEAR is character or numeric

# Group by quantiles and summarize COUNTS
summary_counts <- data_2010 %>%
  group_by(HCDPR_standardized_quantiles) %>%
  summarize(
    Mean_COUNTS = mean(COUNTS, na.rm = TRUE),
    Median_COUNTS = median(COUNTS, na.rm = TRUE),
    Count = n()
  )

print(summary_counts)


summary_AADJ_RATE <- data_2010 %>%
  group_by(HCDPR_standardized_quantiles) %>%
  summarize(
    Mean_AADJ_RATE = mean(AADJ_RATE, na.rm = TRUE),
    Median_AADJ_RATE = median(AADJ_RATE, na.rm = TRUE),
    Count = n()
  )
print (summary_AADJ_RATE)

#===================================COUNTS============================================
ggplot(data_2010, aes(x = factor(HCDPR_standardized_quantiles), y = COUNTS)) +
  geom_boxplot() +
  labs(x = "Quantile", y = "COUNTS", title = "COUNTS by Quantile in 2010") +
  theme_minimal()

ggplot(summary_counts, aes(x = factor(HCDPR_standardized_quantiles), y = Mean_COUNTS, fill = factor(HCDPR_standardized_quantiles))) +
  geom_bar(stat = "identity") +
  labs(x = "Quantile", y = "Mean COUNTS", title = "Mean COUNTS by Quantile in 2010") +
  theme_minimal()

anova_test <- aov(COUNTS ~ factor(HCDPR_standardized_quantiles), data = data_2010)
summary(anova_test)

#===================================AADJ_RATE===================================
ggplot(data_2010, aes(x = factor(HCDPR_standardized_quantiles), y = AADJ_RATE)) +
  geom_boxplot() +
  labs(x = "Quantile", y = "AADJ_RATE", title = "AADJ_RATE by Quantile in 2010") +
  theme_minimal()

ggplot(summary_AADJ_RATE, aes(x = factor(HCDPR_standardized_quantiles), y = Mean_AADJ_RATE, fill = factor(HCDPR_standardized_quantiles))) +
  geom_bar(stat = "identity") +
  labs(x = "Quantile", y = "Mean AADJ_RATE", title = "Mean AADJ_RATE by Quantile in 2010") +
  theme_minimal()

anova_test <- aov(AADJ_RATE ~ factor(HCDPR_standardized_quantiles), data = data_2010)
summary(anova_test)

#===============================================================================

#====================2014===========================================================
# Filter for the year 2014
data_2014 <- merged_data %>% 
  filter(YEAR == "2014" | YEAR == 2014) # Adjust based on whether YEAR is character or numeric

# Group by quantiles and summarize COUNTS
summary_counts <- data_2014 %>%
  group_by(HCDPR_standardized_quantiles) %>%
  summarize(
    Mean_COUNTS = mean(COUNTS, na.rm = TRUE),
    Median_COUNTS = median(COUNTS, na.rm = TRUE),
    Count = n()
  )

print(summary_counts)


summary_AADJ_RATE <- data_2014 %>%
  group_by(HCDPR_standardized_quantiles) %>%
  summarize(
    Mean_AADJ_RATE = mean(AADJ_RATE, na.rm = TRUE),
    Median_AADJ_RATE = median(AADJ_RATE, na.rm = TRUE),
    Count = n()
  )
print (summary_AADJ_RATE)

#===================================COUNTS============================================
ggplot(data_2014, aes(x = factor(HCDPR_standardized_quantiles), y = COUNTS)) +
  geom_boxplot() +
  labs(x = "Quantile", y = "COUNTS", title = "COUNTS by Quantile in 2014") +
  theme_minimal()

ggplot(summary_counts, aes(x = factor(HCDPR_standardized_quantiles), y = Mean_COUNTS, fill = factor(HCDPR_standardized_quantiles))) +
  geom_bar(stat = "identity") +
  labs(x = "Quantile", y = "Mean COUNTS", title = "Mean COUNTS by Quantile in 2014") +
  theme_minimal()

anova_test <- aov(COUNTS ~ factor(HCDPR_standardized_quantiles), data = data_2014)
summary(anova_test)

#===================================AADJ_RATE============================================
ggplot(data_2014, aes(x = factor(HCDPR_standardized_quantiles), y = AADJ_RATE)) +
  geom_boxplot() +
  labs(x = "Quantile", y = "AADJ_RATE", title = "AADJ_RATE by Quantile in 2014") +
  theme_minimal()

ggplot(summary_AADJ_RATE, aes(x = factor(HCDPR_standardized_quantiles), y = Mean_AADJ_RATE, fill = factor(HCDPR_standardized_quantiles))) +
  geom_bar(stat = "identity") +
  labs(x = "Quantile", y = "Mean AADJ_RATE", title = "Mean AADJ_RATE by Quantile in 2014") +
  theme_minimal()

anova_test <- aov(AADJ_RATE ~ factor(HCDPR_standardized_quantiles), data = data_2014)
summary(anova_test)


#====================2016===========================================================
# Filter for the year 2016
data_2016 <- merged_data %>% 
  filter(YEAR == "2016" | YEAR == 2016) # Adjust based on whether YEAR is character or numeric

# Group by quantiles and summarize COUNTS
summary_counts <- data_2016 %>%
  group_by(HCDPR_standardized_quantiles) %>%
  summarize(
    Mean_COUNTS = mean(COUNTS, na.rm = TRUE),
    Median_COUNTS = median(COUNTS, na.rm = TRUE),
    Count = n()
  )

print(summary_counts)


summary_AADJ_RATE <- data_2016 %>%
  group_by(HCDPR_standardized_quantiles) %>%
  summarize(
    Mean_AADJ_RATE = mean(AADJ_RATE, na.rm = TRUE),
    Median_AADJ_RATE = median(AADJ_RATE, na.rm = TRUE),
    Count = n()
  )
print (summary_AADJ_RATE)

#===================================COUNTS============================================
ggplot(data_2016, aes(x = factor(HCDPR_standardized_quantiles), y = COUNTS)) +
  geom_boxplot() +
  labs(x = "Quantile", y = "COUNTS", title = "COUNTS by Quantile in 2016") +
  theme_minimal()

ggplot(summary_counts, aes(x = factor(HCDPR_standardized_quantiles), y = Mean_COUNTS, fill = factor(HCDPR_standardized_quantiles))) +
  geom_bar(stat = "identity") +
  labs(x = "Quantile", y = "Mean COUNTS", title = "Mean COUNTS by Quantile in 2016") +
  theme_minimal()

anova_test <- aov(COUNTS ~ factor(HCDPR_standardized_quantiles), data = data_2016)
summary(anova_test)

#===================================AADJ_RATE============================================
ggplot(data_2016, aes(x = factor(HCDPR_standardized_quantiles), y = AADJ_RATE)) +
  geom_boxplot() +
  labs(x = "Quantile", y = "AADJ_RATE", title = "AADJ_RATE by Quantile in 2016") +
  theme_minimal()

ggplot(summary_AADJ_RATE, aes(x = factor(HCDPR_standardized_quantiles), y = Mean_AADJ_RATE, fill = factor(HCDPR_standardized_quantiles))) +
  geom_bar(stat = "identity") +
  labs(x = "Quantile", y = "Mean AADJ_RATE", title = "Mean AADJ_RATE by Quantile in 2016") +
  theme_minimal()

anova_test <- aov(AADJ_RATE ~ factor(HCDPR_standardized_quantiles), data = data_2016)
summary(anova_test)

#====================2018===========================================================
# Filter for the year 2018
data_2018 <- merged_data %>% 
  filter(YEAR == "2018" | YEAR == 2018) # Adjust based on whether YEAR is character or numeric

# Group by quantiles and summarize COUNTS
summary_counts <- data_2018 %>%
  group_by(HCDPR_standardized_quantiles) %>%
  summarize(
    Mean_COUNTS = mean(COUNTS, na.rm = TRUE),
    Median_COUNTS = median(COUNTS, na.rm = TRUE),
    Count = n()
  )

print(summary_counts)


summary_AADJ_RATE <- data_2018 %>%
  group_by(HCDPR_standardized_quantiles) %>%
  summarize(
    Mean_AADJ_RATE = mean(AADJ_RATE, na.rm = TRUE),
    Median_AADJ_RATE = median(AADJ_RATE, na.rm = TRUE),
    Count = n()
  )
print (summary_AADJ_RATE)

#===================================COUNTS============================================
ggplot(data_2018, aes(x = factor(HCDPR_standardized_quantiles), y = COUNTS)) +
  geom_boxplot() +
  labs(x = "Quantile", y = "COUNTS", title = "COUNTS by Quantile in 2018") +
  theme_minimal()

ggplot(summary_counts, aes(x = factor(HCDPR_standardized_quantiles), y = Mean_COUNTS, fill = factor(HCDPR_standardized_quantiles))) +
  geom_bar(stat = "identity") +
  labs(x = "Quantile", y = "Mean COUNTS", title = "Mean COUNTS by Quantile in 2018") +
  theme_minimal()

anova_test <- aov(COUNTS ~ factor(HCDPR_standardized_quantiles), data = data_2018)
summary(anova_test)

#===================================AADJ_RATE============================================
ggplot(data_2018, aes(x = factor(HCDPR_standardized_quantiles), y = AADJ_RATE)) +
  geom_boxplot() +
  labs(x = "Quantile", y = "AADJ_RATE", title = "AADJ_RATE by Quantile in 2018") +
  theme_minimal()

ggplot(summary_AADJ_RATE, aes(x = factor(HCDPR_standardized_quantiles), y = Mean_AADJ_RATE, fill = factor(HCDPR_standardized_quantiles))) +
  geom_bar(stat = "identity") +
  labs(x = "Quantile", y = "Mean AADJ_RATE", title = "Mean AADJ_RATE by Quantile in 2018") +
  theme_minimal()

anova_test <- aov(AADJ_RATE ~ factor(HCDPR_standardized_quantiles), data = data_2018)
summary(anova_test)


#====================2020===========================================================
# Filter for the year 2020
data_2020 <- merged_data %>% 
  filter(YEAR == "2020" | YEAR == 2020) # Adjust based on whether YEAR is character or numeric

# Group by quantiles and summarize COUNTS
summary_counts <- data_2020 %>%
  group_by(HCDPR_standardized_quantiles) %>%
  summarize(
    Mean_COUNTS = mean(COUNTS, na.rm = TRUE),
    Median_COUNTS = median(COUNTS, na.rm = TRUE),
    Count = n()
  )

print(summary_counts)


summary_AADJ_RATE <- data_2020 %>%
  group_by(HCDPR_standardized_quantiles) %>%
  summarize(
    Mean_AADJ_RATE = mean(AADJ_RATE, na.rm = TRUE),
    Median_AADJ_RATE = median(AADJ_RATE, na.rm = TRUE),
    Count = n()
  )
print (summary_AADJ_RATE)

#===================================COUNTS============================================
ggplot(data_2020, aes(x = factor(HCDPR_standardized_quantiles), y = COUNTS)) +
  geom_boxplot() +
  labs(x = "Quantile", y = "COUNTS", title = "COUNTS by Quantile in 2020") +
  theme_minimal()

ggplot(summary_counts, aes(x = factor(HCDPR_standardized_quantiles), y = Mean_COUNTS, fill = factor(HCDPR_standardized_quantiles))) +
  geom_bar(stat = "identity") +
  labs(x = "Quantile", y = "Mean COUNTS", title = "Mean COUNTS by Quantile in 2020") +
  theme_minimal()

anova_test <- aov(COUNTS ~ factor(HCDPR_standardized_quantiles), data = data_2020)
summary(anova_test)

#===================================AADJ_RATE============================================
ggplot(data_2020, aes(x = factor(HCDPR_standardized_quantiles), y = AADJ_RATE)) +
  geom_boxplot() +
  labs(x = "Quantile", y = "AADJ_RATE", title = "AADJ_RATE by Quantile in 2020") +
  theme_minimal()

ggplot(summary_AADJ_RATE, aes(x = factor(HCDPR_standardized_quantiles), y = Mean_AADJ_RATE, fill = factor(HCDPR_standardized_quantiles))) +
  geom_bar(stat = "identity") +
  labs(x = "Quantile", y = "Mean AADJ_RATE", title = "Mean AADJ_RATE by Quantile in 2020") +
  theme_minimal()

anova_test <- aov(AADJ_RATE ~ factor(HCDPR_standardized_quantiles), data = data_2020)
summary(anova_test)


#====================ALL===========================================================
# Group by quantiles and summarize COUNTS
summary_counts <- merged_data %>%
  group_by(HCDPR_standardized_quantiles) %>%
  summarize(
    Mean_COUNTS = mean(COUNTS, na.rm = TRUE),
    Median_COUNTS = median(COUNTS, na.rm = TRUE),
    Count = n()
  )

print(summary_counts)


summary_AADJ_RATE <- merged_data %>%
  group_by(HCDPR_standardized_quantiles) %>%
  summarize(
    Mean_AADJ_RATE = mean(AADJ_RATE, na.rm = TRUE),
    Median_AADJ_RATE = median(AADJ_RATE, na.rm = TRUE),
    Count = n()
  )
print (summary_AADJ_RATE)

#===================================COUNTS============================================
ggplot(merged_data, aes(x = factor(HCDPR_standardized_quantiles), y = COUNTS)) +
  geom_boxplot() +
  labs(x = "Quantile", y = "COUNTS", title = "COUNTS by Quantile") +
  theme_minimal()

ggplot(summary_counts, aes(x = factor(HCDPR_standardized_quantiles), y = Mean_COUNTS, fill = factor(HCDPR_standardized_quantiles))) +
  geom_bar(stat = "identity") +
  labs(x = "Quantile", y = "Mean COUNTS", title = "Mean COUNTS by Quantile") +
  theme_minimal()

anova_test <- aov(COUNTS ~ factor(HCDPR_standardized_quantiles), data = merged_data)
summary(anova_test)

#===================================AADJ_RATE============================================
ggplot(merged_data, aes(x = factor(HCDPR_standardized_quantiles), y = AADJ_RATE)) +
  geom_boxplot() +
  labs(x = "Quantile", y = "AADJ_RATE", title = "AADJ_RATE by Quantile") +
  theme_minimal()

ggplot(summary_AADJ_RATE, aes(x = factor(HCDPR_standardized_quantiles), y = Mean_AADJ_RATE, fill = factor(HCDPR_standardized_quantiles))) +
  geom_bar(stat = "identity") +
  labs(x = "Quantile", y = "Mean AADJ_RATE", title = "Mean AADJ_RATE by Quantile") +
  theme_minimal()

anova_test <- aov(AADJ_RATE ~ factor(HCDPR_standardized_quantiles), data = merged_data)
summary(anova_test)



