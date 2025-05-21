# ------------------------------------------------------------------------------
# Script: phi_comparison_table_script.R
# Purpose: Compare phi coefficients between two datasets (e.g., FinnGen and UKB / male and female)
#          Setting non-significant values to NA.
#
# Instructions:
#   - Set the significance level and output file name below.
#   - Edit input file names for the two datasets as needed.
#   - This script expects input files with columns: Category1, Category2, Phi, P_value.
#
# Output:
#   - A tab-separated file with columns: Category1, Category2, Phi_data1, Phi_data2
# ------------------------------------------------------------------------------


library(dplyr)
library(tidyr)

# Edit significance level here
significance_level = 0.05

#Edit output file name here
output_filename="phi_comparison_male.txt"

# Edit input names here
data1 <- read.table("male_phi_coefficient_finngen.txt", header = TRUE, sep = "\t")
data2 <- read.table("male_phi_coefficient_ukb.txt", header = TRUE, sep = "\t")

# Join key
data1 <- data1 %>% mutate(CategoryKey = paste(Category1, Category2, sep = "||"))
data2 <- data2 %>%mutate(CategoryKey = paste(Category1, Category2, sep = "||"))

# Select columns
data1_sub <- data1 %>%select(CategoryKey, Phi_data1 = Phi, P_data1 = P_value)

data2_sub <- data2 %>%select(CategoryKey, Phi_data2 = Phi, P_data2 = P_value)

# Full join (all category pairs present in either dataset)
comparison <- full_join(data1_sub, data2_sub, by = "CategoryKey")

# Set Phi to NA where p >= significance level
comparison <- comparison %>%mutate(Phi_data1 = ifelse(P_data1 < significance_level, Phi_data1, NA), Phi_data2 = ifelse(P_data2 < significance_level, Phi_data2, NA))

# Split CategoryKey back to Category1 and Category2
comparison <- comparison %>%separate(CategoryKey, into = c("Category1", "Category2"), sep = "\\|\\|")

# Reorder columns for clarity
comparison <- comparison %>%select(Category1, Category2, Phi_data1, Phi_data2)

write.table(comparison, output_filename, sep="\t", row.names=FALSE, quote=FALSE)
