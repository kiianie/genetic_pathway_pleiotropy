library(psych)
library(corrplot)

# Pairwise phi ----------------------------------------------------------------
data_table <- read.delim('female_common_pathways_finngen.txt')

# Extract category (binary values) columns
category_data <- data_table[, c("cancer", "chf", "copd", "diabetes", "dementia", "mi", "stroke")]

# Initialize data frame
results <- data.frame()

# Loop through all combinations
category_names <- colnames(category_data)
for (i in 1:(ncol(category_data)-1)) {
  for (j in (i+1):ncol(category_data)) {
    var1 <- category_names[i]
    var2 <- category_names[j]
    
    tbl <- table(category_data[[var1]], category_data[[var2]])

    if (length(tbl) < 4) {
      full_tbl <- matrix(0, nrow = 2, ncol = 2)
      rownames(full_tbl) <- colnames(full_tbl) <- c("0", "1")
      full_tbl[rownames(tbl), colnames(tbl)] <- tbl
      tbl <- full_tbl
    }
    
    phi <- phi(tbl)
    # Incase there are zeros in the matrix
    chi <- tryCatch(chisq.test(tbl), error = function(e) list(p.value = NA))

    # Extract counts
    count_00 <- tbl[1, 1]
    count_01 <- tbl[1, 2]
    count_10 <- tbl[2, 1]
    count_11 <- tbl[2, 2]
    
    results <- rbind(results, data.frame(
      Category1 = var1,
      Category2 = var2,
      Phi = round(phi, 3),
      P_value = round(chi$p.value, 100),
      Count_11 = count_11,
      Count_10 = count_10,
      Count_01 = count_01,
      Count_00 = count_00
    ))
  }
}

write.table(results, "male_phi_coefficient.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Corrplot --------------------------------------------------------------

data_table <- read.delim('female_common_pathways_finngen.txt')

disease_data <- data_table[, c("cancer", "chf", "copd", "diabetes", "dementia", "mi", "stroke")] #only binary values

# Filters out any columns with no variation
disease_data <- disease_data[, sapply(disease_data, function(col) length(unique(col)) > 1)]

phi_matrix <- cor(disease_data, method = "pearson")

corrplot(phi_matrix, method = "color", addCoef.col = "black", type = "lower", tl.cex = 1)

write.table(phi_matrix, "male_phi_matrix_finngen.txt", sep = "\t", row.names = FALSE, quote = FALSE)