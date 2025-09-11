library(readr)

A <- data.frame()

# Load all CSV files
files1 <- list.files(pattern = "^.*\\.csv$")

for (i in seq_along(files1)) {
  df <- read_csv(files1[i], show_col_types = FALSE)
  
  # Assuming F_score is a column in each file
  B_f_score <- df$F_score
  
  # Compute confidence interval (percentile method)
  sig <- 0.05
  q_sig1 <- quantile(B_f_score, sig / 2, na.rm = TRUE)
  q_sig2 <- quantile(B_f_score, 1 - sig / 2, na.rm = TRUE)
  
  # Clean filename to use as method name
  NAME <- gsub("^Result_soil_|^One_|\\.csv$", "", files1[i])
  
  # Combine into a row and bind to A
  CI_percentile <- data.frame(Method = NAME, CI_Lower = q_sig1, CI_Upper = q_sig2)
  A <- rbind(A, CI_percentile)
}

# View result
print(A)
write_csv(A,"CI.csv")




