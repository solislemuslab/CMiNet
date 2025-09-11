library(tidyverse)
library(ggplot2)
library(viridis)
library(stringr)
library(patchwork)

# === Load CMiNet multi-score data ===
files1 <- list.files(pattern = "^Result_soil_.*\\.csv$")

df_list1 <- lapply(files1, function(file) {
  df <- read_csv(file, show_col_types = FALSE)
  parts <- str_match(file, "Result_soil_(.*?)_(\\d+)\\.csv")
  df$Group <- parts[2]
  df$Score <- as.integer(parts[3])
  return(df)
})

df1 <- bind_rows(df_list1)

# Beautify group names
df1$Group <- recode(df1$Group,
                    "cor" = "CMiNet_Correlation",
                    "conditional" = "CMiNet_Conditional",
                    "all" = "CMiNet_All")

# Create Method_Label to unify same-score color across CMiNet groups
df1 <- df1 %>%
  mutate(Method_Label = paste0(Group, "_", Score))

# === Load single-score method results ===
files2 <- list.files(pattern = "^One_.*\\.csv$")

df_list2 <- lapply(files2, function(file) {
  df <- read_csv(file, show_col_types = FALSE)
  method <- str_remove(str_remove(file, "One_"), ".csv")
  df$Group <- method
  df$Score <- NA
  df$Method_Label <- method
  return(df)
})

df2 <- bind_rows(df_list2)

# === Combine data ===
df_combined <- bind_rows(df1, df2)

# Define method order
method_order <- c(
  paste0("CMiNet_All_", 0:9),
  paste0("CMiNet_Conditional_", 0:4),
  paste0("CMiNet_Correlation_", 0:4),
  "pearson","spearman","bicor","sparcc","spiecEasi_mb", "spiecEasi_glasso", "cclasso", "CMIMN", "gcoda", "SPRING"
)
df_combined$Method_Label <- factor(df_combined$Method_Label, levels = method_order)

# === Define consistent color palette ===
threshold_colors <- viridis(10, option = "plasma")  # 10 thresholds max
names(threshold_colors) <- as.character(0:9)

# Assign color by threshold for CMiNet scores
df_combined <- df_combined %>%
  mutate(Color = ifelse(!is.na(Score), as.character(Score), Method_Label))

# Define full color map
method_colors <- c(
  threshold_colors,  # for Score 0–9
  "spiecEasi_mb" = "#9491D9",
  "spiecEasi_glasso" = "#3F8C61",
  "cclasso" = "#F24405",
  "CMIMN" = "#8C2E62",
  "gcoda" = "#F2B705",
  "SPRING" = "#11A0D9"
)

# === Add vertical lines to separate method groups ===
vline_positions <- c(
  length(paste0("CMiNet_All_", 0:9)),
  length(paste0("CMiNet_All_", 0:9)) + length(paste0("CMiNet_Conditional_", 0:4)),
  length(paste0("CMiNet_All_", 0:9)) + length(paste0("CMiNet_Conditional_", 0:4)) + length(paste0("CMiNet_Correlation_", 0:4))
)

# === Plot ===
p <- ggplot(df_combined, aes(x = Method_Label, y = F_score, fill = Color)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.9) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5, color = "black") +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0.4, 1), expand = expansion(mult = c(0, 0.05)))+
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13, face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  ) +
  labs(
    y = "F-score (Soil microbiome)"
  ) +
  geom_vline(xintercept = vline_positions + 0.5, linetype = "dashed", color = "black", linewidth = 0.9)

# === Add section labels under the plot ===
section_labels <- tibble(
  x = c(5, 13, 18.5, 25),
  y = rep(1, 4),
  label = c("CMiNet_All", "CMiNet_Conditional", "CMiNet_Correlation", "Single Method")
)

p_final <- p + 
  geom_text(data = section_labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, angle = 0, size = 4.5, fontface = "bold")

# === Save the figure ===
ggsave("fscore_soil_combined_final.png", p_final, width = 14, height = 5, dpi = 300)



library(dplyr)

# Compute FP and FN
df_summary <- df_combined %>%
  mutate(
    FP = number_edges_phy  - TP,
    FN = number_edges_cminet - TP
  )

# Summarize by Method_Label
summary_table <- df_summary %>%
  group_by(Method_Label) %>%
  summarise(
    `Ref. Edges` = number_edges_cminet[1],
    `Mean Bootstrap Edges ± SD` = paste0(round(mean(number_edges_phy), 1), " ± ", round(sd(number_edges_phy), 1)),
    `Mean TP` = round(mean(TP), 1),
    `Mean FP` = round(mean(FP), 1),
    `Mean Precision` = round(mean(TP / (TP + FP), na.rm = TRUE), 2),
    `Mean Recall` = round(mean(TP / (TP + FN), na.rm = TRUE), 2),
    `Mean F1` = round(mean(F_score, na.rm = TRUE), 2),
    .groups = "drop"
  )

# View the result
print(summary_table)

# Optional: save to CSV
write.csv(summary_table, "network_performance_soil.csv", row.names = FALSE)
