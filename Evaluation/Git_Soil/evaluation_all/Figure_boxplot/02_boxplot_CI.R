library(tidyverse)
library(ggplot2)

# Load your CI results
df <- read_csv("CI.csv")

# Add midpoint for plotting
df$CI_Mid <- (df$CI_Lower + df$CI_Upper) / 2

# === Define method order ===
method_order <- c(
  paste0("all_", 0:9),
  paste0("conditional_", 0:4),
  paste0("cor_", 0:4),
  "pearson","spearman","bicor","sparcc",
  "spiecEasi_mb", "spiecEasi_glasso", "cclasso", "CMIMN", "gcoda", "SPRING"
)

df$Method <- factor(df$Method, levels = method_order)

# === Vertical line positions between groups ===
vline_positions <- c(
  length(paste0("CMiNet_All_", 0:9)),
  length(paste0("CMiNet_All_", 0:9)) + length(paste0("CMiNet_Conditional_", 0:4)),
  length(paste0("CMiNet_All_", 0:9)) + length(paste0("CMiNet_Conditional_", 0:4)) + length(paste0("CMiNet_Correlation_", 0:4))
)

# === Plot ===
p <- ggplot(df, aes(x = Method, y = CI_Mid)) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.3, color = "black") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 10),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    y = "Bootstrap Confidence Interval (Soil microbiome)"
  ) +
  geom_vline(xintercept = vline_positions + 0.5, linetype = "dashed", color = "gray30")

# === Section labels under the x-axis ===
section_labels <- tibble(
  x = c(5, 13, 18.5, 25),
  y = rep(0.98, 4),
  label = c("CMiNet_All", "CMiNet_Conditional", "CMiNet_Correlation", "Single Method")
)

p_final <- p + 
  geom_text(data = section_labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, angle = 0, size = 3, fontface = "bold", vjust = -1) +
  scale_y_continuous(limits = c(0.4, 1), expand = expansion(mult = c(0, 0.05)))
# === Save the plot ===
ggsave("Fscore_CI_plot.png", p_final, width = 14, height = 5, dpi = 300)

# Print it
print(p_final)
