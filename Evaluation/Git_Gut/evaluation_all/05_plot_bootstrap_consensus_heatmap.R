library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(grid)   # for unit()

# read the grid summary which resulted by running 04_gridsearch_m0_mstar_cstar.R
gs <- read_csv("grid_summary_m0_mstar_cstar.csv", show_col_types = FALSE)
# % retained and two-line label (first line = joint, second = percent)
gs <- gs %>%
  mutate(
    pct_retained = if_else(method_only > 0, 100 * joint / method_only, NA_real_),
    label = sprintf("%d\n%.1f", joint, pct_retained)
  )

# order axes
gs$m_star     <- factor(gs$m_star,  levels = sort(unique(gs$m_star)))
c_levels      <- sort(unique(gs$c_star))
gs$c_star_lab <- factor(sprintf("%.2f", gs$c_star), levels = sprintf("%.2f", c_levels))
gs$m0_strict  <- factor(gs$m0_strict, levels = sort(unique(gs$m0_strict)))

# dynamic label color for contrast
gs$txt_col <- ifelse(gs$pct_retained >= 75, "white", "black")

# color scale: white at the MIN of your data, dark at 100
lo <- floor(min(gs$pct_retained, na.rm = TRUE) / 5) * 5
hi <- 100
pal <- c("#FFFFFF", "#cfe6f5", "#9ecae1", "#6baed6", "#3182bd", "#08519c")

p <- ggplot(gs, aes(x = c_star_lab, y = m_star, fill = pct_retained)) +
  geom_tile(color = "white", size = 0.4) +
  geom_text(aes(label = label, color = txt_col), size = 3, lineheight = 0.9, fontface = "bold") +
  scale_color_identity() +
  facet_wrap(
    ~ m0_strict,
    nrow = 1,
    labeller = label_bquote(italic(m)[0] == .(as.character(m0_strict)))
  ) +
  scale_fill_gradientn(
    name   = "% retained",
    colors = pal,
    limits = c(lo, hi),   # white starts at the minimum observed
    oob    = squish
  ) +
  labs(
    x = expression(paste(italic(c),"*")),
    y = expression(paste(italic(m),"*")),
    title = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid     = element_blank(),
    strip.text     = element_text(face = "italic"),
    axis.text.x    = element_text(angle = 45, hjust = 1),
    plot.title     = element_text(face = "bold"),
    legend.position = "right",          # move to the right
    legend.direction = "vertical",      # vertical colorbar
    legend.box      = "vertical",
    legend.title    = element_text(vjust = 0.8)
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",           # title above the bar
      barwidth = unit(0.1, "in"),       # narrow width for vertical bar
      barheight = unit(3, "in"),        # tall for vertical layout
      ticks.colour = "grey30"
    )
  )
ggsave("bootstrap_consensus_heatmap_Gut.pdf", p, width = 12, height = 4.8)
