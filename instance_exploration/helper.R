## loading some libraries
library(ggplot2)
library(egg)
source("rremove.R")

## setting default color palette and theme
#cbPalette <- c("#D55E00", "#E69F00", "#F0E442", "#0072B2", "#56B4E9",
#               "#009E73", "#999999", "#CC79A7", "#ff0000", "#00ff00", "#0000ff")
#ggplot <- function(...) ggplot2::ggplot(...) +
#  scale_color_manual(values = cbPalette) +
#  scale_fill_manual(values = cbPalette)

#theme_set(theme_bw())

## log-log plots
log_log_plot <- function() list(scale_x_log10(), scale_y_log10(), annotation_logticks())

