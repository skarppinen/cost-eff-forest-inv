include("fig-helpers.jl");
using RCall
R"""
suppressPackageStartupMessages({
    library(ggplot2)
    library(latex2exp)
})
FIGURE_WIDTH_CM <- $FIGURE_WIDTH_CM
PLOT_OUTPUT_PATH <- $PLOT_OUTPUT_PATH
TEXTSIZE <- unit(10, "pt")
THEME_SIZES <- theme(axis.text = element_text(size = 10),
                     axis.title = element_text(size = 10),
                     strip.text = element_text(size = 10),
                     legend.title = element_text(size = 10))
THEME_TIGHT <- theme(panel.spacing = unit(1, "pt"),
                     plot.margin = unit(rep(2, 4), "pt"))
THEME_TIGHT_TOP_LEGEND <- theme(legend.position = "top", legend.box.margin = unit(c(0, 0, -10, 0), "pt"))
"""
