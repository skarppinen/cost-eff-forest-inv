using JLD2, DataFrames
using Statistics
include("r-fig-helpers.jl");

datapath = joinpath(@__DIR__, "../../../..", "no-vc", "output",
                    "run-random-sweep-new-impl");
filepaths = readdir(datapath, join = true);
procfun = filepath -> begin
    jldopen(filepath, "r") do file
        args = file["args"];
        objvalues = file["objvalues"];
        DataFrame(xI = args["xI"],
                  ninits = args["ninits"],
                  strategy = args["strategy"],
                  index = 1:length(objvalues),
                  obj = objvalues);
    end
end
d = mapreduce(procfun, vcat, filepaths);
minimums_df = combine(
    groupby(d, [:xI]), :obj => minimum => :min_obj
);
d = leftjoin(d, minimums_df, on = :xI);
d_quantile = combine(
    groupby(d, [:xI, :ninits, :strategy]), 
    :obj => (x -> quantile(x, 0.05)) => :q05,
    :obj => (x -> quantile(x, 0.95)) => :q95,
    :obj => (x -> quantile(x, 0.01)) => :q01,
    :obj => (x -> quantile(x, 0.99)) => :q99
);

global_opt_filepath = joinpath(datapath, "..", "run-global-optim", "jobid-1.jld2")
global_opt_data = jldopen(global_opt_filepath, "r") do file
    (objvalue = file["objvalue"], 
     objbound = file["objbound"],
     xI = file["args"]["xI"])
end
d_rel_error = filter(r -> r[:xI] == global_opt_data.xI, d);
d_rel_error[!, :min_obj] .= global_opt_data.objvalue;
d_rel_error[!, :rel_error] = (d_rel_error[!, :obj] - d_rel_error[!, :min_obj]) ./ d_rel_error[!, :min_obj] 

d_rel_err_stat = combine(groupby(d_rel_error, [:ninits, :strategy]),
                         :rel_error => minimum => :min_rel_error,
                         :rel_error => maximum => :max_rel_error,
                         :rel_error => median => :md_rel_error)

R"""
d <- $d_quantile
strategy_names <- c(smallest_negative = "Smallest negative", 
    greatest_negative = "Greatest negative")
d[["strategy_fact"]] <- factor(strategy_names[d[["strategy"]]], levels = unname(strategy_names)) 
inventory_names <- c(`1` = "5 plots", `2` = "10 plots", `3` = "20 plots")
d[["xI_fact"]] = factor(inventory_names[d[["xI"]]], levels = unname(inventory_names))

plt <- ggplot(d, aes(x = factor(ninits),  
       group = strategy_fact, linetype = strategy_fact)) +
       #geom_ribbon(aes(ymin = q01, ymax = q99), color = "black", fill = NA, alpha = 0.3) +
       geom_line(aes(y = q01)) +
       geom_line(aes(y = q99)) +
       #geom_line(aes(y = best_obj)) +
  facet_grid(~ xI_fact) +
  scale_y_continuous(expand = expansion(add = 0.0, mult = 0.01)) +
  scale_x_discrete(expand = expansion(add = 0.0, mult = 0.025)) +
  labs(x = "Number of random initialisations", 
       y = "1% and 99% quantiles of final objective function values",
       linetype = "Descent strategy") +
  theme_bw() +
  THEME_SIZES + 
  THEME_TIGHT + 
  theme(legend.position = "top",
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.box.margin = unit(c(0, 0, -10, 0), "pt"))

fp <- file.path(PLOT_OUTPUT_PATH, "fig-random-sweep.pdf")
ggsave(plot = plt, filename = fp,
       device = "pdf", width = FIGURE_WIDTH_CM, height = FIGURE_WIDTH_CM * 0.66, units = "cm")
cat(paste0("Wrote file ", fp, "\n"))
"""

R"""
d <- $d_rel_error
d <- d[d$ninits > 25, , drop = FALSE]
d <- d[d$strategy == "greatest_negative", , drop = FALSE]
strategy_names <- c(smallest_negative = "Smallest negative", 
    greatest_negative = "Greatest negative")
d[["strategy_fact"]] <- factor(strategy_names[d[["strategy"]]], levels = unname(strategy_names)) 

plt_rel_error <- ggplot(d, aes(x = factor(ninits))) +
geom_violin(aes(y = rel_error)) +
scale_y_continuous(breaks = seq(0.0, 0.5, by = 0.05)) +
labs(x = "Number of random initialisations", y = "Relative error wrt. best found minimum") +
theme_bw() +
THEME_SIZES +
THEME_TIGHT

fp <- file.path(PLOT_OUTPUT_PATH, "fig-rel-error.pdf")
ggsave(plot = plt_rel_error, filename = fp,
       device = "pdf", width = FIGURE_WIDTH_CM, height = 0.5 * FIGURE_WIDTH_CM, units = "cm")
cat(paste0("Wrote file ", fp, "\n"))
"""
