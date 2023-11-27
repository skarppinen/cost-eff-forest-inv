using JLD2, DataFrames, RCall
project_root = joinpath(@__DIR__, "../../../..");
include(joinpath(project_root, "src/julia/lib/model.jl"));
include("r-fig-helpers.jl");


function get_budget_ranges_df(inputfolder) 
    inputfiles = readdir(inputfolder, join = true);
    df = get_planning_problem_evaluations(inputfiles);

    # Find optimal decisions for candidate budgets.
    candidate_budgets = minimum(df[!, :costs]):1000:maximum(df[!, :costs]);
    optimal_decisions = map(candidate_budgets) do budget
	sub = filter(r -> r[:costs] <= budget, df);
	maxpov, maxind = findmax(sub[!, :povs]);
	maxpov, sub[maxind, :costs], sub[maxind, :xI]; 
    end
    stda = jldopen(joinpath(project_root, "data", "planning-problem-data.jld2"), "r") do file
	file["stda"];
    end
    # Construct DataFrame with columns: stand_id, inventory method, budget_lower, budget_upper
    inventory_mat = transpose(hcat(getindex.(optimal_decisions, 3)...));
    inventory_decisions = sort(unique(inventory_mat));
    budget_ranges_df = 
	DataFrame(stand_id = repeat(1:size(inventory_mat, 2), inner = length(candidate_budgets)),
	    budget = repeat(candidate_budgets, size(inventory_mat, 2)),
	    inventory_method = vec(inventory_mat)); 

    # Order stands in order specified by `sort_variable`.
    sort_variable = mapslices(sum, stda.sigma2prior, dims = 2)[:, 1];
    order_df = DataFrame(stand_id = sortperm(sort_variable), order = 1:length(sort_variable));
    budget_ranges_df = leftjoin(budget_ranges_df, order_df, on = [:stand_id]);   
    budget_ranges_df;
end

inputfolder_normal = joinpath(project_root, "no-vc", "output", "run-planning-problem"); 
inputfolder_lognormal = joinpath(project_root, "no-vc", "output", "run-planning-problem-lognormal", "output"); 
budget_ranges_df_normal = get_budget_ranges_df(inputfolder_normal);
budget_ranges_df_lognormal = get_budget_ranges_df(inputfolder_lognormal);
budget_ranges_df_normal[!, :model] .= "normal";
budget_ranges_df_lognormal[!, :model] .= "lognormal";
budget_ranges_df = vcat(budget_ranges_df_normal, budget_ranges_df_lognormal);

R"""
library(ggplot2)
d <- $budget_ranges_df
inventory_methods <- c(`1` = "5 plots", `2` = "10 plots", `3` = "20 plots")
model_names <- c(normal = "Normal model", lognormal = "Lognormal model")
inventory_colors <- structure(c("white", "gray", "black"), .Names = unname(inventory_methods))
d$inventory_method_fact <- factor(inventory_methods[d$inventory_method], levels = unname(inventory_methods))
d$model_fact <- factor(model_names[d$model], levels = unname(model_names))
expansion_obj <- expansion(mult = 0.005, add = 0)

plt <- ggplot(d, aes(x = factor(budget), y = factor(order), fill = inventory_method_fact)) +
    geom_tile(width = 1, height = 1) +
    scale_fill_manual(name = "Optimal inventory method", values = inventory_colors) +
    scale_x_discrete(expand = expansion_obj) + 
    scale_y_discrete(expand = expansion_obj) + 
    facet_wrap(~ model_fact, ncol = 1) + 
    guides(fill = guide_legend(override.aes = list(color = "black"))) +
    labs(y = "Tracts in order of increasing total prior volume variance",
         x = "Inventory budget") + 
    theme_bw() +
    THEME_SIZES + 
    THEME_TIGHT +
    theme(legend.position = "top",
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 0.75)) +
    THEME_TIGHT_TOP_LEGEND

fp <- file.path(PLOT_OUTPUT_PATH, "fig-inventory-per-budget-normal-vs-lognormal.pdf")
suppressMessages(ggsave(plot = plt, filename = fp, 
    width = FIGURE_WIDTH_CM, height = 0.8 * FIGURE_WIDTH_CM, units = "cm"))
cat(paste0("Wrote file ", fp, "\n"))
plt
"""
