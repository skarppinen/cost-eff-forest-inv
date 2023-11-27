using JLD2, DataFrames, RCall
project_root = joinpath(@__DIR__, "../../../..");
include(joinpath(project_root, "src/julia/lib/model.jl"));
include("r-fig-helpers.jl");

df_normal = let 
    inputfolder = joinpath(project_root, "no-vc", "output", "run-planning-problem"); 
    inputfiles = readdir(inputfolder, join = true);
    df = get_planning_problem_evaluations(inputfiles);
    df[!, :model] .= "normal";
    df;
end
df_lognormal = let
    inputfolder = joinpath(project_root, "no-vc", "output", "run-planning-problem-lognormal", "output"); 
    inputfiles = readdir(inputfolder, join = true);
    df = get_planning_problem_evaluations(inputfiles);
    df[!, :model] .= "lognormal";
    df;
end
df = vcat(df_normal, df_lognormal);
candidate_budgets = minimum(df[!, :costs]):1000:maximum(df[!, :costs]);
df[!, :budget] = map(df[!, :costs]) do cost 
    i = findfirst(budget -> cost <= budget, candidate_budgets);
    candidate_budgets[i];
end

gdf = groupby(df, [:model, :budget]) 
budget_vs_pov_df = combine(gdf, :povs => maximum => :max_pov);


R"""
d <- $budget_vs_pov_df
model_names <- c(normal = "Normal model", lognormal = "Lognormal model")
d$model_fact <- factor(model_names[d$model], levels = unname(model_names))
plt <- ggplot(d, aes(x = factor(budget), y = max_pov)) +
    geom_point(aes(shape = model_fact), size = 2) +
    labs(y = "Largest PoV", x = "Budget") +
    guides(shape = guide_legend(override.aes = list(size = 5))) +
    theme_bw() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 30, hjust = 0.75)) +
    THEME_SIZES + 
    THEME_TIGHT +
    THEME_TIGHT_TOP_LEGEND +
    theme(legend.title = element_blank())

fp <- file.path(PLOT_OUTPUT_PATH, "fig-budget-vs-pov.pdf")
suppressMessages(ggsave(plot = plt, filename = fp, 
    width = FIGURE_WIDTH_CM, height = 0.6 * FIGURE_WIDTH_CM, units = "cm"))
cat(paste0("Wrote file ", fp, "\n"))
plt
"""

#normal = select(df_normal, :xI, :costs, :povs => :normal_pov)
#lognormal = select(df_lognormal, :xI, :povs => :lognormal_pov)
#both = sort!(leftjoin(normal, lognormal, on = :xI), :costs)
#both[!, :id] = 1:nrow(both);
#both = stack(both, [:normal_pov, :lognormal_pov]);
#
#R"""
#d <- $both
#ggplot(d, aes(x = id, y = value)) +
#    geom_point(aes(color = variable))
#"""

