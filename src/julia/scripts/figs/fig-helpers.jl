PAPER_WIDTH_CM = 21; 
MARGIN_LEFT_CM = 2; MARGIN_RIGHT_CM = 2;
FIGURE_WIDTH_CM = PAPER_WIDTH_CM - MARGIN_LEFT_CM - MARGIN_RIGHT_CM;
PLOT_OUTPUT_PATH = joinpath(@__DIR__, "../../../..", "plots"); 
mkpath(PLOT_OUTPUT_PATH);
function get_planning_problem_evaluations(filepaths::AbstractVector{String})
    f = filepath -> begin
        jldopen(filepath, "r") do file
            args = file["args"];
            jobids = args["l"]:args["u"];
            xIs = file["xIs"];
            DataFrame(povs = file["povs"],
                      costs = file["costs"],
                      jobid = jobids,
                      xI = xIs);
        end
    end
    mapreduce(f, vcat, filepaths);
end

