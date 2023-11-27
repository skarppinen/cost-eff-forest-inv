include(joinpath(@__DIR__, "../lib/stand-ordering.jl"));

stda = StandData("test-data-4.jld2");
xIs = feasible_set_reduction(stda);
filepath = joinpath(@__DIR__, "../../..", "data", "planning-problem-data.jld2");
jldopen(filepath, "w") do file
    file["stda"] = stda;
    file["xIs"] = xIs;
end
println("Wrote file $filepath");
