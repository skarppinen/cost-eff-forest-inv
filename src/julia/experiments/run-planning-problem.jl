# Script evaluates feasible points of the planning problem.

include("experiment-config.jl");
args = ArgParse.parse_args(ARGUMENT_CONFIG["run-planning-problem"]);

# Parameters.
strategy = Symbol(args["strategy"]);
ninits = args["ninits"];
nmc = args["nmc"];
model = args["model"];
maxsweeps = args["maxsweeps"];
seed = args["seed"];
datafile = args["datafile"];
l = args["l"];
u = args["u"];
jobid = args["jobid"];
verbose = args["verbose"];
outfolder = args["outfolder"];

verbose && println("Got arguments: $args");
@assert l <= u
@assert ninits > 0
@assert model in (1, 2) "model must be 1 (normal) or 2 (lognormal)"
@assert nmc > 0
@assert maxsweeps > 0
@assert strategy in (:smallest_negative, :greatest_negative)

using Random
include(joinpath(@__DIR__, "../lib/model.jl"));
include(joinpath(@__DIR__, "../lib/random-sweep.jl"));
include(joinpath(@__DIR__, "../lib/pov.jl"));
include("experiment-helpers.jl");

function inventory_cost(xI::AbstractVector{<: Integer}, cost::AbstractVector{Float64}) 
    C = 0.0;
    for xi in xI
        C += cost[xi]; 
    end
    return C;
end

function evaluate(xIs, stda, seeds;  
        datafile::AbstractString, maxsweeps::Int, ninits::Int,
        verbose::Bool, strategy::Symbol, model::Int)
    
    (nS, nT, nA, nI) = size(stda);
    vp = VolumePosterior(nS, nA);
    rs = RandomSweepStorage(nS, nT; strategy = strategy);
    n = length(xIs); 
    povs = zeros(Float64, n); 
    costs = zeros(Float64, n);

    t1 = time();
    report_interval = 1; 
    for i in 1:n
        xI = xIs[i];
        xI_cuint = Cuint.(xI) .- Cuint(1); 
        povs[i] = ccall_PoV!(xI_cuint, rs, vp, stda, seeds; 
                             ninits = ninits,
                             maxsweeps = maxsweeps, model = model); 
        costs[i] = inventory_cost(xI, stda.costs);
        if verbose && i % report_interval == 0 
            println("Evaluation $i: got PoV = $(povs[i])");
            println("Evaluation $i: got cost = $(costs[i])");
            report_progress(i, n, time() - t1); 
        end
    end
    povs, costs;
end

# Read u - l + 1 inventory decisions to evaluate.
datapath = joinpath(@__DIR__, "../../../data", datafile); 
xIs, stda = jldopen(datapath, "r") do file
        @assert l < length(file["xIs"])
        @assert u <= length(file["xIs"])
        file["xIs"][l:u], file["stda"]; 
end
rng = Xoshiro(seed);
seeds = zeros(Culong, nmc); 
rand!(rng, seeds);
povs, costs = evaluate(xIs, stda, seeds; datafile = datafile, 
              maxsweeps = maxsweeps, ninits = ninits, verbose = verbose, 
              strategy = strategy, model = model);

outfile = joinpath(outfolder, string("jobid-", jobid, ".jld2")); 
mkpath(dirname(outfile));
jldopen(outfile, "w") do file
	file["args"] = args;
    file["xIs"] = xIs;
    file["povs"] = povs;
    file["costs"] = costs;
end
verbose && println("Wrote file $outfile");

