## Experiment investigates variability of PoV computation given
## number of Monte Carlo samples, number of initialisations and 
## choice of inventory decision. 

include("experiment-config.jl");
args = ArgParse.parse_args(ARGUMENT_CONFIG["run-pov"]);

# Parameters.
strategy = Symbol(args["strategy"]);
ninits = args["ninits"];
nmc = args["nmc"];
model = args["model"];
maxsweeps = args["maxsweeps"];
nruns = args["nruns"];
seed = args["seed"];
dataid = args["dataid"];
xI_choice = args["xI"];
jobid = args["jobid"];
verbose = args["verbose"];
outfolder = args["outfolder"];

@assert ninits > 0
@assert model in (1, 2) "model must be 1 (normal) or 2 (lognormal)"
@assert nmc > 0
@assert maxsweeps > 0
@assert nruns > 0
@assert strategy in (:smallest_negative, :greatest_negative)
@assert dataid in (2, 3, 4) "invalid `datasetid` = $dataid";
@assert xI_choice in (1, 2, 3) "invalid `xI_choice` = $xI"; 

using Random
include(joinpath(@__DIR__, "../lib/model.jl"));
include(joinpath(@__DIR__, "../lib/random-sweep.jl"));
include(joinpath(@__DIR__, "../lib/pov.jl"));
include("experiment-helpers.jl");

verbose && println("Got arguments: $args");
rng = Xoshiro(seed);

function run_pov(rng::AbstractRNG; nmc::Int,
        xI_choice::Int, dataid::Int, maxsweeps::Int, ninits::Int, nruns::Int,
        verbose::Bool, strategy::Symbol, model::Int)
    stda = StandData(id = dataid);
    (nS, nT, nA, nI) = size(stda);
    vp = VolumePosterior(nS, nA);
    rs = RandomSweepStorage(nS, nT; strategy = strategy);
    xI = Cuint.(repeat([xI_choice], nS)) .- Cuint(1); 
    seeds = zeros(Culong, nmc); 
    povs = zeros(Float64, nruns); 
    t1 = time();
    report_interval = max(round(Int, 0.01 * nruns), 1);
    for i in 1:nruns
        rand!(rng, seeds); 
        povs[i] = ccall_PoV!(xI, rs, vp, stda, seeds; 
                             ninits = ninits,
                             maxsweeps = maxsweeps, model = model); 
        if verbose && i % report_interval == 0 
            report_progress(i, nruns, time() - t1); 
        end
    end
    povs;
end

out = run_pov(rng; nmc = nmc, xI_choice = xI_choice, dataid = dataid, 
              maxsweeps = maxsweeps, ninits = ninits, nruns = nruns,
              verbose = verbose, strategy = strategy, model = model);
outfile = joinpath(outfolder, string("jobid-", jobid, ".jld2")); 
mkpath(dirname(outfile));
jldopen(outfile, "w") do file
	file["args"] = args;
    file["povs"] = out;
end
verbose && println("Wrote file $outfile");

