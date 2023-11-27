# Experiment investigates variance of results of random sweep algorithm with different configurations.

using BenchmarkTools, JLD2
include(joinpath(@__DIR__, "experiment-config.jl"));

# Extract arguments.
args = ArgParse.parse_args(ARGUMENT_CONFIG["run-random-sweep"]);
jobid = args["jobid"]; 
dataid = args["dataid"];
algseed = args["algseed"];
dataseed = args["dataseed"];
max_sweeps = args["maxsweeps"];
inits = args["ninits"];
nruns = args["nruns"]; 
xI_choice = args["xI"]; 
verbose = args["verbose"];
strategy = Symbol(args["strategy"]);
outfolder = args["outfolder"];

verbose && println("Got arguments: $args");

include(joinpath(@__DIR__, "../lib/model.jl"));
include(joinpath(@__DIR__, "../lib/random-sweep.jl"));
include("experiment-helpers.jl"); 

function run_random_sweep_true_data(dataid::Integer, algrng::AbstractRNG, datarng::AbstractRNG; 
        max_sweeps::Int, inits::Int, nruns::Int, xI_choice::Int,
        verbose::Bool = true, strategy::Symbol)

	@assert nruns > 0
	@assert xI_choice in (1, 2, 3)
	@assert max_sweeps > 0
	@assert inits > 0
    @assert dataid == 4
    
    stda = StandData(id = dataid);
    (nS, nT, nA, nI) = size(stda);

	xI = repeat([xI_choice], nS);
	vp = random_posterior(stda, xI, datarng);
    prob = build_problem(vp, stda.demands);
    scaling = mean(prob.Q); # To get a more reasonable scale for input data. (converting back below)
	rs = RandomSweepStorage(prob.Q ./ scaling, prob.c ./ scaling, prob.r ./ scaling;
                            strategy = strategy)

	objvalues = zeros(Float64, nruns);	
    t1 = time();
	for i in eachindex(objvalues)
		#objvalues[i] = random_sweep!(rs, algrng; max_sweeps, inits, verbose = false);	
        seed = rand(algrng, Culong);
        objvalues[i] = ccall_random_sweep!(rs, seed; max_sweeps = max_sweeps, inits = inits);
        verbose && i % 100 == 0 && report_progress(i, nruns, time() - t1); 
	end
    objvalues .= scaling .* objvalues;

	# Timing.
	#bout = @benchmark random_sweep!($rs, $algrng; max_sweeps = $max_sweeps, inits = $inits, verbose = false); 

	objvalues #, bout.times;
end

# Do experiment.
algrng = Xoshiro(algseed); 
datarng = Xoshiro(dataseed);
objvalues = run_random_sweep_true_data(dataid, algrng, datarng; 
								max_sweeps, inits, nruns, xI_choice, verbose, strategy)

outfile = joinpath(outfolder, string("jobid-", jobid, ".jld2")); 
mkpath(dirname(outfile));
jldopen(outfile, "w") do file
	file["args"] = args;
    file["objvalues"] = objvalues;
	#file["timings"] = timings;
end
verbose && println("Wrote file $outfile");
				                     


