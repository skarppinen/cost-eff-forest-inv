# Script generates tasklists for running experiments.
# The output tasklist files assume that the environment variables 
# SKOGFORSK_EXPERIMENT_FOLDER and 
# SKOGFORSK_HOME_FOLDER are set.
include(joinpath(@__DIR__, "../lib/build-task-list.jl"));
include(joinpath(@__DIR__, "../experiments/experiment-config.jl"));
using Random, Dates, JLD2
TASKLIST_OUTPATH = joinpath(@__DIR__, "../../bash");
SBATCH_DEFAULTS = Dict(
    "account" => "jkarvane",
    "partition" => "small",
    "mem-per-cpu" => "4000",
    "cpus-per-task" => "1"
);

function sbatch_setting_string(settings::Dict, modules::Vector{String})
    out = raw"#!/bin/bash" * "\n";
    for k in keys(settings)
        out *= string(raw"#SBATCH --", k, "=", settings[k]) * "\n"; 
    end
    out *= "\n";
    for mod in modules
        out *= string("module load ", mod, "\n");
    end
    out;
end

function pad2(x)
    lpad(x, 2, '0');  
end

function hms(t::DateTime)
    ndays = day(t) - 1;
    ndays >= 10 && throw(ArgumentError("runtime over 10 days.."));
    if ndays == 0
        string(pad2(hour(t)), ":", pad2(minute(t)), ":", pad2(second(t))); 
    else
        string(pad2(ndays), "-", pad2(hour(t)), ":", pad2(minute(t)), ":", pad2(second(t))); 
    end
end

# Writes a zip file containing all jobs to `TASKLIST_OUTPATH`.
function compress_taskfiles(experiment::AbstractString, 
        jobstrs::AbstractVector{String})
    tmpoutfolder = mktempdir(prefix = experiment);
    outfiles = Vector{String}(undef, 0);
    for (jobid, jobstr) in enumerate(jobstrs)
        outfile = string(experiment, "-", jobid);
        open(joinpath(tmpoutfolder, outfile), "w") do io
            write(io, jobstr);
        end;   
        push!(outfiles, outfile);
    end
    archive_filepath = joinpath(TASKLIST_OUTPATH, experiment * "-tasklists.zip");
    cmd = Cmd(`zip -r $archive_filepath -q $outfiles`, dir = tmpoutfolder);
    run(pipeline(cmd, stdout = devnull)); 
end

let experiment = "run-random-sweep"
    outfile = experiment * "-tasklist";
    script_name = joinpath(raw"${SKOGFORSK_EXPERIMENT_FOLDER}", experiment);
    prog = raw"srun julia --project=${SKOGFORSK_HOME}";
    fn = "";
    dataid = 4;
    nruns = 10000;
    maxsweeps = 1000_000; # Just to be sure that method converges. Does not really affect run time.
    dataseed = 29032023; # This needs to be constant for this experiment.
    jobid = 0;
    verbose = true;
    jobstrs = Vector{String}(undef, 0);
    modules = ["julia/1.8.5"];
    additional_time = Minute(20);
    sec_per_run_1000_inits = Dict("smallest_negative" => 0.125,
                                  "greatest_negative" => 0.21);

    ninits_choices = [5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 20000];
    xI_choices = [1, 2, 3];
    strategy_choices = ["smallest_negative", "greatest_negative"];
    ntasks = length(ninits_choices) * length(xI_choices) * length(strategy_choices); 
    Random.seed!(28032023);
    algseeds = abs.(rand(Int, ntasks)); 
    @assert length(algseeds) == length(unique(algseeds)) "duplicate algseeds"
    maximal_time = DateTime(0);

    for ninits in ninits_choices, xI in xI_choices, strategy in strategy_choices 
        jobid += 1;
        
        # Estimate runtime based on small tests, 1.4 is a multiplier.
        multiplier = 1.25;
        sec_per_run_1_init = multiplier * sec_per_run_1000_inits[strategy] / 1000;
        est_time_sec = ceil(sec_per_run_1_init * ninits * nruns);
        estimated_time = DateTime(0) + Second(est_time_sec) + additional_time;
        if estimated_time > maximal_time
            maximal_time = estimated_time;
        end

        args = Dict("jobid" => jobid,
                    "algseed" => algseeds[jobid],
                    "dataseed" => dataseed,
                    "dataid" => dataid,
                    "nruns" => nruns,
                    "ninits" => ninits,
                    "strategy" => strategy,
                    "xI" => xI,
                    "maxsweeps" => maxsweeps,
                    "verbose" => verbose);
        sbatch_config = Dict(SBATCH_DEFAULTS...,
                             "job-name" => string(experiment, "-", jobid),
                             "time" => hms(estimated_time),
                             "output" => string(jobid, "-", "output.txt"),
                             "error" => string(jobid, "-", "error.txt"),
                             "mem-per-cpu" => "1500")

        sbatch_str = sbatch_setting_string(sbatch_config, modules);
        fn = generate_call(script_name, args, ARGUMENT_CONFIG[experiment]; prog = prog);
        push!(jobstrs, string(sbatch_str, fn));
    end
    outfolder = joinpath(TASKLIST_OUTPATH);
    mkpath(outfolder);
    compress_taskfiles(experiment, jobstrs);    
    msg = string("Wrote ", length(jobstrs), " jobs",
                 " for experiment ", experiment);
    print(msg);
    println(", maximal jobtime is $(day(maximal_time) - 1)d, $(hour(maximal_time))h, $(minute(maximal_time))m");
end

let experiment = "run-pov"
    outfile = experiment * "-tasklist";
    script_name = joinpath(raw"${SKOGFORSK_EXPERIMENT_FOLDER}", experiment);
    prog = raw"srun julia --project=${SKOGFORSK_HOME}";
    fn = "";
    nruns = 1000;
    dataid = 4;
    maxsweeps = 9999; 
    rng = Xoshiro(02052023); 
    jobid = 0;
    verbose = true;
    jobstrs = Vector{String}(undef, 0);
    modules = ["julia/1.8.5", "gsl"];
    sec_per_unit = Dict("smallest_negative" => 
                        Dict(1 => 1.26 / (100 * 100)),
            #                     2 => 6.5 / (100 * 100)),
                        "greatest_negative" => 
                        Dict(1 => 2.1 / (100 * 100)));
            #                     2 => 7.4 / (100 * 100)));

    ninits_choices = [50, 100, 250, 500];
    nmc_choices = [50, 100, 250, 500];
    xI_choices = [1, 2, 3];
    strategy_choices = ["smallest_negative", "greatest_negative"];
    model_choices = [1]#, 2];
    ntasks = length(ninits_choices) * length(xI_choices) * length(nmc_choices) * 
        length(strategy_choices) * length(model_choices); 
    maximal_time = DateTime(0);
    additional_time = Minute(20);

    for model in model_choices, ninits in ninits_choices, nmc in nmc_choices, 
        xI in xI_choices, strategy in strategy_choices  

        jobid += 1;
        
        # Estimate runtime based on small tests.
        multiplier = 1.25;
        est_time_sec = ceil(multiplier * sec_per_unit[strategy][model] * ninits * nmc * nruns);
        estimated_time = DateTime(0) + Second(est_time_sec) + Minute(additional_time);
        if estimated_time > maximal_time
            maximal_time = estimated_time;
        end

        args = Dict("jobid" => jobid,
                    "seed" => rand(rng, Cuint), 
                    "nruns" => nruns,
                    "ninits" => ninits,
                    "strategy" => strategy,
                    "nmc" => nmc,
                    "model" => model,
                    "xI" => xI,
                    "dataid" => dataid,
                    "maxsweeps" => maxsweeps,
                    "verbose" => verbose);
        sbatch_config = Dict(SBATCH_DEFAULTS...,
                             "job-name" => string(experiment, "-", jobid),
                             "time" => hms(estimated_time),
                             "output" => string(jobid, "-", "output.txt"),
                             "error" => string(jobid, "-", "error.txt"),
                             "mem-per-cpu" => "2000")

        sbatch_str = sbatch_setting_string(sbatch_config, modules);
        fn = generate_call(script_name, args, ARGUMENT_CONFIG[experiment]; prog = prog);
        push!(jobstrs, string(sbatch_str, fn));
    end
    outfolder = joinpath(TASKLIST_OUTPATH);
    mkpath(outfolder);
    compress_taskfiles(experiment, jobstrs);
    msg = string("Wrote ", length(jobstrs), " jobs",
                 " for experiment ", experiment);
    print(msg);
    println(", maximal jobtime is $(day(maximal_time) - 1)d, $(hour(maximal_time))h, $(minute(maximal_time))m");
end

let experiment = "run-planning-problem"
    outfile = experiment * "-tasklist";
    script_name = joinpath(raw"${SKOGFORSK_EXPERIMENT_FOLDER}", experiment);
    prog = raw"srun julia --project=${SKOGFORSK_HOME}";
    fn = "";
    datafile = "planning-problem-data.jld2";
    maxsweeps = 9999; 
    seed = 05062023; 
    verbose = true;
    jobstrs = Vector{String}(undef, 0);
    modules = ["julia/1.8.5", "gsl"];
    sec_per_unit = Dict(1 => 622 * 2.5 / 6,
                        2 => 622 * 2.5);

    evals_per_job = 26;
    njobs_cap = 200;
    total_evals_needed = jldopen(joinpath(@__DIR__, "../../../data", datafile), "r") do file
        length(file["xIs"]);
    end
    if total_evals_needed / evals_per_job > njobs_cap
        throw(ArgumentError("insufficient `njobs_cap = $njobs_cap`, adjust parameters.."));
    end
    if total_evals_needed % evals_per_job == 0
        njobs = convert(Int, total_evals_needed / evals_per_job);
        us = collect(1:njobs) * evals_per_job;
        ls = us .- (evals_per_job - 1);
    else
        njobs = ceil(Int, total_evals_needed / evals_per_job); 
        us = collect(1:(njobs - 1)) * evals_per_job;
        ls = us .- (evals_per_job - 1); 

        # "residual job" (has less than `evals_per_job`).
        push!(ls, us[end] + 1);
        push!(us, total_evals_needed);
    end

    ninits = 2500;
    nmc = 2500; 
    strategy = "greatest_negative"; 
    models = [1, 2];
    maximal_time = DateTime(0);
    additional_time = Minute(20);
    jobid_running = 0;
    
    
    for model in models
        for jobid in 1:njobs 
            jobid_running += 1;
            n_evals = us[jobid] - ls[jobid] + 1;
            # Estimate runtime based on small tests.
            multiplier = 1.5;
            est_time_sec = ceil(multiplier * sec_per_unit[model] * n_evals);
            estimated_time = DateTime(0) + Second(est_time_sec) + Minute(additional_time);
            if estimated_time > maximal_time
                maximal_time = estimated_time;
            end

            args = Dict("jobid" => jobid_running,
                        "seed" => seed, 
                        "ninits" => ninits,
                        "strategy" => strategy,
                        "nmc" => nmc,
                        "model" => model,
                        "datafile" => datafile,
                        "l" => ls[jobid],
                        "u" => us[jobid],
                        "maxsweeps" => maxsweeps,
                        "verbose" => verbose);
            sbatch_config = Dict(SBATCH_DEFAULTS...,
                                 "job-name" => string(experiment, "-", jobid_running),
                                 "time" => hms(estimated_time),
                                 "output" => string(jobid_running, "-", "output.txt"),
                                 "error" => string(jobid_running, "-", "error.txt"),
                                 "mem-per-cpu" => "2000")

            sbatch_str = sbatch_setting_string(sbatch_config, modules);
            fn = generate_call(script_name, args, ARGUMENT_CONFIG[experiment]; prog = prog);
            push!(jobstrs, string(sbatch_str, fn));
        end
    end
    outfolder = joinpath(TASKLIST_OUTPATH);
    mkpath(outfolder);
    
    compress_taskfiles(experiment, jobstrs);
    msg = string("Wrote ", length(jobstrs), " jobs",
                 " for experiment ", experiment);
    print(msg);
    println(", maximal jobtime is $(day(maximal_time) - 1)d, $(hour(maximal_time))h, $(minute(maximal_time))m");
end


