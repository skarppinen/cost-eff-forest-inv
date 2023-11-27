import ArgParse

ARGUMENT_CONFIG = Dict();
let
    s = ArgParse.ArgParseSettings();
    ArgParse.@add_arg_table! s begin
        "--jobid"
        help = "An integer used to distinguish separate jobs, used in naming output file"
        arg_type = Int
        required = true

        "--algseed"
        help = "Random seed used in algorithm during experiment"
        arg_type = Int
        required = true

        "--dataseed"
        help = "Random seed used to simulate data"
        arg_type = Int
        default = 29032023

        "--dataid"
        help = "Stand dataset to use, must be 4 (checked)"
        arg_type = Int
        required = true

        "--nruns"
        help = "Times algorithm is run"
        arg_type = Int
        required = true

        "--maxsweeps"
        help = "Maximum number of sweeps to do in algorithm"
        arg_type = Int
        default = 9999

        "--ninits"
        help = "Number of initialisations of algorithm to run per each run of the algorithm"
        arg_type = Int
        required = true

        "--xI"
        help = "Inventory method choice for each stand. should take value in {1, 2, 3} (checked)."
        arg_type = Int
        required = true
        
        "--strategy"
        help = "strategy used in descent: smallest_negative, random_negative or greatest_negative"
        arg_type = String
        required = true
        
        "--outfolder", "-o"
        help = """Folder to place results relative to path from which script is invoked. An absolute \
        path may also be provided. The folder is created if it does not exists.

        """
        arg_type = String
        default = "output"

        "--verbose"
        help = "display progress reports?" 
        action = :store_true
    end;
    global ARGUMENT_CONFIG["run-random-sweep"] = s;
end

let
    s = ArgParse.ArgParseSettings();
    ArgParse.@add_arg_table! s begin
        "--jobid"
        help = "An integer used to distinguish separate jobs, used in naming output file"
        arg_type = Int
        required = true

        "--seed"
        help = "Random seed"
        arg_type = Cuint
        required = true

        "--model"
        help = "Integer specifying the model to use. must be 1 (normal) or 2 (lognormal). checked"
        arg_type = Int
        required = true

        "--nmc"
        help = "Number of Monte Carlo samples"
        arg_type = Int
        required = true

        "--maxsweeps"
        help = "Maximum number of sweeps to do in random sweeping algorithm"
        arg_type = Int
        default = 9999

        "--ninits"
        help = "Number of initialisations of random sweeping algorithm to run per each Monte Carlo sample"
        arg_type = Int
        required = true
        
        "--nruns"
        help = "Number of computations of PoV (with `nmc` and `ninits`)"
        arg_type = Int
        required = true

        "--xI"
        help = "Inventory method choice for each stand. should take value in {1, 2, 3} (checked)."
        arg_type = Int
        required = true
        
        "--strategy"
        help = "strategy used in descent: smallest_negative or greatest_negative"
        arg_type = String
        required = true

        "--dataid"
        help = "dataset to use, must be 4 (checked)"
        arg_type = Int
        required = true
        
        "--outfolder", "-o"
        help = """Folder to place results relative to path from which script is invoked. An absolute \
        path may also be provided. The folder is created if it does not exists.

        """
        arg_type = String
        default = "output"

        "--verbose"
        help = "display progress reports?" 
        action = :store_true
    end;
    global ARGUMENT_CONFIG["run-pov"] = s;
end

let
    s = ArgParse.ArgParseSettings();
    ArgParse.@add_arg_table! s begin
        "--jobid"
        help = "An integer used to distinguish separate jobs, used in naming output file"
        arg_type = Int
        required = true

        "--seed"
        help = "Random seed"
        arg_type = Cuint
        required = true

        "--model"
        help = "Integer specifying the model to use. must be 1 (normal) or 2 (lognormal). checked"
        arg_type = Int
        required = true

        "--nmc"
        help = "Number of Monte Carlo samples"
        arg_type = Int
        required = true

        "--maxsweeps"
        help = "Maximum number of sweeps to do in random sweeping algorithm"
        arg_type = Int
        default = 9999

        "--ninits"
        help = "Number of initialisations of random sweeping algorithm to run per each Monte Carlo sample"
        arg_type = Int
        required = true
        
        "--strategy"
        help = "strategy used in descent: smallest_negative or greatest_negative"
        arg_type = String
        required = true

        "--datafile"
        help = """dataset to use, looked up from the data folder at project root. \
        the dataset should contain both the stand data and the inventory decisions \
        to consider, see arguments `l` and `u`."""
        arg_type = String
        required = true

        "--l"
        help = """lower index marking the index of the first inventory decision to lookup from `datafile`. \
        must be less than the maximal number of inventory decisions in `datafile` (checked). used \
        together with argument `l` to select `u - l + 1` inventory decisions to evaluate from `datafile`. \
        `l` must be less than `u` (checked)"""
        arg_type = Int
        required = true

        "--u"
        help = """upper index marking the index of the last inventory decision to lookup from `datafile`. \
        must be less than the maximal number of inventory decisions in `datafile` (checked). used \
        together with argument `l` to select `u - l + 1` inventory decisions to evaluate from `datafile`. \
        `l` must be less than `u` (checked)"""
        arg_type = Int
        required = true

        
        "--outfolder", "-o"
        help = """Folder to place results relative to path from which script is invoked. An absolute \
        path may also be provided. The folder is created if it does not exists.
        """
        arg_type = String
        default = "output"

        "--verbose"
        help = "display progress reports?" 
        action = :store_true
    end;
    global ARGUMENT_CONFIG["run-planning-problem"] = s;
end

let
    s = ArgParse.ArgParseSettings();
    ArgParse.@add_arg_table! s begin
        "--jobid"
        help = "An integer used to distinguish separate jobs, used in naming output file"
        arg_type = Int
        required = true

        "--dataid"
        help = "index of stand data to use, must be 4 (checked)"
        arg_type = Int
        required = true 

        "--dataseed"
        help = "Random seed for data generation"
        arg_type = Int
        default = 29032023
        
        "--timelimit"
        help = "Maximum time in seconds to spend in finding global minimum"
        arg_type = Float64
        default = 12.0 * 60 * 60

        "--xI"
        help = "Inventory decision (for all stands), must be in (1, 2, 3) (checked)"
        arg_type = Int
        required = true

        "--outfolder", "-o"
        help = """Folder to place results relative to path from which script is invoked. An absolute \
        path may also be provided. The folder is created if it does not exists.
        """
        arg_type = String
        default = "output"

        "--verbose"
        help = "display progress reports?" 
        action = :store_true
    end;
    global ARGUMENT_CONFIG["run-global-optim"] = s;
end
