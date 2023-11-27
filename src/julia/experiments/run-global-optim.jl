# Experiment investigates true objective function value for particular
# inventory decision and sampled dataset.

include(joinpath(@__DIR__, "experiment-config.jl"));
args = ArgParse.parse_args(ARGUMENT_CONFIG["run-global-optim"]);
jobid = args["jobid"];
timelimit = args["timelimit"];
dataseed = args["dataseed"];
dataid = args["dataid"];
xI_choice = args["xI"]; 
verbose = args["verbose"];
outfolder = args["outfolder"];
verbose && println("Got arguments: $args");
@assert xI_choice in (1, 2, 3)
@assert dataid == 4
@assert timelimit > 0.0

using JLD2
using JuMP
using Gurobi
include(joinpath(@__DIR__, "../lib/model.jl"));
include(joinpath(@__DIR__, "../lib/random-sweep.jl"));

stda = StandData(id = dataid);
nS, nT, nA, nI = size(stda);
xI = repeat([xI_choice], nS);
rng = Xoshiro(dataseed);
vp = random_posterior(stda, xI, rng); 
prob = build_problem(vp.muplus, vp.sigma2plus, stda.demands);

#A = sparse_hcat([sparse(Matrix{Float64}(I, nS, nS)) for i in 1:nT]...);
A = zeros(nS, nS * nT);
b = repeat([1], nS); 
Q = zeros(nS * nT, nS * nT);
for i in 1:nT
    l = (i - 1) * nS + 1;
    u = i * nS;
    Q[l:u, l:u] .= prob.Q;
    A[1:nS, l:u] .= Matrix{Float64}(I, nS, nS);
end
r = prob.r;
c = prob.c;
const RS = RandomSweepStorage(prob.Q, prob.c, prob.r, strategy = :greatest_negative)

function random_timing!(x)
    sumx = sum(x);
    if rand() < (1.0 - sumx)
        x .= 0.0;
    else
        i = wsample(x);
        x .= 0.0;
        x[i] = 1.0;
    end
    nothing;
end

function build_model(Q, c, r, A, b, solstart; timelimit)
    N = size(A, 2);
    model = direct_model(Gurobi.Optimizer());
    @variable(model, x[1:N], Bin);
    @constraint(model, cons, A * x .<= b) 
    @objective(model, Min, 
               0.5 * x' * Q * x + c' * x + r 
    )
    for i in 1:length(model[:x])
        MOI.set(model, Gurobi.VariableAttribute("Start"), model[:x][i], solstart[i]) 
    end
    set_optimizer_attribute(model, "MIPFocus", 3);
    set_optimizer_attribute(model, "Cuts", 2);
    set_optimizer_attribute(model, "TimeLimit", timelimit);

    cb = let model = model, N = N 
        function f(cb_data, cb_where::Cint)
            # This is run when Gurobi finds new feasible solution.
            if cb_where == GRB_CB_MIPSOL
                Gurobi.load_callback_variable_primal(cb_data, cb_where);
                for i in 1:N
                    RS.X[i] = callback_value(cb_data, model[:x][i]);
                end
                obj = random_sweep!(RS, Random.GLOBAL_RNG; max_sweeps = 999999, inits = 0);  
                MOI.submit(model, MOI.HeuristicSolution(cb_data), model[:x], vec(RS.X))
            #elseif cb_where == GRB_CB_MIPNODE 
            #    resultP = Ref{Cint}();
            #    GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
            #    if resultP[] != GRB_OPTIMAL
            #        return  # Solution is something other than optimal.
            #    end
            #    # Now safe to do this.
            #    Gurobi.load_callback_variable_primal(cb_data, cb_where);
            #    for i in 1:N
            #        RS.X[i] = callback_value(cb_data, model[:x][i]);
            #    end
            #    for s in 1:nS
            #        random_timing!(view(RS.X, s, :));  
            #    end
            #    random_sweep!(RS, Random.GLOBAL_RNG; max_sweeps = 999999, inits = 0);  
            #    MOI.submit(model, MOI.HeuristicSolution(cb_data), model[:x], vec(RS.X))
            end
            return;
        end
    end
    MOI.set(model, Gurobi.CallbackFunction(), cb)
    model;
end

# Get initial feasible solution using heuristic.
rs = RandomSweepStorage(prob.Q, prob.c, prob.r, strategy = :greatest_negative)
algrng = Random.GLOBAL_RNG;
objvalue = random_sweep!(rs, algrng, max_sweeps = 9999999, inits = 30000);
model = build_model(Q, vec(c), r, A, b, vec(rs.Xopt); timelimit = timelimit);
optimize!(model);

outfile = joinpath(outfolder, string("jobid-", jobid, ".jld2")); 
mkpath(dirname(outfile));
jldopen(outfile, "w") do file
	file["args"] = args;
    file["objvalue"] = objective_value(model);
    file["objbound"] = objective_bound(model);
end
verbose && println("Wrote file $outfile");

