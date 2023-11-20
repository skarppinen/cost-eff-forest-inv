using Random
using LinearAlgebra

const DESCENT_CHOICES = (:smallest_negative, :random_negative, :greatest_negative);
const DESCENT_ID = (smallest_negative = 0, 
                    greatest_negative = 1,
                    random_negative = 2);

struct DescentStrategy{S}
    tmp::Vector{Int}
    function DescentStrategy(strategy::Symbol, n::Integer) 
        @assert strategy in DESCENT_CHOICES "`strategy` must be in $DESCENT_CHOICES"; 
        tmp = zeros(Int, n + 1);
        new{strategy}(tmp);
    end
end

get_strategy(ds::DescentStrategy{S}) where S = S;

struct RandomSweepStorage{S}
    ds::DescentStrategy{S}
    X::Matrix{Float64}
    Xopt::Matrix{Float64}
    C::Matrix{Float64}
    Qb::Matrix{Float64}
    QbX::Matrix{Float64}
    Qbx::Vector{Float64}
    order::Vector{Int}
    status_lu::Vector{Int}
    increments::Vector{Float64}
    r::Base.RefValue{Float64}
    function RandomSweepStorage(Q::AbstractMatrix{Float64}, C::AbstractMatrix{Float64}, r::AbstractFloat;
            strategy::Symbol = :smallest_negative)
        @assert issymmetric(Q) "`Q` must be a symmetric matrix."
        nS = size(Q, 1);
        nT = size(C, 2);
        @assert size(C, 1) == nS "invalid first dimension of `C` matrix, got $(size(C, 1)), should be $nS";
        X = zeros(Float64, nS, nT);
        Xopt = zeros(Float64, size(X));
        QbX = zeros(Float64, nS, nT);
        Qbx = zeros(Float64, nS);
        order = collect(1:nS);
        status_lu = zeros(Int, nS); 
        increments = zeros(Float64, nT + 1); 
        ds = DescentStrategy(strategy, nT);
        rref = Ref(r);
        new{strategy}(ds, copy(X), Xopt, copy(C), copy(Q), 
                      QbX, Qbx, order, status_lu, increments, rref);
    end
end

# Easily initialise storage.
function RandomSweepStorage(nS::Integer, nT::Integer; strategy::Symbol = :smallest_negative)
    RandomSweepStorage(zeros(nS, nS), zeros(nS, nT), 0.0; strategy = strategy); 
end

struct RandomSweepStorageC
    X::Ptr{Cdouble}
    Xopt::Ptr{Cdouble}
    C::Ptr{Cdouble}
    Qb::Ptr{Cdouble}
    QbX::Ptr{Cdouble}
    Qbx::Ptr{Cdouble}
    order::Ptr{Cuint}
    status_lu::Ptr{Cuint}
    increments::Ptr{Cdouble}
    r::Cdouble
    nS::Cuint
    nT::Cuint
    strategy::Cint
    function RandomSweepStorageC(rs::RandomSweepStorage)
        @assert rs.ds != :random_negative "not implemented!";
        did = getfield(DESCENT_ID, get_strategy(rs.ds)); 
        nS, nT = size(rs.X);
        new(pointer(rs.X), pointer(rs.Xopt), pointer(rs.C), 
            pointer(rs.Qb), pointer(rs.QbX), pointer(rs.Qbx), pointer(rs.order), 
            pointer(rs.status_lu), pointer(rs.increments), rs.r[], 
            Cuint(nS), Cuint(nT), did);
    end
end

function random_feasible_sol!(X::AbstractMatrix{Float64}, rng::AbstractRNG = Random.GLOBAL_RNG)
    (nS, nT) = size(X);
    X .= zero(eltype(X));
    for i in 1:nS
        j = floor(Int, rand(rng) * (nT + 1)); # + 1 for case of zero.
        j == 0 && (continue;);
        @inbounds X[i, j] = 1.0; 
    end
    nothing;
end

function build_sol(status_lu::AbstractVector{Int}, nT::Int)
    nS = length(status_lu);
    X = zeros(Float64, nS, nT);
    for (i, s) in enumerate(status_lu)
        s == 0 && (continue;);
        X[i, s] = 1.0;
    end
    return X;
end

#"""
#Get increment when stand `s` transitions from status `from` to status `to`.
#Note that these statuses are not indexes, since they include the status 0.
#"""
#function _get_increment(s::Integer, from::Integer, to::Integer, X, cmat, Qblock) 
#    from == to && (return 0.0); 
#    @inbounds qss = Qblock[s, s];
#    vQ = view(Qblock, :, s);
#    if from == 0
#        # NOTE: `to` is non-zero here.
#        return @inbounds cmat[s, to] + dot(vQ, view(X, :, to)) + 0.5 * qss; 
#    end
#    if to == 0
#        # NOTE: `from` is non-zero here.
#        return @inbounds -cmat[s, from] - dot(vQ, view(X, :, from)) + 0.5 * qss;
#    end
#    # NOTE: `from` and `to` must be non-zero here.
#    return @inbounds -cmat[s, from] + cmat[s, to] +
#           dot(vQ, view(X, :, to)) - dot(vQ, view(X, :, from)) + qss;  
#end

"""
Get increment when stand `s` transitions from status `from` to status `to`.
Note that these statuses are not indexes, since they include the status 0.
"""
function _get_increment(s::Integer, from::Integer, to::Integer, C::AbstractMatrix{Float64}, 
                        QbX::AbstractMatrix{Float64}, qss::Float64) 
    from == to && (return 0.0); 
    if from == 0
        # NOTE: `to` is non-zero here.
        return @inbounds C[s, to] + QbX[s, to] + 0.5 * qss; 
    end
    if to == 0
        # NOTE: `from` is non-zero here.
        return @inbounds -C[s, from] - QbX[s, from] + 0.5 * qss;
    end
    # NOTE: `from` and `to` must be non-zero here.
    return @inbounds -C[s, from] + C[s, to] + QbX[s, to] - QbX[s, from] + qss;  
end

@inline function fill_increments!(rs::RandomSweepStorage, stand::Int, status::Int)
    for j in eachindex(rs.increments)
        #@inbounds rs.increments[j] = _get_increment(stand, status, j - 1, rs.X, rs.C, rs.Qb);  
        @inbounds rs.increments[j] = _get_increment(stand, status, j - 1, rs.C, rs.QbX, rs.Qb[stand, stand]);  
    end
    nothing;
end

"""
Compute objective sum_{t = 1}^{n_T} 0.5 * x_t^' * Qb * x_t + sum_{t=1}^{n_T} c_t^' * x_t + r 
where x_t and c_t are columns of `X` and `C`
"""
function objective(X, C, Qb, r)
    nT = size(X, 2);
    obj = 0.0;
    for t in 1:nT
        vx = view(X, :, t);
        obj += 0.5 * dot(vx, Qb, vx);     
        obj += dot(view(C, :, t), vx);
    end
    obj += r;
    obj;
end

function init_statuses!(status_lu::AbstractVector{Int}, X::AbstractMatrix{Float64}) 
    @assert length(status_lu) == size(X, 1) "length of `status_lu` must match number of rows of `X`"; 
    for i in eachindex(status_lu)
        status = findfirst(x -> isapprox(1.0, x), view(X, i, :));
        if status == nothing
            @inbounds status_lu[i] = 0;
        else
            @inbounds status_lu[i] = status; 
        end
    end
    nothing;
end



@inline function random_feasible_sol!(rs::RandomSweepStorage, rng::AbstractRNG = Random.GLOBAL_RNG)
    random_feasible_sol!(rs.X, rng);
end

@inline function init_statuses!(rs::RandomSweepStorage)
    init_statuses!(rs.status_lu, rs.X); 
end

@inline function objective(rs::RandomSweepStorage)
    objective(rs.X, rs.C, rs.Qb, rs.r[]); 
end

function is_feasible(X::AbstractMatrix{Float64})
    all(mapslices(sum, X, dims = 2)[:, 1] .<= 1.0);
end

function update_solution!(rs::RandomSweepStorage, stand::Int, status::Int, newstatus::Int)
    status == newstatus && (return nothing;);
    if status == 0
        # NOTE: `newstatus` can't be zero here.
        @inbounds rs.X[stand, newstatus] = 1.0;
        return nothing;
    end
    if newstatus == 0
        # Note: `status` can't be zero here.
        @inbounds rs.X[stand, status] = 0.0;
        return nothing;
    end
    @inbounds rs.X[stand, newstatus] = 1.0;
    @inbounds rs.X[stand, status] = 0.0;
    nothing;
end

"""
Function updates the matrix `rs.QbX` i.e `rs.Qb * rs.X`, knowing that stand `stand` 
transitions from `status` to `newstatus`. This implies moving a one in the matrix `X` and has very little
implications for `rs.QbX`, therefore this function exists.
"""
function update_dotprod_cache!(rs::RandomSweepStorage, stand::Integer, status::Integer,
                               newstatus::Integer)
    status == newstatus && (return nothing;);
    vQb = view(rs.Qb, :, stand);
    if status == 0
        # Note `newstatus` is nonzero here.
        @views @inbounds rs.QbX[:, newstatus] .+= vQb;
        return nothing;
    end
    if newstatus == 0
        # Note `status` is nonzero here.
        @views @inbounds rs.QbX[:, status] .-= vQb;
        return nothing;
    end
    # Both nonzero and different.
    @views @inbounds rs.QbX[:, newstatus] .+= vQb;
    @views @inbounds rs.QbX[:, status] .-= vQb;
    nothing;
end


function random_sweep(X::AbstractMatrix{Float64}, C::AbstractMatrix{Float64}, 
							 Qb::AbstractMatrix{Float64}, r::Float64, 
							 rng::AbstractRNG = Random.GLOBAL_RNG;
							 max_sweeps::Int = 10, inits::Int = 1, 
							 verbose::Bool = false)
    @assert max_sweeps > 0 "`max_sweeps` should be > 0."
    rs = RandomSweepStorage(X, C, Qb; r = r);
    random_sweep!(rs, rng; max_sweeps = max_sweeps, inits = inits, verbose = verbose);
end

function random_sweep!(rs::RandomSweepStorage, rng::AbstractRNG = Random.GLOBAL_RNG; 
                       max_sweeps::Int = 10, inits::Int = 1, verbose::Bool = false)
    nS, nT = size(rs.X);
    best_objvalue = Inf;
    rs.order .= 1:nS; # Not strictly necessary but required to ensure reproducibility (wrt. rng).

    # `inits <= 0` activates mode where current solution in `rs` is used.
    if inits <= 0
        initialise = false;
        inits = 1; # So that following loop runs once.
    else
        initialise = true;
    end
    for k in Base.OneTo(inits)
        if initialise 
            random_feasible_sol!(rs, rng);
        end
        mul!(rs.QbX, rs.Qb, rs.X); # Initialise dot product cache.
        objvalue = objective(rs); 
        objvalue < best_objvalue && (best_objvalue = objvalue;);  
        init_statuses!(rs);
        verbose && println("Initial objective on initialisation $k is: $objvalue");
        for i in Base.OneTo(max_sweeps)
            shuffle!(rng, rs.order); # Process stands in random order.
            changed = false; # Whether a stand was updated during the sweep.
            for stand in rs.order
                @inbounds status = rs.status_lu[stand]; # Get current status of stand.

                # Get all plausible objective function increments in response to changes to current stand, 
                # given current solution.
                fill_increments!(rs, stand, status); 
                # Get objective function increment for stand given descent strategy. 
                incr, ind = find_increment(rs.ds, status, rs.increments, rng);
                newstatus = ind - 1; # Get new status of stand.
                if newstatus != status 
                    changed = true;
                    update_solution!(rs, stand, status, newstatus); 
                    @inbounds rs.status_lu[stand] = newstatus; # Update current status of stand.
                    update_dotprod_cache!(rs, stand, status, newstatus); # Update QbX accordingly. 
                    objvalue += incr; # Modify objective in response to changed status of stand.
                    if objvalue < best_objvalue 
                        best_objvalue = objvalue;  
                        rs.Xopt .= rs.X;
                    end
                end
                # Debugging:
                #@assert is_feasible(rs.X) "X not feasible"
                #@assert incr <= 0.0 "increrror"
                #@assert isapprox(objvalue, objective_value(rs)) "objerror" 
                #@assert isapprox(build_sol(rs.status_lu, size(rs.X, 2)), rs.X)  "solerror" 
            end
            verbose && println("Sweep $i of initialisation $k found objective $objvalue");
            if !changed # early exit if sweep did not do anything (algorithm converged)
                verbose && println("Exiting early, no improvement.");
                break
            end
        end
    end
    return best_objvalue;
end

function ccall_random_sweep!(rs::RandomSweepStorage, seed::Integer; 
                             max_sweeps::Integer = 100000, inits::Integer = 100)

    out = Ref(0.0);
    inits = Cuint(inits);
    seed = Culong(seed);
    max_sweeps = Cuint(max_sweeps);
    rsC = RandomSweepStorageC(rs);
    ret = ccall((:random_sweep_w_seed, "libpov"), Cint, 
          (Ref{Cdouble}, Ref{RandomSweepStorageC}, Cuint, Cuint, Culong),
          out, rsC, max_sweeps, inits, seed);
    if ret != 0
        throw(ArgumentError("random sweep algorithm did not converge."));
    end
    out[];
end

@inline function find_increment(ds::DescentStrategy{:smallest_negative}, status::Integer, 
                                increments::AbstractVector{Float64}, 
                                rng::AbstractRNG = Random.GLOBAL_RNG)
    findmin(increments);
end
@inline function find_increment(ds::DescentStrategy{:random_negative}, status::Integer, 
                                increments::AbstractVector{Float64}, 
                                rng::AbstractRNG = Random.GLOBAL_RNG)
    (incr, ind) = rand_among_negative(increments, ds.tmp, rng);
    ind < 0 && (return (0.0, status + 1););
    (incr, ind);
end
@inline function find_increment(ds::DescentStrategy{:greatest_negative}, status::Integer,
                                increments::AbstractVector{Float64}, 
                                rng::AbstractRNG = Random.GLOBAL_RNG)
    (incr, ind) = greatest_negative(increments);  
    ind < 0 && (return (0.0, status + 1););
    (incr, ind);
end

"""
    `greatest_negative(x)`

Find the value and index of the greatest strictly negative value in the vector `x`.
If `x` contains no negative values, the tuple `(-Inf, -1)` is returned.
"""
function greatest_negative(x::AbstractVector{<: Real})
    greatest = -Inf;
    greatest_i = -1;
    for i in eachindex(x)
        if @inbounds x[i] < 0.0 && x[i] > greatest
            @inbounds greatest = x[i];
            greatest_i = i;
        end
    end
    greatest, greatest_i
end

"""
    `rand_among_negative(x, tmp, rng)`

Draw a random negative value and its associated index from the vector `x`.
If `x` contains no negative values, the tuple `(-Inf, - 1)` is returned.
`tmp` should be an integer vector of at least the length of `x`.
`rng` is the random number generator used in the sampling.
"""
function rand_among_negative(x::AbstractVector{<: Real}, tmp::AbstractVector{Int}, 
                             rng::AbstractRNG = Random.GLOBAL_RNG)
    n = 0;
    for i in eachindex(x)
        if @inbounds x[i] < zero(eltype(x))
            n += 1;
            @inbounds tmp[n] = i;  
        end
    end
    n == 0 && (return (-Inf, -1));
    ind = rand(rng, view(tmp, 1:n));
    @inbounds (x[ind], ind);
end
