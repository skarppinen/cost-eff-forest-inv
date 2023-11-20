include(joinpath(@__DIR__, "data-functions.jl"));
using Random, Distributions, LinearAlgebra



struct StandData
    stand_ids::Vector{String}
    muprior::Matrix{Float64}
    sigma2prior::Matrix{Float64}
    sigma2meas::Array{Float64, 3}    
    demands::Matrix{Float64}
    costs::Vector{Float64}
    areas::Vector{Float64}
    coordinates::Matrix{Float64}
    function StandData(stand_ids::AbstractVector{String},
                       muprior::AbstractMatrix{Float64}, sigma2prior::AbstractMatrix{Float64},
                       sigma2meas::AbstractArray{Float64}, demands::AbstractMatrix{Float64},
                       costs::AbstractVector{Float64}, areas::AbstractVector{Float64},
                       coordinates::AbstractMatrix{Float64})
        @assert size(muprior) == size(sigma2prior) "`muprior` and `sigma2prior` sizes must match.";
        (nS, nA) = size(muprior);
        nT = size(demands, 1);
        @assert nA == size(demands, 2); "invalid second dimension size in `demands`.";
        nI = size(sigma2meas, 1);
        @assert nS == size(sigma2meas, 2) "invalid second dimension size in `sigma2meas`";
        @assert nA == size(sigma2meas, 3) "invalid third dimension size in `sigma2meas`";
        @assert all(sigma2meas .> 0.0) "invalid negative measurement variances";
        @assert all(sigma2prior .> 0.0) "invalid negative prior variances";
        @assert all(demands .> 0.0) "invalid negative demands";
        @assert length(costs) == nI "invalid length of `costs`";
        @assert size(coordinates) == (nS, 2) "invalid size of `coordinates`";
        @assert length(areas) == nS "invalid length of `areas`";
        @assert length(stand_ids) == nS "invalid length of `stand_ids`";

        new(stand_ids, muprior, sigma2prior, sigma2meas, demands, costs, areas, coordinates);
    end
end

function StandData(filename::AbstractString)
    stda = get_stand_data(filename); 
    StandData(stda.stand_ids, stda.muprior, stda.varprior, stda.measvars, 
              stda.demands, stda.costs, stda.areas, stda.coordinates);
end

function StandData(; id::Integer)
    filename = "test-data-$id.jld2";
    StandData(filename); 
end

import Base.size;
"""
Return the dimensions of the `StandData` object,
i.e (# of stands, # time points, # assortments, # inventory methods).
"""
function size(stda::StandData)
    nS, nA = size(stda.muprior); 
    nI = size(stda.sigma2meas, 1);
    nT = size(stda.demands, 1);
    nS, nT, nA, nI;
end

struct StandDataC
    nS::Cuint
    nT::Cuint
    nA::Cuint
    nI::Cuint
    muprior::Ptr{Cdouble}
    sigma2prior::Ptr{Cdouble}
    sigma2meas::Ptr{Cdouble}
    demands::Ptr{Cdouble}
    function StandDataC(stda::StandData) 
        (nI, nS, nA) = size(stda.sigma2meas);
        nT = size(stda.demands, 1);
        new(Cuint(nS), Cuint(nT), Cuint(nA), Cuint(nI), 
            pointer(stda.muprior), pointer(stda.sigma2prior),
            pointer(stda.sigma2meas), pointer(stda.demands)); 
    end
end

struct VolumePosterior
    muplus::Matrix{Float64}
    sigma2plus::Matrix{Float64}
    y::Matrix{Float64}
    function VolumePosterior(muplus::AbstractMatrix{Float64}, 
                             sigma2plus::AbstractMatrix{Float64},
                             y::AbstractMatrix{Float64} = zeros(Float64, size(muplus)))
        @assert size(muplus) == size(sigma2plus) == size(y) "sizes must match.";
        @assert all(sigma2plus .> 0.0) "negative variances"; 
        new(muplus, sigma2plus, y);
    end
end
function VolumePosterior(nS::Integer, nA::Integer)
    VolumePosterior(zeros(Float64, nS, nA), ones(Float64, nS, nA));
end
function VolumePosterior(stda::StandData) 
    nS, nT, nA, nI = size(stda);
    VolumePosterior(nS, nA);
end

struct VolumePosteriorC
    nS::Cuint
    nA::Cuint
    muplus::Ptr{Cdouble}
    sigma2plus::Ptr{Cdouble}
    y::Ptr{Cdouble}
    function VolumePosteriorC(vp::VolumePosterior)
        (nS, nA) = size(vp.muplus);
        new(Cuint(nS), Cuint(nA), pointer(vp.muplus), pointer(vp.sigma2plus),
            pointer(vp.y));
    end
end
struct SeedVectorC
    arr::Ptr{Culong}
    n::Cuint
end

"""
Compute the posterior mean and sd of the normal distribution `m | y`, when 
`m ~ N(mean0, sd0)` and `y ~ N(m, meas_sd)`.
"""
function normal_posterior_mean_sd(y, mean0, sd0, meas_sd)
	var0 = sd0 * sd0;
	meas_var = meas_sd * meas_sd;
	post_var = inv(inv(var0) + inv(meas_var)); 
	post_mean = post_var * (mean0 * inv(var0) + y * inv(meas_var));
	post_mean, sqrt(post_var);
end


function simulate_normal_stand_data!(ys, xI, mu0s, sigma0s, meas_sds, rng) 
    (nS, nA) = size(mu0s);
    @assert size(ys) == size(mu0s)
    @assert length(xI) == nS
    @assert size(mu0s) == size(sigma0s)
    @assert size(meas_sds) == (3, nS, nA)
    for a in 1:nA
        for (s, xi) in enumerate(xI)
            mu = mu0s[s, a];
            sd = sqrt(sigma0s[s, a] * sigma0s[s, a] + 
                      meas_sds[xi, s, a] * meas_sds[xi, s, a]);
            ys[s, a] = rand(rng, Normal(mu, sd));
        end
    end
    nothing;
end

"""
Function returns posterior means, variances and assortment demands given an inventory
decision `xI`, assumed to be a vector of length 200, with each element taking a value in
{1, 2, 3}.
The posterior means and variances are computed by sampling from the marginal distribution
of the observations with the specified `rng`, and then using the normal update formulas to
compute the posterior mean and variance.
"""
function generate_normal_test_data(xI::AbstractVector{<: Integer}, rng::AbstractRNG = Random.GLOBAL_RNG) 
    mu0s, sigma0s, demands = get_prior_data_and_demands();
    nS, nA = size(mu0s);
    @assert length(xI) == nS "`xI` must equal $nS in length";
    nT = size(demands, 1);
    meas_sds = get_meas_sds();

    # Simulate from p(y | x^{(I)}).
    ys = zeros(nS, nA);
    simulate_normal_stand_data!(ys, xI, mu0s, sigma0s, meas_sds, rng); 

    # Compute p(V | y).
    muplus = zeros(nS, nA);
    sigmaplus = zeros(nS, nA);
    for a in 1:nA
        for (s, xi) in enumerate(xI)
            muplus[s, a], sigmaplus[s, a] = normal_posterior_mean_sd(ys[s, a], mu0s[s, a], 
											 sigma0s[s, a], meas_sds[xi, s, a]); 
        end
    end
    sigma2plus = sigmaplus .^ 2.0;
    muplus, sigma2plus, demands;
end

function generate_normal_test_data(stda::StandData, xI::AbstractVector{<: Integer},
        rng::AbstractRNG = Random.GLOBAL_RNG)
    (nS, nT, nA, nI) = size(stda);
    @assert length(xI) == nS "`xI` must equal $nS in length (is $(length(xI)))";
    
    # Simulate from p(y | x^{(I)}).
    ys = zeros(nS, nA);
    sigma0s = sqrt.(stda.sigma2prior); 
    meas_sds = sqrt.(stda.sigma2meas);
    simulate_normal_stand_data!(ys, xI, stda.muprior, sigma0s, meas_sds, rng); 

    # Compute p(V | y).
    muplus = zeros(nS, nA);
    sigmaplus = zeros(nS, nA);
    for a in 1:nA
        for (s, xi) in enumerate(xI)
            muplus[s, a], sigmaplus[s, a] = normal_posterior_mean_sd(ys[s, a], stda.muprior[s, a], 
											 sigma0s[s, a], meas_sds[xi, s, a]); 
        end
    end
    sigma2plus = sigmaplus .^ 2.0;
    muplus, sigma2plus, ys;
end

function random_posterior(stda::StandData, xI::AbstractVector{<: Integer}, 
        rng::AbstractRNG = Random.GLOBAL_RNG; model::Integer = 1)
    @assert model in (1, 2) "invalid model choice $model";  
    if model == 1
        # Normal model. 
        muplus, sigma2plus, ys = generate_normal_test_data(stda, xI, rng); 
    else
        throw(ArgumentError("not implemented"));
    end
    VolumePosterior(muplus, sigma2plus, ys);
end

"""
Function constructs matrices and vectors of the central quadratic programming problem
min 0.5 * x'Qx + c'x + r, wrt. Ax = 1.
These quantities are all functions of the inputs.
Q is a block diagonal matrix, whose repeating block is given by the output `Qblock`,
(or as vector `vQblock`, in column major order)
"""
function build_iqp_system(muplus::AbstractMatrix{<: AbstractFloat}, 
                          sigma2plus::AbstractMatrix{<: AbstractFloat},
                          demands::AbstractMatrix{<: AbstractFloat})
    if size(muplus) != size(sigma2plus) 
        msg = "dimensions of `muplus` and `sigma2plus` must match.";
        throw(ArgumentError(msg));
    end
    nS, nA = size(muplus);
    if size(demands, 2) != nA
        msg = "invalid second dimension of `demands`, should equal nA = $nA";
        throw(ArgumentError(msg));
    end
    nT = size(demands, 1);
    musum = map(1:size(muplus, 2)) do i
        v = view(muplus, :, i);
        v * transpose(v);
    end |> (x -> reduce(+, x));
    Sigma2sum = map(1:size(sigma2plus, 2)) do i
        diagm(view(sigma2plus, :, i));
    end |> (x -> reduce(+ , x));
    Qblock = 2.0 * (Sigma2sum + musum);
    c = -2.0 * vec(muplus * transpose(demands));
    demandsqsum = sum(demands .* demands);

    Q = zeros(nS * nT, nS * nT);
    for i in 1:nT
        l = (i - 1) * nS + 1;
        u = i * nS;
        Q[l:u, l:u] .= Qblock;
    end
    A = zeros(nS, nS * nT);
    for i in 1:nT
        A[1:nS, ((i - 1) * nS + 1):(i * nS)] .= I[1:nS, 1:nS];
    end
    b = repeat([1.0], nS);
    (Qblock = Qblock, vQblock = vec(Qblock), c = c, nS = nS, nT = nT, nA = nA,
     r = demandsqsum, Q = Q, A = A, b = b);
end

"""
Function constructs the ingredients Q, c and r of the inner 
minimisation problem based on posterior means and variances and the demands. 

The inner minimisation problem is of the form:
min 0.5\\sum_{t = 1}^{nT} x_t' Q x_t + \\sum_{t=1}^{nT} c_t'x_t + r, 
where `c_t` is the `t`th column of the output matrix `c`.
"""
function build_problem(muplus::AbstractMatrix{<: AbstractFloat}, 
                       sigma2plus::AbstractMatrix{<: AbstractFloat},
                       demands::AbstractMatrix{<: AbstractFloat})
    if size(muplus) != size(sigma2plus) 
        msg = "dimensions of `muplus` and `sigma2plus` must match.";
        throw(ArgumentError(msg));
    end
    @assert all(sigma2plus .> 0.0) "expected non-positive variances!";
    @assert all(demands .>= 0.0) "expected non-negative demands!";

    nS, nA = size(muplus);
    if size(demands, 2) != nA
        msg = "invalid second dimension of `demands`, should equal nA = $nA";
        throw(ArgumentError(msg));
    end
    nT = size(demands, 1);
    @assert nT > 0 "invalid demand matrix dimensions (nT, nA) = ($nT, $nA)";

    musum = map(1:size(muplus, 2)) do i
        v = view(muplus, :, i);
        v * transpose(v);
    end |> (x -> reduce(+, x));
    Sigma2sum = map(1:size(sigma2plus, 2)) do i
        diagm(view(sigma2plus, :, i));
    end |> (x -> reduce(+ , x));
    Q = 2.0 * (Sigma2sum + musum);
    c = -2.0 * muplus * transpose(demands);
    demandsqsum = sum(demands .* demands);
    (Q = Q, c = c, r = demandsqsum);
end

function build_problem(vp::VolumePosterior, demands::AbstractMatrix{<: AbstractFloat})
    build_problem(vp.muplus, vp.sigma2plus, demands);
end
