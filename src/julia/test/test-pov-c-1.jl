include("../lib/random-sweep.jl");
include("../lib/model.jl");
include("../lib/data-functions.jl");
using Random
using Test

function ccall_build_problem!(rs::RandomSweepStorage, 
                              vp::VolumePosterior,
                              stda::StandData)
    rsC = RandomSweepStorageC(rs);
    vpC = VolumePosteriorC(vp);
    stdaC = StandDataC(stda);

    ccall((:build_problem, "libpov"), Cvoid,
          (Ref{RandomSweepStorageC}, Ref{VolumePosteriorC}, Ref{Cdouble}),
          rsC, vpC, stda.demands);
end

function ccall_compute_posterior!(vp::VolumePosterior, xIC::Vector{Cuint}, stda::StandData)
    vpC = VolumePosteriorC(vp);
    stdaC = StandDataC(stda);
    ccall((:normal_model_compute_posterior, "libpov"), Cvoid,
          (Ref{VolumePosteriorC}, Ref{Cuint}, Ref{StandDataC}),
          vpC, xIC, stdaC);
end

xI_choice = 3;
xI = repeat([xI_choice], 200);
datarng = Xoshiro(20042023);
muplus, sigma2plus, demands = generate_normal_test_data(xI, datarng);
prob = build_problem(muplus, sigma2plus, demands);
nS = size(muplus, 1); 
nT = size(demands, 1); 
nA = size(muplus, 2);
rs = RandomSweepStorage(prob.Q, prob.c, prob.r);
rs.C .= 0.0;
rs.Qb .= 0.0;
rs2 = RandomSweepStorage(prob.Q, prob.c, prob.r); 
vp = VolumePosterior(copy(muplus), copy(sigma2plus));
muprior, sdprior, demands = get_prior_data_and_demands(); 
sigma2prior = sdprior .* sdprior;
sigmameas = get_meas_sds();
sigma2meas = sigmameas .* sigmameas;

stda = StandData(repeat([""], nS), copy(muprior), copy(sigma2prior), copy(sigma2meas), copy(demands),
                 [0.0, 0.0, 0.0], repeat([0.0], nS), zeros(nS, 2));

ccall_build_problem!(rs, vp, stda);
@testset "check that C implementation builds problem right" begin
    @test isapprox(rs.C, rs2.C)
    @test isapprox(rs.Qb, rs2.Qb)
end

ys = zeros(nS, nA);
datarng = Xoshiro(20042023);
simulate_normal_stand_data!(ys, xI, muprior, sdprior, get_meas_sds(), datarng); 
vp.y .= ys;
vp.muplus .= 0.0;
vp.sigma2plus .= 0.0;
xIC = convert.(Cuint, xI) .- Cuint(1); # For 0 based indexing.

ccall_compute_posterior!(vp, xIC, stda);
@testset "check that C implementation builds normal posterior right given data" begin
    @test isapprox(muplus, vp.muplus)
    @test isapprox(sigma2plus, vp.sigma2plus)
end







