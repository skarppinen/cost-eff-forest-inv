include("../lib/pov.jl");
include("../lib/r-data-functions.jl");
using Test

stda = StandData(id = 3);
nS, nT, nA, nI = size(stda);
vp = VolumePosterior(nS, nA);
rs = RandomSweepStorage(nS, nT);

seeds = Culong.(1:10);
xI = Cuint.(repeat([1], nS));
ninits = 50;

for model in [1]
    # PoV with 10 samples.
    rs.r[] = 0.0;
    pov_all = ccall_PoV!(xI, rs, vp, stda, seeds; ninits = ninits, model = model);

    # PoV with r modified on Julia side. 
    rs.r[] = sum(stda.demands .* stda.demands); 
    pov_mod_r = ccall_PoV!(xI, rs, vp, stda, seeds; ninits = ninits, model = model);

    # Computing separately.
    pov_separate = zeros(length(seeds));
    for i in eachindex(seeds)
        oneseed = [seeds[i]];
        pov_separate[i] = ccall_PoV!(xI, rs, vp, stda, oneseed; ninits = ninits, model = model);
    end

    # Test that gives the same.
    msg = """
    Model $model, test that separately computing PoV for each sample and taking mean is same as \
    computing PoV approximation for all samples with single call.
    """
    @testset "$msg" begin
        @test isapprox(pov_all, mean(pov_separate));
    end

    @testset "Model $model, check that modifying r in Julia has no effect to what C implementation computes" begin
        @test isapprox(pov_all, pov_mod_r);
    end
end

function ccall_random_posterior!(xI::Vector{Cuint}, vp::VolumePosterior, 
        stda::StandData; seed::Integer, model::Integer) 

    @assert model in (1, 2)
    model = convert(Cint, model) - Cint(1);
    seed = convert(Culong, seed);
    vpC = VolumePosteriorC(vp);
    stdaC = StandDataC(stda);
    ccall((:random_posterior_w_seed, "libpov"), Cvoid,
          (Ref{VolumePosteriorC}, Ref{Cuint}, Ref{StandDataC}, Cint, Culong),
          vpC, xI, stdaC, model, seed);
end

#@testset "different models yield different random posteriors with same seed" begin
#
#vp_normal = deepcopy(vp);
##vp_lognormal = deepcopy(vp);
#ccall_random_posterior!(xI, vp_normal, stda; seed = seeds[1], model = 1);
##ccall_random_posterior!(xI, vp_lognormal, stda; seed = seeds[1], model = 2);
#
#@test !isapprox(mean(vp_normal.muplus), mean(vp_lognormal.muplus))
##@test !isapprox(mean(vp_normal.sigma2plus), mean(vp_lognormal.sigma2plus))
#
#end


