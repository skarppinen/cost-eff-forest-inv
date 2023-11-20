include("../lib/random-sweep.jl");
include("../lib/model.jl");
using Test

seed = 24032023;
Random.seed!(seed);
xI = repeat([2], 200);
muplus, sigma2plus, demands = generate_normal_test_data(xI, Random.GLOBAL_RNG);
prob = build_problem(muplus, sigma2plus, demands);
nS = size(muplus, 1); 
nT = size(demands, 1);
X = zeros(Float64, nS, nT);
Xcopy = copy(X); Xcopy2 = copy(X);
random_feasible_sol!(X);
cmat = rand(Float64, nS, nT);
#cmat = -2.0 * muplus * transpose(demands);
Qblock = prob.Q;
r = prob.r;

Random.seed!();
@testset "Test increment formulas associated with the random sweep algorithm" begin
    ## Test changing 1 to 0:
    stand = 2; time = 5;
    @assert X[stand, time] == 1 "something wrong with setup"
    cur_with_1 = objective(X, cmat, Qblock, r);
    Xcopy .= X;
    Xcopy[stand, time] = 0.0;
    cur_with_0 = objective(Xcopy, cmat, Qblock, r);
    increment_1_to_0 = (-cmat[stand, time] - dot(Qblock[:, stand], X[:, time]) + 0.5 * Qblock[stand, stand]);
    cur_with_0_up = cur_with_1 + increment_1_to_0; 
    @testset "Test 1 => 0 rule" begin
        @test isapprox(cur_with_0_up, cur_with_0)
        @test isapprox(cur_with_1, cur_with_0_up - increment_1_to_0);
        @test isapprox(increment_1_to_0, _get_increment(stand, time, 0, cmat, Qblock * X, Qblock[stand, stand]));
    end

    ## Test changing 0 to 1 rule:
    stand = 193;
    @assert isapprox(sum(X[stand, :]), 0.0) "something wrong with setup"
    cur_with_0 = objective(X, cmat, Qblock, r);
    Xcopy .= X;
    time = rand(1:nT);
    Xcopy[stand, time] = 1.0; 
    cur_with_1 = objective(Xcopy, cmat, Qblock, r);

    increment_0_to_1 = cmat[stand, time] + dot(Qblock[:, stand], X[:, time]) + 0.5 * Qblock[stand, stand];
    cur_with_1_up = cur_with_0 + increment_0_to_1; 
    @testset "Test 0 => 1 rule" begin
        @test isapprox(cur_with_1_up, cur_with_1); 
        @test isapprox(increment_0_to_1, _get_increment(stand, 0, time, cmat, Qblock * X, Qblock[stand, stand]));
    end

    ## Test changing 1 to 1 rule:
    stand = 2; time = 5;
    @assert X[stand, time] == 1 "something wrong with setup"
    cur_with_1_init = objective(X, cmat, Qblock, r);
    Xcopy .= X;
    times = collect(1:nT); splice!(times, time);
    totime = rand(times); 
    Xcopy[stand, time] = 0.0; 
    Xcopy2 .= Xcopy;
    Xcopy[stand, totime] = 1.0;
    cur_with_1_changed = objective(Xcopy, cmat, Qblock, r);
    increment_1_to_0 = (-cmat[stand, time] - dot(Qblock[:, stand], X[:, time]) + 0.5 * Qblock[stand, stand]);
    increment_0_to_1 = -(-cmat[stand, totime] - dot(Qblock[:, stand], X[:, totime]) - 0.5 * Qblock[stand, stand]);
    increment_1_to_1 = increment_1_to_0 + increment_0_to_1;
    increment_1_to_1_alt = -cmat[stand, time]  + cmat[stand, totime] - dot(Qblock[:, stand], X[:, time]) +
    dot(Qblock[:, stand], X[:, totime]) + Qblock[stand, stand];
    @testset "Test 1 => 1 rule" begin
        @test isapprox(cur_with_1_changed, cur_with_1_init + increment_1_to_1)
        @test isapprox(cur_with_1_changed, cur_with_1_init + increment_1_to_1_alt)
        @test isapprox(increment_1_to_1, _get_increment(stand, time, totime, cmat, Qblock * X, Qblock[stand, stand]));
        @test isapprox(increment_1_to_1_alt, _get_increment(stand, time, totime, cmat, Qblock * X, Qblock[stand, stand]));
    end
end

@testset "Test RandomSweepStorage" begin
    Random.seed!();
    rs = RandomSweepStorage(Qblock, cmat, r);
    random_feasible_sol!(rs);
    mul!(rs.QbX, rs.Qb, rs.X); # Initialise QbX.
    init_statuses!(rs);
    objvalue = objective(rs);
    for stand in 1:size(rs.X, 1)
        status = rs.status_lu[stand];
        row = rs.X[stand, :];
        fill_increments!(rs, stand, status);
        if status == 0
            @test isapprox(rs.increments[1], 0.0) 
        else
            @test isapprox(rs.increments[status + 1], 0.0); 
        end
        for j in eachindex(rs.increments)
            update_solution!(rs, stand, status, j - 1);  
            @test is_feasible(rs.X)
            if j == 1
                @test isapprox(sum(rs.X[stand, :]), 0.0)
            else
                @test isapprox(rs.X[stand, j - 1], 1.0)
            end
            @test isapprox(objective(rs), objvalue + rs.increments[j]); 
            rs.X[stand, :] .= row;
        end
    end
end


