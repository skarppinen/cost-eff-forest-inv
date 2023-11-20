using Test
include(joinpath(@__DIR__, "../lib/random-sweep.jl"));

@testset "test `greatest_negative`" begin
    # Case of zero length.
    @test greatest_negative(Float64[]) == (-Inf, -1);    

    # No negatives.
    @test greatest_negative([0.0, 1.0]) == (-Inf, -1);
    @test greatest_negative([5.0, 1.0]) == (-Inf, -1);

    # Has negatives.
    @test greatest_negative([-1.0, 1.0]) == (-1.0, 1);
    @test greatest_negative([-1.0, -1.0, 1.0]) == (-1.0, 1);
    @test greatest_negative([-1.0, -1.0, 1.0, -0.5]) == (-0.5, 4);
    @test greatest_negative([1.0, -2.0, -1.0, 0.0]) == (-1.0, 3);
end

@testset "test `rand_among_negative`" begin
    tmp = zeros(Int, 4);
    @test rand_among_negative(Float64[], tmp) == (-Inf, -1);
    @test rand_among_negative([1.0, 1.0, 0.0, 0.0], tmp) == (-Inf, -1); 
    @test rand_among_negative([1.0, -1.0, 0.0, 0.0], tmp) == (-1.0, 2); 
    @test rand_among_negative([1.0, -1.0, -1.0, 0.0], tmp)[1] == -1.0; 
    @test rand_among_negative([1.0, -1.0, -2.0, 0.0], tmp)[1] in (-2.0, -1.0); 
    @test rand_among_negative([1.0, -1.0, -2.0, 0.0], tmp)[2] in (2, 3); 
end
