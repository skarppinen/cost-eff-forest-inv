include(joinpath(@__DIR__, "../lib", "stand-ordering.jl"));
using Test

@testset "test that `ordered_choices` works" begin
    @test ordered_choices(1, 3, 1) == [Int[1], Int[2], Int[3]];
    @test ordered_choices(1, 2, 2) == [Int[1, 1], Int[1, 2], Int[2, 2]];

    @test ordered_choices(1, 2, 0) == Vector{Int}[];
    @test ordered_choices(1, 2, -1) == Vector{Int}[];
    @test ordered_choices(0, 0, -1) == Vector{Int}[]; 
    @test ordered_choices(2, 1, -1) == Vector{Int}[];
    @test ordered_choices(-2, -56, -1) == Vector{Int}[]; 

    @test_throws AssertionError ordered_choices(2, 1, 10);

    for upper in [2, 3, 4, 5]
        for n in 1:20
            out = ordered_choices(1, upper, n);
            # Test against analytical formula for number of possibilities.
            @test length(out) == binomial(n + upper - 1, n);

            # All are sorted and no duplicates.
            @test length(out) == length(unique(out))
            @test all(issorted.(out))
        end
    end
end
