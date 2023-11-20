include("model.jl");

function _do_ordered_choices(first::Integer, upper::Integer, n::Integer)
    if n <= 0
        return Vector{typeof(first)}[];
    end
    if first == upper
        return [repeat([first], n)];
    end
    out = Vector{Vector{typeof(first)}}(undef, 0);
    for i in first:upper
        firstvec = [i];
        rest = _do_ordered_choices(i, upper, n - 1);
        if isempty(rest)
            push!(out, firstvec);
        else
            for r in rest
                push!(out, vcat(firstvec, r));
            end
        end
    end
    out;
end

"""
Return a vector containing all possible vectors of length `n` that are ordered
and consist of elements in the set { `lower`, `lower + 1`, ..., `upper` }. 
"""
function ordered_choices(lower::Integer, upper::Integer, n::Integer)
    if n < 0
        return Vector{typeof(lower)}[];
    end
    @assert lower > 0
    @assert upper >= lower
    _do_ordered_choices(lower, upper, n);
end

""" 
Return the smaller feasible set, that is, all possible inventory decision vectors that 
satisfy the property that for any two stands with total prior assortment variances `x` and 
`y`, with `x < y`, a more expensive inventory is never carried out for the stand with 
the total prior assortment variance `x`. 
The output vectors are in the order of the stand in the input data.
"""
function feasible_set_reduction(stda::StandData) 
    (nS, nT, nA, nI) = size(stda);
    vartotal = mapslices(sum, stda.sigma2prior, dims = 2)[:, 1];

    # varorder = "Stand indices such that `vartotal[varorder]` is increasing".
    # i.e first index contains stand index that has smallest total variance,
    # second index contains stand index that has second smallest total variance, and so on..
    varorder = sortperm(vartotal)
    
    # Get all possible ordered inventory decision vectors.
    ordered_decisions = ordered_choices(1, nI, nS); 
    
    # Construct all inventory decisions to consider, that take into account order of total variance.
    inventory_decisions = [zeros(Int, nS) for i in 1:length(ordered_decisions)]; 
    for i in 1:length(inventory_decisions)
        inventory_decisions[i][varorder] .= ordered_decisions[i];
    end
    inventory_decisions;
end
