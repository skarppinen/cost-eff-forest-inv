using JLD2

function get_stand_data(filename::AbstractString)
    stand_ids, muprior, varprior, measvars, demands, costs, areas, coordinates = 
    jldopen(joinpath(@__DIR__, "../../../data/$filename"), "r") do file
        file["stand_ids"],
        file["muprior"], 
        file["sdprior"] .^ 2.0, 
        file["meas_sds"] .^ 2.0,
        file["demands"],
        file["costs"],
        file["areas"],
        file["coordinates"];
    end;
    (; stand_ids, muprior, varprior, measvars, demands, costs, areas, coordinates);
end

function get_prior_data_and_demands()
    stda = get_stand_data("test-data-2.jld2");
    stda.muprior, sqrt.(stda.varprior), stda.demands;
end

function get_meas_sds()
    stda = get_stand_data("test-data-2.jld2");
    sqrt.(stda.measvars); 
end
