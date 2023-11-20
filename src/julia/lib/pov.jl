include("model.jl");
include("random-sweep.jl");

function ccall_PoV!(xI::Vector{Cuint},
                    rs::RandomSweepStorage, 
                    vp::VolumePosterior,
                    stda::StandData,
                    seeds::Vector{Culong};
                    maxsweeps::Integer = Cuint(9999),
                    ninits::Integer,
                    model::Integer)
    @assert model in (1, 2, 3) "invalid model choice";
    model = convert(Cint, model) - Cint(1);
    rsC = RandomSweepStorageC(rs);
    vpC = VolumePosteriorC(vp);
    stdaC = StandDataC(stda);
    svC = SeedVectorC(pointer(seeds), Cuint(length(seeds)));
    maxsweeps = convert(Cuint, maxsweeps);
    ninits = convert(Cuint, ninits);
    ccall((:PoV, "libpov"), Cdouble,
          (Ref{Cuint}, Ref{RandomSweepStorageC}, Ref{VolumePosteriorC}, 
           Ref{StandDataC}, Ref{SeedVectorC}, Cuint, Cuint, Cint),
          xI, rsC, vpC, stdaC, svC, maxsweeps, ninits, model); 
end

function ccall_PoV_w_algseed!(xI::Vector{Cuint},
                    rs::RandomSweepStorage, 
                    vp::VolumePosterior,
                    stda::StandData,
                    seeds::Vector{Culong};
                    algseed::Culong,
                    maxsweeps::Integer = Cuint(9999),
                    ninits::Integer)
    
    rsC = RandomSweepStorageC(rs);
    vpC = VolumePosteriorC(vp);
    stdaC = StandDataC(stda);
    svC = SeedVectorC(pointer(seeds), Cuint(length(seeds)));
    maxsweeps = convert(Cuint, maxsweeps);
    ninits = convert(Cuint, ninits);
    ccall((:PoV_w_algseed, "libpov"), Cdouble,
          (Ref{Cuint}, Ref{RandomSweepStorageC}, Ref{VolumePosteriorC}, 
           Ref{StandDataC}, Ref{SeedVectorC}, Culong, Cuint, Cuint),
          xI, rsC, vpC, stdaC, svC, algseed, maxsweeps, ninits);
end
