module DataLoader
using JSON
include("json2type.jl")

function LoadForcing(ULC_lat::Int, ULC_lon::Int; dirpath="input")
        #load forcing data -------------------------------------------------------------
        forcing_path = "$(dirpath)/FORCING_JSONfiles/FORCING_ULC_" * string(ULC_lon) * "_" * string(ULC_lat) * ".json"
        jsontxt = open(forcing_path,"r") do file
                dicttxt = read(file, String)
        end
        FORCING = json2type.typenarrow!(JSON.parse(jsontxt))
        #in the case nan values exisit in Ptot
        for i = 1:length(FORCING["data"]["Ptot"])
                if FORCING["data"]["Ptot"][i]==nothing
                        FORCING["data"]["Ptot"][i] = 0.0
                end
        end
        FORCING["data"]["Ptot"] = convert(Array{Float64,1},FORCING["data"]["Ptot"])
        #snow fall in winter (snow fall is in m/day - snow water equivalent - SWE!!!)
        FORCING["data"]["snowfall"] = (FORCING["data"]["Tair"].<=0.0).*FORCING["data"]["Ptot"]/1000.0;
        return FORCING
end

function LoadParams(ULC_lat::Int, ULC_lon::Int; dirpath="input")
        #load parameter ----------------------------------------------------------------
        para_path = "$(dirpath)/PARA_JSONfiles/PARA_ULC_" * string(ULC_lon) * "_" * string(ULC_lat) * ".json"
        jsontxt = open(para_path,"r") do file
            dicttxt = read(file, String)
        end
        PARA = json2type.typenarrow!(JSON.parse(jsontxt))
end
end
