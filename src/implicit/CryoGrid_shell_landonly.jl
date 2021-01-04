#Start script for CryoGrid
#Created on Fri Jul 17 19:06:49 2020
#author: Moritz Langer (moritz.langer@awi.de)
#runs with Julia Version 0.6.2 (2017-12-13 18:08 UTC)
#https://julialang.org/downloads/oldreleases/
#note that the script might require to install additional packages e.g. Pkg.add("Interpolations")
# ========================================================================
#ArcticLakes shell script to run the CryoGridImplicit Model
using MAT
using Glob
using JLD
using JSON
using Dates

include("CryoGridImplicit.jl")
include("matlab.jl")
include("ArcticLakesSetup.jl")
include("CryoGridTyps.jl")
include("json2type.jl")

# ==============================================================================
#Set site coordinates as upper left corner (grid: 0-360E, 55-90N)
ULC_lat = 72
ULC_lon = 126
# ==============================================================================

#load forcing data -------------------------------------------------------------
forcing_path = "input/FORCING_JSONfiles/FORCING_ULC_" * string(ULC_lon) * "_" * string(ULC_lat) * ".json"
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

#load parameter ----------------------------------------------------------------
para_path = "input/PARA_JSONfiles/PARA_ULC_" * string(ULC_lon) * "_" * string(ULC_lat) * ".json"
jsontxt = open(para_path,"r") do file
    dicttxt = read(file, String)
end
PARA = json2type.typenarrow!(JSON.parse(jsontxt))

GRID = PARA["GRID"] #soil grid with depth
PARA = PARA["PARA"] #soil and snow parameters

#Set to land only (no lakes)
tile_number = 1 #use only one tile for reference run
LAKESTAT = Dict()
LAKESTAT["LandMaskArea"]=[1.0]
LAKESTAT["LakeDistance"]=[1.0]
LAKESTAT["LakeNumber"]=[0.0]
LAKESTAT["LakeDepth"]=[0.0]
LAKESTAT["SnowDepth"]=[NaN]
LAKESTAT["VorArea"]=[1.0]
LAKESTAT["LandMaskArea"]=[1.0]
LAKESTAT["LakeArea"]=[0.0]
LAKESTAT["LakePerimeter"]=[NaN]

#BUILD input structures --------------------------------------------------------
FOR, PAR, TEM, GRI, STRA, STAT, OUT = ArcticLakesSetup.SetUpInputStructs(FORCING, GRID, PARA, LAKESTAT, tile_number)

#RUN model ---------------------------------------------------------------------
Implicit.Model(FOR, GRI, STRA, PAR, STAT, TEM, OUT; start=DateTime(1900,1,1))

#REDUCE output precision -------------------------------------------------------
SAVE = CryoGridTyps.save(OUT.Date, OUT.H_av, OUT.T_av, OUT.T_min, OUT.T_max, OUT.W_av, OUT.W_min, OUT.W_max, OUT.Q_lat, OUT.FDD, OUT.TDD, OUT.FrostDays, OUT.SnowDepth_av, OUT.SnowDepth_max, OUT.SnowDays)

#SAVE result as JSON for each year ---------------------------------------------
subfolder = "RESULT_JSONfiles"
mkpath(subfolder);
subsubfolder = "ULC_" * string(ULC_lon) * "_" * string(ULC_lat)
mkpath(subfolder * "/" * subsubfolder);
for idx = 1:length(SAVE.Date)
        out_dict = Dict()
        year = SAVE.Date[idx]
        FieldsInStruct = fieldnames(typeof(SAVE))
        for i = 1:length(FieldsInStruct)
                #Check fields
                Value = getfield(SAVE, FieldsInStruct[i])
                if size(Value)[2]>=idx
                        value = Value[:,idx,:]
                else
                        value = Value[:,end,:]
                end
                name = string(FieldsInStruct[i])
                out_dict[name] = value
        end
        #save as json uncompressed
        saveloc = subfolder * "/" * subsubfolder * "/" * "RESULT_" * string(year) * ".zjson"
        open(saveloc,"w") do file
                write(file,JSON.json(out_dict))
        end
end
