# SURFEX paramterization
Base.@kwdef struct SURFEX{Tρs,Tρo,Tpor} <: SoilParameterization
    ρ_soc::Tρs
    ρ_org::Tρo
    por_org::Tpor
end
porosity(soil::Soil{<:SURFEX}) = nothing