"""
    HomogeneousMixture{Tpor,Tsat,Torg} <: SoilParameterization

Represents a simple, uniform organic/mineral soil mixutre in terms of its characteristic fractions:
i.e. natural porosity, saturation, and organic solid fraction.
"""
Base.@kwdef struct HomogeneousMixture{Tpor,Tsat,Torg} <: SoilParameterization
    por::Tpor = 0.5 # natural porosity
    sat::Tsat = 1.0 # saturation
    org::Torg = 0.0 # organic fraction of solid; mineral fraction is 1-org
end
CryoGrid.parameterize(para::HomogeneousMixture) = HomogeneousMixture(
    por = CryoGrid.parameterize(para.por, domain=0..1),
    sat = CryoGrid.parameterize(para.sat, domain=0..1),
    org = CryoGrid.parameterize(para.org, domain=0..1),
)
# Helper functions for obtaining soil compositions from characteristic fractions.
soilcomponent(::Val{var}, para::HomogeneousMixture) where var = soilcomponent(Val{var}(), para.por, para.sat, para.org)
soilcomponent(::Val{:θp}, ϕ, θ, ω) = ϕ
soilcomponent(::Val{:θwi}, ϕ, θ, ω) = ϕ*θ
soilcomponent(::Val{:θm}, ϕ, θ, ω) = (1-ϕ)*(1-ω)
soilcomponent(::Val{:θo}, ϕ, θ, ω) = (1-ϕ)*ω
"""
    from_components(θo, por, sat=1.0)

Constructs a `HomogeneousMixture` soil parameterization from the component volumetric fractions of
organic content and natural porosity `por`. Saturation `sat` can also optionally be specfied.
"""
function from_components(θo, sat=1.0)
    @assert zero(θo) <= θo <= one(θo)
    @assert zero(por) <= por <= one(por)
    @assert zero(sat) <= por <= one(sat)
    @assert zero(por) <= θo + por <= one(por)
    θm = 1 - por - θo
    org = θo / (θm + θo)
    return HomogeneousMixture(;por, sat, org)
end
"""
    MineralSediment{Tpor} <: SoilParameterization

Represents a simple, mineral sediment with (possibly) spatially variable porosity.
"""
Base.@kwdef struct MineralSediment{Tpor} <: SoilParameterization
    por::Tpor = 0.5
end
CryoGrid.parameterize(para::MineralSediment) = MineralSediment(
    por = CryoGrid.parameterize(para.por, domain=0..1),
)
