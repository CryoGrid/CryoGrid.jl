"""
    HomogeneousMixture{Tpor,Tsat,Torg} <: SoilParameterization

Represents a simple, uniform organic/mineral soil mixutre in terms of its characteristic fractions:
i.e. natural porosity, saturation, organic solid fraction, and excess ice fraction.
"""
Base.@kwdef struct HomogeneousMixture{Tpor,Tsat,Torg,Txic} <: SoilParameterization
    por::Tpor = 0.5 # natural porosity
    sat::Tsat = 1.0 # saturation
    org::Torg = 0.0 # organic fraction of solid; mineral fraction is 1-org
    xic::Txic = 0.0 # excess ice fraction of total volume
end
CryoGrid.parameterize(para::HomogeneousMixture) = HomogeneousMixture(
    por = CryoGrid.parameterize(para.por, domain=0..1),
    sat = CryoGrid.parameterize(para.sat, domain=0..1),
    org = CryoGrid.parameterize(para.org, domain=0..1),
    xic = CryoGrid.parameterize(para.xic, domain=0..1),
)
# Helper functions for obtaining soil compositions from characteristic fractions.
soilcomponent(::Val{var}, para::HomogeneousMixture) where var = soilcomponent(Val{var}(), para.por, para.sat, para.org, para.xic)
soilcomponent(::Val{:θwi}, ϕ, θ, ω, χ) = χ + ϕ*θ
soilcomponent(::Val{:θp}, ϕ, θ, ω, χ) = (1-χ)*ϕ
soilcomponent(::Val{:θa}, ϕ, θ, ω, χ) = (1-χ)*ϕ*(1-θ)
soilcomponent(::Val{:θm}, ϕ, θ, ω, χ) = (1-χ)*(1-ϕ)*(1-ω)
soilcomponent(::Val{:θo}, ϕ, θ, ω, χ) = (1-χ)*(1-ϕ)*ω
"""
    from_components(θo, por, sat=1.0)

Constructs a `HomogeneousMixture` soil parameterization from the component volumetric fractions of
organic content and natural porosity `por`. Saturation `sat` and excess ice `xic` can also optionally be specfied.
"""
function from_components(θo, por, sat=1.0, xic=0.0)
    @assert zero(θo) <= θo <= one(θo)
    @assert zero(por) <= por <= one(por)
    @assert zero(sat) <= por <= one(sat)
    @assert zero(xic) <= xic <= one(xic)
    @assert zero(por) <= θo + max(por, xic) <= one(por)
    θm = 1 - max(por, xic) - θo
    org = θo / (θm + θo)
    return HomogeneousMixture(;por, sat, org, xic)
end
"""
    MineralSediment{Tpor,Tsat} <: SoilParameterization

Represents a simple, mineral sediment with (possibly) spatially variable porosity.
"""
Base.@kwdef struct MineralSediment{Tpor,Tsat} <: SoilParameterization
    por::Tpor = 0.5
    sat::Tsat = 1.0
end
CryoGrid.parameterize(para::MineralSediment) = MineralSediment(
    por = CryoGrid.parameterize(para.por, domain=0..1),
    sat = CryoGrid.parameterize(para.sat, domain=0..1),
)
