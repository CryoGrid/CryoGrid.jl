abstract type SnowAblationScheme end
struct ConstantAblation <: SnowAblationScheme
    rate::Trate = 2.0u"mm/hr"
end

abstract type SnowAccumulationScheme end
struct LinearAccumulation <: SnowAccumulationScheme end

abstract type SnowDensityScheme end
Base.@kwdef struct ConstantDensity{Tρsn}
    ρsn::Tρsn = 250.0u"kg/m^3" # snow density [kg m^-3]
end

abstract type SnowMassParameterization end
Base.@kwdef struct Prescribed{Tswe,Tρsn} <: SnowMassParameterization
    swe::Tswe = 0.0u"m" # depth snow water equivalent [m]
    ρsn::Tρsn = 250.0u"kg/m^3" # snow density [kg m^-3]
end
Base.@kwdef struct Dynamic{TAcc,TMelt,TDensity} <: SnowMassParameterization
    accumulation::TAcc = LinearAccumulation()
    ablation::TMelt = InstantaneousRunoff()
    density::TDensity = ConstantDensity()
end

Base.@kwdef struct SnowMassBalance{Tpara} <: CryoGrid.SubSurfaceProcess
    para::Tpara = Dynamic()
end

snowvariables(::Snowpack, ::SnowMassBalance) = (
    Diagnostic(:dsn, Scalar, u"m"),
    Diagnostic(:T_ub, Scalar, u"°C"),
)

snow_water_equivalent(::Snowpack, ::SnowMassBalance, state) = state.swe
snow_water_equivalent(::Snowpack, smb::SnowMassBalance{<:Prescribed}, state) = smb.para.swe
snow_water_equivalent(::Snowpack, smb::SnowMassBalance{<:Prescribed{<:Forcing}}, state) = smb.para.swe(state.t)
snowdensity(::Snowpack, ::SnowMassBalance, state) = state.ρsn
snowdensity(::Snowpack, smb::SnowMassBalance{<:Prescribed}, state) = smb.para.ρsn
snowdensity(::Snowpack, smb::SnowMassBalance{<:Prescribed{Tswe,<:Forcing}}, state) where {Tswe} = smb.para.ρsn(state.t)

# Boundary conditions

struct Snowfall{Tsn<:Forcing} <: BoundaryProcess{SnowMassBalance}
    snowfall::Tsn
end
CryoGrid.BoundaryStyle(::Snowfall) = Neumann()
@inline boundaryvalue(bc::Snowfall, ::Top, ::SnowMassBalance, ::Snowpack, s1, s2) = bc.snowfall(s1.t)

# Implementations

# for prescribed snow depth/density, the mass balance is given so we do not need to do anything here
CryoGrid.prognosticstep!(::Snowpack, ::SnowMassBalance{<:Prescribed}, ssnow) = nothing

include("snow_bulk.jl")
