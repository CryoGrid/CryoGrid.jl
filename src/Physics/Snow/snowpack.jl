"""
    SnowMassBalance{TAcc,TAbl} <: CryoGrid.SubSurfaceProcess

Subsurface process for snow layers governing how snow is accumulated and ablated.
"""
Base.@kwdef struct SnowMassBalance{TAcc,TAbl,TAux,TDt} <: CryoGrid.SubSurfaceProcess
    accumulation::TAcc = LinearAccumulation()
    ablation::TAbl = DegreeDayMelt()
    dtlim::TDt = CryoGrid.MaxDelta(0.01)
    aux::TAux = nothing
end

"""
    SnowBC

Type alias for any `BoundaryProcess` compatible with `SnowMassBalance`.
"""
const SnowBC = BoundaryProcess{T} where {SnowMassBalance<:T<:SubSurfaceProcess}

"""
    SnowThermalProperties{Tcond<:SnowThermalConductivity,Tprop}

Specifies the thermal properties of the snowpack.
"""
Base.@kwdef struct SnowThermalProperties{Tcond<:SnowThermalConductivity,Tprop}
    cond::Tcond = SturmQuadratic()
    prop::Tprop = ThermalProperties()
end

"""
    Snowpack{Tpara<:SnowpackParameterization,Tmass<:SnowMassBalance,Twater<:WaterBalance,Theat<:HeatBalance,Taux} <: CryoGrid.SubSurface

Generic representation of a snowpack "subsurface" layer.
"""
Base.@kwdef struct Snowpack{Tpara<:SnowpackParameterization,Tmass<:SnowMassBalance,Twater<:WaterBalance,Theat<:HeatBalance,Taux} <: CryoGrid.SubSurface
    para::Tpara = Bulk()
    mass::Tmass = SnowMassBalance()
    water::Twater = WaterBalance()
    heat::Theat = HeatBalance()
    aux::Taux = nothing
end

# Processes type aliases
const CoupledSnowWaterHeat{Tmass,Twater,Theat} = Coupled(SnowMassBalance, WaterBalance, HeatBalance)

"""
    Snowpack(para::SnowpackParameterization; kwargs...)

Convenience constructor that accepts the parameterization as a positional argument.
"""
Snowpack(para::SnowpackParameterization; kwargs...) = Snowpack(; para, kwargs...)

### Default implmentations of other module methods common to all Snowpack types ###

# Hydrology methods;
Hydrology.hydraulicproperties(snow::Snowpack) = snow.para.water

# for snow, use constant hydraulic conductivity
Hydrology.hydraulicconductivity(snow::Snowpack, water::WaterBalance, θw, θwi, θsat) = Hydrology.kwsat(snow, water)

# max (fully saturated) water content
Hydrology.maxwater(::Snowpack, ::WaterBalance, state) = 1.0

# Heat methods:

Heat.thermalproperties(snow::Snowpack) = snow.para.heat.prop

# Default implementations of CryoGrid methods for Snowpack
CryoGrid.processes(snow::Snowpack) = Coupled(snow.mass, snow.water, snow.heat)

CryoGrid.isactive(snow::Snowpack, state) = CryoGrid.thickness(snow, state) > threshold(snow)

CryoGrid.Volume(::Type{<:Snowpack{T,<:SnowMassBalance}}) where {T} = CryoGrid.DiagnosticVolume()

# volumetric fractions for snowpack
function CryoGrid.volumetricfractions(::Snowpack, state, i)
    @inbounds let θwi = state.θwi[i],
        θw = state.θw[i],
        θa = 1.0 - θwi,
        θi = θwi - θw;
        return (θw, θi, θa)
    end
end

# handle snow depth initializers; since swe is the actual prognostic state variable,
# we need to make sure this is computed as a function of the user initialized snow depth
function CryoGrid.initialcondition!(init!::VarInitializer{:dsn}, snow::Snowpack, state)
    # set snow depth variable
    init!(snow, state)
    # compute snow density
    snowdensity!(snow, snow.mass, state)
    # compute snow water equivalent
    ρsn = state.ρsn
    ρw = waterdensity(snow)
    @. state.swe = state.dsn * ρsn / ρw
end

function CryoGrid.computediagnostic!(
    snow::Snowpack,
    state,
)
    mass, water, heat = processes(snow)
    computediagnostic!(snow, mass, state)
    computediagnostic!(snow, water, state)
    computediagnostic!(snow, heat, state)
end

function CryoGrid.computediagnostic!(
    snow::Snowpack,
    mass::SnowMassBalance,
    state,
)
    # update snow density
    snowdensity!(snow, mass, state)
    # update snow depth;
    # by default, we just use the current layer thickness
    @setscalar state.dsn = getscalar(state.Δz)
end

# Special overrides for heat timestep control on snow layer

CryoGrid.timestep(::Snowpack, heat::HeatBalance{<:FreeWater,THeatOp,<:CryoGrid.CFL}, state) where {THeatOp} = error("CFL is not supported on snow layer")
function CryoGrid.timestep(snow::Snowpack, heat::HeatBalance{<:FreeWater,THeatOp,<:CryoGrid.MaxDelta}, state) where {THeatOp}
    Δx = Δ(state.grid)
    dtmax = Inf
    if isactive(snow, state)
        @inbounds for i in eachindex(Δx)
            dtmax = min(dtmax, heat.dtlim(state.dH[i], state.H[i], state.t))
        end
        dtmax = isfinite(dtmax) && dtmax > 0 ? dtmax : Inf
    end
    return dtmax
end
