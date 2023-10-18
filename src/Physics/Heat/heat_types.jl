"""
    HeatOperator{progvar}

Base type for different numerical formulations of heat transfer.
"""
abstract type HeatOperator{progvar} end
"""
Type alias for `HeatOperator{:H}`, i.e. enthalpy-based heat transfer operators.
"""
const EnthalpyBased = HeatOperator{:H}
"""
Type alias for `HeatOperator{:T}`, i.e. temperature-based heat transfer operators.
"""
const TemperatureBased = HeatOperator{:T}

# default step limiters
default_dtlim(::TemperatureBased) = CryoGrid.CFL(maxdelta=CryoGrid.MaxDelta(Inf))
default_dtlim(::EnthalpyBased) = CryoGrid.MaxDelta(50u"kJ")
default_dtlim(::HeatOperator) = nothing

# Heat Balance type
"""
    HeatBalance{Tfc<:FreezeCurve,THeatOp<:HeatOperator,Tdt,Tprop} <: SubSurfaceProcess

Represents subsurface heat transfer processes. The formulation of heat transfer is governed by
the `HeatOperator`, `op`. 
"""
Base.@kwdef struct HeatBalance{Tfc<:FreezeCurve,THeatOp<:HeatOperator,Tdt,Tprop} <: SubSurfaceProcess
    freezecurve::Tfc = FreeWater()
    op::THeatOp = Diffusion1D(:H)
    prop::Tprop = HeatBalanceProperties()
    dtlim::Tdt = default_dtlim(op)  # timestep limiter
    advection::Bool = true # whether or not to include advective fluxes when coupled with WaterBalance
    function HeatBalance(freezecurve, op, prop, dtlim, advection)
        # check that heat configuration is valid
        _validate_heat_config(freezecurve, op)
        return new{typeof(freezecurve),typeof(op),typeof(dtlim),typeof(prop)}(freezecurve, op, prop, dtlim, advection)
    end
end
# convenience constructors for HeatBalance
HeatBalance(var::Symbol; kwargs...) = HeatBalance(Val{var}(); kwargs...)
HeatBalance(::Val{:H}; freezecurve::FreezeCurve=FreeWater(), kwargs...) = HeatBalance(; op=Diffusion1D(:H), freezecurve, kwargs...)
HeatBalance(::Val{:T}; freezecurve::FreezeCurve, kwargs...) = HeatBalance(; op=Diffusion1D(:T), freezecurve, kwargs...)
HeatBalance(op::HeatOperator; kwargs...) = HeatBalance(; op, kwargs...)

# validation of HeatBalance freezecurve/operator configuration
_validate_heat_config(::FreezeCurve, ::HeatOperator) = nothing # do nothing when valid
_validate_heat_config(::FreeWater, ::TemperatureBased) = error("Invalid heat balance configuration; temperature formulations of the heat operator are not compatible with the free water freeze curve.")

# Heat operators
"""
    Diffusion1D{progvar,Tcond,Thc} <: HeatOperator{progvar}

Represents a standard method-of-lines (MOL) forward diffusion operator in 1 dimension.
"""
struct Diffusion1D{progvar,Tcond,Thc} <: HeatOperator{progvar}
    cond::Tcond
    hc::Thc
    Diffusion1D(progvar::Symbol=:H, cond=quadratic_parallel_conductivity, hc=weighted_average_heatcapacity) = new{progvar,typeof(cond),typeof(hc)}(cond, hc)
end
ConstructionBase.constructorof(::Type{<:Diffusion1D{progvar}}) where {progvar} = (cond, hc) -> Diffusion1D(progvar, cond, hc)

"""
    EnthalpyImplicit{Tcond,Thc} <: HeatOperator{:H}

Implicit enthalpy formulation of Swaminathan and Voller (1992) and Langer et al. (2022). Note that this
heat operator formulation does not compute a divergence `dH` but only computes the necessary diffusion
coefficients for use by an appropriate solver. See the `LiteImplicit` module for the appropriate
solver algorithms.
"""
struct EnthalpyImplicit{Tcond,Thc} <: HeatOperator{:H}
    cond::Tcond
    hc::Thc
    EnthalpyImplicit(cond=quadratic_parallel_conductivity, hc=weighted_average_heatcapacity) = new{typeof(cond),typeof(hc)}(cond, hc)
end

"""
Numerical constants for pararameterizing heat processes.
"""
Utils.@properties HeatBalanceProperties(
    ρw = CryoGrid.Constants.ρw,
    Lsl = CryoGrid.Constants.Lsl,
    L = ρw*Lsl,
)
# do not parameterize heat properties
CryoGrid.parameterize(prop::HeatBalanceProperties) = prop
