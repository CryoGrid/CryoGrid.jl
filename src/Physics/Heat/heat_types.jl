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
    HeatBalance{THeatOp<:HeatOperator,Tdt,Tprop} <: SubSurfaceProcess

Represents subsurface heat transfer processes. The formulation of heat transfer is governed by
the `HeatOperator`, `op`. 
"""
Base.@kwdef struct HeatBalance{THeatOp<:HeatOperator,Tdt,Tprop} <: SubSurfaceProcess
    op::THeatOp = Diffusion1D(:H)
    prop::Tprop = HeatBalanceProperties()
    dtlim::Tdt = default_dtlim(op)  # timestep limiter
    advection::Bool = true # whether or not to include advective fluxes when coupled with WaterBalance
end

# convenience constructors for HeatBalance
HeatBalance(var::Symbol; kwargs...) = HeatBalance(Val{var}(); kwargs...)
HeatBalance(::Val{:H}; kwargs...) = HeatBalance(; op=Diffusion1D(:H), kwargs...)
HeatBalance(::Val{:T}; kwargs...) = HeatBalance(; op=Diffusion1D(:T), kwargs...)
HeatBalance(op::HeatOperator; kwargs...) = HeatBalance(; op, kwargs...)

# Heat operators
"""
    Diffusion1D{progvar} <: HeatOperator{progvar}

Represents a standard method-of-lines (MOL) forward diffusion operator in 1 dimension.
"""
struct Diffusion1D{progvar} <: HeatOperator{progvar}
    Diffusion1D(progvar::Symbol=:H) = new{progvar}()
end
ConstructionBase.constructorof(::Type{<:Diffusion1D{progvar}}) where {progvar} = () -> Diffusion1D(progvar)

"""
    EnthalpyImplicit <: HeatOperator{:H}

Implicit enthalpy formulation of Swaminathan and Voller (1992) and Langer et al. (2022). Note that this
heat operator formulation does not compute a divergence `dH` but only computes the necessary diffusion
coefficients for use by an appropriate solver. See the `LiteImplicit` module for the appropriate
solver algorithms.
"""
struct EnthalpyImplicit <: HeatOperator{:H} end

"""
Numerical constants for pararameterizing heat processes.
"""
Utils.@properties HeatBalanceProperties(
    ρw = CryoGrid.Constants.ρw,
    Lsl = CryoGrid.Constants.Lsl,
    L = ρw*Lsl,
)
