"""
    HeatOperator{progvar}

Base type for different numerical formulations of heat transfer.
"""
abstract type HeatOperator{progvar} end
"""
Type alias for `HeatOperator{:H}`, i.e. enthalpy-based heat transfer operators.
"""
const Enthalpy = HeatOperator{:H}
"""
Type alias for `HeatOperator{:T}`, i.e. temperature-based heat transfer operators.
"""
const Temperature = HeatOperator{:T}

# default step limiters
default_dtlim(::Temperature) = CryoGrid.CFL(maxdelta=CryoGrid.MaxDelta(Inf))
default_dtlim(::Enthalpy) = CryoGrid.MaxDelta(50u"kJ")
default_dtlim(::HeatOperator) = nothing

# default freezecurve solvers
default_fcsolver(::FreeWater) = nothing
default_fcsolver(::SFCC) = SFCCPreSolver()

# Heat Balance type
"""
    HeatBalance{Tfc<:FreezeCurve,THeatOp<:HeatOperator,Tdt,Tprop} <: SubSurfaceProcess

Represents subsurface heat transfer processes. The formulation of heat transfer is governed by
the `HeatOperator`, `op`. 
"""
Base.@kwdef struct HeatBalance{Tfc<:FreezeCurve,THeatOp<:HeatOperator,Tdt,Tprop} <: SubSurfaceProcess
    freezecurve::Tfc = FreeWater()
    op::THeatOp = EnthalpyForm(default_fcsolver(freezecurve))
    prop::Tprop = HeatBalanceProperties()
    dtlim::Tdt = default_dtlim(op)  # timestep limiter
    function HeatBalance(freezecurve, op, prop, dtlim)
        # check that heat configuration is valid
        _validate_heat_config(freezecurve, op)
        op = deepcopy(op) # make a deep copy in case of freeze curve solver caches
        return new{typeof(freezecurve), typeof(op), typeof(dtlim), typeof(prop)}(freezecurve, op, prop, dtlim)
    end
end
# convenience constructors for HeatBalance
HeatBalance(var::Symbol; kwargs...) = HeatBalance(Val{var}(); kwargs...)
HeatBalance(::Val{:H}; freezecurve::FreezeCurve=FreeWater, fcsolver=default_fcsolver(freezecurve), kwargs...) = HeatBalance(; op=EnthalpyForm(fcsolver), freezecurve, kwargs...)
HeatBalance(::Val{:T}; freezecurve::FreezeCurve, kwargs...) = HeatBalance(; op=TemperatureForm(), freezecurve, kwargs...)
HeatBalance(op::HeatOperator; kwargs...) = HeatBalance(; op, kwargs...)
# validation of HeatBalance freezecurve/operator configuration
_validate_heat_config(::FreezeCurve, ::HeatOperator) = nothing # do nothing when valid
_validate_heat_config(::FreeWater, ::Temperature) = error("Invalid heat balance configuration; temperature formulations of the heat operator are not compatible with the free water freeze curve.")
# Heat operators
"""
    TemperatureForm{Tcond,Thc} <: HeatOperator{:T}

Represents a standard method-of-lines (MOL) forward diffusion operator for heat conduction with
temperature `T` as the prognostic variable. The time derivative is scaled by the reciprocal of
the apparent heat capacity `dH/dT` to account for latent heat effects due to phase change.
"""
struct TemperatureForm{Tcond,Thc} <: HeatOperator{:T}
    cond::Tcond
    hc::Thc
    TemperatureForm(cond=quadratic_parallel_conductivity, hc=weighted_average_heatcapacity) = new{typeof(cond),typeof(hc)}(cond, hc)
end
"""
    EnthalpyForm{Tsolver,Tcond,Thc} <: HeatOperator{:H}

Represents a standard method-of-lines (MOL) forward diffusion operator for heat conduction with
enthalpy `H` as the prognostic variable and a nonlinear solver for resolving the inverse
enthalpy -> temperature mapping when applicable. This formulation should generally be preferred
over `TemperatureForm` since it is energy-conserving and embeds the latent heat storage directly
in the prognostic state.
"""
struct EnthalpyForm{Tsolver,Tcond,Thc} <: HeatOperator{:H}
    fcsolver::Tsolver
    cond::Tcond
    hc::Thc
    EnthalpyForm(fcsolver=nothing, cond=quadratic_parallel_conductivity, hc=weighted_average_heatcapacity) = new{typeof(fcsolver),typeof(cond),typeof(hc)}(fcsolver, cond, hc)
end
"""
    EnthalpyImplicit <: HeatOperator{:H}

Implicit enthalpy formulation of Swaminathan and Voller (1992) and Langer et al. (2022). Note that this
heat operator formulation does not compute a divergence `∂H∂t` but only computes the necessary diffusion
coefficients for use by an appropriate solver. See the `LiteImplicit` module for the appropriate
solver algorithms.
"""
struct EnthalpyImplicit{Tsolver,Tcond,Thc} <: HeatOperator{:H}
    fcsolver::Tsolver
    cond::Tcond
    hc::Thc
    EnthalpyImplicit(fcsolver=nothing, cond=quadratic_parallel_conductivity, hc=weighted_average_heatcapacity) = new{typeof(fcsolver),typeof(cond),typeof(hc)}(fcsolver, cond, hc)
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

"""
    ThermalProperties

Basic material thermal properties.
"""
Utils.@properties ThermalProperties(
    kh_w = 0.57u"J/s/m/K", # thermal conductivity of water [Hillel (1982)]
    kh_i = 2.2u"J/s/m/K", # thermal conductivity of ice [Hillel (1982)]
    kh_a = 0.025u"J/s/m/K", # thermal conductivity of air [Hillel (1982)]
    ch_w = 4.2e6u"J/K/m^3", # heat capacity of water
    ch_i = 1.9e6u"J/K/m^3", # heat capacity of ice
    ch_a = 0.00125e6u"J/K/m^3", # heat capacity of air
)
function CryoGrid.parameterize(props::ThermalProperties)
    @set! props.kh_w = CryoGrid.parameterize(props.kh_w)
    @set! props.kh_i = CryoGrid.parameterize(props.kh_i)
    @set! props.kh_a = CryoGrid.parameterize(props.kh_a)
    @set! props.ch_w = CryoGrid.parameterize(props.ch_w)
    @set! props.ch_i = CryoGrid.parameterize(props.ch_i)
    @set! props.ch_a = CryoGrid.parameterize(props.ch_a)
    return props
end
