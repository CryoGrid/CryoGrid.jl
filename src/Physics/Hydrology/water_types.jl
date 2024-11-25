"""
    WaterBalanceProperties

Numerical constants shared across water balance implementations.
"""
Base.@kwdef struct WaterBalanceProperties{Tρw,TLsg,Trfs}
    ρw::Tρw = CryoGrid.Constants.ρw
    Lsg::TLsg = CryoGrid.Constants.Lsg
    rf_smoothness::Trfs = 0.3
end

"""
    HydraulicProperties

Default material hydraulic properties.
"""
Utils.@properties HydraulicProperties(
    kw_sat = Param(1e-5, units=u"m/s"),
)

"""
    WaterFlow

Base type for different formulations of water flow in `WaterBalance`.
"""
abstract type WaterFlow end

"""
    Evapotranspiration

Base type for parameterizations of evapotranspiration (ET).
"""
abstract type Evapotranspiration end

"""
    WaterBalance{TFlow<:WaterFlow,TET<:Union{Nothing,Evapotranspiration},Tdt,Taux,TProp} <: CryoGrid.SubSurfaceProcess

Represents subsurface water transport processes.
"""
struct WaterBalance{TFlow<:WaterFlow,TET<:Union{Nothing,Evapotranspiration},Tdt,Taux,TProp<:WaterBalanceProperties} <: CryoGrid.SubSurfaceProcess
    flow::TFlow # vertical flow scheme
    et::TET # evapotranspiration scheme
    prop::TProp # hydraulic parameters/constants
    dtlim::Tdt # dtlim
    aux::Taux # user-defined specialization
end

"""
    NoFlow <: WaterFlow

Represents a zero flow scheme 
"""
struct NoFlow <: WaterFlow end

"""
    BucketScheme <: WaterFlow

"Bucket" water scheme for downward advective flow due to gravity.
"""
Base.@kwdef struct BucketScheme <: WaterFlow end

# default dt limiters
default_dtlim(::BucketScheme) = CryoGrid.MaxDelta(0.01)
default_dtlim(::WaterFlow) = CryoGrid.MaxDelta(Inf)
# default ET scheme
default_ET(::BucketScheme) = DampedET()
default_ET(::WaterFlow) = nothing

WaterBalance(flow::WaterFlow = NoFlow(), et=default_ET(flow); prop = WaterBalanceProperties(), dtlim = default_dtlim(flow), aux = nothing) = WaterBalance(flow, et, prop, dtlim, aux)
