"""
    WaterBalanceProperties

Numerical constants shared across water balance implementations.
"""
Base.@kwdef struct WaterBalanceProperties{Tρw,TLsg,Trb,Trc}
    ρw::Tρw = CryoGrid.Constants.ρw
    Lsg::TLsg = CryoGrid.Constants.Lsg
    r_β::Trb = 1e3 # reduction factor scale parameter
    r_c::Trc = 0.96325 # reduction factor shift parameter
end
CryoGrid.parameterize(prop::WaterBalanceProperties) = prop
"""
    HydraulicProperties

Default material hydraulic properties.
"""
Utils.@properties HydraulicProperties(
    kw_sat = 1e-5u"m/s",
)
hydraulicproperties(::SubSurface) = HydraulicProperties()
function CryoGrid.parameterize(prop::HydraulicProperties)
    return HydraulicProperties(
        map(values(prop)) do val
            # this currently assumes that all properties have a strictly positive domain!
            CryoGrid.parameterize(val, domain=StrictlyPositive)
        end
    )
end
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
    WaterBalance{TFlow<:WaterFlow,TET<:Union{Nothing,Evapotranspiration},Tdt,Tsp,TProp} <: CryoGrid.SubSurfaceProcess

Represents subsurface water transport processes.
"""
struct WaterBalance{TFlow<:WaterFlow,TET<:Union{Nothing,Evapotranspiration},Tdt,Tsp,TProp<:WaterBalanceProperties} <: CryoGrid.SubSurfaceProcess
    flow::TFlow # vertical flow scheme
    et::TET # evapotranspiration scheme
    prop::TProp # hydraulic parameters/constants
    dtlim::Tdt # dtlim
    sp::Tsp # user-defined specialization
end
"""
    NoFlow <: WaterFlow

Represents a zero flow scheme 
"""
struct NoFlow <: WaterFlow end
"""
    BucketScheme{Tfc} <: WaterFlow

"Bucket" water scheme for downward advective flow due to gravity.
"""
Base.@kwdef struct BucketScheme{Tfc} <: WaterFlow
    fieldcap::Tfc = 0.2
end
CryoGrid.parameterize(flow::BucketScheme) = BucketScheme(
    fieldcap = CryoGrid.parameterize(flow.fieldcap, domain=0..1, desc="Minimum saturation level, a.k.a 'field capacity'."),
)
# default dt limiters
default_dtlim(::BucketScheme) = CryoGrid.MaxDelta(0.1)
default_dtlim(::WaterFlow) = CryoGrid.MaxDelta(Inf)
# default ET scheme
default_ET(::BucketScheme) = DampedET()
default_ET(::WaterFlow) = nothing

WaterBalance(flow::WaterFlow = BucketScheme(), et=default_ET(flow); prop = WaterBalanceProperties(), dtlim = default_dtlim(flow), sp = nothing) = WaterBalance(flow, et, prop, dtlim, sp)
