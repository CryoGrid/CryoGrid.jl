"""
    SoilParameterization

Base type for parameterizations of soil consituents.
"""
abstract type SoilParameterization <: GroundParameterization end

"""
    Soil{Tpara,Theat,Twater}

Type alias for any `AbstractGround` layer with a `SoilParameterization`.
"""
const Soil{Tpara,Theat,Twater} = AbstractGround{Tpara,Theat,Twater} where {Tpara<:SoilParameterization,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance}}

"""
    Heterogeneous{Tpara,Taux} <: SoilParameterization

Special `SoilParameterization` which wraps another soil parameterization type
to indicate that it should be heterogeneous with over depth. Parameterizations
that support such configurations should provide dispatches for `Heterogeneous{...}`
that instantiate the relevant soil properties as on-grid state variables.
"""
Base.@kwdef struct Heterogeneous{Tpara,Taux} <: SoilParameterization
    para::Tpara
    aux::Taux = nothing
    Heterogeneous(para::SoilParameterization, aux=nothing) = new{typeof(para),typeof(aux)}(para, aux)
end

# forward getproperty to nested parameterization
Base.propertynames(para::Heterogeneous) = (:para, Base.propertynames(getfield(para, :para))...)
Base.getproperty(para::Heterogeneous, name::Symbol) = name == :para ? getfield(para, :para) : getproperty(getfield(para, :para), name)

checksolver!(::Heterogeneous, ::SFCCPreSolver) = error("SFCCPreSolver requires homogeneous soil properties in each layer.")
