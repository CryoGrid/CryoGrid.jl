"""
    SWRCFunction

Base type for soil water retention curve SWRC function implementations.
"""
abstract type SWRCFunction end
"""
    SWRC{F,∇F}

Soil water retention curve with function type `F` and derivative `∇F`.
"""
struct SWRC{F,∇F}
    f::F # soil water retention curve function f: (ψ,...) -> θ
    ∇f::∇F # derivative of SWRC function
    SWRC(f::F,∇f::∇F) where {F<:SWRCFunction,∇F<:Function} = new{F,∇F}(f,∇f)
end
"""
    SWRC(f::SWRCFunction)

Convenience constructor for SWRC that automatically generates an analytical derivative of the given
retention curve function `f` using ModelingToolkit/Symbolics.jl.
"""
function SWRC(f::SWRCFunction; dvar=:ψ, choosefn=first, context_module=Numerics)
    ∇f = ∇(f, dvar; choosefn=choosefn, context_module=context_module)
    # we wrap ∇f with Base.splat here to avoid a weird issue with in-place splatting causing allocations
    # when applied to runtime generated functions.
    SWRC(f, Base.splat(∇f), s)
end
"""
    VanGenuchten <: SWRCFunction

van Genuchten MT, 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
    Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj 1980.03615995004400050002x.
"""
struct VanGenuchten <: SWRCFunction end
function (f::VanGenuchten)(ψ,θres,θsat,α,n)
    let m = 1-1/n;
        IfElse.ifelse(ψ <= zero(ψ), θres + (θsat - θres)*(1 + (-α*ψ)^n)^(-m), θsat)
    end
end
