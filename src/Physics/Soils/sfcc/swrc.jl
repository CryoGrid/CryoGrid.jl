"""
    SWRCFunction

Base type for soil water retention curve SWRC function implementations.
"""
abstract type SWRCFunction end
Base.inv(f::SWRCFunction) = (args...) -> f(inv, args...)
"""
    SWRC{F}

Soil water retention curve with function type `F`.
"""
struct SWRC{F}
    f::F # soil water retention curve function f: (ψ,...) -> θ
    SWRC(f::F) where {F<:SWRCFunction} = new{F}(f)
end
"""
    VanGenuchten <: SWRCFunction

van Genuchten MT, 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
    Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj 1980.03615995004400050002x.
"""
Base.@kwdef struct VanGenuchten{Tα,Tn} <: SWRCFunction
    α::Tα = Param(1.0, units=u"1/m")
    n::Tn = Param(2.0)
end
function (f::VanGenuchten)(ψ, θres, θsat, α=f.α, n=f.n)
    let m = 1-1/n;
        IfElse.ifelse(ψ <= zero(ψ), θres + (θsat - θres)*(1 + abs(-α*ψ)^n)^(-m), θsat)
    end
end
function (f::VanGenuchten)(::typeof(inv), θ, θres, θsat, α=f.α, n=f.n)
    let m = 1-1/n;
        IfElse.ifelse(θ < θsat, -1/α*(((θ-θres)/(θsat-θres))^(-1/m)-1.0)^(1/n), zero(1/α))
    end
end
