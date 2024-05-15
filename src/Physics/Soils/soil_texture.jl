"""
    SoilTexture{Tsand,Tclay,Tsilt}

Represents soil "texture" as a simple mixture of sand, silt, and clay.
"""
Base.@kwdef struct SoilTexture{Tsand,Tclay,Tsilt}
    sand::Tsand = 1.0
    clay::Tclay = 0.0
    silt::Tsilt = 1 - sand - clay
    function SoilTexture(sand, silt, clay)
        @assert sand + silt + clay ≈ 1.0 "sand, silt, and clay fractions must sum to unity"
        return new{typeof(sand),typeof(silt),typeof(clay)}(sand, silt, clay)
    end
end

# Presets for common soil textures
SoilTexture(name::Symbol) = SoilTexture(Val{name}())
SoilTexture(::Val{:sand}) = SoilTexture(sand=1.0, silt=0.0, clay=0.0)
SoilTexture(::Val{:silt}) = SoilTexture(sand=0.0, silt=1.0, clay=0.0)
SoilTexture(::Val{:clay}) = SoilTexture(sand=0.0, silt=0.0, clay=1.0)
SoilTexture(::Val{:sandyclay}) = SoilTexture(sand=0.50, silt=0.0, clay=0.50)
SoilTexture(::Val{:siltyclay}) = SoilTexture(sand=0.0, silt=0.50, clay=0.50)
SoilTexture(::Val{:loam}) = SoilTexture(sand=0.40, silt=0.40, clay=0.20)
SoilTexture(::Val{:sandyloam}) = SoilTexture(sand=0.80, silt=0.10, clay=0.10)
SoilTexture(::Val{:siltyloam}) = SoilTexture(sand=0.10, silt=0.80, clay=0.10)
SoilTexture(::Val{:clayloam}) = SoilTexture(sand=0.30, silt=0.30, clay=0.40)

function Heat.freezecurve(texture::SoilTexture)
    return if texture.sand >= 0.90
        PainterKarra(swrc=VanGenuchten(:sand))
    elseif texture.sand >= 0.50 && texture.clay < 0.20
        PainterKarra(swrc=VanGenuchten(:sandyloam))
    elseif texture.sand >= 0.50
        PainterKarra(swrc=VanGenuchten(:sandyclay))
    elseif texture.silt >= 0.90
        PainterKarra(swrc=VanGenuchten(:silt))
    elseif texture.silt >= 0.75
        PainterKarra(swrc=VanGenuchten(:siltloam))
    elseif texture.silt >= 0.40 && texture.clay >= 0.45
        PainterKarra(swrc=VanGenuchten(:siltyclay))
    elseif texture.clay >= 0.55
        PainterKarra(swrc=VanGenuchten(:clay))
    else
        # van Genuchten parameters for loam
        # Carsel and Parrish (1988)
        PainterKarra(swrc=VanGenuchten(α=0.036u"cm^-1", n=1.56))
    end
end
