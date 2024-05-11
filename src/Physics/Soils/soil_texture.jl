"""
    SoilTexture{Tsand,Tsilt,Tclay}

Represents soil "texture" as a simple mixture of sand, silt, and clay.
"""
Base.@kwdef struct SoilTexture{Tsand,Tsilt,Tclay}
    sand::Tsand = 1.0
    silt::Tsilt = 0.0
    clay::Tclay = 0.0
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

