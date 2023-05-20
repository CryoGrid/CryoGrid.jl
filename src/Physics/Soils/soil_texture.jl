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
# Constants for common soil textures
const Sand = SoilTexture(sand=1.0, silt=0.0, clay=0.0)
const Silt = SoilTexture(sand=0.0, silt=1.0, clay=0.0)
const Clay = SoilTexture(sand=0.0, silt=0.0, clay=1.0)
const SandyClay = SoilTexture(sand=0.50, silt=0.0, clay=0.50)
const SiltyClay = SoilTexture(sand=0.0, silt=0.50, clay=0.50)
const Loam = SoilTexture(sand=0.40, silt=0.40, clay=0.20)
const SandyLoam = SoilTexture(sand=0.80, silt=0.10, clay=0.10)
const SiltyLoam = SoilTexture(sand=0.10, silt=0.80, clay=0.10)
const ClayLoam = SoilTexture(sand=0.30, silt=0.30, clay=0.40)
