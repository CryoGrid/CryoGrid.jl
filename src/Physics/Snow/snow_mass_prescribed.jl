abstract type PrescribedSnowMassBalance <: SnowMassBalance end

"""
    PrescribedSWE{Tswe} <: SnowMassBalance

"Prescribed" snow mass balance, i.e. where the snow water equivalent is given as a constant or forcing.
"""
Base.@kwdef struct PrescribedSWE{Tswe<:Forcing{u"m"}} <: PrescribedSnowMassBalance
    swe::Tswe = ConstantForcing(0.0u"m") # depth snow water equivalent [m]
end

"""
    PrescribedSnowDepth{Tswe} <: SnowMassBalance

"Prescribed" snow mass balance, i.e. where the snow water equivalent is given as a constant or forcing.
"""
Base.@kwdef struct PrescribedSnowDepth{Tdsn<:Forcing{u"m"}} <: PrescribedSnowMassBalance
    dsn::Tdsn = ConstantForcing(0.0u"m") # real depth [m]
end

snowwater(::Snowpack, smb::PrescribedSWE{<:Forcing{u"m"}}, state) = smb.para.snowwater(state.t)

snowdepth(::Snowpack, smb::PrescribedSnowDepth{<:Forcing{u"m"}}, state) = smb.para.dsn(state.t)

CryoGrid.Volume(::Type{<:Snowpack{T,<:PrescribedSnowMassBalance}}) where {T} = CryoGrid.DiagnosticVolume()
