"""
    SnowpackParameterization

Base type for snowpack paramterization schemes.
"""
abstract type SnowpackParameterization <: CryoGrid.Parameterization end

"""
    SnowAblationScheme

Base type for different snow ablation (i.e. melting or redistribution) schemes.
"""
abstract type SnowAblationScheme end

"""
    SnowAccumulationScheme

Base type for different snow accumulation schemes.
"""
abstract type SnowAccumulationScheme end

"""
    SnowDensityScheme

Base type for different snow density schemes.
"""
abstract type SnowDensityScheme end


"""
    SnowThermalConductivity

Base type for snow thermal conductivity parameterizations.
"""
abstract type SnowThermalConductivity end
