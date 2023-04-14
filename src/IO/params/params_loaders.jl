abstract type ParameterFormat end
"""
JSON parameter input format (from CryoGridLite). Not yet implemented.
"""
struct ParamsJSON{Version} <: ParameterFormat end
"""
YAML parameter input format matching that of the CryoGrid community model. Not yet implemented.
"""
struct ParamsYAML{Version} <: ParameterFormat end

filesuffix(::ParamsJSON) = "json"
filesuffix(::ParamsYAML) = "yml"
