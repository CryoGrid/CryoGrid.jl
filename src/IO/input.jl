"""
    Input{name,attrsType<:NamedTuple}

The `Input` type represents a placeholder for (typically time-varying)
inputs that will be replaced with the numeric value at runtime. The
`name` of the `Input` will be used to match it to the corresponding
input function.
"""
struct Input{name,attrsType<:NamedTuple}
    attrs::attrsType
    Input(name::Symbol, attrs::NamedTuple) = new{name,typeof(attrs)}(attrs)
    """
        Input(name::Symbol; kwargs...)

    Constructs a `Input` with the given identifier `name` and optional
    attribtues given as keyword arguments.
    """
    function Input(name::Symbol; kwargs...)
        attrs = (; kwargs...)
        new{name,typeof(attrs)}(attrs)
    end
end

ConstructionBase.constructorof(::Type{<:Input{name}}) where {name} = attrs -> Input(name, attrs)

"""
    nameof(::Input{name}) where {name}

Extract and return the `name` of the given `Input`.
"""
Base.nameof(::Input{name}) where {name} = name

function Base.getproperty(f::Input, name::Symbol)
    attrs = getfield(f, :attrs)
    if name ∈ keys(attrs)
        return getproperty(attrs, name)
    else
        return getfield(f, name)
    end
end

"""
    InputProvider

Base type for `Input` providers that realize one or more `Input`s.
"""
abstract type InputProvider end

"""
    inputmap(provider::InputProvider)

Returns a mapping-like object (e.g. a named tuple) with key/value pairs
corresponding to names and realized inputs.
"""
inputmap(::InputProvider) = error("not implemented")

Base.propertynames(provider::InputProvider) = propertynames(inputmap(provider))
Base.getproperty(provider::InputProvider, name::Symbol) = getproperty(inputmap(provider), name)
Base.keys(provider::InputProvider) = propertynames(provider)
Base.values(provider::InputProvider) = values(inputmap(provider))

"""
    InputFunctionProvider

Generic input provider that simply wraps a `NamedTuple` of functions.
"""
struct InputFunctionProvider{TF<:NamedTuple} <: InputProvider
    inputs::TF
end

InputFunctionProvider(; input_functions...) = InputFunctionProvider((; input_functions...))

function (ifp::InputFunctionProvider)(::Input{name}, args...)
    f = ifp.inputs[name]
    return f(args...)
end

inputmap(ifp::InputFunctionProvider) = getfield(ifp, :inputs)

"""
    inputs(; named_inputs...)

Alias for the constructor of `InputFunctionProvider`.
"""
inputs(; named_inputs...) = InputFunctionProvider(; named_inputs...)
