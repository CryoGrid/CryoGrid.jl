"""
    param_ids(ps::CryoGridParams)

Constructs semi-unique identifiers for all parameters in `ps`. Ids are constructed from the field and stratigraphy layer names.
"""
function param_ids(ps::CryoGridParams)
    namestr(layer, name) = "$layer.$name"
    namestr(layer, parafield, name) = "$layer.$parafield.$name"
    condense(str) = join(filter(!=("nothing"), split(str, ".")), ".")
    pnames = if haskey(ps.table, :parafield)
        map(condense ∘ namestr, ps[:layer], ps[:parafield], ps[:fieldname])
    else
        map(condense ∘ namestr, ps[:layer], ps[:fieldname])
    end
    counts = counter(pnames) # count instances
    acc = counter(eltype(pnames)) # accumulator for labeling duplicated names
    return map(pnames) do name
        if counts[name] > 1
            inc!(acc, name) # increment accumulator
            Symbol(name, "[$(acc[name])]")
        else
            Symbol(name)
        end
    end
end
