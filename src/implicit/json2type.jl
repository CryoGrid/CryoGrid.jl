module json2type

function typenarrow!(d::Dict)
    for k in keys(d)
        if d[k] isa Array{Any,1}
            d[k] = typenarrow(d[k])
        elseif d[k] isa Dict
            typenarrow!(d[k])
        end
    end
    return d
end

function typenarrow(v::Array{Any,1})
    for T in [String,Float64,Int64,Bool,Vector{Float64}]
        try
            #return( convert(Vector{T},v) )
            return collect(Iterators.flatten(Iterators.flatten(convert(Vector{T},v))))
        catch; end
    end
    return v
end

end
