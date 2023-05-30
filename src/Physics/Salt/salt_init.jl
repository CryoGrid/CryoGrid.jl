Base.@kwdef struct SedimentCompactionInitializer <: VarInitializer{:por}
    porosityZero = 0.4
    porosityMin = 0.03
    # TODO: maybe consider moving some of the hardcoded constants below here?
end

function compaction!(por::AbstractVector, init::SedimentCompactionInitializer, grid::Grid)
    z = cells(grid)
    # shift to make z semipositive
    z = z .- minimum(z)
    porosityZero = init.porosityZero;
    bulkDensityZero = 1.0 ./ ((porosityZero .+ 0.6845) ./ 1.8)
    bulkDensity = bulkDensityZero .+ 0.0037 .* z .^ 0.766
    @. por = min(1.80 * bulkDensity ^ (-1.0) - 0.6845, init.porosityMin)
end

function CryoGrid.initialcondition!(::Soil, state, init::SedimentCompactionInitializer)
    compaction!(state.por, init, grid)
    return nothing
end
