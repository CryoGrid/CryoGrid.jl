Base.@kwdef struct SedimentCompactionInitializer{TZ,TM} <: VarInitializer{:por}
    porosityZero::TZ = 0.4
    porosityMin::TM = 0.03
    # TODO: maybe consider moving some of the hardcoded constants below here?
end

(init::SedimentCompactionInitializer)(::Layer, state) = compaction!(state.por, state.grid, init.porosityZero, init.porosityMin)

function compaction!(por::AbstractVector, grid::Grid, porosityZero, porosityMin)
    z = cells(grid)
    # shift to make z semipositive
    z = z .- minimum(z)
    porosityZero = porosityZero;
    bulkDensityZero = 1.0 ./ ((porosityZero .+ 0.6845) ./ 1.8)
    bulkDensity = bulkDensityZero .+ 0.0037 .* z .^ 0.766
    @. por = max(1.80 * bulkDensity ^ (-1.0) - 0.6845, porosityMin)
end
