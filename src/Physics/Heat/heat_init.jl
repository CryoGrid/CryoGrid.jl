"""
    LinearTwoPhaseInitialTempProfile{TT,Tz} <: VarInitializer{:T}

Simple, piecewise linear temprature initializer that uses two temperature values
and three characteristic depths to initialize the temperature profile. `T0` is the
initial surface temperature, `T1` is the permafrost temperature at `z_deep`,
`z_thaw` and `z_base` are the top and bottom freezing fronts which are both assumed
to be have temperature equal to `Tmelt`.
"""
Base.@kwdef struct LinearTwoPhaseInitialTempProfile{TT,Tz} <: VarInitializer{:T}
    T0::TT = 1.0u"°C"
    T1::TT = -10.0u"°C"
    Tmelt::TT = 0.0u"°C"
    z_top::Tz = 0.0u"m"
    z_deep::Tz = 20.0u"m"
    z_thaw::Tz = 0.5u"m"
    z_base::Tz = 500.0u"m"
end

# make initializer callable
"""
    (init::LinearTwoPhaseInitialTempProfile)(T::AbstractVector, grid::Grid)

Evaluate the initializer on the given temperature vector `T` which should match the length of
`cells(grid)`.
"""
function (init::LinearTwoPhaseInitialTempProfile)(T::AbstractVector, grid::Grid)
    @assert length(T) == length(cells(grid))
    T0 = ustrip(init.T0)
    T1 = ustrip(init.T1)
    Tm = ustrip(init.Tmelt)
    z_top = ustrip(init.z_top)
    z_deep = ustrip(init.z_deep)
    z_thaw = ustrip(init.z_thaw)
    z_base = ustrip(init.z_base)
    Ts = [T0, Tm, T1, Tm]
    zs = [z_top, z_thaw, z_deep, z_base]
    f = Interp.extrapolate(Interp.interpolate((zs,), Ts, Interp.Gridded(Interp.Linear())), Interp.Flat())
    # evaluate interpolant at grid cell midpoints
    T .= f.(cells(grid))
end

# implement initialcondition! for initializer type
function CryoGrid.initialcondition!(init!::LinearTwoPhaseInitialTempProfile, ::SubSurface, ::HeatBalance, state)
    init!(state.T, state.grid)
    return nothing
end
