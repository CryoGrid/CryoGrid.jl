"""
    LinearTwoPhaseTempProfile{TT1,TT2,TTb,TTm,Tz1,Tz2,Tz3,Tz4} <: VarInitializer{:T}

Simple, piecewise linear temprature initializer that uses three temperature values
and four characteristic depths to initialize the temperature profile. `T0` is the
initial surface temperature, `T1` is the permafrost temperature at `z_deep`,
`z_thaw` and `z_base` are the top and bottom freezing fronts which are both assumed
to be have temperature equal to `Tm`.
"""
Base.@kwdef struct LinearTwoPhaseTempProfile{TT1,TT2,TTb,TTm,Tz1,Tz2,Tz3,Tz4} <: VarInitializer{:T}
    T0::TT1 = 1.0u"°C"
    T1::TT2 = -10.0u"°C"
    Tm::TTm = 0.0u"°C"
    Tbot::TTb = 10.0u"°C"
    z_top::Tz1 = 0.0u"m"
    z_deep::Tz2 = 20.0u"m"
    z_thaw::Tz3 = 0.5u"m"
    z_base::Tz4 = 500.0u"m"
end

# make initializer callable
"""
    (init::LinearTwoPhaseTempProfile)(T::AbstractVector, grid::Grid)

Evaluate the initializer on the given temperature vector `T` which should match the length of
`cells(grid)`.
"""
function (init::LinearTwoPhaseTempProfile)(T::AbstractVector, grid::Grid)
    @assert length(T) == length(cells(grid))
    T0 = ustrip(init.T0)
    T1 = ustrip(init.T1)
    Tbot = ustrip(init.Tbot)
    Tm = ustrip(init.Tm)
    z_top = ustrip(init.z_top)
    z_deep = ustrip(init.z_deep)
    z_thaw = ustrip(init.z_thaw)
    z_base = ustrip(init.z_base)
    z_bot = parent(grid)[end]
    Ts = [T0, Tm, T1, Tm, Tbot]
    zs = [z_top, z_thaw, z_deep, z_base, z_bot]
    f = Interp.extrapolate(Interp.interpolate((zs,), Ts, Interp.Gridded(Interp.Linear())), Interp.Flat())
    # evaluate interpolant at grid cell midpoints
    T .= f.(cells(grid))
end

# implement initialcondition! for initializer type
function CryoGrid.initialcondition!(init!::LinearTwoPhaseInitialTempProfile, ::SubSurface, ::HeatBalance, state)
    init!(state.T, state.grid)
    return nothing
end
