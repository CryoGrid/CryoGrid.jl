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

(init::LinearTwoPhaseTempProfile)(::SubSurface, state) = init(state.T, state.grid)

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

Base.@kwdef struct ThermalSteadyStateInit{TT,TQ} <: CryoGrid.VarInitializer{:T}
    T0::TT = 0.0u"°C"
    Qgeo::TQ = 0.053u"W/m^2"
    maxiters::Int = 100
    thresh::Float64 = 1e-6
end

(init::ThermalSteadyStateInit)(::SubSurface, state) = nothing

function CryoGrid.initialcondition!(init::ThermalSteadyStateInit, top::Top, sub::SubSurface, stop, ssub)
    steadystate!(sub, ssub, init.T0, init.Qgeo, init.maxiters, init.thresh)
end

function CryoGrid.initialcondition!(init::ThermalSteadyStateInit, sub1::SubSurface, sub2::SubSurface, s1, s2)
    T_upper = s1.T[end] # bottom temperature from upper layer
    steadystate!(sub2, s2, T_upper, init.Qgeo, init.maxiters, init.thresh)
end

function steadystate!(sub::SubSurface, state, T0, Qgeo, maxiters::Int, convergence_thresh::Float64)
    i = 1
    z = cells(state.grid)
    ΔT = Inf
    while i < maxiters && ΔT > convergence_thresh
        T_new = T0 .+ z*Qgeo ./ state.kc
        ΔT = maximum(abs.(T_new .- state.T))
        state.T .= T_new
        # recompute initial condition and diagnostic state
        initialcondition!(sub, state)
        computediagnostic!(sub, state)
        i += 1
    end
    # return temperature profile
    return state.T
end
