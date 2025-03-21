"""
    PermafrostTemperatureInit{TT1,TT2,TTb,TTm,Tz1,Tz2,Tz3,Tz4} <: VarInitializer{:T}

Simple, piecewise linear temprature initializer that uses three temperature values
and four characteristic depths to initialize the temperature profile. `T0` is the
initial surface temperature, `Tpf` is the permafrost temperature at `z_deep`,
`z_thaw` and `z_base` are the top and bottom freezing fronts which are both assumed
to be have temperature equal to `Tm`.
"""
Base.@kwdef struct PermafrostTemperatureInit{TT1,TT2,TTb,TTm,Tz1,Tz2,Tz3,Tz4} <: VarInitializer{:T}
    T0::TT1 = 1.0u"°C" # surface temperature
    Tpf::TT2 = -10.0u"°C" # permafrost temperature
    Tm::TTm = 0.0u"°C" # melting point
    Tbot::TTb = 10.0u"°C" # bottom temperature
    z0::Tz1 = 0.0u"m" # surface depth
    thaw_depth_offset::Tz3 = 0.5u"m" # thaw depth relative to z0
    perm_depth_offset::Tz2 = 20.0u"m" # depth offset of peramfrost temperature Tpf, relative to z_thaw
    base_depth_offset::Tz4 = 480.0u"m" # depth offset of permafrost base, relative to z_pf
end

(init::PermafrostTemperatureInit)(::SubSurface, state) = init(state.T, state.grid)

"""
    (init::PermafrostTemperatureInit)(T::AbstractVector, grid::Grid)

Evaluate the initializer on the given temperature vector `T` which should match the length of
`cells(grid)`.
"""
function (init::PermafrostTemperatureInit)(T::AbstractVector, grid::Grid)
    @assert length(T) == length(cells(grid))
    T0 = ustrip(init.T0)
    Tpf = ustrip(init.Tpf)
    Tbot = ustrip(init.Tbot)
    Tm = ustrip(init.Tm)
    z0 = ustrip(init.z0)
    z_thaw = ustrip(init.thaw_depth_offset)
    z_deep = z_thaw + ustrip(init.perm_depth_offset)
    z_base = z_deep + ustrip(init.base_depth_offset)
    z_bot = parent(grid)[end]
    Ts = [T0, Tm, Tpf, Tm, Tbot]
    zs = [z0, z_thaw, z_deep, z_base, z_bot]
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
        T_new = adstrip.(T0 .+ z*Qgeo ./ state.kc)
        ΔT = maximum(abs.(T_new .- state.T))
        state.T .= one(eltype(state.T))*T_new
        # recompute initial condition and diagnostic state variables
        initialcondition!(sub, state)
        computediagnostic!(sub, state)
        i += 1
    end
    # compute final temperature profile
    @. state.T = T0 + z*Qgeo / state.kc
    # recompute initial condition and diagnostic state variables
    initialcondition!(sub, state)
    computediagnostic!(sub, state)
    # return temperature profile
    return state.T
end
