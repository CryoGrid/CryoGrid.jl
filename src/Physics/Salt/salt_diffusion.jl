function resetfluxes!(::MarineSediment, salt::SaltMassBalance, state)
    state.dc_F .= zero(eltype(state.∂c∂t))
    state.jc .= zero(eltype(state.∂c∂t))
end

# Heat methods

function Heat.freezethaw!(sediment::MarineSediment, ps::Coupled2{<:SaltMassBalance, THeat}, state) where {THeat<:HeatBalance{<:DallAmicoSalt, Temperature}}
    ∇(f, x) = ∇(typeof(x), f, x)
    function ∇(::Type{T}, f, x) where {T}
        dθw = ForwardDiff.gradient(f,x)#!(ForwardDiff.DiffResult(zero(T), x), f, x)
        return dθw
    end

    salt, heat = ps
    sfcc = freezecurve(heat)
    let L = heat.prop.L,
        f = sfcc.f;
        @inbounds @fastmath for i in 1:length(state.T)
            θfracs = volumetricfractions(sediment, state, i)
            T = state.T[i]
            c = state.c[i]
            x = @SVector[T, c]
            θw = f(x[1]; θsat=state.θp[i], θtot=state.θwi[i], saltconc=x[2])
            ∇θw = ∇(x -> f(x[1]; θsat=state.θp[i], θtot=state.θwi[i], saltconc=x[2]), x)
            state.θw[i] = θw
            state.∂θw∂T[i] = ∇θw[1]
            state.∂θw∂c[i] = ∇θw[2]
            state.C[i] = C = heatcapacity(sediment, heat, θfracs...)
            state.∂H∂T[i] = Heat.C_eff(T, C, L, ∇θw[1], heat.prop.cw, heat.prop.ci)
            state.H[i] = enthalpy(T, C, L, θw)
            state.dₛ[i] = salt.prop.dₛ₀ * θw / salt.prop.τ
        end
    end
end

# CryoGrid methods

CryoGrid.processes(sediment::MarineSediment) = Coupled(sediment.heat, sediment.salt)
CryoGrid.processes(sediment::MarineSediment{Tpara,Theat,Tsalt,<:WaterBalance}) where {Tpara,Theat,Tsalt} = Coupled(sediment.water, sediment.salt, sediment.heat)

CryoGrid.variables(sediment::MarineSediment, ps::Coupled2{<:SaltMassBalance, THeat}) where {THeat<:HeatBalance{<:SFCC, Temperature}} = (
    variables(sediment, ps[2])...,  # all the variables from heat
    Prognostic(:c, Numerics.OnGrid(Cells), u"mol/m^3"),
    Diagnostic(:jc, Numerics.OnGrid(Edges), u"mol/m^3/s"),
    Diagnostic(:Tmelt, Numerics.OnGrid(Cells), u"°C"),
    Diagnostic(:dₛ, Numerics.OnGrid(Cells), u"m^2/s"),
    Diagnostic(:∂θw∂c, Numerics.OnGrid(Cells)),
    Diagnostic(:dc_F, Numerics.OnGrid(Cells))
)

CryoGrid.initialcondition!(sediment::MarineSediment, ps::Coupled2{<:SaltMassBalance, THeat}, state) where {THeat<:HeatBalance{<:SFCC, Temperature}} = CryoGrid.diagnosticstep!(sediment, ps, state)

function CryoGrid.updatestate!(sediment::MarineSediment, ps::Coupled2{<:SaltMassBalance, THeat}, state) where {THeat<:HeatBalance{<:SFCC, Temperature}}
    salt, heat = ps
    Heat.resetfluxes!(sediment, heat, state)
    resetfluxes!(sediment, salt, state)
    # Evaluate freeze/thaw processes
    freezethaw!(sediment, ps, state)
    # Update thermal conductivity
    thermalconductivity!(sediment, heat, state)
    # thermal conductivity at boundaries
    # assumes boundary conductivities = cell conductivities
    @inbounds state.k[1] = state.kc[1]
    @inbounds state.k[end] = state.kc[end]
    # Harmonic mean of inner conductivities
    @inbounds let k = (@view state.k[2:end-1]),
        Δk = Δ(state.grids.k);
        Numerics.harmonicmean!(k, state.kc, Δk)
    end
    return nothing # ensure no allocation
end

function CryoGrid.interact!(sediment1::MarineSediment, ::SaltMassBalance, sediment2::MarineSediment, ::SaltMassBalance, state1, state2)
    z₁ = CryoGrid.midpoint(sediment1, state1, last)
    z₂ = CryoGrid.midpoint(sediment2, state2, first)
    δ = z₂ - z₁

    Δ₁ = CryoGrid.thickness(sediment1, state1, last)
    Δ₂ = CryoGrid.thickness(sediment2, state2, first)

    flux = state1.dₛ[end] * (state2.c[1] - state1.c[end]) / abs(δ)

    state1.∂c∂t[end] = + flux / Δ₁
    state2.∂c∂t[1] = - flux / Δ₂
end

function CryoGrid.timestep(::MarineSediment, salt::SaltMassBalance{T, <:CryoGrid.CFL}, state) where {T}
    dtmax = Inf
    @inbounds for i in 1:length(state.sat)
        dt = water.dtlim(state.∂c∂t[i], state.c[i], state.t)
        dt = isfinite(dt) && dt > zero(dt) ? dt : Inf # make sure it's +Inf
        dtmax = min(dtmax, dt)
    end
    return dtmax
end

function CryoGrid.computefluxes!(::MarineSediment, ps::Coupled2{<:SaltMassBalance, THeat}, state) where {THeat<:HeatBalance{<:SFCC, Temperature}}
    salt, heat = ps

    #read current state
    T = state.T;
    c = state.c;

    #read constants and grid;
    k = state.k
    L = heat.prop.L
    dₛ = state.dₛ

    layerThick = Δ(state.grids.k)
    midptThick = Δ(state.grids.c)

    #heat flux divergence: dT_F = dH
    Numerics.nonlineardiffusion!(state.∂H∂t, state.jH, T, midptThick, k, layerThick)

    #ion flux divergence: dc_F
    Numerics.nonlineardiffusion!(state.dc_F, state.jc, c, midptThick, dₛ, layerThick)

    #Derivative of water content
    dθwdT = state.∂θw∂T
    ∂θw∂c = state.∂θw∂c
    θw = state.θw
    
    #put everything together
    A = state.∂H∂T
    B = L * ∂θw∂c
    D = state.∂H∂t
    E = (θw .+ c .* ∂θw∂c)
    F = c .* dθwdT
    G = state.dc_F

    @. state.∂T∂t = (-B * G + D * E) / (A * E - B * F)
    @. state.∂c∂t = (-F * D + A * G) / (A * E - B * F)
end
