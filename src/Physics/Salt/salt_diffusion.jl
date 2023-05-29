function resetfluxes!(::MarineSediment, salt::SaltMassBalance, state)
    state.dc_F .= zero(eltype(state.∂c∂t))
    state.jc .= zero(eltype(state.∂c∂t))
end

# Heat methods

function Heat.freezethaw!(
    sediment::MarineSediment,
    ps::CoupledHeatSalt{THeat},
    state
) where {THeat<:HeatBalance{<:DallAmicoSalt,<:Temperature}}
    ∇(f, x) = ∇(typeof(x), f, x)
    function ∇(::Type{T}, f, x) where {T}
        dθw = Numerics.ForwardDiff.gradient(f,x)
        return dθw
    end

    salt, heat = ps
    sfcc = heat.freezecurve
    thermalprops = Heat.thermalproperties(sediment)
    @unpack ch_w, ch_i = thermalprops
    let L = heat.prop.L;
        @inbounds @fastmath for i in 1:length(state.T)
            θfracs = volumetricfractions(sediment, state, i)
            T = state.T[i]
            c = state.c[i]
            θsat = Soils.porosity(sediment, state, i)
            sat = Soils.saturation(sediment, state, i)
            x = @SVector[T, c]
            θw = sfcc(x[1], sat; θsat=θsat, saltconc=x[2])
            ∇θw = ∇(x -> sfcc(x[1], sat; θsat=θsat, saltconc=x[2]), x)
            state.θw[i] = θw
            state.∂θw∂T[i] = ∇θw[1]
            state.∂θw∂c[i] = ∇θw[2]
            state.C[i] = C = heatcapacity(sediment, heat, θfracs...)
            state.∂H∂T[i] = Heat.dHdT(T, C, L, ∇θw[1], ch_w, ch_i)
            state.H[i] = enthalpy(T, C, L, θw)
            state.dₛ[i] = salt.prop.dₛ₀ * θw / salt.prop.τ
        end
    end
end

# CryoGrid methods

CryoGrid.variables(::SaltMassBalance) = (
    Prognostic(:c, OnGrid(Cells), u"mol/m^3"),
    Diagnostic(:jc, OnGrid(Edges), u"mol/m^3/s"),
    Diagnostic(:Tmelt, OnGrid(Cells), u"°C"),
    Diagnostic(:dₛ, OnGrid(Cells), u"m^2/s"),
    Diagnostic(:∂θw∂c, OnGrid(Cells)),
    Diagnostic(:dc_F, OnGrid(Cells)),
)

function CryoGrid.initialcondition!(sediment::MarineSediment, ps::CoupledHeatSalt, state)
    CryoGrid.updatestate!(sediment, ps, state)
end

function CryoGrid.updatestate!(
    sediment::MarineSediment,
    ps::CoupledHeatSalt{THeat},
    state
) where {THeat<:HeatBalance{<:SFCC,<:Temperature}}
    salt, heat = ps
    Heat.resetfluxes!(sediment, heat, state)
    resetfluxes!(sediment, salt, state)
    # Evaluate freeze/thaw processes
    freezethaw!(sediment, Coupled(salt, heat), state)
    # Update thermal conductivity
    thermalconductivity!(sediment, heat, state)
    # thermal conductivity at boundaries
    # assumes boundary conductivities = cell conductivities
    @inbounds state.k[1] = state.kc[1]
    @inbounds state.k[end] = state.kc[end]
    # Harmonic mean of inner conductivities
    @inbounds let k = (@view state.k[2:end-1]),
        Δk = Δ(state.grid);
        Numerics.harmonicmean!(k, state.kc, Δk)
    end
    return nothing # ensure no allocation
end

# interaction for two MarineSediment layers
function CryoGrid.interact!(sediment1::MarineSediment, sediment2::MarineSediment, state1, state2)
    # water interaction
    interact!(sediment1, sediment1.water, sediment2, sediment2.water, state1, state2)
    # heat interaction
    interact!(sediment1, sediment1.heat, sediment2, sediment2.heat, state1, state2)
    # salt interaction
    interact!(sediment1, sediment1.salt, sediment2, sediment2.salt, state1, state2)
end

# interaction for salt mass balance
function CryoGrid.interact!(sediment1::MarineSediment, ::SaltMassBalance, sediment2::MarineSediment, ::SaltMassBalance, state1, state2)
    z₁ = CryoGrid.midpoint(sediment1, state1, last)
    z₂ = CryoGrid.midpoint(sediment2, state2, first)
    δ = z₂ - z₁

    # flux positive downward
    flux = -state1.dₛ[end] * (state2.c[1] - state1.c[end]) / abs(δ)

    state1.jc[end] = state2.jc[1] = flux
end

function CryoGrid.timestep(::MarineSediment, salt::SaltMassBalance{T,<:CryoGrid.CFL}, state) where {T}
    dtmax = Inf
    @inbounds for i in 1:length(state.sat)
        dt = water.dtlim(state.∂c∂t[i], state.c[i], state.t)
        dt = isfinite(dt) && dt > zero(dt) ? dt : Inf # make sure it's +Inf
        dtmax = min(dtmax, dt)
    end
    return dtmax
end

function CryoGrid.computefluxes!(
    ::MarineSediment,
    ps::CoupledHeatSalt{THeat},
    state
) where {THeat<:HeatBalance{<:SFCC,<:Temperature}}
    salt, heat = ps

    #read current state
    T = state.T;
    c = state.c;

    #read constants and grid;
    k = state.k
    L = heat.prop.L
    dₛ = state.dₛ

    layerThick = Δ(state.grid)
    midptThick = Δ(cells(state.grid))

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
