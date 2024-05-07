# Heat methods

function Heat.freezethaw!(
    soil::SalineGround,
    ps::CoupledHeatSalt{THeat},
    state
) where {THeat<:HeatBalance{<:DallAmicoSalt,<:TemperatureBased}}
    salt, heat = ps
    sfcc = heat.freezecurve
    let L = heat.prop.L;
        @inbounds @fastmath for i in eachindex(state.T)
            @unpack ch_w, ch_i = Heat.thermalproperties(soil, state, i)
            T = state.T[i]
            c = state.c[i]
            θsat = Soils.porosity(soil, state, i)
            sat = Soils.saturation(soil, state, i)
            x = @SVector[T, c]
            x_dual = Numerics.dual(x, typeof(sfcc))
            res_dual = sfcc(x_dual[1], sat, Val{:all}(); θsat=θsat, saltconc=x_dual[2])
            state.θw[i] = θw = ForwardDiff.value(res_dual.θw)
            state.∂θw∂T[i] = ∂θw∂T = ForwardDiff.partials(res_dual.θw)[1]
            state.∂θw∂c[i] = ForwardDiff.partials(res_dual.θw)[2]
            state.C[i] = C = heatcapacity(soil, state, i)
            state.∂H∂T[i] =  C + L*∂θw∂T
            state.H[i] = enthalpy(T, C, L, θw)
            state.dₛ_mid[i] = salt.prop.dₛ₀ * θw / salt.prop.τ
        end
    end
    return nothing
end

# CryoGrid methods

CryoGrid.variables(::SaltMassBalance) = (
    Prognostic(:c, OnGrid(Cells), u"mol/m^3"),
    Diagnostic(:jc, OnGrid(Edges), u"mol/m^3/s"),
    Diagnostic(:Tmelt, OnGrid(Cells), u"°C"),
    Diagnostic(:dₛ, OnGrid(Edges), u"m^2/s"),
    Diagnostic(:dₛ_mid, OnGrid(Cells), u"m^2/s"),
    Diagnostic(:∂θw∂c, OnGrid(Cells)),
    Diagnostic(:dc_F, OnGrid(Cells)),
    Diagnostic(:ctmp_B, OnGrid(Cells)),
    Diagnostic(:ctmp_E, OnGrid(Cells)),
    Diagnostic(:ctmp_F, OnGrid(Cells)),
)

function CryoGrid.initialcondition!(soil::SalineGround, ps::CoupledHeatSalt, state)
    CryoGrid.computediagnostic!(soil, ps, state)
end

function CryoGrid.computediagnostic!(
    soil::SalineGround,
    ps::CoupledHeatSalt{THeat},
    state
) where {THeat<:HeatBalance{<:SFCC,<:TemperatureBased}}
    # Reset energy flux to zero; this is redundant when H is the prognostic variable
    # but necessary when it is not.
    salt, heat = ps
    # Evaluate freeze/thaw processes
    freezethaw!(soil, ps, state)
    # Update thermal conductivity
    thermalconductivity!(soil, heat, state)
    # thermal conductivity at boundaries
    # assumes boundary conductivities = cell conductivities
    @inbounds state.k[1] = state.kc[1]
    @inbounds state.k[end] = state.kc[end]
    # Harmonic mean of inner conductivities
    @inbounds let k = (@view state.k[2:end-1]),
        Δk = Δ(state.grid);
        Numerics.harmonicmean!(k, state.kc, Δk)
    end
    # salt diffusivity at boundaries
    state.dₛ[1] = state.dₛ_mid[1]
    state.dₛ[end] = state.dₛ_mid[end]
    # Harmonic mean of inner conductivities
    @inbounds let dₛ = (@view state.dₛ[2:end-1]),
        Δdₛ = Δ(state.grid);
        Numerics.harmonicmean!(dₛ, state.dₛ_mid, Δdₛ)
    end
    return nothing # ensure no allocation
end

# interaction for two SalineGround layers
function CryoGrid.interact!(sediment1::SalineGround, sediment2::SalineGround, state1, state2)
    # water interaction
    interact!(sediment1, sediment1.water, sediment2, sediment2.water, state1, state2)
    # heat interaction
    interact!(sediment1, sediment1.heat, sediment2, sediment2.heat, state1, state2)
    # salt interaction
    interact!(sediment1, sediment1.salt, sediment2, sediment2.salt, state1, state2)
end

# interaction for salt mass balance
function CryoGrid.interact!(sediment1::SalineGround, ::SaltMassBalance, sediment2::SalineGround, ::SaltMassBalance, state1, state2)
    thick1 = CryoGrid.thickness(sediment1, state1, last)
    thick2 = CryoGrid.thickness(sediment2, state2, first)
    z1 = last(cells(state1.grid))
    z2 = first(cells(state2.grid))
    Δz = abs(z2 - z1)

    # calculate dₛ on shared edge
    state1.dₛ[end] = state2.dₛ[1] = Numerics.harmonicmean(state1.dₛ_mid[end], state2.dₛ_mid[1], thick1, thick2)

    # flux positive downward
    flux = -state1.dₛ[end] * (state2.c[1] - state1.c[end]) / Δz
    state1.jc[end] = state2.jc[1] = flux
    return nothing
end

function CryoGrid.timestep(::SalineGround, salt::SaltMassBalance{T,<:CryoGrid.CFL}, state) where {T}
    dtlim = salt.dtlim
    Δx = Δ(state.grid)
    Δc_max = dtlim.maxdelta.Δmax
    dtmax = Inf
    @inbounds for i in eachindex(Δx)
        dtmax = let v = 1 / state.dₛ[i],
            Δx = Δx[i],
            courant_number = dtlim.courant_number,
            Δt = courant_number*v*Δx^2;
            # minimum of CFL and maxium saltConc change per timestep
            min(min(dtmax, Δt), Δc_max / abs(state.dc[i]))
        end
    end
    return dtmax
end

function CryoGrid.computefluxes!(
    ::SalineGround,
    ps::CoupledHeatSalt{THeat},
    state
) where {THeat<:HeatBalance{<:SFCC,<:TemperatureBased}}
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
    Numerics.nonlineardiffusion!(state.dH, state.jH, T, midptThick, k, layerThick)

    #ion flux divergence: dc_F
    Numerics.nonlineardiffusion!(state.dc_F, state.jc, c, midptThick, dₛ, layerThick)

    #Derivative of water content
    ∂θw∂T = state.∂θw∂T
    ∂θw∂c = state.∂θw∂c
    θw = state.θw

    #put everything together
    A = state.∂H∂T
    B = state.ctmp_B .= L * ∂θw∂c
    D = state.dH
    E = state.ctmp_E .= (θw .+ c .* ∂θw∂c)
    F = state.ctmp_F .= c .* ∂θw∂T
    G = state.dc_F

    @. state.dT = (-B * G + D * E) / (A * E - B * F)
    @. state.dc = (-F * D + A * G) / (A * E - B * F)
    return nothing
end

function CryoGrid.resetfluxes!(::SalineGround, salt::SaltMassBalance, state)
    state.dc_F .= zero(eltype(state.dc))
    state.jc .= zero(eltype(state.dc))
    return nothing
end
