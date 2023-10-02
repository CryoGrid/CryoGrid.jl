CryoGrid.variables(::Top, ::SurfaceEnergyBalance) = (
    Diagnostic(:Sout, Scalar, u"W/(m^2)"),    # outgoing shortwave radiation [J/(s*m^2)]
    Diagnostic(:Lout, Scalar, u"W/(m^2)"),    # outgoing longwave radiation [J/(s*m^2)]
    Diagnostic(:Qnet, Scalar, u"W/(m^2)"),    # net radiation budget at surface [J/(s*m^2)]
    Diagnostic(:Qh, Scalar, u"W/(m^2)"),      # sensible heat flux [J/(s*m^2)]
    Diagnostic(:Qe, Scalar, u"W/(m^2)"),      # latent heat flux [J/(s*m^2)]
    Diagnostic(:Qg, Scalar, u"W/(m^2)"),      # ground heat flux [J/(s*m^2)]
    Diagnostic(:Lstar, Scalar, u"m"),         # Obukhov length [m]
    Diagnostic(:ustar, Scalar, u"m/s"),       # friction velocity [m/s]
    Diagnostic(:T_ub, Scalar, u"°C"),         # air temperature
)

function CryoGrid.initialcondition!(top::Top, seb::SurfaceEnergyBalance, state)
    state.Sout .= zero(state.Sout);
    state.Lout .= zero(state.Lout);
    state.Qnet .= zero(state.Qnet);
    state.Qh .= zero(state.Qh);
    state.Qe .= zero(state.Qe);
    state.Qg .= zero(state.Qg);
    state.Lstar .= -1e5*oneunit(eltype(state.Lstar));
    state.ustar .= 10.0*oneunit(eltype(state.ustar));
    state.T_ub .= seb.forcings.Tair(state.t)
end

CryoGrid.BCKind(::Type{<:SurfaceEnergyBalance}) = CryoGrid.Neumann()

CryoGrid.boundaryvalue(::SurfaceEnergyBalance, state) = getscalar(state.Qg)

# interact! with soil layer
function CryoGrid.interact!(::Top, seb::SurfaceEnergyBalance, soil::Soil, ::HeatBalance, stop, ssoil)
    seb_output = updateseb!(seb, soil, stop, ssoil)
    ssoil.jH[1] += seb_output.Qg
    return nothing
end
# interact! with snowpack
function CryoGrid.interact!(::Top, seb::SurfaceEnergyBalance, snow::Snowpack, ::HeatBalance, stop, ssnow)
    seb_output = updateseb!(seb, snow, stop, ssnow)
    ssnow.jH[1] += seb_output.Qg
    @setscalar ssnow.T_ub = getscalar(stop.T_ub)
    return nothing
end

function updateseb!(seb::SurfaceEnergyBalance, sub::SubSurface, stop, ssub)
    state = SEBState(seb, sub, stop, ssub)
    seb_output = seb(state)
    # copy outputs into state variables
    @setscalar stop.Sout = seb_output.Sout
    @setscalar stop.Lout = seb_output.Lout
    @setscalar stop.Qnet = seb_output.Qnet
    @setscalar stop.Qg = seb_output.Qg
    @setscalar stop.Qh = seb_output.state.Qh
    @setscalar stop.Qe = seb_output.state.Qe
    @setscalar stop.Lstar = seb_output.state.Lstar
    @setscalar stop.ustar = seb_output.state.ustar
    # TODO: in the future, consider near surface air convection?
    @setscalar stop.T_ub = state.inputs.Tair
    return seb_output
end
function updateseb!(seb::SurfaceEnergyBalance{<:Numerical}, sub::SubSurface, stop, ssub)
    initialstate = SEBState(seb, sub, stop, ssub)
    seb_output = solve(seb, initialstate)
    # copy outputs into state variables
    @setscalar stop.Sout = seb_output.Sout
    @setscalar stop.Lout = seb_output.Lout
    @setscalar stop.Qnet = seb_output.Qnet
    @setscalar stop.Qg = seb_output.Qg
    @setscalar stop.Qh = seb_output.state.Qh
    @setscalar stop.Qe = seb_output.state.Qe
    @setscalar stop.Lstar = seb_output.state.Lstar
    @setscalar stop.ustar = seb_output.state.ustar
    # TODO: in the future, consider near surface air convection?
    @setscalar stop.T_ub = initialstate.inputs.Tair
    return seb_output
end

function (seb::SurfaceEnergyBalance)(state::SEBState)
    # 1. calculate radiation budget
    # outgoing shortwave radiation as reflected
    Sout = let α = state.inputs.surf.α,
        Sin = state.inputs.Sin;
        -α * Sin  # Eq. (2) in Westermann et al. (2016)
    end
    # outgoing longwave radiation composed of emitted and reflected radiation
    Lout = let ϵ = state.inputs.surf.ϵ,
        σ = seb.para.σ,
        T₀ = state.inputs.Ts,
        Lin = state.inputs.Lin;
        -ϵ * σ * normalize_temperature(T₀)^4 - (1 - ϵ) * Lin # Eq. (3) in Westermann et al. (2016)
    end
    # net radiation budget
    Qnet = let Sin = state.inputs.Sin,
        Lin = state.inputs.Lin;
        Sin + Sout + Lin + Lout
    end

    # 2. calcuate turbulent heat flux budget
    Qh, Qe, Lstar, ustar = turbulent_fluxes(seb, state)

    # 3. determine ground heat flux as the residual of the radiative and turbulent fluxes
    Qg = Qnet - Qh - Qe # essentially Eq. (1) in Westermann et al. (2016)
    newstate = SEBState(; Qh, Qe, Lstar, ustar, inputs=state.inputs)
    return SEBOutputs(; state=newstate, Qg, Qnet, Sout, Lout)
end

function turbulent_fluxes(seb::SurfaceEnergyBalance, state::SEBState)
    # determine atmospheric stability conditions
    Lstar = L_star(seb, state);
    ustar = u_star(seb, state);
    # sensible heat flux
    Qh = Q_H(seb, state);
    # latent heat flux
    Qe = Q_E(seb, state);
    return Qh, Qe, Lstar, ustar
end

function turbulent_fluxes(seb::SurfaceEnergyBalance{<:Iterative}, state::SEBState)
    # determine atmospheric stability conditions
    Lstar = L_star(seb, state);
    # update state for iterative scheme
    @set! state.Lstar = Lstar
    ustar = u_star(seb, state);
    @set! state.ustar = ustar
    # sensible heat flux
    Qh = Q_H(seb, state);
    # latent heat flux
    Qe = Q_E(seb, state);
    return Qh, Qe, Lstar, ustar
end

"""
Density of air at given tempeature and pressure
"""
density_air(seb::SurfaceEnergyBalance, Tair, pr) = pr / (normalize_temperature(Tair) * seb.para.Rₐ);

"""
Saturation pressure of water/ice according to the empirical August-Roche-Magnus formula
"""
estar(T) = (T > 0) ? 611.2 * exp(17.62 * T / (243.12 + T)) : 611.2 * exp(22.46 * T / (272.62 + T)); # Eq. (B3) in Westermann et al. (2016)

"""
Latent heat of evaporation/condensation of water in [J/kg]
according to https://en.wikipedia.org/wiki/Latent_heat#cite_note-RYfit-11
"""
L_lg(T) = 1000 * (2500.8 - 2.36 * T + 0.0016 * T^2 - 0.00006 * T^3);

"""
Latent heat of sublimation/resublimation of water in [J/kg]
accodring to https://en.wikipedia.org/wiki/Latent_heat#cite_note-RYfit-11
"""
L_sg(T) = 1000 * (2834.1 - 0.29 * T - 0.004 * T^2);

"""
Friction velocity according to Monin-Obukhov theory
"""
function u_star(seb::SurfaceEnergyBalance, state::SEBState)
    let κ = seb.para.κ,
        uz = state.inputs.wind, # wind speed at height z
        z = state.inputs.z, # height z of wind forcing
        z₀ = state.inputs.surf.z₀, # aerodynamic roughness length [m]
        Lstar = state.Lstar;
        κ * uz ./ (log(z / z₀) - Ψ_M(seb, z / Lstar, z₀ / Lstar)) # Eq. (7) in Westermann et al. (2016)
    end
end

"""
Obukhov length according to Monin-Obukhov theory, iterative determination as in CryoGrid3 / Westermann et al. 2016
- uses the turubulent fluxes Qe and Qh as well as the friction velocity of the previous time step
"""
function L_star(seb::SurfaceEnergyBalance, state::SEBState)
    res = let κ = seb.para.κ,
        g = seb.para.g,
        Rₐ = seb.para.Rₐ,
        cₚ = seb.para.cₐ / seb.para.ρₐ, # specific heat capacity of air at constant pressure
        Tair = state.inputs.Tair,
        Tₕ = normalize_temperature(Tair), # air temperature at height z over surface
        Tₛ = state.inputs.Ts,
        pr = state.inputs.pr, # atmospheric pressure at surface
        ustar = state.ustar,
        Qₑ = state.Qe,
        Qₕ = state.Qh,
        Llg = L_lg(ustrip(Tₛ)),
        ρₐ = density_air(seb, Tair, pr); # density of air at surface air temperature and surface pressure [kg/m^3]
        -ρₐ * cₚ * Tₕ / (κ * g) * ustar^3 / (Qₕ + 0.61 * cₚ / Llg * Tₕ * Qₑ) # Eq. (8) in Westermann et al. (2016)
    end
    # upper and lower limits for Lstar
    res = (abs(res) < 1e-7) ? sign(res) * 1e-7 : res
    res = (abs(res) > 1e+7) ? sign(res) * 1e+7 : res
end

"""
Obukhov length, analytical solution of Monin-Obukhov theory according to Byun 1990
- implicitly assumes the Businger 1971 stability functions
- only uses the current surface temperature as input
"""
function L_star(seb::SurfaceEnergyBalance{Analytical}, state::SEBState)
    res = let g = seb.para.g,
            Rₐ = seb.para.Rₐ,
            cₚ = seb.para.cₐ / seb.para.ρₐ, # specific heat capacity of air at constant pressure
            Tₕ = normalize_temperature(state.inputs.Tair), # air temperature at height z over surface
            T₀ = normalize_temperature(state.inputs.Ts), # surface temperature
            p = state.inputs.pr, # atmospheric pressure at surface (height z)
            p₀ = state.inputs.pr, # normal pressure (for now assumed to be equal to p)
            uz = state.inputs.wind, # wind speed at height z
            z = state.inputs.z, # height z of wind forcing
            z₀ = state.inputs.surf.z₀, # aerodynamic roughness length [m]
            Pr₀ = seb.para.Pr₀, # turbulent Prandtl number
            γₕ = seb.para.γₕ,
            γₘ = seb.para.γₘ,
            βₕ = seb.para.βₕ,
            βₘ = seb.para.βₘ;

        Θₕ = Tₕ * (p₀/p)^(Rₐ/cₚ); # potential temperature (for now identical to actual temperature)
        Θ₀ = T₀ * (p₀/p)^(Rₐ/cₚ);

        # calcuate bulk Richardson number
        Ri_b = g / Θ₀ * (Θₕ - Θ₀) * (z - z₀) / uz^2; # eq. (9) in Byun 1990

        # calulate ζ
        a = (z / (z-z₀)) * log(z/z₀);
        if Ri_b>0 # stable conditions
            b = 2 * βₕ * (βₘ * Ri_b - 1);
            c = -(2 * βₕ * Ri_b - 1) - (1 + (4 * (βₕ - βₘ) * Ri_b) / Pr₀ )^0.5
            ζ = a / b * c;                                                      # eq. (19) in Byun 1990
        else # unstable conditions
            s_b = Ri_b / Pr₀;                                                   # eq. (30) in Byun 1990
            Q_b = 1/9 * ( 1/γₘ^2 + 3 * γₕ/γₘ * s_b^2 );                         # eq. (31) in Byun 1990
            P_b = 1/54 * ( -2/γₘ^3 + 9/γₘ * (3 - γₕ/γₘ) * s_b^2 );              # eq. (32) in Byun 1990
            θ_b = acos( P_b / Q_b^(3/2) );                                      # eq. (33) in Byun 1990
            T_b = ( (P_b^2 - Q_b^3)^0.5 + abs(P_b) )^(1/3);                     # eq. (34) in Byun 1990
            if Q_b^3-P_b^2 >= 0
                ζ = a * (-2 * Q_b^0.5 * cos(θ_b/3) + 1/(3*γₘ) )                 # eq. (28) in Byun 1990
            else
                ζ = a * (-(T_b + Q_b/T_b ) + 1/(3*γₘ) );                        # eq. (29) in Byun 1990
            end
        end

        # calculate Lstar
        z / ζ;

    end
    # upper and lower limits for Lstar
    res = (abs(res) < 1e-7) ? sign(res) * 1e-7 : res
    res = (abs(res) > 1e+7) ? sign(res) * 1e+7 : res
end

"""
Sensible heat flux, defined as positive if it is a flux towards the surface
"""
function Q_H(seb::SurfaceEnergyBalance, state::SEBState)
    let κ = seb.para.κ,
        Rₐ = seb.para.Rₐ,
        Tₕ = state.inputs.Tair, # air temperature
        T₀ = state.inputs.Ts, # surface temperature
        cₚ = seb.para.cₐ / seb.para.ρₐ, # specific heat capacity of air at constant pressure
        z = state.inputs.z, # height at which forcing data are provided
        Lstar = state.Lstar,
        ustar = state.ustar,
        pr = state.inputs.pr,
        z₀ = state.inputs.surf.z₀,
        ρₐ = density_air(seb, Tₕ, pr); # density of air at surface air temperature and surface pressure [kg/m^3]

        rₐᴴ = (κ * ustar)^-1 * (log(z / z₀) - Ψ_HW(seb, z / Lstar, z₀ / Lstar)) # Eq. (6) in Westermann et al. (2016)

        # calculate Q_H; Eq. (4) in Westermann et al. (2016)
        -ρₐ * cₚ * (Tₕ - T₀) / rₐᴴ
    end
end

"""
Latent heat flux, defined as positive if it is a flux towards the surface.
Represents evapo(transpi)ration/condensation at positive surface temperatures and
sublimation/resublimation at negative surface temperatures
"""
function Q_E(seb::SurfaceEnergyBalance, state::SEBState)
    let κ = seb.para.κ,
        γ = seb.para.γ,
        Rₐ = seb.para.Rₐ,
        Tₕ = state.inputs.Tair, # air temperature at height z over surface
        T₀ = state.inputs.Ts, # surface temperature
        p = state.inputs.pr, # atmospheric pressure at surface
        qₕ = state.inputs.qh, # specific humidity at height h over surface
        z = state.inputs.z, # height at which forcing data are provided
        rₛ = state.inputs.surf.rₛ, # surface resistance against evapotranspiration / sublimation [1/m]
        Lstar = state.Lstar,
        ustar = state.ustar,
        Llg = L_lg(state.inputs.Ts),
        Lsg = L_sg(state.inputs.Ts),
        z₀ = state.inputs.surf.z₀, # aerodynamic roughness length [m]
        ρₐ = density_air(seb, Tₕ, state.inputs.pr); # density of air at surface air temperature and surface pressure [kg/m^3]

        q₀ = γ * estar(T₀) / p # saturation pressure of water/ice at the surface; Eq. (B1) in Westermann et al (2016)
        rₐᵂ = (κ * ustar)^-1 * (log(z / z₀) - Ψ_HW(seb, z / Lstar, z₀ / Lstar)) # aerodynamic resistance Eq. (6) in Westermann et al. (2016)
        L = (T₀ <= 0.0) ? Lsg : Llg # latent heat of sublimation/resublimation or evaporation/condensation [J/kg]

        # calculate Q_E
        res = if qₕ > q₀
            -ρₐ * L * (qₕ - q₀) / (rₐᵂ) # Eq. (5) in Westermann et al. (2016) # condensation / deposition (no aerodynamics resistance)
        else
            -ρₐ * L * (qₕ - q₀) / (rₐᵂ + rₛ); # evaporation / sublimation (account for surface resistance against evapotranspiration/sublimation)
        end
        
        # for now: set sublimation and deposition to zero.
        res = (T₀ <= 0.0) ? zero(res) : res;
    end
end

"""
Integrated stability function for heat/water transport
    Högström, 1988
    SHEBA, Uttal et al., 2002, Grachev et al. 2007
"""
function Ψ_HW(seb::SurfaceEnergyBalance{T,HøgstrømSHEBA}, ζ₁, ζ₂) where T
    if ζ₁ <= 0 # neutral and unstable conditions (according to Høgstrøm, 1988)
        # computed using WolframAlpha command "Integrate[ (1-(0.95*(1-11.6x)^(-1/2)))/x ] assuming x<0"
        res = real( log(Complex(ζ₁)) + 1.9 * atanh(Complex((1 - 11.6 * ζ₁)^0.5)) -
                   (log(Complex(ζ₂)) + 1.9 * atanh(Complex((1 - 11.6 * ζ₂)^0.5)) ) )
        # numerical integration using the QuadGK package
        # res, err = quadgk(x -> (1-(0.95*(1-11.6*x)^(-1/2)))/x, ζ₂, ζ₁, rtol=1e-6);
        return res
    else     # stable stratification (according to Grachev et al. 2007)
        # computed using WolframAlpha command "Integrate[ (1-(1+(5x*(1+x))/(1+3x+x^2)))/x ] assuming x>0"
        res = real( 0.5 * ((-5 + 5^0.5) * log(Complex(-3 + 5^0.5 - 2 * ζ₁)) - (5 + 5^0.5) * log(Complex(3 + 5^0.5 + 2 * ζ₁))) -
                    0.5 * ((-5 + 5^0.5) * log(Complex(-3 + 5^0.5 - 2 * ζ₂)) - (5 + 5^0.5) * log(Complex(3 + 5^0.5 + 2 * ζ₂)))  )
        # numerical integration using the QuadGK package
        # res, err = quadgk(x -> (1-(1+(5*x*(1+x))/(1+3*x+x^2)))/x, ζ₂, ζ₁, rtol=1e-6);
        return res
    end
end

"""
Integrated stability function for momentum transport
    Högström, 1988
    SHEBA, Uttal et al., 2002, Grachev et al. 2007
"""
function Ψ_M(seb::SurfaceEnergyBalance{T,HøgstrømSHEBA}, ζ₁, ζ₂) where T
    if ζ₁ <= 0 # neutral and unstable conditions (according to Høgstrøm, 1988)
        # computed using WolframAlpha command "Integrate[ (1-(1-19.3x)^(-1/4))/x ] assuming x<0"
        # log(ζ₁) - 2*atan((1-19.3*ζ₁)^(1/4)) + 2*atanh((1-19.3*ζ₁)^(1/4)) -
        # (log(ζ₂) - 2*atan((1-19.3*ζ₂)^(1/4)) + 2*atanh((1-19.3*ζ₂)^(1/4)) )
        # copied from CryoGrid MATLAB code (Note: Høgstrøm (1988) suggests phi_M=(1-19.3x)^(-1/4) while here (1-19.0x)^(-1/4) is used.)
        res = real( -2 * atan((1 - 19.3 * ζ₁)^(1 / 4)) + 2 * log(Complex(1 + (1 - 19.3 * ζ₁)^(1 / 4))) + log(Complex(1 + (1 - 19.3 * ζ₁)^0.5)) -
                   (-2 * atan((1 - 19.3 * ζ₂)^(1 / 4)) + 2 * log(Complex(1 + (1 - 19.3 * ζ₂)^(1 / 4))) + log(Complex(1 + (1 - 19.3 * ζ₂)^0.5))) )
        # numerical integration using the QuadGK package
        # res, err = quadgk(x -> (1-(1-19.3*x)^(-1/4))/x, ζ₂, ζ₁, rtol=1e-6);
        return res
    else     # stable stratification (according to Grachev et al. 2007)
        # computed using WolframAlpha command "Integrate[ (1-(1+6.5x*(1+x)^(1/3)/(1.3+x)))/x ] assuming x>0"
        # -19.5*(1 + ζ₁)^(1/3) - 7.5367*atan(0.57735 - 1.72489*(1+ζ₁)^(1/3)) + 4.35131*log(1.44225 + 2.15443*(1+ζ₁)^(1/3)) - 2.17566*log(2.08008 - 3.10723*(1+ζ₁)^(1/3) + 4.64159*(1+ζ₁)^(2/3)) -
        # (-19.5*(1 + ζ₂)^(1/3) - 7.5367*atan(0.57735 - 1.72489*(1+ζ₂)^(1/3)) + 4.35131*log(1.44225 + 2.15443*(1+ζ₂)^(1/3)) - 2.17566*log(2.08008 - 3.10723*(1+ζ₂)^(1/3) + 4.64159*(1+ζ₂)^(2/3)) )
        # copied from CryoGrid MATLAB code
        res = real( -19.5 * (1 + ζ₁)^(1 / 3) - 7.5367 * atan(0.57735 - 1.72489 * (1 + ζ₁)^(1 / 3)) + 4.35131 * log(3       + 4.4814 * (1 + ζ₁)^(1 / 3)) - 2.17566 * log(3       - 4.4814 * (1 + ζ₁)^(1 / 3) + 6.69433 * (1 + ζ₁)^(2 / 3)) -
                   (-19.5 * (1 + ζ₂)^(1 / 3) - 7.5367 * atan(0.57735 - 1.72489 * (1 + ζ₂)^(1 / 3)) + 4.35131 * log(3       + 4.4814 * (1 + ζ₂)^(1 / 3)) - 2.17566 * log(3       - 4.4814 * (1 + ζ₂)^(1 / 3) + 6.69433 * (1 + ζ₂)^(2 / 3))) )
        # numerical integration using the QuadGK package
        # res, err = quadgk(x -> (1-(1+6.5*x*(1+x)^(1/3)/(1.3+x)))/x, ζ₂, ζ₁, rtol=1e-6);
        return res
    end
end

"""
Integrated stability function for heat/water transport
    Businger 1971
"""
function Ψ_HW(seb::SurfaceEnergyBalance{T,Businger}, ζ₁::Float64, ζ₂::Float64) where T
    let γₕ = seb.para.γₕ,
        βₕ = seb.para.βₕ;
        if ζ₁ <= 0 # neutral and unstable conditions (according to Businger, 1971)
            return 2*log( ( (1-γₕ*ζ₁)^0.5 + 1 ) / ( (1-γₕ*ζ₂)^0.5 + 1 ) );      # eq. (15) in Byun 1990
        else     # stable stratification (according to Busigner, 1971)
            return -βₕ*(ζ₁-ζ₂);                                                 # eq. (13) in Byun 1990
        end
    end
end

"""
Integrated stability function for momentum transport
    Busigner 1971
"""
function Ψ_M(seb::SurfaceEnergyBalance{T,Businger}, ζ₁::Float64, ζ₂::Float64) where T
    let γₘ = seb.para.γₘ,
        βₘ = seb.para.βₘ;
        if ζ₁ <= 0 # neutral and unstable conditions (according to Businger, 1971)
            x=(1-γₘ*ζ₁)^0.25;
            x₀=(1-γₘ*ζ₂)^0.25;
            return 2*log( (1+x)/(1+x₀) ) + log( (1+x^2) / (1+x₀^2) ) -
                   2*atan(x) + 2*atan(x) ;                                      # eq. (14) in Byun 1990

        else     # stable stratification (according to Busigner, 1971)
            return -βₘ*(ζ₁-ζ₂);                                                 # eq. (12) in Byun 1990
        end
    end
end
