variables(top::Top, seb::SurfaceEnergyBalance) = (
    Diagnostic(:Sout, Scalar, u"W/(m^2)"),    # outgoing shortwave radiation [J/(s*m^2)]
    Diagnostic(:Lout, Scalar, u"W/(m^2)"),    # outgoing longwave radiation [J/(s*m^2)]
    Diagnostic(:Qnet, Scalar, u"W/(m^2)"),    # net radiation budget at surface [J/(s*m^2)]
    Diagnostic(:Qh, Scalar, u"W/(m^2)"),      # sensible heat flux [J/(s*m^2)]
    Diagnostic(:Qe, Scalar, u"W/(m^2)"),      # latent heat flux [J/(s*m^2)]
    Diagnostic(:Qg, Scalar, u"W/(m^2)"),      # ground heat flux [J/(s*m^2)]
    Diagnostic(:Lstar, Scalar, u"m"),         # Obukhov length [m]
    Diagnostic(:ustar, Scalar, u"m/s"),       # friction velocity [m/s]
)

function initialcondition!(top::Top, seb::SurfaceEnergyBalance, state)
    @setscalar state.Sout = 0.;
    @setscalar state.Lout = 0.;
    @setscalar state.Qnet = 0.;
    @setscalar state.Qh = 0.;
    @setscalar state.Qe = 0.;
    @setscalar state.Qg = 0.;
    @setscalar state.Lstar = -1e5;
    @setscalar state.ustar = 10.;
end

BoundaryStyle(::Type{<:SurfaceEnergyBalance}) = Neumann()

"""
Top interaction, ground heat flux from surface energy balance. (no snow, no water body, no infiltration)
"""
function boundaryvalue(seb::SurfaceEnergyBalance, ::Top, ::Heat, ::Soil, stop, ssoil)

    # TODO (optimize): pre-compute all forcings at time t here, then pass to functions

    # 1. calculate radiation budget
    # outgoing shortwave radiation as reflected
    @setscalar stop.Sout = let α = seb.sebparams.α,
        Sin = seb.forcings.Sin(stop.t);
        -α * Sin  # Eq. (2) in Westermann et al. (2016)
    end

    # outgoing longwave radiation composed of emitted and reflected radiation
    @setscalar stop.Lout = let ϵ = seb.sebparams.ϵ,
        σ = seb.sebparams.σ,
        T₀ = ssoil.T[1],
        Lin = seb.forcings.Lin(stop.t);
        -ϵ * σ * normalize_temperature(T₀)^4 - (1 - ϵ) * Lin # Eq. (3) in Westermann et al. (2016)
    end

    # net radiation budget
    @setscalar stop.Qnet = let Sin = seb.forcings.Sin(stop.t),
        Lin = seb.forcings.Lin(stop.t),
        Sout = getscalar(stop.Sout),
        Lout = getscalar(stop.Lout);
        Sin + Sout + Lin + Lout
    end

    # 2. calcuate turbulent heat flux budget
    # determine atmospheric stability conditions
    @setscalar stop.Lstar = Lstar(seb, stop, ssoil);
    @setscalar stop.ustar = ustar(seb, stop);

    # sensible heat flux
    @setscalar stop.Qh = Q_H(seb, stop, ssoil);

    # latent heat flux
    @setscalar stop.Qe = Q_E(seb, stop, ssoil);

    # 3. determine ground heat flux as the residual of the radiative and turbulent fluxes
    @setscalar stop.Qg = let Qnet = getscalar(stop.Qnet), Qₕ = getscalar(stop.Qh), Qₑ = getscalar(stop.Qe);
        Qnet - Qₕ - Qₑ # essentially Eq. (1) in Westermann et al. (2016)
    end

    # 4. return the ground heat flux to the uppermost soil grid cell
    return stop.Qg |> getscalar
end

"""
Density of air at given tempeature and pressure
"""
density_air(seb::SurfaceEnergyBalance,T, p) = p / (normalize_temperature(T) * seb.sebparams.Rₐ);

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
function ustar(seb::SurfaceEnergyBalance, stop)
    let κ = seb.sebparams.κ,
        uz = seb.forcings.wind(stop.t),                                           # wind speed at height z
        z = seb.forcings.z,                                                       # height z of wind forcing
        z₀ = seb.sebparams.z₀,                                                   # aerodynamic roughness length [m]
        Lstar = stop.Lstar |> getscalar;
        κ * uz ./ (log(z / z₀) - Ψ_M(seb, z / Lstar, z₀ / Lstar))                # Eq. (7) in Westermann et al. (2016)
    end
end

"""
Obukhov length according to Monin-Obukhov theory, iterative determination as in CryoGrid3 / Westermann et al. 2016
- uses the turubulent fluxes Qe and Qh as well as the friction velocity of the previous time step
"""
function Lstar(seb::SurfaceEnergyBalance{Iterative}, stop, ssoil)
    res = let κ = seb.sebparams.κ,
        g = seb.sebparams.g,
        Rₐ = seb.sebparams.Rₐ,
        cₚ = seb.sebparams.cₐ / seb.sebparams.ρₐ, # specific heat capacity of air at constant pressure
        Tair = seb.forcings.Tair(stop.t),
        Tₕ = normalize_temperature(Tair), # air temperature at height z over surface
        p = seb.forcings.p(stop.t), # atmospheric pressure at surface
        ustar = stop.ustar |> getscalar,
        Qₑ = stop.Qe |> getscalar,
        Qₕ = stop.Qh |> getscalar,
        Llg = L_lg(ssoil.T[1]),
        ρₐ = density_air(seb, Tair, p); # density of air at surface air temperature and surface pressure [kg/m^3]
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
function Lstar(seb::SurfaceEnergyBalance{Analytical}, stop, ssoil)
    res = let g = seb.sebparams.g,
            Rₐ = seb.sebparams.Rₐ,
            cₚ = seb.sebparams.cₐ / seb.sebparams.ρₐ,                             # specific heat capacity of air at constant pressure
            Tₕ = normalize_temperature(seb.forcings.Tair(stop.t)),                                        # air temperature at height z over surface
            T₀ = normalize_temperature(ssoil.T[1]),                                                      # surface temperature
            p = seb.forcings.p(stop.t),                                            # atmospheric pressure at surface (height z)
            p₀ = seb.forcings.p(stop.t),                                           # normal pressure (for now assumed to be equal to p)
            uz = seb.forcings.wind(stop.t),                                        # wind speed at height z
            z = seb.forcings.z,                                                    # height z of wind forcing
            z₀ = seb.sebparams.z₀,                                                # aerodynamic roughness length [m]
            Pr₀ = seb.sebparams.Pr₀,                                              # turbulent Prandtl number
            γₕ = seb.sebparams.γₕ,
            γₘ = seb.sebparams.γₘ,
            βₕ = seb.sebparams.βₕ,
            βₘ = seb.sebparams.βₘ;

        Θₕ = Tₕ * (p₀/p)^(Rₐ/cₚ);                                                 # potential temperature (for now identical to actual temperature)
        Θ₀ = T₀ * (p₀/p)^(Rₐ/cₚ);

        # calcuate bulk Richardson number
        Ri_b = g / Θ₀ * (Θₕ - Θ₀) * (z - z₀) / uz^2;                                # bulk Richardson number, eq. (9) in Byun 1990

        # calulate ζ
        a = (z / (z-z₀)) * log(z/z₀);
        if Ri_b>0 # stable conditions
            b = 2 * βₕ * (βₘ * Ri_b - 1);
            c = -(2 * βₕ * Ri_b - 1) - (1 + (4 * (βₕ - βₘ) * Ri_b) / Pr₀ )^0.5
            ζ = a / b * c;                                                        # eq. (19) in Byun 1990
        else # unstable conditions
            s_b = Ri_b / Pr₀;                                                     # eq. (30) in Byun 1990
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
function Q_H(seb::SurfaceEnergyBalance, stop, ssoil)
    let κ = seb.sebparams.κ,
        Rₐ = seb.sebparams.Rₐ,
        Tₕ = seb.forcings.Tair(stop.t),                                                # air temperature
        T₀ = ssoil.T[1],                                                              # surface temperature
        cₚ = seb.sebparams.cₐ / seb.sebparams.ρₐ,                                     # specific heat capacity of air at constant pressure
        z = seb.forcings.z,                                                            # height at which forcing data are provided
        Lstar = stop.Lstar |> getscalar,
        ustar = stop.ustar |> getscalar,
        p = seb.forcings.p(stop.t),
        z₀ = seb.sebparams.z₀,
        ρₐ = density_air(seb, Tₕ, p); # density of air at surface air temperature and surface pressure [kg/m^3]

        rₐᴴ = (κ * ustar)^-1 * (log(z / z₀) - Ψ_HW(seb, z / Lstar, z₀ / Lstar))            # Eq. (6) in Westermann et al. (2016)

        # calculate Q_H
        -ρₐ * cₚ * (Tₕ - T₀) / rₐᴴ                                                    # Eq. (4) in Westermann et al. (2016)
    end
end

"""
Latent heat flux, defined as positive if it is a flux towards the surface.
Represents evapo(transpi)ration/condensation at positive surface temperatures and
sublimation/resublimation at negative surface temperatures
"""
function Q_E(seb::SurfaceEnergyBalance, stop, ssoil)
    let κ = seb.sebparams.κ,
        γ = seb.sebparams.γ,
        Rₐ = seb.sebparams.Rₐ,
        Tₕ = seb.forcings.Tair(stop.t),                                                # air temperature at height z over surface
        T₀ = ssoil.T[1],                                                              # surface temperature
        p = seb.forcings.p(stop.t),                                                    # atmospheric pressure at surface
        qₕ = seb.forcings.q(stop.t),                                                   # specific humidity at height h over surface
        z = seb.forcings.z,                                                            # height at which forcing data are provided
        rₛ = seb.sebparams.rₛ,                                                        # surface resistance against evapotranspiration / sublimation [1/m]
        Lstar = stop.Lstar |> getscalar,
        ustar = stop.ustar |> getscalar,
        Llg = L_lg(ssoil.T[1]),
        Lsg = L_sg(ssoil.T[1]),
        z₀ = seb.sebparams.z₀,                                                        # aerodynamic roughness length [m]
        ρₐ = density_air(seb, Tₕ, seb.forcings.p(stop.t));       # density of air at surface air temperature and surface pressure [kg/m^3]

        q₀ = γ * estar(T₀) / p                                                        # saturation pressure of water/ice at the surface; Eq. (B1) in Westermann et al (2016)
        rₐᵂ = (κ * ustar)^-1 * (log(z / z₀) - Ψ_HW(seb, z / Lstar, z₀ / Lstar))       # aerodynamic resistance Eq. (6) in Westermann et al. (2016)
        L = (T₀ <= 0.0) ? Lsg : Llg                                                   # latent heat of sublimation/resublimation or evaporation/condensation [J/kg]

        # calculate Q_E
        res = (qₕ > q₀) ? -ρₐ * L * (qₕ - q₀) / (rₐᵂ)  :                              # Eq. (5) in Westermann et al. (2016) # condensation / deposition (no aerodynamics resistance)
                        -ρₐ * L * (qₕ - q₀) / (rₐᵂ + rₛ) ;                            # evaporation / sublimation (account for surface resistance against evapotranspiration/sublimation)

        res = (T₀ <= 0.0) ? 0 : res ;                                                 # for now: set sublimation and deposition to zero.
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
function Ψ_HW(seb::SurfaceEnergyBalance{T,Businger},ζ₁::Float64, ζ₂::Float64) where T
    let γₕ = seb.sebparams.γₕ,
        βₕ = seb.sebparams.βₕ;
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
function Ψ_M(seb::SurfaceEnergyBalance{T,Businger},ζ₁::Float64, ζ₂::Float64) where T
    let γₘ = seb.sebparams.γₘ,
        βₘ = seb.sebparams.βₘ;
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
