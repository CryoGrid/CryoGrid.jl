variables(top::Top, seb::SurfaceEnergyBalance) = (
    Diagnostic(:Sout, Float"W/(m^2)", Scalar),    # outgoing shortwave radiation [J/(s*m^2)]
    Diagnostic(:Lout, Float"W/(m^2)", Scalar),    # outgoing longwave radiation [J/(s*m^2)]
    Diagnostic(:Qnet, Float"W/(m^2)", Scalar),    # net radiation budget at surface [J/(s*m^2)]
    Diagnostic(:Qh, Float"W/(m^2)", Scalar),      # sensible heat flux [J/(s*m^2)]
    Diagnostic(:Qe, Float"W/(m^2)", Scalar),      # latent heat flux [J/(s*m^2)]
    Diagnostic(:Qg, Float"W/(m^2)", Scalar),      # ground heat flux [J/(s*m^2)]
    Diagnostic(:Lstar, Float"m", Scalar),         # Obukhov length [m]
    Diagnostic(:ustar, Float"m/s", Scalar),       # friction velocity [m/s]
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

"""
Top interaction, ground heat flux from surface energy balance. (no snow, no water body, no infiltration)
"""
function (seb::SurfaceEnergyBalance)(top::Top, soil::Soil, heat::Heat, stop, ssoil)

    # TODO (optimize): pre-compute all forcings at time t here, then pass to functions

    # 1. calculate radiation budget
    # outgoing shortwave radiation as reflected
    @setscalar stop.Sout = let α = seb.sebparams.α, Sin = seb.forcing.Sin(stop.t);
        -α * Sin                                                        # Eq. (2) in Westermann et al. (2016)
    end

    # outgoing longwave radiation composed of emitted and reflected radiation
    @setscalar stop.Lout = let ϵ = seb.sebparams.ϵ, σ = seb.sebparams.σ, T₀ = ssoil.T[1], Lin = seb.forcing.Lin(stop.t);
        -ϵ * σ * T₀^4 - (1 - ϵ) * Lin # Eq. (3) in Westermann et al. (2016)
    end

    # net radiation budget
    @setscalar stop.Qnet = let Sin = seb.forcing.Sin(stop.t), Lin = seb.forcing.Lin(stop.t), Sout = getscalar(stop.Sout), Lout = getscalar(stop.Lout);
        Sin + Sout + Lin + Lout
    end

    # 2. calcuate turbulent heat flux budget
    # determine atmospheric stability conditions
    @setscalar stop.ustar = ustar(seb, stop);
    @setscalar stop.Lstar = Lstar(seb, stop, ssoil);

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
density_air(seb::SurfaceEnergyBalance,T::Real"K",p::Real"Pa") = p / (T * seb.sebparams.Rₐ);

"""
Saturation pressure of water/ice according to the empirical August-Roche-Magnus formula
Note: T is passed [K] and converted to [°C]
"""
estar(T::Real"K") = ( (T > 0) ? 611.2 * exp(17.62 * (T - 273.15) / (243.12 + (T - 273.15)))# Eq. (B3) in Westermann et al. (2016)
: 611.2 * exp(22.46 * (T - 273.15) / (272.62 + (T - 273.15))) ;
)

"""
Latent heat of evaporation/condensation of water in [J/kg]
according to https://en.wikipedia.org/wiki/Latent_heat#cite_note-RYfit-11
"""
L_lg(T::Real"K") = 1000 * (2500.8 - 2.36 * (T - 273.15) + 0.0016 * (T - 273.15)^2 - 0.00006 * (T - 273.15)^3);

"""
Latent heat of sublimation/resublimation of water in [J/kg]
accodring to https://en.wikipedia.org/wiki/Latent_heat#cite_note-RYfit-11
"""
L_sg(T::Real"K") = 1000 * (2834.1 - 0.29 * (T - 273.15) - 0.004 * (T - 273.15)^2);

"""
Friction velocity according to Monin-Obukhov theory
"""
function ustar(seb::SurfaceEnergyBalance, stop)
    let κ = seb.sebparams.κ,
        uz = seb.forcing.wind(stop.t),                                           # wind speed at height z
        z = seb.forcing.z,                                                       # height z of wind forcing
        z₀ = seb.sebparams.z₀,                                                   # aerodynamic roughness length [m]
        Lstar = stop.Lstar |> getscalar;
        κ * uz ./ (log(z / z₀) - Ψ_M(z / Lstar, z₀ / Lstar))                                 # Eq. (7) in Westermann et al. (2016)
    end
end

"""
Obukhov length according to Monin-Obukhov theory
"""
function Lstar(seb::SurfaceEnergyBalance, stop, ssoil)
    res = let κ = seb.sebparams.κ,
            g = seb.sebparams.g,
            Rₐ = seb.sebparams.Rₐ,
            cₚ = seb.sebparams.cₐ / seb.sebparams.ρₐ,                             # specific heat capacity of air at constant pressure
            Tₕ = seb.forcing.Tair(stop.t),                                        # air temperature at height z over surface
            p = seb.forcing.p(stop.t),                                            # atmospheric pressure at surface
            ustar = stop.ustar |> getscalar,
            Qₑ = stop.Qe |> getscalar,
            Qₕ = stop.Qh |> getscalar,
            Llg = L_lg(ssoil.T[1]),
            ρₐ = density_air(seb, Tₕ, p); # density of air at surface air temperature and surface pressure [kg/m^3]
        -ρₐ * cₚ * Tₕ / (κ * g) * ustar^3 / (Qₕ + 0.61 * cₚ / Llg * Tₕ * Qₑ)        # Eq. (8) in Westermann et al. (2016)
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
        Tₕ = seb.forcing.Tair(stop.t),                                                # air temperature
        T₀ = ssoil.T[1],                                                              # surface temperature
        cₚ = seb.sebparams.cₐ / seb.sebparams.ρₐ,                                     # specific heat capacity of air at constant pressure
        z = seb.forcing.z,                                                            # height at which forcing data are provided
        Lstar = stop.Lstar |> getscalar,
        ustar = stop.ustar |> getscalar,
        p = seb.forcing.p(stop.t),
        z₀ = seb.sebparams.z₀,
        ρₐ = density_air(seb, Tₕ, p); # density of air at surface air temperature and surface pressure [kg/m^3]

        rₐᴴ = (κ * ustar)^-1 * (log(z / z₀) - Ψ_HW(z / Lstar, z₀ / Lstar))            # Eq. (6) in Westermann et al. (2016)

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
        Tₕ = seb.forcing.Tair(stop.t),                                                # air temperature at height z over surface
        T₀ = ssoil.T[1],                                                              # surface temperature
        p = seb.forcing.p(stop.t),                                                    # atmospheric pressure at surface
        qₕ = seb.forcing.q(stop.t),                                                   # specific humidity at height h over surface
        z = seb.forcing.z,                                                            # height at which forcing data are provided
        rₛ = seb.sebparams.rₛ,                                                        # surface resistance against evapotranspiration / sublimation [1/m]
        Lstar = stop.Lstar |> getscalar,
        ustar = stop.ustar |> getscalar,
        Llg = L_lg(ssoil.T[1]),
        Lsg = L_sg(ssoil.T[1]),
        z₀ = seb.sebparams.z₀,                                                        # aerodynamic roughness length [m]
        ρₐ = density_air(seb, seb.forcing.Tair(stop.t), seb.forcing.p(stop.t));       # density of air at surface air temperature and surface pressure [kg/m^3]

        q₀ = γ * estar(T₀) / p                                                        # saturation pressure of water/ice at the surface; Eq. (B1) in Westermann et al (2016)
        rₐᵂ = (κ * ustar)^-1 * (log(z / z₀) - Ψ_HW(z / Lstar, z₀ / Lstar))            # aerodynamic resistance Eq. (6) in Westermann et al. (2016)
        L = (T₀ <= 273.15) ? Lsg : Llg                                                # latent heat of sublimation/resublimation or evaporation/condensation [J/kg]

        # calculate Q_E
        res = (qₕ > q₀) ? -ρₐ * L * (qₕ - q₀) / (rₐᵂ)  :                              # Eq. (5) in Westermann et al. (2016) # condensation / deposition (no aerodynamics resistance)
                        -ρₐ * L * (qₕ - q₀) / (rₐᵂ + rₛ) ;                            # evaporation / sublimation (account for surface resistance against evapotranspiration/sublimation)

        res = (T₀ <= 273.15) ? 0 : res ;                                              # for now: set sublimation and deposition to zero.
    end
end

"""
Integrated stability function for heat/water transport
    Högström, 1988
    SHEBA, Uttal et al., 2002, Grachev et al. 2007
"""
function Ψ_HW(ζ₁::Real, ζ₂::Real)
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
"""
function Ψ_M(ζ₁::Real, ζ₂::Real)
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

export SurfaceEnergyBalance
