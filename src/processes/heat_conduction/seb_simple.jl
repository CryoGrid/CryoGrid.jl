@with_kw struct SEBParams{} <: Params
    # surface properties --> should be associated with the Stratigraphy and maybe made state variables
    α::Float"1" = 0.2xu"1"                          # surface albedo [-]
    ϵ::Float"1" = 0.97xu"1"                         # surface emissivity [-]
    z₀::Float"m" = 1e-3xu"m"                        # surface roughness length [m]
    rₛ::Float"1/m" = 50.0xu"1/m"                    # surface resistance against evapotranspiration and sublimation [1/m]

    # "natural" constant
    σ::Float"J/(s*m^2*K^4)" = 5.6704e-8xu"J/(s*m^2*K^4)"   # Stefan-Boltzmann constant
    κ::Float"1" = 0.4xu"1"                          # von Kármán constant [-]
    γ::Float"1" = 0.622xu"1"                        # Psychrometric constant [-]
    Rₐ::Float"J/(kg*K)" = 287.058xu"J/(kg*K)"       # specific gas constant of air [J/(kg*K)]
    g::Float"m/s^2" = 9.81xu"m/s^2"                 # gravitational acceleration [m/s^2]

    # material properties (assumed to be constant)
    ρₐ::Float"kg/m^3" = 1.293xu"kg/m^3"             # density of air at standard pressure and 0°C [kg/m^3]
    cₐ::Float"J/(m^3*K)"= 1005.7xu"J/(kg*K)"*ρₐ     # volumetric heat capacity of dry air at standard pressure and 0°C [J/(m^3*K)]
end

struct SurfaceEnergyBalance{F} <: BoundaryProcess{Heat}
    forcing::F
    sebparams::SEBParams
    SurfaceEnergyBalance(   Tair::TimeSeriesForcing, p::TimeSeriesForcing, q::TimeSeriesForcing,
                            wind::TimeSeriesForcing, Lin::TimeSeriesForcing,
                            Sin::TimeSeriesForcing, z::Float"m") =
                            begin
                                sebparams = SEBParams()
                                forcing = (Tair=Tair,p=p,q=q,wind=wind,Lin=Lin,Sin=Sin,z=z)
                                new{typeof(forcing)}(forcing)
                            end

end

BoundaryStyle(::Type{<:SurfaceEnergyBalance}) = Neumann()

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
    state.Sout[1] = 0.;
    state.Lout[1] = 0.;
    state.Qnet[1] = 0.;
    state.Qh[1] = 0.;
    state.Qe[1] = 0.;
    state.Qg[1] = 0.;
    state.Lstar[1] = -1e5;
    state.ustar[1] = 10.;
end

"""
Top interaction, ground heat flux from surface energy balance. (no snow, no water body, no infiltration)
"""
function (seb::SurfaceEnergyBalance)(top::Top, soil::Soil, heat::Heat, stop, ssoil)

    # 1. calculate radiation budget
    # outgoing shortwave radiation as reflected
    stop.Sout[1] = let α=seb.sebparams.α, Sin=seb.forcing.Sin(stop.t);
                -α * Sin                                                    # Eq. (2) in Westermann et al. (2016)
    end

    # outgoing longwave radiation composed of emitted and reflected radiation
    stop.Lout[1] = let ϵ=seb.sebparams.ϵ, σ=seb.sebparams.σ, T₀=ssoil.T[1], Lin=seb.forcing.Lin(stop.t);
                -ϵ * σ * T₀^4 - (1-ϵ) * Lin                                  # Eq. (3) in Westermann et al. (2016)
    end

    # net radiation budget
    stop.Qnet[1] = let Sin=seb.forcing.Sin(stop.t), Lin=seb.forcing.Lin(stop.t), Sout=stop.Sout[1], Lout=stop.Lout[1];
                Sin + Sout + Lin + Lout
    end

    # 2. calcuate turbulent heat flux budget
    # determine atmospheric stability conditions
    stop.ustar[1] = ustar(seb, stop);
    stop.Lstar[1] = Lstar(seb, stop, ssoil);

    # sensible heat flux
    stop.Qh[1] = Q_H(seb, stop, ssoil);

    # latent heat flux
    stop.Qe[1] = Q_E(seb, stop, ssoil);

    # 3. determine ground heat flux as the residual of the radiative and turbulent fluxes
    stop.Qg[1] = let Qnet=stop.Qnet[1], Qₕ=stop.Qh[1], Qₑ=stop.Qe[1];
                    Qnet - Qₕ - Qₑ                                                  # essentially Eq. (1) in Westermann et al. (2016)
    end

    # 4. return the ground heat flux to the uppermost soil grid cell
    return stop.Qg[1]
end

"""
Density of air at giben tempeature and pressure
"""
density_air(seb::SurfaceEnergyBalance,T::Float"K",p::Float"Pa") = p/(T*seb.sebparams.Rₐ);

"""
Saturation pressure of water/ice according to the empirical August-Roche-Magnus formula
Note: T is passed [K] and converted to [°C]
"""
estar(T::Float"K") = ( (T>0) ? 611.2 * exp(17.62*(T-273.15)/(243.12-273.15+T))     # Eq. (B3) in Westermann et al. (2016)
                             : 611.2 * exp(22.46*(T-273.15)/(272.62-273.15+T)) ;
)

"""
Latent heat of evaporation/condensation of water according to https://en.wikipedia.org/wiki/Latent_heat#cite_note-RYfit-11
"""
Llg(T::Float"K") = 1000 * (2500.8 - 2.36*T + 0.0016*T^2 - 0.00006*T^3);

"""
Latent heat of sublimation/resublimation of water accodring to https://en.wikipedia.org/wiki/Latent_heat#cite_note-RYfit-11
"""
Lsg(T::Float"K") = 1000 * (2834.1 - 0.29*T - 0.004*T^2);

"""
Friction velocity according to Monin-Obukhov theory
"""
function ustar(seb::SurfaceEnergyBalance, stop)
    ustar = let κ = seb.sebparams.κ,
             uz = seb.forcing.wind(stop.t),                                                  # wind speed at height z
             z = seb.forcing.z,                                                              # height z of wind forcing
             z₀ = seb.sebparams.z₀,                                                          # aerodynamic roughness length [m]
             Lstar = stop.Lstar[1];
        κ * uz * (log(z/z₀) - Ψ_M(z/Lstar,z₀/Lstar))^-1                                      #Eq. (7) in Westermann et al. (2016)
    end
end

"""
Obukhov length according to Monin-Obukhov theory
"""
function Lstar(seb::SurfaceEnergyBalance, stop, ssoil)
    Lstar = let κ = seb.sebparams.κ,
                g = seb.sebparams.g,
                Rₐ = seb.sebparams.Rₐ,
                cₚ = seb.sebparams.cₐ / seb.sebparams.ρₐ,                                    # specific heat capacity of air at constant pressure
                Tₕ = seb.forcing.Tair(stop.t),                                               # air temperature at height z over surface
                p = seb.forcing.p(stop.t),                                            # atmospheric pressure at surface
                ustar = stop.ustar[1],
                Qₑ = stop.Qe[1],
                Qₕ = stop.Qh[1],
                Llg = Llg(ssoil.T[1]),
                ρₐ = density_air(seb,seb.forcing.Tair(stop.t),seb.forcing.p(stop.t));  # density of air at surface air temperature and surface pressure [kg/m^3]
            -ρₐ * cₚ * Tₕ / (κ * g) * ustar^3 / (Qₕ + 0.61*cₚ / Llg * Tₕ * Qₑ)                # Eq. (8) in Westermann et al. (2016)
    end
    # upper and lower limits for Lstar
    Lstar = (abs(Lstar)<1e-7) && sign(Lstar)*1e-7
    Lstar = (abs(Lstar)>1e+7) && sign(Lstar)*1e+7
end

"""
Sensible heat flux, defined as positive if it is a flux towards the surface
"""
function Q_H(seb::SurfaceEnergyBalance, stop, ssoil)
    let κ = seb.sebparams.κ,
        Rₐ = seb.sebparams.Rₐ,
        Tₕ = seb.forcing.Tair(stop.t),                                       # air temperature
        T₀ = ssoil.T[1],                                                         # surface temperature
        cₚ = seb.sebparams.cₐ / seb.sebparams.ρₐ,                                                       # specific heat capacity of air at constant pressure
        z = seb.forcing.z,                                                          # height at which forcing data are provided
        Lstar = stop.Lstar[1],
        ustar = stop.ustar[1],
        p = seb.forcing.p(stop.t),
        z₀ = seb.sebparams.z₀,
        ρₐ = density_air(seb,seb.forcing.Tair(stop.t),seb.forcing.p(stop.t));# density of air at surface air temperature and surface pressure [kg/m^3]

        rₐᴴ = (κ * ustar)^-1 * (log(z/z₀) - Ψ_HW(z/Lstar,z₀/Lstar))                          # Eq. (6) in Westermann et al. (2016)
        Q_H = -ρₐ * cₚ * (Tₕ-T₀) / rₐᴴ                                              # Eq. (4) in Westermann et al. (2016)
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
        Tₕ = seb.forcing.Tair(stop.t),                                                   # air temperature at height z over surface
        T₀ = ssoil.T[1],                                                         # surface temperature
        p = seb.forcing.p(stop.t),                                                              # atmospheric pressure at surface
        qₕ = seb.forcing.q(stop.t),                                                             # specific humidity at height h over surface
        z = seb.forcing.z(stop.t),                                                              # height at which forcing data are provided
        rₛ = seb.sebparams.rₛ,                                                               # surface resistance against evapotranspiration / sublimation [1/m]
        Lstar = stop.Lstar[1],
        ustar = stop.ustar[1],
        Llg = Llg(ssoil.T[1]),
        Lsg = Lsg(ssoil.T[1]),
        z₀ = seb.sebparams.z₀,                                                               # aerodynamic roughness length [m]
        ρₐ = density_air(seb,seb.forcing.Tair(stop.t),seb.forcing.p(stop.t));                        # density of air at surface air temperature and surface pressure [kg/m^3]

        q₀ = γ*estar(T₀)/p                                                              # saturation pressure of water/ice at the surface; Eq. (B1) in Westermann et al (2016)
        rₐᵂ = (κ * ustar)^-1 * (log(z/z₀) - Ψ_HW(z/Lstar,z₀/Lstar))                     # aerodynamic resistance Eq. (6) in Westermann et al. (2016)
        (T₀<=273.15) ? L = Lsg : L = Llg                                                # latent heat of sublimation/resublimation or evaporation/condensation [J/kg]

        # calculate Q_E
        (qₕ>q₀) ? Q_E = -ρₐ * L * (qₕ-q₀) / (rₐᵂ)   :                                   # Eq. (5) in Westermann et al. (2016) # condensation / resublimation (no aerodynamics resistance)
                  Q_E = -ρₐ * L * (qₕ-q₀) / (rₐᵂ+rₛ)                                    # evaporation / sublimation (account for surface resistance against evapotranspiration/sublimation)
    end
end

"""
Integrated stability function for heat/water transport
"""
function Ψ_HW(ζ₁::Float64, ζ₂::Float64)
    if ζ₁<=0 # neutral and unstable conditions
        1.9 * atanh((1 - 11.6 * ζ₁)^0.5) + log(ζ₁) - (1.9 * atanh((1 - 11.6 * ζ₂)^0.5) + log(ζ₂))
    else     # stable stratification
        0.5*((-5 + 5^0.5) * log(-3 + 5^0.5- 2*ζ₁) - (5 + 5^0.5) * log(3 + 5^0.5 + 2*ζ₁)) -
        0.5*((-5 + 5^0.5) * log(-3 + 5^0.5- 2*ζ₂) - (5 + 5^0.5) * log(3 + 5^0.5 + 2*ζ₂))
    end
end

"""
Integrated stability function for momentum transport
"""
function Ψ_M(ζ₁::Float64, ζ₂::Float64)
    if ζ₁<=0 # neutral and unstable conditions
         -2*atan((1 - 19*ζ₁)^(1/4)) + 2*log(1 + (1 - 19*ζ₁)^(1/4)) + log(1 + (1 - 19*ζ₁)^0.5) -
        (-2*atan((1 - 19*ζ₂)^(1/4)) + 2*log(1 + (1 - 19*ζ₂)^(1/4)) + log(1 + (1 - 19*ζ₂)^0.5))
    else     # stable stratification
         -19.5*(1 + ζ₁)^(1/3) - 7.5367*atan(0.57735 - 1.72489*(1 + ζ₁)^(1/3)) + 4.35131*log(3+4.4814*(1+ζ₁)^(1/3)) - 2.17566*log(3 - 4.4814*(1 + ζ₁)^(1/3) + 6.69433*(1 + ζ₁)^(2/3)) -
        (-19.5*(1 + ζ₂)^(1/3) - 7.5367*atan(0.57735 - 1.72489*(1 + ζ₂)^(1/3)) + 4.35131*log(3+4.4814*(1+ζ₂)^(1/3)) - 2.17566*log(3 - 4.4814*(1 + ζ₂)^(1/3) + 6.69433*(1 + ζ₂)^(2/3)))
    end
end

export SurfaceEnergyBalance
