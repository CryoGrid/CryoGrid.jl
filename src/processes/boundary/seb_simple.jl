variables(soil::Soil, heat::Heat{UT"J"}) = (
    Diagnostic(:Sout, Float"W/(m^2)", Scalar()),    # outgoing shortwave radiation [J/(s*m^2)]
    Diagnostic(:Lout, Float"W/(m^2)", Scalar()),    # outgoing longwave radiation [J/(s*m^2)]
    Diagnostic(:Qnet, Float"W/(m^2)", Scalar()),    # net radiation budget at surface [J/(s*m^2)]
    Diagnostic(:Qh, Float"W/(m^2)", Scalar()),      # sensible heat flux [J/(s*m^2)]
    Diagnostic(:Qe, Float"W/(m^2)", Scalar()),      # latent heat flux [J/(s*m^2)]
    Diagnostic(:Qg, Float"W/(m^2)", Scalar()),      # ground heat flux [J/(s*m^2)]
    Diagnostic(:L∗, Float"m", Scalar()),            # Obukhov length [m]
    Diagnostic(:u∗, Float"m/s", Scalar()),          # friction velocity [m/s]
)

initialcondition!(soil::Soil, heat::Heat{UT"J"}, state) = (
    @setscalar state.S↑ = 0.;
    @setscalar state.L↑ = 0.;
    @setscalar state.Qnet = 0.;
    @setscalar state.Qh = 0.;
    @setscalar state.Qe = 0.;
    @setscalar state.Qg = 0.;
    @setscalar state.L∗ = -1e5;
    @setscalar state.u∗ = 10.;
)

@with_kw struct SEBParams{T} <: Params
    # surface properties --> should be associated with the Stratigraphy
    @setscalar α::Float"1" = 0.2xu"1"                          # surface albedo [-]
    @setscalar ϵ::Float"1" = 0.97xu"1"                         # surface emissivity [-]
    @setscalar z₀::Float"m" = 1e-3xu"m"                        # surface roughness length [m]
    @setscalar rₛ::Float"1/m" = 50.0xu"1/m"                      # surface resistance against evapotranspiration and sublimation [1/m]

    # "natural" constant
    @setscalar σ::Float64 = 5.6704e-8                          # Stefan-Boltzmann constant
    @setscalar κ::Float"1" = 0.4xu"1"                          # von Kármán constant [-]
    @setscalar γ::Float"1" = 0.622xu"1"                        # Psychrometric constant [-]
    @setscalar Rₐ::Float"J/(kg*K)" = 287.058xu"J/(kg*K)"       # specific gas constant of air [J/(kg*K)]
    @setscalar g::Float"m/s^2" = 9.81xu"m/s^2"                 # gravitational acceleration [m/s^2]

    # material properties (assumed to be constant)
    @setscalar ρₐ::Float"kg/m^3" = 1.293xu"kg/m^3"             # density of air at standard pressure and 0°C [kg/m^3]
    @setscalar cₐ::Float"J/(m^3*K)"= 1005.7xu"J/(kg*K)"*ρₐ     # volumetric heat capacity of dry air at standard pressure and 0°C [J/(m^3*K)]
end

@with_kw struct SurfaceEnergyBalance{P,T} <: BoundaryProcess{P}
    forcing::T
    sebparams::SEBParams
    SurfaceEnergyBalance{P}(forcing::T) where {P<:Heat,T} = new{P,T}(forcing)
end

BoundaryStyle(::Type{<:SurfaceEnergyBalance}) = Neumann()

"""
Top interaction, ground heat flux from surface energy balance. (no snow, no water body, no infiltration)
"""
function (seb::SurfaceEnergyBalance)(top::Top, soil::Soil, heat::Heat, stop, ssoil)

    # 1. calculate radiation budget
    # outgoing shortwave radiation as reflected
    @setscalar stop.S↑ = let α=seb.sebparams.α, S↓=seb.forcing.Sin(stop.t)
                -α * S↓                                                    # Eq. (2) in Westermann et al. (2016)
    end

    # outgoing longwave radiation composed of emitted and reflected radiation
    @setscalar stop.L↑ = let ϵ=seb.sebparams.ϵ, σ=seb.sebparams.σ, T₀=ssoil.T[1], L↓=seb.forcing.Lin(stop.t);
                -ϵ * σ * T₀^4 - (1-ϵ) * L↓                                  # Eq. (3) in Westermann et al. (2016)
    end

    # net radiation budget
    @setscalar stop.Qnet = let S↓=seb.forcing.Sin(stop.t), L↓=seb.forcing.Lin(stop.t), Sout=stop.S↑, Lout=stop.L↑;
                S↓ + S↑ + L↓ + L↑
    end

    # 2. calcuate turbulent heat flux budget
    # determine atmospheric stability conditions
    @setscalar stop.u∗ = u∗(seb, stop)
    @setscalar stop.L∗ = L∗(seb, stop)

    # sensible heat flux
    @setscalar stop.Qh = Q_H(seb, stop, ssoil)

    # latent heat flux
    @setscalar stop.Qe = Q_E(seb, stop, ssoil)

    # 3. determine ground heat flux as the residual of the radiative and turbulent fluxes
    @setscalar stop.Qg = let Qnet=stop.Qnet, Qₕ=stop.Qh, Qₑ=stop.Qe;
                    Qnet - Qₕ - Qₑ                                                  # essentially Eq. (1) in Westermann et al. (2016)

    # 4. return the ground heat flux to the uppermost soil grid cell
    return stop.Qg
end

"""
Density of air at giben tempeature and pressure
"""
density_air(seb::SurfaceEnergyBalance,T::Float"K",p::Float"Pa") = p/(T*seb.sebparams.Rₐ);

"""
Saturation pressure of water/ice according to the empirical August-Roche-Magnus formula
Note: T is passed [K] and converted to [°C]
"""
e∗(T::Float"K") = ( (T>0) ? 611.2 * exp(17.62*(T-273.15)/(243.12-273.15+T))     # Eq. (B3) in Westermann et al. (2016)
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
function u∗(seb::SurfaceEnergyBalance, stop)
    u∗ = let κ = seb.sebparams.κ,
             uz = seb.forcing.wind(stop.t),                                                  # wind speed at height z
             z = seb.forcing.z,                                                              # height z of wind forcing
             z₀ = seb.sebparams.z₀,                                                          # aerodynamic roughness length [m]
             L∗ = stop.L∗;
        κ * uz * (log(z/z₀) - Ψ_M(z/L∗,z₀/L∗))^-1                                            #Eq. (7) in Westermann et al. (2016)
    end
end

"""
Obukhov length according to Monin-Obukhov theory
"""
function L∗(seb::SurfaceEnergyBalance, stop, ssoil)
    L∗ = let κ = seb.sebparams.κ,
                g = seb.sebparams.g,
                Rₐ = seb.sebparams.Rₐ,
                cₚ = seb.sebparams.cₐ / seb.sebparams.ρₐ,                                                       # specific heat capacity of air at constant pressure
                Tₕ = seb.forcing.Tair(stop.t),                                                   # air temperature at height z over surface
                p = seb.forcing.p(stop.t),                                                              # atmospheric pressure at surface
                u∗ = stop.u∗,
                Qₑ = stop.Qe,
                Llg = Llg(ssoil.T[1]),
                ρₐ = density_air(seb,seb.forcing.Tair(stop.t),seb.forcing.p(stop.t));                        # density of air at surface air temperature and surface pressure [kg/m^3]
            -ρₐ * cₚ * Tₕ / (κ * g) * u∗^3 / (Qₕ + 0.61*cₚ / Llg * Tₕ * Qₑ))                # Eq. (8) in Westermann et al. (2016)
    end
    # upper and lower limits for L∗
    (abs(L∗)<1e-7) ? L∗=sign(L∗)*1e-7 : ;
    (abs(L∗)>1e+7) ? L∗=sign(L∗)*1e7 : ;

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
        L∗ = stop.L∗,
        u∗ = stop.u∗,
        p = seb.forcing.p(stop.t),
        z₀ = seb.sebparams.z₀,
        ρₐ = density_air(seb,seb.forcing.Tair(stop.t),seb.forcing.p(stop.t));# density of air at surface air temperature and surface pressure [kg/m^3]

        rₐᴴ = (κ * u∗)^-1 * (log(z/z₀) - Ψ_HW(z/L∗,z₀/L∗))                          # Eq. (6) in Westermann et al. (2016)
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
        L∗ = stop.L∗,
        u∗ = stop.u∗,
        Llg = Llg(ssoil.T[1]),
        Lsg = Lsg(ssoil.T[1]),
        z₀ = seb.sebparams.z₀,                                                               # aerodynamic roughness length [m]
        ρₐ = density_air(seb,seb.forcing.Tair(stop.t),seb.forcing.p(stop.t));                        # density of air at surface air temperature and surface pressure [kg/m^3]

        q₀ = γ*e∗(T₀)/p                                                          # saturation pressure of water/ice at the surface; Eq. (B1) in Westermann et al (2016)
        rₐᵂ = (κ * u∗)^-1 * (log(z/z₀) - Ψ_HW(z/L∗,z₀/L∗))                          # aerodynamic resistance Eq. (6) in Westermann et al. (2016)
        (T₀<=273.15) ? L = Lsg : L = Llg                                            # latent heat of sublimation/resublimation or evaporation/condensation [J/kg]

        # calculate Q_E
        (qₕ>q₀) ? Q_E = -ρₐ * L * (qₕ-q₀) / (rₐᵂ)                                    # Eq. (5) in Westermann et al. (2016) # condensation / resublimation (no aerodynamics resistance)
                : Q_E = -ρₐ * L * (qₕ-q₀) / (rₐᵂ+rₛ)                                 # evaporation / sublimation (account for surface resistance against evapotranspiration/sublimation)
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
