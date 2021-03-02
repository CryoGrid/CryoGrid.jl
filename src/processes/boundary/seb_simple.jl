function variables(soil::Soil, heat::Heat{UT"J"}) = (
    Diagnostic(:Sout, Float"W/(m^2)", Scalar()),    # outgoing shortwave radiation [J/(s*m^2)]
    Diagnostic(:Lout, Float"W/(m^2)", Scalar()),    # outgoing longwave radiation [J/(s*m^2)]
    Diagnostic(:Qnet, Float"W/(m^2)", Scalar()),    # net radiation budget at surface [J/(s*m^2)]
    Diagnostic(:Qh, Float"W/(m^2)", Scalar()),      # sensible heat flux [J/(s*m^2)]
    Diagnostic(:Qe, Float"W/(m^2)", Scalar()),      # latent heat flux [J/(s*m^2)]
    Diagnostic(:Qg, Float"W/(m^2)", Scalar()),      # ground heat flux [J/(s*m^2)]
    Diagnostic(:L∗, Float"m", Scalar()),            # Obukhov length [m]
    Diagnostic(:u∗, Float"m/s", Scalar()),          # friction velocity [m/s]
)

function initialcondition!(soil::Soil, heat::Heat{UT"J"}, state)
    state.Sout = 0.;
    state.Lout = 0.;
    state.Qnet = 0.;
    state.Qh = 0.;
    state.Qe = 0.;
    state.Qg = 0.;
    state.L∗ = -1e5;
    state.u∗ = 10.;
end

@with_kw struct SurfaceEnergyBalance{P,T} <: BoundaryProcess{P}
    forcing::T
    SurfaceEnergyBalance{P}(forcing::T) where {P<:Heat,T} = new{P,T}(forcing)
end

@with_kw struct SEBParams{T} <: Params
    # "natural" constant
    σ::Float64 = 5.6704e-8                          # Stefan-Boltzmann constant
    κ::Float"1" = 0.4xu"1"                          # von Kármán constant [-]
    γ::Float"1" = 0.622xu"1"                        # Psychrometric constant [-]
    Rₐ::Float"J/(kg*K)" = 287.058xu"J/(kg*K)"       # specific gas constant of air [J/(kg*K)]
    g::Float"m/s^2" = 9.81xu"m/s^2"                 # gravitational acceleration

    # material properties (assumed to be constant)
    ρₐ::Float"kg/m^3" = 1.293xu"kg/m^3"             # density of air at standard pressure and 0°C [kg/m^3]
    cₐ::Float"J/(m^3*K)"= 1005.7xu"J/(kg*K)"*ρₐ     # volumetric heat capacity of dry air at standard pressure and 0°C [J/(m^3*K)]
end

BoundaryStyle(::Type{<:SurfaceEnergyBalance}) = Neumann()

"""
Top interaction, ground heat flux from surface energy balance. (no snow, no water body, no infiltration)
"""
function interact!(top::Top, seb::SurfaceEnergyBalance, soil::Soil, heat::Heat, stop, ssoil)

    # 1. calculate radiation budget
    # outgoing shortwave radiation as reflected
    stop.Sout = let α=stop.α, Sin=forcing.Sin;
                -α * Sin                                                    # Eq. (2) in Westermann et al. (2016)
    end

    # outgoing longwave radiation composed of emitted and reflected radiation
    stop.Lout = let ϵ=stop.ϵ, σ=seb.σ, Tₛ=(ssoil.T[1]+273.15), Lin=forcing.Lin;
                -ϵ * σ * T^4 - (1-ϵ) * Lin                                  # Eq. (3) in Westermann et al. (2016)
    end

    # net radiation budget
    stop.Qnet = let Sin=forcing.Sin, Lin=forcing.Lin, Sout=seb.Sout, Lout=seb.Lout;
                Sin + Sout + Lin + Lout
    end

    # 2. calcuate turbulent heat flux budget
    # determine atmospheric stability conditions
    stop.u∗ = ustar(forcing, seb, stop)
    stop.L∗ = Lstar(forcing, seb, stop)

    # sensible heat flux
    stop.Qh = Q_H(forcing, seb, stop)

    # latent heat flux
    stop.Qe = Q_E(forcing, seb, stop)

    # 3. determine ground heat flux as the residual of the radiative and turbulent fluxes
    stop.Qg = let Qnet=stop.Qnet, Qₕ=stop.Qh, Qₑ=stop.Qe;
                    Qnet - Qₕ - Qₑ                                                  # essentially Eq. (1) in Westermann et al. (2016)

    # 4. apply ground heat flux to the uppermost soil grid cell taking into account its size
    Δk = Δ(ssoil.grids.k)
    ssoil.dH[1] += let Qg=stop.Qg , a=(Δk[1]);
        Qg/a
    end
    return nothing # ensure no allocation
end

"""
Density of air at gibven tempeature and pressure
"""
density_air(seb,T,p) = p/(T*seb.Rₐ);

"""
Friction velocity according to Monin-Obukhov theory
"""
function ustar(forcing, seb, stop)
    ustar = let κ = seb.κ,
             uz = forcing.wind,                                                          # wind speed at height z
             z = forcing.z,                                                              # height z of wind forcing
             z₀ = stop.z0,                                                               # aerodynamic roughness length [m]
             L∗ = seb.L∗;
        κ * uz * (log(z/z₀) - Ψ_M(z/L∗,z₀/L∗))^-1                              #Eq. (7) in Westermann et al. (2016)
    end
end

"""
Obukhov length according to Monin-Obukhov theory
"""
function Lstar(forcing, seb)
    Lstar = let κ = seb.κ,
                g = seb.g,
                Rₐ = seb.Rₐ,
                cₚ = seb.cₐ / seb.ρₐ,                                                       # specific heat capacity of air at constant pressure
                Tₕ = forcing.Tair+273.15,                                                   # air temperature at height z over surface
                p = forcing.p,                                                              # atmospheric pressure at surface
                u∗ = seb.u∗,
                Qₕ = seb.Qh,
                Qₑ = seb.Qe,
                Llg = 1000 * (2500.8 - 2.36*stop.T + 0.0016*stop.T^2 - 0.00006*stop.T^3),   #[J/kg] latent heat of evaporation/condensation of water according to https://en.wikipedia.org/wiki/Latent_heat#cite_note-RYfit-11
                ρₐ = density_air(seb,forcing.Tair+273.15,forcing.p);                        # density of air at surface air temperature and surface pressure [kg/m^3]
            -ρₐ * cₚ * Tₕ / (κ * g) * u∗^3 / (Qₕ + 0.61*cₚ / Llg * Tₕ * Qₑ))                # Eq. (8) in Westermann et al. (2016)
    end
    # upper and lower limits for Lstar
    (abs(Lstar)<1e-7) ? Lstar=sign(Lstar)*1e-7 : ;
    (abs(Lstar)>1e+7) ? Lstar=sign(Lstar)*1e7 : ;

end

"""
Sensible heat flux, defined as positive if it is a flux towards the surface
"""
function Q_H(forcing, seb, stop)
    let κ = seb.κ,
        Rₐ = seb.Rₐ,
        Tₕ = forcing.Tair+273.15,                                                   # air temperature
        T₀ = stop.T+273.15,                                                         # surface temperature
        cₚ = seb.cₐ / seb.ρₐ,                                                       # specific heat capacity of air at constant pressure
        z = forcing.z,                                                              # height at which forcing data are provided
        L∗ = seb.L∗,
        u∗ = seb.u∗,
        p = forcing.p,
        z₀ = stop.z0,
        ρₐ = density_air(seb,forcing.Tair+273.15,forcing.p);                        # density of air at surface air temperature and surface pressure [kg/m^3]

        rₐᴴ = (κ * u∗)^-1 * (log(z/z₀) - Ψ_HW(z/L∗,z₀/L∗))                          # Eq. (6) in Westermann et al. (2016)
        Q_H = -ρₐ * cₚ * (Tₕ-T₀) / rₐᴴ                                               # Eq. (4) in Westermann et al. (2016)
    end
end

"""
Latent heat flux, defined as positive if it is a flux towards the surface.
Represents evapo(transpi)ration/condensation at positive surface temperatures and
sublimation/resublimation at negative surface temperatures
"""
function Q_E(forcing, seb, stop)
    let κ = seb.κ,
        γ = seb.γ,
        Rₐ = seb.Rₐ,
        Tₕ = forcing.Tair+273.15,                                                   # air temperature at height z over surface
        T₀ = stop.T+273.15,                                                         # surface temperature
        p = forcing.p,                                                              # atmospheric pressure at surface
        qₕ = forcing.q,                                                             # specific humidity at height h over surface
        z = forcing.z,                                                              # height at which forcing data are provided
        rₛ = stop.rₛ,                                                               # surface resistance against evapotranspiration / sublimation [1/m]
        L∗ = seb.L∗,
        u∗ = seb.u∗,
        Llg = 1000 * (2500.8 - 2.36*stop.T + 0.0016*stop.T^2 - 0.00006*stop.T^3),   #[J/kg] latent heat of evaporation/condensation of water according to https://en.wikipedia.org/wiki/Latent_heat#cite_note-RYfit-11
        Lsg = 1000 * (2834.1 - 0.29*stop.T - 0.004*stop.T^2),                       #[J/kg] latent heat of sublimation/resublimation of water accodring to https://en.wikipedia.org/wiki/Latent_heat#cite_note-RYfit-11
        z₀ = stop.z0,                                                               # aerodynamic roughness length [m]
        ρₐ = density_air(seb,forcing.Tair+273.15,forcing.p);                        # density of air at surface air temperature and surface pressure [kg/m^3]

        q₀ = γ*estar(T₀)/p                                                          # saturation pressure of water/ice at the surface; Eq. (B1) in Westermann et al (2016)
        rₐᵂ = (κ * u∗)^-1 * (log(z/z₀) - Ψ_HW(z/L∗,z₀/L∗))                          # aerodynamic resistance Eq. (6) in Westermann et al. (2016)
        (T₀<=273.15) ? L = Lsg : L = Llg                                            # latent heat of sublimation/resublimation or evaporation/condensation [J/kg]

        # calculate Q_E
        (qₕ>q₀) ? Q_E = -ρₐ * L * (qₕ-q₀) / (rₐᵂ)                                    # Eq. (5) in Westermann et al. (2016) # condensation / resublimation (no aerodynamics resistance)
                : Q_E = -ρₐ * L * (qₕ-q₀) / (rₐᵂ+rₛ)                                 # evaporation / sublimation (account for surface resistance against evapotranspiration/sublimation)
    end
end

"""
Saturation pressure of water/ice according to the empirical August-Roche-Magnus formula
Note: T is in [K]
"""
function estar(T::Float"K")
    if T>0
        611.2 * exp(17.62*(T-273.15)/(243.12-273.15+T))     # Eq. (B3) in Westermann et al. (2016)
    else
        611.2 * exp(22.46*(T-273.15)/(272.62-273.15+T))     # Eq. (B3) in Westermann et al. (2016)
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
