########
# Simple heat conduction model (no forcings)
# Uses ModelingToolkit
# Currently does not work because the MOL discretization strategy does not
# support variables without time derivatives, which is required for this model.
########
module HeatConduction

using DifferentialEquations, ModelingToolkit, DiffEqBase

export model, thermalConductivity, heatCapacity

function thermalConductivity(waterIce, water, mineral, organic)
    ka = 0.025;       #air [Hillel(1982)]
    kw = 0.57;        #water [Hillel(1982)]
    ko = 0.25;        #organic [Hillel(1982)]
    km = 3.8;         #mineral [Hillel(1982)]
    ki = 2.2;         #ice [Hillel(1982)]
    ice = waterIce - water;
    air = 1.0 - waterIce - mineral - organic;
    (water*kw^0.5 + ice*ki^0.5 + mineral*km^0.5 + organic*ko^0.5 + air*ka^0.5)^2
end

function heatCapacity(waterIce, water, mineral, organic)
    cs = 2000; #[J/kgK]  heat capacity solid
    cl = 2500; #[J/kgK]  heat capacity liquid
    cm = 1000; #[J/kgK]  heat capacity matrix
    cp = 1.0;  #[J/kgK]  heat capacity pore space

    cw = 4.2*10^6; #[J/m^3K] heat capacity water
    co = 2.5*10^6; #[J/m^3K]  heat capacity organic
    cm = 2*10^6; #[J/m^3K]  heat capacity mineral
    ca = 0.00125*10^6;#[J/m^3K]  heat capacity pore space
    ci = 1.9*10^6;#[J/m^3K]  heat capacity ice

    air = 1.0 - waterIce - mineral - organic
    ice = waterIce - water
    water*cw + ice*ci + mineral*cm + organic*co + air*ca
end

function enthalpyInv(H, θ, dHdT)
    ρ = 1000.; #[kg/m^3]
    Lsl = 334000.; #[J/kg]
    L = ρ*Lsl;#[J/m^3]
    if θ <= 0.0
        θ = 1e-8
    end
    if H > 0.0 && H <= L*θ
        0.0
    elseif H > L*θ
        (H - L*θ) / dHdT
    elseif H < 0.0
        H / dHdT
    end
end

function model()
    ivars = @parameters t z
    dvars = @variables H(..) T(..) k(..) dHdT(..) Wₗ(..) W(..) M(..) O(..)
    @derivatives Dt'~t
    @derivatives Dz'~z
    @register thermalConductivity(waterIce, water, mineral, organic)
    @register heatCapacity(waterIce, water, mineral, organic)
    @register enthalpyInv(H,W,dHdT)

    # Space and time domains
    domains = [t ∈ IntervalDomain(0.0,1.0),
               z ∈ IntervalDomain(-5.0,100.0)]
    eqs  = [
        W(z) ~ 0.5,
        Wₗ(z) ~ 0.5,
        M(z) ~ 0.3,
        O(z) ~ 0.2,
        k(z) ~ thermalConductivity(W(z),Wₗ(z),M(z),O(z)),
        dHdT(z) ~ heatCapacity(W(z),Wₗ(z),M(z),O(z)),
        T(z,t) ~ enthalpyInv(H(z,t),W(z),dHdT(z)),
        Dt(H(z,t)) ~ Dz(k(z)*Dz(T(z,t)))
    ]
    T₀(z) = -0.1*z
    bcs = [
        T(z,0) ~ T₀(z),
        T(-5.0,t) ~ 1.0,
        T(100,t) ~ 10.0
    ]
    PDESystem(eqs, bcs, domains, ivars, dvars)
end

function model2()
    W,Wₗ,M,O = 0.5,0.5,0.3,0.2
    # dHdT(W,Wₗ,M,O) = heatCapacity(W,Wₗ,M,O)
    # k(W,Wₗ,M,O) = thermalConductivity(W,Wₗ,M,O)
    # T(H,W,dHdT) = enthalpyInv(H,W,dHdT)
    ivars = @parameters t z # W Wₗ M O
    dvars = @variables H(..)
    @derivatives Dt'~t
    @derivatives Dz'~z
    #@register thermalConductivity(waterIce, water, mineral, organic)
    #@register heatCapacity(waterIce, water, mineral, organic)
    @register enthalpyInv(H,W,dHdT)

    # Space and time domains
    domains = [t ∈ IntervalDomain(0.0,1.0),
               z ∈ IntervalDomain(-5.0,100.0)]
    eq = Dt(H(z,t)) ~ Dz(thermalConductivity(W,Wₗ,M,O)*Dz(enthalpyInv(H(z,t),W,heatCapacity(W,Wₗ,M,O))))
    bcs = [
        H(z,0) ~ 1.0,
        H(-5.0,t) ~ 1.0,
        H(100,t) ~ 10.0
    ]
    PDESystem(eq, bcs, domains, ivars, dvars)
end

end
