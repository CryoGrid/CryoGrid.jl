"""
    SEBInputs{TT1,TT2,TR,TP,TQ,TW,TZ,Tsurf}

Non-prognostic input variables to the SEB that are assumed to be known a priori.
"""
Base.@kwdef struct SEBInputs{TT1,TT2,TR,TP,TQ,TW,TZ,Tsurf}
    Ts::TT1
    Tair::TT2
    Lin::TR
    Sin::TR
    pr::TP
    qh::TQ
    wind::TW
    z::TZ
    # surface properties
    surf::Tsurf = SurfaceProperties()
end
function SEBInputs(seb::SurfaceEnergyBalance, sub::SubSurface, stop, ssub)
    Ts = ssub.T[1] # soil temperature at surface
    Tair = seb.forcings.Tair(stop.t)
    Lin = seb.forcings.Lin(stop.t)
    Sin = seb.forcings.Sin(stop.t)
    pr = seb.forcings.pr(stop.t)
    qh = seb.forcings.qh(stop.t)
    wind = seb.forcings.wind(stop.t)
    z = seb.forcings.z
    surf = surfaceproperties(seb, sub)
    return SEBInputs(; Ts, Tair, Lin, Sin, pr, qh, wind, z, surf)
end

"""
    SEBState{TQ,TL,TU,TIn}

Unknown state variables that constitute a candidate solution to the SEB.
"""
Base.@kwdef struct SEBState{TQ,TL,TU,TIn}
    Qh::TQ
    Qe::TQ
    Lstar::TL
    ustar::TU
    inputs::TIn
end
function SEBState(seb::SurfaceEnergyBalance, sub::SubSurface, stop, ssub)
    Qh = getscalar(stop.Qh)
    Qe = getscalar(stop.Qe)
    Lstar = getscalar(stop.Lstar)
    ustar = getscalar(stop.ustar)
    inputs = SEBInputs(seb, sub, stop, ssub)
    return SEBState(; Qh, Qe, Lstar, ustar, inputs)
end


"""
    SEBOutputs{TQs,TQg,TQnet,TSout,TLout,TL,TU}

Outputs of surface energy balance computation that include both the updated `state` as well as
the net radiation `Qnet` and outgoing radiation components `Sout` and `Lout`.
"""
Base.@kwdef struct SEBOutputs{TQs,TQg,TQnet,TSout,TLout,TL,TU}
    state::SEBState{TQs,TL,TU}
    Qg::TQg
    Qnet::TQnet
    Sout::TSout
    Lout::TLout
end
