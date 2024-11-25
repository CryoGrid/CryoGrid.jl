"""
    SEBInputs{TT1,TT2,TR,TP,TQ,TW,Tsurf}

Non-prognostic input variables to the SEB that are assumed to be given.
"""
Base.@kwdef struct SEBInputs{TT1,TT2,TR,TP,TQ,TW,Tsurf}
    Ts::TT1 = 0.0u"°C"
    Tair::TT2 = 0.0u"°C"
    Lin::TR = 0.0u"W/m^2"
    Sin::TR = 0.0u"W/m^2"
    pr::TP = 0.0u"m/s"
    qh::TQ = 0.0
    wind::TW = 1e-8u"m/s"
    # surface properties
    surf::Tsurf = SurfaceProperties()
end

function SEBInputs(seb::SurfaceEnergyBalance, sub::SubSurface, stop, ssub)
    Ts = ssub.T[1] # soil temperature at surface
    Tair = seb.inputs.Tair
    Lin = seb.inputs.Lin
    Sin = seb.inputs.Sin
    pr = seb.inputs.pr
    qh = seb.inputs.qh
    wind = seb.inputs.wind
    surf = surfaceproperties(seb, sub)
    return SEBInputs(; Ts, Tair, Lin, Sin, pr, qh, wind, surf)
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
