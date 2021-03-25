# Fallback (error) implementation for freeze curve
(fc::FreezeCurve)(layer::SubSurface, heat::Heat, state) =
    error("freeze curve $(typeof(fc)) not implemented for layer $(typeof(layer))")

"""
Implementation of "free water" freeze curve for any subsurface layer. Assumes that
'state' contains at least temperature (T), enthalpy (H), heat capacity (C),
total water content (θw), and liquid water content (θl).
"""
@inline function (fc::FreeWater)(layer::SubSurface, heat::Heat{u"J"}, state)
    @inline function enthalpyinv(H, C, L, θtot)
        let θtot = max(1.0e-8,θtot),
            Lθ = L*θtot,
            I_t = H > Lθ,
            I_f = H <= 0.0;
            T = (I_t*(H-Lθ) + I_f*H)/C + 273.15
        end
    end
    @inline function freezethaw(H, C, L, θtot)
        let θtot = max(1.0e-8,θtot),
            Lθ = L*θtot,
            I_t = H > Lθ,
            I_c = (H > 0.0) && (H <= Lθ);
            liquidfraction = I_c*(H/Lθ) + I_t
        end
    end
    L = heat.params.L
    @. state.T = enthalpyinv(state.H, state.C, L, state.θw)
    @. state.θl = freezethaw(state.H, state.C, L, state.θw)*state.θw
end
