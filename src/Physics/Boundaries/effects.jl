abstract type BoundaryEffect end
"""
    Damping{D,K} <: BoundaryEffect

Generic implementation of bulk conductive damping at the boundary.
"""
Base.@kwdef struct Damping{D,K} <: BoundaryEffect
    depth::D = t -> 0.0 # function of t -> damping depth; defaults to zero function
    k::K = Param(1.0, bounds=(0.0,Inf)) # conductivity of medium
end
function (damp::Damping)(u_top, u_sub, k_sub, Δsub, t)
    let d_med = damp.depth(t), # length of medium
        d_sub = Δsub / 2, # half length of upper grid cell
        d = d_med + d_sub, # total flux distance (from grid center to "surface")
        k_med = damp.k, # conductivity of damping medium (assumed uniform)
        k_sub = k_sub, # conductivity at grid center
        k = Numerics.harmonicmean(k_med, k_sub, d_med, d_sub);
        -k*(u_sub - u_top)/d/Δsub # flux per unit volume
    end
end
