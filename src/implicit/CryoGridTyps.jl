module CryoGridTyps
export stratigraphy, temporary, forcing, statvar, grid, out, para
struct stratigraphy
    Water::Array{Float64,2}
    Mineral::Array{Float64,2}
    Organic::Array{Float64,2}
    WaterIce::Array{Float64,2}
end
struct temporary
    lat_flux::Array{Float64,2}
    SnowDepth::Array{Float64,2}
    kp::Array{Float64,2}
    cp::Array{Float64,2}
    kn::Array{Float64,1}
    ks::Array{Float64,1}
end
struct forcing
    t_span::Array{Float64,1}
    Tair::Array{Float64,1}
    snowfall::Array{Float64,1}
end
struct statvar
    T::Array{Float64,2}
    H::Array{Float64,2}
end
struct grid
    Zp::Array{Float64,1}
    Zn::Array{Float64,1}
    Zs::Array{Float64,1}
    dxo::Array{Float64,2}
    dxp::Array{Float64,1}
    dxn::Array{Float64,1}
    dxs::Array{Float64,1}
    An::Array{Float64,2}
    As::Array{Float64,2}
    Ao::Array{Float64,2}
    Vp::Array{Float64,2}
end
struct out
    #output time stamp
    Date::Array{Int64,2}

    #Temparatur and water content on output time step
    T::Array{Float64,3}
    Water::Array{Float64,3}
    WaterIce::Array{Float64,3}

    #average Energy profile [J/m³]
    H_av::Array{Float64,3}

    #average and min max soil temperature profiles [°C]
    T_av::Array{Float64,3}
    T_min::Array{Float64,3}
    T_max::Array{Float64,3}

    #average and min max liquid water content profiles [vol fraction]
    W_av::Array{Float64,3}
    W_min::Array{Float64,3}
    W_max::Array{Float64,3}

    #annual sum of lateral heat fluxes of tiles [W d]
    Q_lat::Array{Float64,3}

    #freezing and thawing degree days
    FDD::Array{Float64,3}
    TDD::Array{Float64,3}
    FrostDays::Array{Int64,3}

    SnowDepth_av::Array{Float64,3}
    SnowDepth_max::Array{Float64,3}
    SnowDays::Array{Int64,3}
end
struct para
    SnowDensityZero::Array{Float64,1}
    SnowDensityMax::Array{Float64,1}
    SnowDensityK1::Array{Float64,1}
    SnowDensityK2::Array{Float64,1}
    SnowCoverMax::Array{Float64,2}
    WaterDensity::Array{Float64,1}
    WaterDepth::Array{Float64,2}
    Qgeo::Array{Float64,1}
end
struct save
    #output time stamp
    Date::Array{Int16,2}

    #Temparatur and water content on output time step
    #T::Array{Float16,3}
    #Water::Array{Float16,3}
    #WaterIce::Array{Float16,3}

    #average Energy profile [J/m³]
    H_av::Array{Float64,3}

    #average and min max soil temperature profiles [°C]
    T_av::Array{Float16,3}
    T_min::Array{Float16,3}
    T_max::Array{Float16,3}

    #average and min max liquid water content profiles [vol fraction]
    W_av::Array{Float16,3}
    W_min::Array{Float16,3}
    W_max::Array{Float16,3}

    #annual sum of lateral heat fluxes of tiles [W d]
    Q_lat::Array{Float32,3}

    #freezing and thawing degree days
    FDD::Array{Float32,3}
    TDD::Array{Float32,3}
    FrostDays::Array{Int64,3}

    SnowDepth_av::Array{Float32,3}
    SnowDepth_max::Array{Float32,3}
    SnowDays::Array{Int64,3}
end

"""
Convenience wrapper for truncating numerics in 'save' type.
"""
function save(
    Date::Array{Int64,2},
    H_av::Array{Float64,3},
    T_av::Array{Float64,3},
    T_min::Array{Float64,3},
    T_max::Array{Float64,3},
    W_av::Array{Float64,3},
    W_min::Array{Float64,3},
    W_max::Array{Float64,3},
    Q_lat::Array{Float64,3},
    FDD::Array{Float64,3},
    TDD::Array{Float64,3},
    FrostDays::Array{Int64,3},
    SnowDepth_av::Array{Float64,3},
    SnowDepth_max::Array{Float64,3},
    SnowDays::Array{Int64,3})
    save(convert(Array{Int16,2}, Date),
         H_av,
         convert(Array{Float16,3}, T_av),
         convert(Array{Float16,3}, T_min),
         convert(Array{Float16,3}, T_max),
         convert(Array{Float16,3}, W_av),
         convert(Array{Float16,3}, W_min),
         convert(Array{Float16,3}, W_max),
         convert(Array{Float32,3}, Q_lat),
         convert(Array{Float32,3}, FDD),
         convert(Array{Float32,3}, TDD),
         FrostDays,
         convert(Array{Float32,3}, SnowDepth_av),
         convert(Array{Float32,3}, SnowDepth_max),
         SnowDays)
end

end
