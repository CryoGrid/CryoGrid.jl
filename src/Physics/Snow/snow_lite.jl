"""
    LiteGridded{T1,T2,T3,T4} <: SnowpackParameterization

Simple bulk density snow scheme minimally adapted from `CryoGridLite`.
"""
Base.@kwdef struct LiteGridded{T1,T2,T3,T4,T5,Theat,Twater} <: SnowpackParameterization
    dsn_max::T1 = 2.0u"m"
    ρsn_0::T2 = 250.0u"kg/m^3"
    ρsn_max::T3 = 500.0u"kg/m^3"
    ρsn_k1::T4 = 0.003
    ρsn_k2::T5 = 0.005
    # heat and water flow properties
    heat::Theat = SnowThermalProperties() # thermal properties
    water::Twater = HydraulicProperties(kw_sat=1e-4) # hydraulic properties
end

const LiteSnowpack = Snowpack{Tpara} where {Tpara<:LiteGridded}

# TODO: consider using the common interface for ablation, compaction, and accumulation?

snowdensity!(::LiteSnowpack, ::SnowMassBalance, state) = nothing

ablation!(::Top, ::SnowBC, ::LiteSnowpack, ::SnowMassBalance, stop, ssnow) = nothing

compaction!(::Top, ::SnowBC, ::LiteSnowpack, ::SnowMassBalance, stop, ssnow) = nothing

accumulation!(::Top, ::SnowBC, ::LiteSnowpack, ::SnowMassBalance, stop, ssnow) = nothing

# Declare snow mass balance variables for Lite scheme
CryoGrid.variables(::LiteSnowpack, ::SnowMassBalance) = (
    Prognostic(:swe, OnGrid(Cells), u"m"),
    Prognostic(:dsn, Scalar, u"m"),
    Diagnostic(:ρsn, OnGrid(Cells), u"kg/m^3"),
    Diagnostic(:snowfall, Scalar, u"m/s"),
    Diagnostic(:T_ub, Scalar, u"°C"),
    Diagnostic(:θw, OnGrid(Cells), domain=0..1),
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
    Diagnostic(:ubc_idx, Scalar, NoUnits, Int),
)

CryoGrid.computediagnostic!(::LiteSnowpack, ::SnowMassBalance, state) = nothing

CryoGrid.computefluxes!(::LiteSnowpack, ::SnowMassBalance, state) = nothing

function CryoGrid.interact!(
    top::Top,
    sbc::SnowBC,
    snow::LiteSnowpack,
    mass::SnowMassBalance{<:LinearAccumulation},
    stop,
    ssnow
)
    # get scaling factor(s)
    rate_scale = mass.accumulation.rate_scale
    # get snowfall from top layer flux
    @. ssnow.snowfall = rate_scale*stop.jw_snow
    return nothing
end

# interact! for implicit heat scheme

function CryoGrid.interact!(top::Top, bc::HeatBC, snowpack::Snowpack, heat::HeatBalanceImplicit, stop, ssnow)
    ssnow.T_ub .= T_ub = getscalar(stop.T_ub)
    ubc_idx = Int(getscalar(ssnow.ubc_idx))
    # get variables
    T = ssnow.T
    an = ssnow.DT_an
    as = ssnow.DT_as
    ap = ssnow.DT_ap
    bp = ssnow.DT_bp
    k = ssnow.k
    ∂H∂T = ssnow.∂H∂T
    dx = Δ(cells(ssnow.grid))
    dxp = Δ(ssnow.grid)
    # topmost cell
    bp[1] = T_ub*k[1] / (dxp[1]/2) / dxp[1]
    ap[1] = k[1] / (dxp[1]/2) / dxp[1]
    ∂H∂T[1] = zero(eltype(∂H∂T))
    T[1] = T_ub
    if ubc_idx > 1
        as[1] = an[1] = zero(eltype(as))
    end
    # inner cells
    @inbounds for i in 2:ubc_idx-1
        bp[i] = T_ub*k[i] / dx[i-1] / dxp[i]
        ap[i] = an[i]
        as[i] = zero(eltype(as))
        an[i] = zero(eltype(an))
        ∂H∂T[i] = zero(eltype(∂H∂T))
        T[i] = T_ub
    end
    return nothing
end
function CryoGrid.interact!(snowpack::Snowpack, ::HeatBalanceImplicit, sub::SubSurface, ::HeatBalanceImplicit, ssnow, ssub)
    Δk₁ = CryoGrid.thickness(snowpack, ssnow, last)
    Δk₂ = CryoGrid.thickness(sub, ssub, first)
    Δz = CryoGrid.midpoint(sub, ssub, first) - CryoGrid.midpoint(snowpack, ssnow, last)
    ubc_idx = Int(getscalar(ssnow.ubc_idx))
    if ubc_idx <= length(ssnow.H)
        k = ssnow.k[end] = ssub.k[1] =
            @inbounds let k₁ = ssnow.kc[end],
                k₂ = ssub.kc[1],
                Δ₁ = Δk₁[end],
                Δ₂ = Δk₂[1];
                harmonicmean(k₁, k₂, Δ₁, Δ₂)
            end
    else
        k = ssub.k[1]
    end
    ssnow.DT_ap[end] += ssnow.DT_as[end] = (ubc_idx <= length(ssnow.H))*k / Δz / Δk₁
    ssub.DT_ap[1] += ssub.DT_an[1] = k / Δz / Δk₂
    return nothing
end

###################################

function CryoGrid.diagnosticstep!(snowpack::LiteSnowpack, state)
    para = snowpack.para
    grid = state.grid
    dz = Δ(grid)
    dsn = getscalar(state.dsn)
    ρw = waterdensity(snowpack)
    sf = getscalar(state.snowfall)*24*3600
    date = convert_t(state.t)
    state.ρsn .= para.ρsn_0
    new_dsn, new_ubc_idx = lite_snow!(state.θw, state.θwi, grid, dsn, para.dsn_max, sf, state.ρsn, para.ρsn_max, para.ρsn_k1, para.ρsn_k2, ρw, date)
    state.dsn .= new_dsn
    state.ubc_idx .= new_ubc_idx
    # calculate swe from updated water/ice contents
    @inbounds for i in eachindex(state.θwi)
        state.swe[i] = (state.θwi[i] - state.θw[i])*min(dz[i], abs(dsn - grid[i+1]))
    end
    @. state.ρsn = ρw*state.swe / dz
    return false
end

function CryoGrid.initialcondition!(snowpack::LiteSnowpack, ::SnowMassBalance, state)
    @. state.ρsn = snowpack.para.ρsn_0
    @. state.swe = 0.0
    @. state.dsn = 0.0
    @. state.ubc_idx = length(state.T)
end

function Hydrology.watercontent!(snow::LiteSnowpack, ::WaterBalance, state)
    ρw = waterdensity(snow)
    # snow water/ice = swe scaled by grid cell thickness
    dz = Δ(state.grid)
    @. state.θwi = state.swe / dz
    # pore space
    @. state.θsat = 1 - state.θwi
    # add pore water (saturation) to total water content
    @. state.θwi += state.θsat*state.sat
    return nothing
end

function lite_snow!(Water, WaterIce, grid, snow_depth, snow_max, snow_fall, SnowDensity, SnowDensityMax, SnowDensityK1, SnowDensityK2, WaterDensity, current_date::DateTime)
    # last cell index is always just the last local grid cell
    snow_lower_idx = length(grid)-1
    # find cell above last snow cell
    # snow_upper_idx = findfirst(>(0), WaterIce)
    snow_upper_idx = max(snow_lower_idx-sum((WaterIce.>0.0) .& (grid[2:end].<=0.0)),1); #cell above last snow cell
    ρsn_0 = SnowDensity[1]

    # snow ablation
    snow_depth, snow_upper_idx = lite_snow_ablation!(Water, WaterIce, grid, SnowDensity, WaterDensity, snow_depth, snow_upper_idx, snow_lower_idx)

    # snow compaction
    snow_depth, snow_upper_idx = lite_snow_compaction!(Water, WaterIce, grid, snow_depth, snow_upper_idx, snow_lower_idx, snow_max, SnowDensity, SnowDensityMax, SnowDensityK1, SnowDensityK2, WaterDensity, current_date::DateTime)

    # snow accumulation
    snow_depth, snow_upper_idx = lite_snow_accumulation!(Water, WaterIce, grid, snow_depth, snow_upper_idx, snow_max, snow_fall, ρsn_0, WaterDensity)

    # get new air-ground interface index
    idx = snow_upper_idx+1;

    return snow_depth, idx
end

function lite_snow_ablation!(Water, WaterIce, grid, SnowDensity, WaterDensity, snow_depth, snow_upper_idx, snow_lower_idx)
    old_snow_upper_idx = snow_upper_idx
    delta_snow_depth = snow_depth - max(-grid[snow_upper_idx], 0.0) #the difference between real and grided snow depth
    #snow depletion
    if snow_upper_idx < snow_lower_idx #snow cover exists
        #route melt water downward
        #water_holding_cap = 0.1;
        meltwater = 0.0; #initlize melt trem
        snow_rm = 0.65; #together with snow grid cell size this value scales the positive degree day factor of snow melt (see e.g. VICTOR C. TSAI & XIAOZHOU RUAN 2018)
        water_hold = 0.025;
        for i=snow_upper_idx+1:snow_lower_idx-1
            water_holding_cap = water_hold*(WaterIce[i]-Water[i]); #2.5 percent of the water equivalent of the snowpack existing in the ice phase (U.S. Army, 1956; Leaf, 1966)
            #route melt water downward
            if Water[i] > water_holding_cap
                meltwater = max(0.0, Water[i]-water_holding_cap);
                WaterIce[i] = WaterIce[i] - meltwater;
                Water[i] = Water[i] - meltwater;
                WaterIce[i+1] = WaterIce[i+1] + meltwater;
                Water[i+1] = Water[i+1] + meltwater;

                #remove snow cell if it contains less ice than 50% of inituial SWE
                if (WaterIce[i]-Water[i])<=(snow_rm*SnowDensity[i]/WaterDensity[1])
                    #treat rest as additional meltwater
                    meltwater = WaterIce[i];
                    WaterIce[i] = 0.0;
                    Water[i] = 0.0;
                    WaterIce[i+1] = WaterIce[i+1] + meltwater;
                    Water[i+1] = Water[i+1] + meltwater;

                    #update upper snow indx and snow depth
                    snow_upper_idx = snow_upper_idx+1;
                    snow_depth = max(-grid[snow_upper_idx],0.0) + delta_snow_depth;
                end
            end
        end

        #last snow cell
        i = snow_lower_idx;
        #remove last snow cell and/or determine routetable water
        if (WaterIce[i]-Water[i])<=(snow_rm*SnowDensity[i]/WaterDensity[1])
            meltwater = WaterIce[i];
            WaterIce[i] = 0.0;
            Water[i] = 0.0;
            #update upper snow indx and snow depth
            snow_upper_idx = snow_upper_idx+1;
            snow_depth = max(-grid[snow_upper_idx],0.0) + delta_snow_depth;
            if snow_upper_idx == snow_lower_idx
                #if the snow cover is completely melted set snow_depth to zero
                snow_depth = 0.0
            end
        else
            meltwater = Water[i];
            #consider this water as free routetable water for upward routing
            Water[i] = Water[i]-meltwater;
            WaterIce[i] = WaterIce[i]-meltwater;
        end

        #route metlwater upward
        if (meltwater>0.0) & (snow_upper_idx<snow_lower_idx)
            #copy water and ice arrays
            WaterIce_cp = WaterIce;
            Water_cp = Water;
            j = snow_lower_idx;
            for k = snow_lower_idx:-1:snow_upper_idx

                #step over previously removed snow cells
                while (WaterIce_cp[j]==0.0) & (j>=old_snow_upper_idx+1)
                    j = j-1;
                end

                #melt water routeting
                if j>0 #only for safety reason
                    WaterIce[k] = WaterIce_cp[j] + meltwater;
                    Water[k] = Water_cp[j] + meltwater;
                    meltwater = 0.0;

                    #further route rest if cell is filled until an uptake
                    #maximum
                    UTmax = 0.9; #1.0 would be completely filled
                    if WaterIce[k]>UTmax
                        meltwater = WaterIce[k]-UTmax;
                        WaterIce[k] = UTmax;
                        Water[k] = Water[k]-meltwater;
                    end
                else #only for safety reason
                    WaterIce[k] = 0.0;
                    Water[k] = 0.0;
                end
                j = j-1;
            end
        end
    end
    #remove all water and ice rests above actual snow cover
    WaterIce[1:snow_upper_idx] .= 0.0;
    Water[1:snow_upper_idx] .= 0.0;
    return snow_depth, snow_upper_idx
end

function lite_snow_compaction!(Water, WaterIce, grid, snow_depth, snow_upper_idx, snow_lower_idx, snow_max, SnowDensity, SnowDensityMax, SnowDensityK1, SnowDensityK2, WaterDensity, current_date::DateTime)
    #snow compaction module according to Sturm et al.(2010)
    if snow_upper_idx < snow_lower_idx #snow cover exists
        #current_date = DateTime(2019,1,30)
        current_year = Dates.year(current_date)
        if current_date < DateTime(current_year,10,1)
            delta_doy = 0
        else
            delta_doy = Dates.dayofyear(DateTime(current_year,10,1))+92
        end
        ##=
        DOY = Dates.dayofyear(current_date)
        DOY_snowseason = DOY-delta_doy
        DOY_snowseason = max(-92,min(181,DOY_snowseason)) #limints according to Sturm et al. 2010
        k1 = SnowDensityK1[1]
        k2 = SnowDensityK2[1]
        SnowDensity_0 = SnowDensity[1]
        SnowDensity_max = SnowDensityMax[1]
        BulkSnowDensity_new = (SnowDensity_max - SnowDensity_0).*(1.0 - exp(-k1*snow_depth - k2*DOY_snowseason)) + SnowDensity_0
        #current snow bulk density according to gridded ice content (only takes Ice into account)
        BulkSnowDensity_current = max(SnowDensity_0, WaterDensity[1]*mean(WaterIce[snow_upper_idx+1:snow_lower_idx]))

        if (current_date >= DateTime(current_year,8,1)) && (current_date < DateTime(current_year,9,1))
            snow_depth = 0.0 #set snow depth to zero in the period from August 1 to August 31 to prevent glaciers
        else
            snow_depth = max(min(snow_depth, snow_max), 0.0); #avoid values smaller 0.0 and larger max snow depth
        end
        #new snow depth according to new snow bulk snow density this emulates compaction only
        snow_depth_compact = snow_depth*BulkSnowDensity_current/BulkSnowDensity_new
        snow_depth_compact = max(min(snow_depth_compact, snow_max), 0.0); #avoid values smaller 0.0 and larger max snow depth
        new_snow_upper_idx = argmin(abs.(snow_depth_compact.-(0. .- grid[2:end])));
        
        #remove snow cover if snow depth is lower than gridded
        if new_snow_upper_idx>snow_upper_idx
            WaterIce[snow_upper_idx:new_snow_upper_idx] .= 0.0;
            Water[snow_upper_idx:new_snow_upper_idx] .= 0.0;
            snow_upper_idx = new_snow_upper_idx;
            snow_depth = max(-grid[snow_upper_idx],0.0)
            WaterIce_current = WaterIce[snow_upper_idx+1:snow_lower_idx]
            WaterIce_new = WaterIce[snow_upper_idx+1:snow_lower_idx]*(BulkSnowDensity_new/BulkSnowDensity_current)
            WaterIce_current = WaterIce_new
            WaterIce_current[WaterIce_current.>0.9] .= 0.9
            WaterIce[snow_upper_idx+1:snow_lower_idx] = WaterIce_current
        end
    end

    return snow_depth, snow_upper_idx
end

function lite_snow_accumulation!(Water, WaterIce, grid, snow_depth, snow_upper_idx, snow_max, snow_fall, ρsn_0, ρw)
    snow_depth = snow_depth + snow_fall*ρw/ρsn_0; #increses snow depth by snow fall
    snow_depth = max(min(snow_depth, snow_max), 0.0); #avoid values smaller 0.0 and larger max snow depth
    new_snow_upper_idx = argmin(abs.(snow_depth.- (0. .-grid)));
    #build a new snow cell on top if snow depth fills a complet new cell
    if new_snow_upper_idx < snow_upper_idx
        WaterIce[new_snow_upper_idx+1:snow_upper_idx] .= ρsn_0/ρw;
        Water[new_snow_upper_idx+1:snow_upper_idx] .= 0.0;
        snow_upper_idx = new_snow_upper_idx;
    end
    return snow_depth, snow_upper_idx
end
