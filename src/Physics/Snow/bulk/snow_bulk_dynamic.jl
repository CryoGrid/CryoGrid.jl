# ==== Dynamic bulk snow scheme ==== #

CryoGrid.variables(snow::BulkSnowpack, ::DynamicSnowMassBalance) = (
    Prognostic(:swe, Scalar, u"m", domain=0..Inf),
    Diagnostic(:ρsn, Scalar, u"kg/m^3", domain=0..Inf),
    Diagnostic(:por, OnGrid(Cells), domain=0..1),
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
    snowvariables(snow)...,
)

# implement ablation! for DegreeDayMelt
function ablation!(
    ::Top,
    ::SnowBC,
    snow::BulkSnowpack,
    mass::DynamicSnowMassBalance{TAcc,<:DegreeDayMelt},
    stop,
    ssnow,
) where {TAcc}
    if isactive(snow, ssnow)
        T_ub = getscalar(ssnow.T_ub) # upper boundary temperature
        dmelt = calculate_degree_day_snow_melt(mass.ablation, T_ub)
        dmelt = min(dmelt, getscalar(ssnow.swe))
        # swe flux
        @. ssnow.dswe -= dmelt
        # thickness flux
        por = getscalar(ssnow.por)
        θis = 1 - por # solid ice
        Δdsn = -dmelt / θis
        @. ssnow.dΔz += Δdsn
        # add water flux due to melt
        sat = getscalar(ssnow.sat)
        ssnow.jw[1] += dmelt - Δdsn*por*sat
    end
end

# simple linear accumulation scheme for bulk snow
function accumulation!(
    ::Top,
    snowbc::SnowBC,
    snowpack::BulkSnowpack,
    mass::DynamicSnowMassBalance{<:LinearAccumulation},
    stop,
    ssnow,
)
    rate_scale = mass.accumulation.rate_scale
    jw_snow = snowfall(snowbc, stop)
    Δswe = rate_scale*jw_snow
    @. ssnow.dswe += Δswe
    por = getscalar(ssnow.por)
    θis = 1 - por # solid ice
    Δdsn = Δswe/ θis
    @. ssnow.dΔz += Δdsn
end
