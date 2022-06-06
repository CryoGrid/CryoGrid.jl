Base.@kwdef struct BucketScheme{Tfc,Tkwsat} <: WaterFlowParameterization
    θfc::Tfc = Param(0.5, domain=0..1)
    kw_sat::Tkwsat = Param(1e-5, domain=0..Inf, units=u"m/s")
end
fieldcapacity(::SubSurface, water::WaterFlow{<:BucketScheme}) = water.para.θfc
hydraulicconductivity(::SubSurface, water::WaterFlow{<:BucketScheme}) = water.para.kw_sat

CryoGrid.variables(::WaterFlow{<:BucketScheme}) = (
    Prognostic(:θwi, OnGrid(Cells)),
    Diagnostic(:θw, OnGrid(Cells)),
    Diagnostic(:θw_sat, OnGrid(Cells)),
    Diagnostic(:jw, OnGrid(Edges), u"1/s"),
    Diagnostic(:kw, OnGrid(Edges), u"m/s"),
    Diagnostic(:kwc, OnGrid(Cells), u"m/s"),
)
function CryoGrid.diagnosticstep!(sub::SubSurface, water::WaterFlow{<:BucketScheme}, state)
    kw_sat = hydraulicconductivity(sub, water)
    por = porosity(sub, water, state)
    @. state.θw_sat = state.θw / por
    @. state.kwc = state.θw_sat*kw_sat
    harmonicmean!(state.kw[2:end-1], state.kwc, Δ(state.grids.kw))
    state.kw[1] = state.kwc[1] 
    state.kw[end] = state.kwc[1]
end
function CryoGrid.prognosticstep!(sub::SubSurface, water::WaterFlow{<:BucketScheme}, state)
    # set downward fluxes according to kw and field capacity
    θfc = fieldcapacity(sub, water)
    @inbounds for i in 2:length(state.kw)-1 # note that the index is over grid *edges*
        let θwi_up = state.θwi[i-1], # cell above edge i
            θw_up = state.θw[i-1],
            θi_up = θwi_up - θw_up,
            θwi_lo = state.θwi[i], # cell below edge i
            por_lo = porosity(sub, water, state, i);
            # downward advective flux due to gravity
            jw = -state.kw[i]*(θw_up >= θfc)*10^(-7*θi_up)
            # limit flux based on i) available water in cell above and ii) free pore space in cell below
            state.jw[i] = min(min(jw, θw_up), min(por_lo - θwi_lo, zero(θwi_lo)))
        end
    end
    @. state.dθwi += (state.jw[2:end] - state.jw[1:end-1]) / Δ(state.grids.jw)
end
