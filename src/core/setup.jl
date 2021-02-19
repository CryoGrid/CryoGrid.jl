struct LayerData{TPvars,TDvars,TDstate,TGrids}
    pvarnames::TPvars
    dvarnames::TDvars
    dstate::TDstate
    grids::TGrids
end

struct CryoGridSetup{S,G,U,D}
    strat::S
    grid::G
    uproto::U
    layerdata::D
    CryoGridSetup(strat::S,grid::G,uproto::U,layerdata::D) where {S<:Stratigraphy,G<:Grid{Edges},U,D} =
        new{S,G,U,D}(strat,grid,uproto,layerdata)
end

"""
Generated step function (i.e. du/dt) for any arbitrary CryoGridSetup. Specialized code is generated and compiled
on the fly via the @generated macro to ensure type stability. The generated code updates each layer in the stratigraphy
in sequence, i.e for each layer 1 < i < N:

diagnostic_step!(layer i, ...)
interact!(layer i-1, ...)
prognostic_step!(layer i, ...)

Note for developers: All sections of code wrapped in quote..end blocks are generated. Code outside of quote blocks
is only executed during compilation and will not appear in the compiled version.
"""
@generated function (setup::CryoGridSetup{TStrat})(du,u,p,t) where {TStrat}
    nodetyps = nodetypes(TStrat)
    N = length(nodetyps)
    expr = Expr(:block)
    # Build ComponentArrays directly from pre-allocated axes (should be fast and allocation free)
    # TODO also need to handle autodiff types for u and du
    @>> quote
    u = ComponentArray(u,getaxes(setup.uproto))
    du = ComponentArray(du,getaxes(setup.uproto))
    end push!(expr.args) # add to expression block
    # Declare variables
    @>> quote
    strat = setup.strat
    layerdata = setup.layerdata
    end push!(expr.args)
    # Iterate over layers
    for i in 1:N-1
        n1,n2 = name(nodetyps[i]), name(nodetyps[i+1])
        # create variable names
        n1state, n2state = Symbol(n1,:state), Symbol(n2,:state)
        n1layer, n2layer = Symbol(n1,:layer), Symbol(n2,:layer)
        n1process, n2process = Symbol(n1,:process), Symbol(n2,:process)
        # generated code for layer updates
        @>> quote
        $n1state = state(layerdata.$n1, u.$n1, du.$n1, p.$n1)
        $n2state = state(layerdata.$n2, u.$n2, du.$n2, p.$n2)
        $n1layer = strat.nodes[i].layer
        $n2layer = strat.nodes[i+1].layer
        $n1process = strat.nodes[i].process
        $n2process = strat.nodes[i+1].process
        diagnosticstep!($n1layer,$n1process,$n1state)
        interact!($n1layer,$n1process,$n2layer,$n2process,$n1state,$n2state)
        prognosticstep!($n1layer,$n1process,$n1state)
        end push!(expr.args)
    end
    # Call diagnostic/prognostic update for bottom; we can assume from the preceding
    # block that "bottom" variables have already been declared.
    @>> quote
    diagnosticstep!(bottomlayer, bottomprocess, bottomstate)
    prognosticstep!(bottomlayer, bottomprocess, bottomstate)
    # all updates are done in-place; return nothing
    return nothing
    end push!(expr.args)
    # final expression block for generated function
    return expr
end

@generated function initialcondition(setup::CryoGridSetup{TStrat}) where {TStrat}
    nodetyps = nodetypes(TStrat)
    N = length(nodetyps)
    expr = Expr(:block)
    # Declare variables
    @>> quote
    u = similar(setup.uproto)
    du = similar(setup.uproto)
    end push!(expr.args) # add to expression block
    @>> quote
    strat = setup.strat
    layerdata = setup.layerdata
    end push!(expr.args)
    # Iterate over layers
    for i in 1:N
        n = name(nodetyps[i])
        # create variable names
        nstate, nlayer, nprocess = Symbol(n,:state), Symbol(n,:layer), Symbol(n,:process)
        # generated code for layer updates
        @>> quote
        $nstate = state(layerdata.$n, u.$n, du.$n, [])
        $nlayer = strat.nodes[i].layer
        $nprocess = strat.nodes[i].process
        initialcondition!($nlayer,$nstate)
        initialcondition!($nlayer,$nprocess,$nstate)
        end push!(expr.args)
    end
    @>> quote
    # call getdata to fetch underlying arrays; this avoids the performance penalty of using a custom array type
    # as the state variable in integration.
    return getdata(u), getdata(du)
    end push!(expr.args)
end

@inline @generated function state(ldata::L, u, du, params) where {TPvars,TDvars,L<:LayerData{TPvars,TDvars}}
    value(val::Val{T}) where T = T
    pvals = tuple(TPvars.parameters...,)
    dvals = tuple(TDvars.parameters...,)
    # extract symbol from Val types
    pnames = map(v -> value(v), pvals)
    dnames = map(v -> value(v), dvals)
    # construct symbols for derivative variables; assumes no existing conflicts
    dpnames = map(n -> Symbol(:d,n), pnames)
    # generate accessor expressions
    paccessors = tuple((:(u[$p]) for p in pvals)...,)
    dpaccessors = tuple((:(du[$p]) for p in pvals)...,)
    daccessors = tuple((:(ldata.dstate[$d]) for d in dvals)...,)
    # build state named tuple
    quote
    NamedTuple{tuple($(pnames...),$(dpnames...),$(dnames...),:params)}(
        tuple($(paccessors...),$(dpaccessors...),$(daccessors...),params)
    )
    end
end

function CryoGridSetup(strat::Stratigraphy, grid::Grid{Edges}, arrayproto::A=zeros()) where {A<:AbstractArray}
    layer_data = OrderedDict()
    for (i,node) in enumerate(strat)
        # determine subgrid for layer
        lo = strat.boundaries[i]
        hi = i < length(strat) ? strat.boundaries[i+1] : grid[end]
        # build subgrid using half open interval to avoid overlapping grid edges
        subgrid = grid[Interval{:closed,:open}(lo,hi)]
        # build layer metadata
        layer_data[name(node)] = buildlayer(node,subgrid,arrayproto)
    end
    # unzip list of tuples into three separate lists
    pvar_arrays, data = zip(values(layer_data)...)
    # construct named tuples containing data for each layer
    nt_parr = NamedTuple{keys(layer_data)...}(pvar_arrays...)
    nt_data = NamedTuple{keys(layer_data)...}(data...)
    # construct prototype of u (prognostic state) array
    uproto = ComponentArray(nt_parr)
    CryoGridSetup(strat,grid,uproto,nt_data)
end

function buildlayer(node::StratNode, grid::Grid{Edges}, arrayproto::A) where {A<:AbstractArray}
    nodename = name(node)
    layer, process = node.layer, node.process
    layer_vars = variables(layer)
    @assert all([isdiagnostic(var) for var in layer_vars]) "Layer variables must be diagnostic."
    process_vars = variables(layer, process)
    all_vars = tuple(layer_vars...,process_vars...)
    # check for (permissible) duplicates between variables
    groups = groupby(var -> name(var), all_vars)
    for (name,gvars) in filter(g -> length(g.second) > 1, groups)
        # if any duplicate variable deifnitions do not match, raise an error
        if !reduce(==, gvars)
            error("Found one or more conflicting definitions of $name in $gvars")
        end
    end
    diag_vars = filter(x -> isdiagnostic(x), all_vars)
    prog_vars = filter(x -> isprognostic(x), all_vars)
    # filter duplicates by converting to Sets
    diag_vars = Set(diag_vars)
    prog_vars = Set(prog_vars)
    # check for re-definition of diagnostic variables as prognostic
    diag_prog = diag_vars ∩ prog_vars
    if !isempty(diag_prog)
        @warn "Variables $(tuple(diag_prog...,)) declared as both prognostic and diagnostic. In-place modifications outside of callbacks may cause integration errors."
    end
    # prognostic takes precedence, so we remove duplicated variables from the diagnostic variable set
    setdiff!(diag_vars, diag_prog)
    # convert back to tuples
    diag_vars, prog_vars = tuple(diag_vars...), tuple(prog_vars...)
    # build prognostic and diagnostic component arrays and corresponding grid tuples
    function buildcomponent(vars::Tuple{<:Var}, grid::Grid{Edges})
        names = @>> vars map(var -> name(var))
        vars_with_names = NamedTuple{names...}(vars...)
        grids = @>> vars_with_names map(var -> var(grid))
        carray = ComponentArray(map(grid -> similar(arrayproto, size(grid))))
        return carray, grids
    end
    diag_carr, diag_grids = buildcomponent(diag_vars, grid)
    prog_carr, prog_grids = buildcomponent(prog_vars, grid)
    # merge grid tuples
    grids = merge(diag_grids, prog_grids)
    # get variable names for diagnostic and prognostic
    dvarnames = map(var->name(var),diag_vars)
    pvarnames = map(var->name(var),prog_vars)
    # return prognostic variable component array for top-level composition;
    # return LayerData with type-stable indices for all variables, diagnostic state array, and layer grids.
    return prog_carr, LayerData(fastindices(pvarnames...),fastindices(dvarnames...),diag_carr,grids)
end

export CryoGridSetup
