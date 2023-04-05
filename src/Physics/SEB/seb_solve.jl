function NonlinearSolve.solve(seb::SurfaceEnergyBalance{<:Numerical}, initialstate::SEBState)
    function resid(u, p)
        state = SEBState(;
            Qh=u[1],
            Qe=u[2],
            Lstar=u[3],
            ustar=u[4],
            inputs=initialstate.inputs
        )
        res = seb(state)
        Qh = res.state.Qh
        Qe = res.state.Qe
        Lstar = res.state.Lstar
        ustar = res.state.ustar
        # Calculate and return residuals
        return @SVector[
            Qh - state.Qh,
            Qe - state.Qe,
            Lstar - state.Lstar,
            ustar - state.ustar,
        ]
    end
    u0 = @SVector[
        initialstate.Qh,
        initialstate.Qe,
        initialstate.Lstar,
        initialstate.ustar,
    ]
    prob = NonlinearProblem{false}(resid, u0)
    sol = solve(prob, seb.solscheme.solver)
    solstate = SEBState(;
        Qh=sol.u[1],
        Qe=sol.u[2],
        Lstar=sol.u[3],
        ustar=sol.u[4],
        inputs=initialstate.inputs
    )
    seb_output = seb(solstate)
    return seb_output
end
