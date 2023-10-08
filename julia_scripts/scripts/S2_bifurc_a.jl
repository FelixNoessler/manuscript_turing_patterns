using OrdinaryDiffEq
using JLD2
using UnPack
using Statistics
import ComponentArrays as ca
import Tables

function local_cb(;local_threshold)
    some_negative = function (u,t,integrator)
        return any(u .< local_threshold .&& u .!= 0.0)
    end

    cb_affect! = function (integrator)
        integrator.u[integrator.u .< local_threshold] .= 0.0
    end

    return DiscreteCallback(
        some_negative,
        cb_affect!;
        save_positions=(false, false))
end

function local_system!(du,u,p,t)
    ### states
    @unpack N, A, H = u

    ### parameters
    @unpack S, D, Nh, rmax, h, e, a = p

    ### first system
    r = rmax*N/(Nh+N)
    g = a*A/(1+a*h*A)

    du.N = D*S -D*N -r*A
    du.A = r*A -g*H -D*A
    du.H = e*g*H -D*H
end


function run_bifurc_a(; nvals)
    ###################### time span of the simulations
    tspan = (0.0, 100_000.0)
    tsave = 95_000.0:100_000.0
    ######################

    ###################### parameters that do not change
    p_prep = (
        S = 4.8,
        D = 0.3,
        Nh = 1.5,
        rmax = 0.7,
        h = 0.53,
        e = 0.33
    )
    ######################

    ###################### inital conditions
    u0 = ca.ComponentArray(
        N = 7.0,
        A = 1.0,
        H = 0.2)
    ######################

    ###################### Callback
    cb = local_cb(; local_threshold=1e-30)
    ######################

    ###################### prepare output
    minresult = Array{Float64, 2}(undef, nvals, 3)
    maxresult = Array{Float64, 2}(undef, nvals, 3)
    ######################

    a_vals = LinRange(0.0, 2.0, nvals)

    for i in eachindex(a_vals)
        p = (p_prep..., a = a_vals[i])


        prob = ODEProblem(
            local_system!,
            u0,
            tspan,
            p;
            callback=cb)

        sol = solve(prob, Vern9();
            saveat=tsave,
            reltol=1e-12,
            abstol=1e-12,
            maxiters=1e20
        )

        ###################### set new start densities
        u0 = sol[end] .+ 0.01

        ###################### store the result
        minresult[i, :] = minimum(sol, dims=2)
        maxresult[i, :] = maximum(sol, dims=2)
    end

    return (;
        a_vals, minresult, maxresult
    )
end

@info "Run simulations for figure S2"
sim_result = run_bifurc_a(; nvals)
jldsave("simulation_results/S2_bifurc_a.jld2"; sim_result...)
