using OrdinaryDiffEq
using JLD2
using UnPack
using Statistics
using Accessors
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


function final_densities(u)
    S = u.Sx[end] + u.Sy[end]
    I = u.Ix[end] + u.Iy[end]
    return S, I
end

function calc_cv(u)
    Sx_cv = std(u.Sx) / mean(u.Sx)
    Sy_cv = std(u.Sy) / mean(u.Sy)
    Ix_cv = std(u.Ix) / mean(u.Ix)
    Iy_cv = std(u.Iy) / mean(u.Iy)
    cv = (Sx_cv + Sy_cv + Ix_cv + Iy_cv) / 4
    return cv
end

function autotroph_diff(u)
    diff = abs.(u.Ax_star .- u.Ay_star)
    return mean(diff)
end


function only_inferior!(du,u,p,t)
    ### states
    @unpack Nx_star, Ny_star, Ax_star, Ay_star, Ix_star, Iy_star = u

    ### parameters
    @unpack S, D, r, h, e, aI, dN, dA, dI = p

    ### first system
    gx = aI*Ax_star/(1+aI*h*Ax_star)
    gy = aI*Ay_star/(1+aI*h*Ay_star)

    du.Nx_star = D*S -D*Nx_star -r*Ax_star*Nx_star +dN*(Ny_star-Nx_star)
    du.Ny_star = D*S -D*Ny_star -r*Ay_star*Ny_star +dN*(Nx_star-Ny_star)
    du.Ax_star = r*Ax_star*Nx_star -gx*Ix_star -D*Ax_star +dA*(Ay_star-Ax_star)
    du.Ay_star = r*Ay_star*Ny_star -gy*Iy_star -D*Ay_star +dA*(Ax_star-Ay_star)
    du.Ix_star = e*gx*Ix_star -D*Ix_star +dI*(Iy_star-Ix_star)
    du.Iy_star = e*gy*Iy_star -D*Iy_star +dI*(Ix_star-Iy_star)
end


function both_competitors!(du,u,p,t)
    ### states
    @unpack Nx, Ny, Ax, Ay, Sx, Sy, Ix, Iy = u

    ### parameters
    @unpack S, D, r, h, e, aS, aI, dN, dA, dS, dI = p

    ### second system
    gxS = aS*Ax/(1+aS*h*Ax)
    gyS = aS*Ay/(1+aS*h*Ay)
    gxI = aI*Ax/(1+aI*h*Ax)
    gyI = aI*Ay/(1+aI*h*Ay)

    du.Nx = D*S -r*Ax*Nx -D*Nx +dN*(Ny-Nx)
    du.Ny = D*S -r*Ay*Ny -D*Ny +dN*(Nx-Ny)
    du.Ax = r*Ax*Nx -gxS*Sx -gxI*Ix -D*Ax +dA*(Ay-Ax)
    du.Ay = r*Ay*Ny -gyS*Sy -gyI*Iy -D*Ay +dA*(Ax-Ay)
    du.Sx = e*gxS*Sx -D*Sx +dS*(Sy-Sx)
    du.Sy = e*gyS*Sy -D*Sy +dS*(Sx-Sy)
    du.Ix = e*gxI*Ix -D*Ix +dI*(Iy-Ix)
    du.Iy = e*gyI*Iy -D*Iy +dI*(Ix-Iy)
end



function run_sim(; nvals)
    ###################### parameters of scenarios
    dI_vals = 10 .^ LinRange(-2.5, 0.5, nvals)
    dS_vals = 10 .^ LinRange(-2.5, 0.5, nvals)
    ######################

    ###################### outputs
    H_density = Array{Float64}(undef, nvals, nvals, 2)
    cvs = Array{Float64}(undef, nvals, nvals)
    auto_diff = Array{Float64}(undef, nvals)
    ######################

    ###################### time span of the simulations
    tspan_onlyI = (0.0, 10_000.0)
    tsave_onlyI = 9_500.0:10_000.0
    tspan = (0.0, 100_000.0)
    tsave = 95_000.0:100_000.0
    ######################

    ###################### parameters that do not change
    ### these will change: dS, dI
    p = (
        S = 5.0,
        D = 0.3,
        r = 0.5,
        h = 0.5,
        e = 0.33,
        aS = 1.3,
        aI = 1.0,
        dN = 4.0,
        dA = 0.004,
        dS = nothing,
        dI = nothing
    )
    p_onlyI = deepcopy(p)
    ######################

    ###################### inital conditions
    u0_onlyI = ca.ComponentArray(
        Nx_star=7.0, Ny_star=2.0,
        Ax_star=1.0, Ay_star=1.0,
        Ix_star=0.2, Iy_star=0.1
    )
    ######################

    ###################### Callback
    cb = local_cb(; local_threshold=1e-30)
    ######################

    progress = 0
    @time for i in eachindex(dI_vals)
        progress += 1
        actual_progress = round(progress / nvals; digits=2)
        println(actual_progress)

        @reset p_onlyI.dI = dI_vals[i]

        prob_onlyI = ODEProblem(
            only_inferior!,
            u0_onlyI,
            tspan_onlyI,
            p_onlyI;)
        sol_onlyI = solve(prob_onlyI, Vern9();
            saveat=tsave_onlyI,
            reltol=1e-12,
            abstol=1e-12,
            maxiters=1e20);

        ############### store outputs
        u_onlyI = Tables.columntable(sol_onlyI.u);
        auto_diff[i] = autotroph_diff(u_onlyI)

        ############### prepare inital conditions
        u0_final = ca.ComponentArray(;
            Nx=u_onlyI.Nx_star[end],
            Ny=u_onlyI.Ny_star[end],
            Ax=u_onlyI.Ax_star[end],
            Ay=u_onlyI.Ay_star[end],
            Ix=u_onlyI.Ix_star[end],
            Iy=u_onlyI.Iy_star[end],
            Sx=0.001, Sy=0.0001
        )

        @Threads.threads for s in eachindex(dS_vals)
            p_final = (;
                p...,
                dI=dI_vals[i],
                dS=dS_vals[s])

            prob = ODEProblem(
                both_competitors!,
                u0_final,
                tspan,
                p_final;
                callback=cb)
            sol = solve(prob, Vern9();
                saveat=tsave,
                reltol=1e-12,
                abstol=1e-12,
                maxiters=1e20);

            ############### store outputs
            u = Tables.columntable(sol.u);
            H_density[i, s, :] .= final_densities(u)
            cvs[i, s] = calc_cv(u)
        end
    end

    return (;
        H_density, cvs, auto_diff,
        dI_vals, dS_vals, nvals)
end

@info "Run simulations for figure S8"
sim_result = run_sim(; nvals)
jldsave("simulation_results/S8_S_invader.jld2"; sim_result...)
