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


function final_densities(u)
    HS = u.HSx[end] + u.HSy[end]
    HI = u.HIx[end] + u.HIy[end]
    return HS, HI
end

function calc_cv(u)
    HSx_cv = std(u.HSx) / mean(u.HSx)
    HSy_cv = std(u.HSy) / mean(u.HSy)
    HIx_cv = std(u.HIx) / mean(u.HIx)
    HIy_cv = std(u.HIy) / mean(u.HIy)
    cv = (HSx_cv + HSy_cv + HIx_cv + HIy_cv) / 4
    return cv
end

function autotroph_diff1(u)
    diff = abs.(u.Ax .- u.Ay)
    return mean(diff)
end


function only_inferior!(du,u,p,t)
    ### states
    @unpack Nx, Ny, Ax, Ay, HIx, HIy = u

    ### parameters
    @unpack S, D, Nh, rmax, h, e, dN, dA = p
    @unpack kI, dmaxI, aI, AcritI = p

    ### system
    rx = rmax*Nx/(Nh+Nx)
    ry = rmax*Ny/(Nh+Ny)
    gx = aI*Ax/(1+aI*h*Ax)
    gy = aI*Ay/(1+aI*h*Ay)
    dHIx = dmaxI/(1+exp(kI*(Ax-AcritI)))
    dHIy = dmaxI/(1+exp(kI*(Ay-AcritI)))

    du.Nx = D*S -D*Nx -rx*Ax +dN*(Ny-Nx)
    du.Ny = D*S -D*Ny -ry*Ay +dN*(Nx-Ny)
    du.Ax = rx*Ax -gx*HIx -D*Ax +dA*(Ay-Ax)
    du.Ay = ry*Ay -gy*HIy -D*Ay +dA*(Ax-Ay)
    du.HIx = e*gx*HIx -D*HIx -dHIx*HIx +dHIy*HIy
    du.HIy = e*gy*HIy -D*HIy -dHIy*HIy +dHIx*HIx
end


function both_competitors!(du,u,p,t)
    ### states
    @unpack Nx, Ny, Ax, Ay, HSx, HSy, HIx, HIy = u

    ### parameters
    @unpack S, D, Nh, rmax, h, e, dN, dA = p
    @unpack kS, kI, dmaxS, dmaxI, aS, aI, AcritS, AcritI = p

    ### second system
    rx = rmax*Nx/(Nh+Nx)
    ry = rmax*Ny/(Nh+Ny)
    gxS = aS*Ax/(1+aS*h*Ax)
    gyS = aS*Ay/(1+aS*h*Ay)
    gxI = aI*Ax/(1+aI*h*Ax)
    gyI = aI*Ay/(1+aI*h*Ay)
    AcritS = D/(aS*(e-h*D))
    AcritI = D/(aI*(e-h*D))
    dHSx = dmaxS/(1+exp(kS*(Ax-AcritS)))
    dHSy = dmaxS/(1+exp(kS*(Ay-AcritS)))
    dHIx = dmaxI/(1+exp(kI*(Ax-AcritI)))
    dHIy = dmaxI/(1+exp(kI*(Ay-AcritI)))

    du.Nx = D*S -rx*Ax -D*Nx +dN*(Ny-Nx)
    du.Ny = D*S -ry*Ay -D*Ny +dN*(Nx-Ny)
    du.Ax = rx*Ax -gxS*HSx -gxI*HIx -D*Ax +dA*(Ay-Ax)
    du.Ay = ry*Ay -gyS*HSy -gyI*HIy -D*Ay +dA*(Ax-Ay)
    du.HSx = e*gxS*HSx -D*HSx -dHSx*HSx +dHSy*HSy
    du.HSy = e*gyS*HSy -D*HSy -dHSy*HSy +dHSx*HSx
    du.HIx = e*gxI*HIx -D*HIx -dHIx*HIx +dHIy*HIy
    du.HIy = e*gyI*HIy -D*HIy -dHIy*HIy +dHIx*HIx
end



function run_sim_S10_patterns(; nvals)
    ###################### parameters of scenarios
    dmaxI_vals = 10 .^ LinRange(-3, 1, nvals)
    dmaxS_vals = 10 .^ LinRange(-3, 1, nvals)
    kI_vals = [0.0, 2.0]
    nscenarios = length(kI_vals)
    ######################

    ###################### outputs
    H_density = Array{Float64}(undef, nvals, nvals, nscenarios, 2)
    cvs = Array{Float64}(undef, nvals, nvals, nscenarios)
    auto_diff = Array{Float64}(undef, nvals, nscenarios)
    ######################

    ###################### time span of the simulations
    tspan_onlyI = (0.0, 10_000.0)
    tsave_onlyI = 9_500.0:10_000.0
    tspan = (0.0, 100_000.0)
    tsave = 95_000.0:100_000.0
    ######################

    ###################### parameters that do not change
    ### these will change: kI, dmaxS, dmaxI
    p = (
        S = 4.8,
        D = 0.3,
        Nh = 1.5,
        rmax = 0.7,
        h = 0.53,
        e = 0.33,
        dN = 1.0,
        dA = 0.001,
        kS = 0.0,
        aS = 1.2,
        aI = 1.0,
        AcritS = 0.3/(1.2*(0.33-0.53*0.3)),
        AcritI = 0.3/(1.0*(0.33-0.53*0.3))
    )
    ######################

    ###################### inital conditions
    u0_onlyI = ca.ComponentArray(
        Nx=7.0, Ny=2.0,
        Ax=1.0, Ay=1.0,
        HIx=0.2, HIy=0.1,
    )
    ######################

    ###################### Callback
    cb = local_cb(; local_threshold=1e-30)
    ######################

    for sce in 1:nscenarios
        progress = 0
        @info "Scenario: $sce"
        kI = kI_vals[sce]

        @time for i in eachindex(dmaxI_vals)
            progress += 1
            actual_progress = round(progress / nvals; digits=2)
            println(actual_progress)

            p_onlyI = (;
                p...,
                dmaxI=dmaxI_vals[i],
                kI)
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
            auto_diff[i, sce] = autotroph_diff1(u_onlyI)

            ############### prepare inital conditions
            u0_final = ca.ComponentArray(;
                Nx=u_onlyI.Nx[end],
                Ny=u_onlyI.Ny[end],
                Ax=u_onlyI.Ax[end],
                Ay=u_onlyI.Ay[end],
                HIx=u_onlyI.HIx[end],
                HIy=u_onlyI.HIy[end],
                HSx=0.001, HSy=0.0001
            )

            @Threads.threads for s in eachindex(dmaxS_vals)
                p_final = (;
                    p...,
                    dmaxI=dmaxI_vals[i],
                    dmaxS=dmaxS_vals[s],
                    kI)
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
                H_density[i, s, sce, :] .= final_densities(u)
                cvs[i, s, sce] = calc_cv(u)
            end
        end
    end

    return (;
        H_density, cvs, auto_diff,
        dmaxI_vals, dmaxS_vals, nvals,
        kI_vals, p)
end

@info "Run simulations for figure S10"
sim_result = run_sim_S10_patterns(; nvals)
jldsave("simulation_results/S10_pattern.jld2"; sim_result...)
