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

function autotroph_diff(u)
    diff = abs.(u.Ax_star .- u.Ay_star)
    return mean(diff)
end


function only_superior!(du,u,p,t)
    ### states
    @unpack Nx_star, Ny_star, Ax_star, Ay_star, HSx_star, HSy_star = u

    ### parameters
    @unpack S, D, Nh, rmax, h, e, dN, dA = p
    @unpack kS, dmaxS, aS, AcritS = p

    ### first system
    rx = rmax*Nx_star/(Nh+Nx_star)
    ry = rmax*Ny_star/(Nh+Ny_star)
    gx = aS*Ax_star/(1+aS*h*Ax_star)
    gy = aS*Ay_star/(1+aS*h*Ay_star)
    dHx = dmaxS/(1+exp(kS*(Ax_star-AcritS)))
    dHy = dmaxS/(1+exp(kS*(Ay_star-AcritS)))

    du.Nx_star = D*S -D*Nx_star -rx*Ax_star +dN*(Ny_star-Nx_star)
    du.Ny_star = D*S -D*Ny_star -ry*Ay_star +dN*(Nx_star-Ny_star)
    du.Ax_star = rx*Ax_star -gx*HSx_star -D*Ax_star +dA*(Ay_star-Ax_star)
    du.Ay_star = ry*Ay_star -gy*HSy_star -D*Ay_star +dA*(Ax_star-Ay_star)
    du.HSx_star = e*gx*HSx_star -D*HSx_star -dHx*HSx_star +dHy*HSy_star
    du.HSy_star = e*gy*HSy_star -D*HSy_star -dHy*HSy_star +dHx*HSx_star
end


function env_heterogeneity!(du,u,p,t)
    ### states
    @unpack Nx_star, Ny_star, Ax_star, Ay_star, HSx_star, HSy_star = u
    @unpack Nx, Ny, Ax, Ay, HSx, HSy, HIx, HIy = u

    ### parameters
    @unpack S, D, Nh, rmax, h, e, dN, dA = p
    @unpack kS, kI, dmaxS, dmaxI, aS, aI, AcritS, AcritI = p

    ### first system
    rx = rmax*Nx_star/(Nh+Nx_star)
    ry = rmax*Ny_star/(Nh+Ny_star)
    gx = aS*Ax_star/(1+aS*h*Ax_star)
    gy = aS*Ay_star/(1+aS*h*Ay_star)
    dHx = dmaxS/(1+exp(kS*(Ax_star-AcritS)))
    dHy = dmaxS/(1+exp(kS*(Ay_star-AcritS)))

    du.Nx_star = D*S -D*Nx_star -rx*Ax_star +dN*(Ny_star-Nx_star)
    du.Ny_star = D*S -D*Ny_star -ry*Ay_star +dN*(Nx_star-Ny_star)
    du.Ax_star = rx*Ax_star -gx*HSx_star -D*Ax_star +dA*(Ay_star-Ax_star)
    du.Ay_star = ry*Ay_star -gy*HSy_star -D*Ay_star +dA*(Ax_star-Ay_star)
    du.HSx_star = e*gx*HSx_star -D*HSx_star -dHx*HSx_star +dHy*HSy_star
    du.HSy_star = e*gy*HSy_star -D*HSy_star -dHy*HSy_star +dHx*HSx_star

    ### calculation of the supply concentration
    Seff_x = S +dN*(Ny_star-Nx_star)/D
    Seff_y = S +dN*(Nx_star-Ny_star)/D

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

    du.Nx = D*Seff_x -rx*Ax -D*Nx
    du.Ny = D*Seff_y -ry*Ay -D*Ny
    du.Ax = rx*Ax -gxS*HSx -gxI*HIx -D*Ax +dA*(Ay-Ax)
    du.Ay = ry*Ay -gyS*HSy -gyI*HIy -D*Ay +dA*(Ax-Ay)
    du.HSx = e*gxS*HSx -D*HSx -dHSx*HSx +dHSy*HSy
    du.HSy = e*gyS*HSy -D*HSy -dHSy*HSy +dHSx*HSx
    du.HIx = e*gxI*HIx -D*HIx -dHIx*HIx +dHIy*HIy
    du.HIy = e*gyI*HIy -D*HIy -dHIy*HIy +dHIx*HIx
end



function run_sim(; nvals)
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
    tspan_onlyS = (0.0, 10_000.0)
    tsave_onlyS = 9_500.0:10_000.0
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
        dA = 0.001,
        dN = 1.0,
        kS = 0.0,
        aS = 1.2,
        aI = 1.0,
        AcritS = 0.3/(1.2*(0.33-0.53*0.3)), # D/(aS*(e-h*D))
        AcritI = 0.3/(1.0*(0.33-0.53*0.3))
    )
    ######################

    ###################### inital conditions
    u0_onlyS = ca.ComponentArray(
        Nx_star=7.0, Ny_star=2.0,
        Ax_star=1.0, Ay_star=1.0,
        HSx_star=0.2, HSy_star=0.1,
    )
    ######################

    ###################### Callback
    cb = local_cb(; local_threshold=1e-30)
    ######################

    for sce in 1:nscenarios
        progress = 0
        @info "Scenario: $sce"

        kI = kI_vals[sce]

        @time for s in eachindex(dmaxS_vals)
            progress += 1
            actual_progress = round(progress / nvals; digits=2)
            println(actual_progress)

            p_onlyS = (;
                p...,
                dmaxS=dmaxS_vals[s])
            prob_onlyS = ODEProblem(
                only_superior!,
                u0_onlyS,
                tspan_onlyS,
                p_onlyS;)
            sol_onlyS = solve(prob_onlyS, Vern9();
                saveat=tsave_onlyS,
                reltol=1e-12,
                abstol=1e-12,
                maxiters=1e20);

            ############### store outputs
            u_onlyS = Tables.columntable(sol_onlyS.u);
            auto_diff[s, sce] = autotroph_diff(u_onlyS)

            ############### prepare inital conditions
            u0_final = ca.ComponentArray(;
                Nx_star=u_onlyS.Nx_star[end],
                Ny_star=u_onlyS.Ny_star[end],
                Ax_star=u_onlyS.Ax_star[end],
                Ay_star=u_onlyS.Ay_star[end],
                HSx_star=u_onlyS.HSx_star[end],
                HSy_star=u_onlyS.HSy_star[end],
                Nx=u_onlyS.Nx_star[end],
                Ny=u_onlyS.Ny_star[end],
                Ax=u_onlyS.Ax_star[end],
                Ay=u_onlyS.Ay_star[end],
                HSx=u_onlyS.HSx_star[end],
                HSy=u_onlyS.HSy_star[end],
                HIx=0.001, HIy=0.0001
            )

            @Threads.threads for i in eachindex(dmaxI_vals)
                p_final = (;
                    p...,
                    dmaxI=dmaxI_vals[i],
                    dmaxS=dmaxS_vals[s],
                    kI)
                prob = ODEProblem(
                    env_heterogeneity!,
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
        kI_vals)
end


@info "Run simulations for figure 2, part 1"
sim_result = run_sim(; nvals)
jldsave("simulation_results/02_env_het.jld2"; sim_result...)
