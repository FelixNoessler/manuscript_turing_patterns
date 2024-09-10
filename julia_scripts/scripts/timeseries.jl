using OrdinaryDiffEq
using UnPack
using Statistics
using CairoMakie
import Tables
import ComponentArrays as ca

function only_superior!(du, u, p, t)
    ### states
    @unpack Nx_star, Ny_star, Ax_star, Ay_star, Sx_star, Sy_star = u

    ### parameters
    @unpack S, D, r, h, e, aS, dN, dA, dS = p

    ### first system
    gx = aS*Ax_star/(1+aS*h*Ax_star)
    gy = aS*Ay_star/(1+aS*h*Ay_star)

    du.Nx_star = D*S -D*Nx_star -r*Ax_star*Nx_star +dN*(Ny_star-Nx_star)
    du.Ny_star = D*S -D*Ny_star -r*Ay_star*Ny_star +dN*(Nx_star-Ny_star)
    du.Ax_star = r*Ax_star*Nx_star -gx*Sx_star -D*Ax_star +dA*(Ay_star-Ax_star)
    du.Ay_star = r*Ay_star*Ny_star -gy*Sy_star -D*Ay_star +dA*(Ax_star-Ay_star)
    du.Sx_star = e*gx*Sx_star -D*Sx_star +dS*(Sy_star-Sx_star)
    du.Sy_star = e*gy*Sy_star -D*Sy_star +dS*(Sx_star-Sy_star)
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

let
    ######################
    tspan_onlyS = (0.0, 10_000.0)
    tsave_onlyS = 9_500.0:10_000.0
    tspan_onlyI = (0.0, 10_000.0)
    tsave_onlyI = 9_500.0:10_000.0
    tsave = 99_000.0:100_000.0    # change!
    tspan = (0.0, maximum(tsave))
    S_invader = false             # change!

    ######################
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
        dS = 1.0,  # change!
        dI = 0.15 # change!
    )

    #########################################
    u0_final = nothing
    if S_invader
        u0_onlyI = ca.ComponentArray(
            Nx_star=7.0, Ny_star=2.0,
            Ax_star=1.0, Ay_star=1.0,
            Ix_star=0.2, Iy_star=0.1)
        prob_onlyI = ODEProblem(only_inferior!, u0_onlyI, tspan_onlyI, p;)
        sol_onlyI = solve(prob_onlyI, Vern9(); saveat=tsave_onlyI,
                        reltol=1e-12, abstol=1e-12, maxiters=1e20);
        u_onlyI = Tables.columntable(sol_onlyI.u);
        u0_final = ca.ComponentArray(;
            Nx=u_onlyI.Nx_star[end],
            Ny=u_onlyI.Ny_star[end],
            Ax=u_onlyI.Ax_star[end],
            Ay=u_onlyI.Ay_star[end],
            Ix=u_onlyI.Ix_star[end],
            Iy=u_onlyI.Iy_star[end],
            Sx=0.001, Sy=0.0001)
    else
        u0_onlyS = ca.ComponentArray(
            Nx_star=7.0, Ny_star=2.0,
            Ax_star=1.0, Ay_star=1.0,
            Sx_star=0.2, Sy_star=0.1)
        prob_onlyS = ODEProblem(only_superior!, u0_onlyS, tspan_onlyS, p;)
        sol_onlyS = solve(prob_onlyS, Vern9(); saveat=tsave_onlyS,
                        reltol=1e-12, abstol=1e-12, maxiters=1e20);
        u_onlyS = Tables.columntable(sol_onlyS.u);
        u0_final = ca.ComponentArray(;
            Nx=u_onlyS.Nx_star[end],
            Ny=u_onlyS.Ny_star[end],
            Ax=u_onlyS.Ax_star[end],
            Ay=u_onlyS.Ay_star[end],
            Sx=u_onlyS.Sx_star[end],
            Sy=u_onlyS.Sy_star[end],
            Ix=0.001, Iy=0.0001)
    end
    #########################################

    ###################### Callback
    cb = local_cb(; local_threshold=1e-30)
    ######################

    prob = ODEProblem(both_competitors!, u0_final, tspan, p; )
    sol = solve(prob, Vern9(); saveat=tsave, reltol=1e-12, abstol=1e-12, maxiters=1e20);
    u = Tables.columntable(sol.u);


    begin
        pattern = Dict(
            :Nx => :solid,
            :Ny => :dash,
            :Ax => :solid,
            :Ay => :dash,
            :Sx => :solid,
            :Sy => :dash,
            :Ix => :solid,
            :Iy => :dash
        )

        color = Dict(
            :Nx => :black,
            :Ny => :black,
            :Ax => :green,
            :Ay => :green,
            :Sx => :red,
            :Sy => :red,
            :Ix => :orange,
            :Iy => :orange
        )

        fig = Figure()
        ax1 = Axis(fig[1,1];
                   xlabel = "time", ylabel = "density",
                   title = "dS = $(p[:dS]), dI = $(p[:dI])",
                   topspinevisible = false, rightspinevisible = false,
                   xgridvisible = false, ygridvisible = false,
                   yscale = log10)

        for k in keys(u)
            u[k][u[k] .< 1e-30] .= NaN

            lines!(sol.t, u[k]; linestyle = pattern[k], color = color[k],
                   label = string(k))
        end
        Legend(fig[1, 2], ax1; framevisible = false)


        save("figures/timeseries.png", fig)
        fig
    end
end
