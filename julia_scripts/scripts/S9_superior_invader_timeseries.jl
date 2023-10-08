using OrdinaryDiffEq
using JLD2
using UnPack
using Statistics
using CairoMakie
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


function only_superior!(du,u,p,t)
    ### states
    @unpack Nx, Ny, Ax, Ay, HSx, HSy = u

    ### parameters
    @unpack S, D, Nh, rmax, h, e, dN, dA = p
    @unpack kS, dmaxS, aS, AcritS = p

    ### system
    rx = rmax*Nx/(Nh+Nx)
    ry = rmax*Ny/(Nh+Ny)
    gx = aS*Ax/(1+aS*h*Ax)
    gy = aS*Ay/(1+aS*h*Ay)
    dHSx = dmaxS/(1+exp(kS*(Ax-AcritS)))
    dHSy = dmaxS/(1+exp(kS*(Ay-AcritS)))

    du.Nx = D*S -D*Nx -rx*Ax +dN*(Ny-Nx)
    du.Ny = D*S -D*Ny -ry*Ay +dN*(Nx-Ny)
    du.Ax = rx*Ax -gx*HSx -D*Ax +dA*(Ay-Ax)
    du.Ay = ry*Ay -gy*HSy -D*Ay +dA*(Ax-Ay)
    du.HSx = e*gx*HSx -D*HSx -dHSx*HSx +dHSy*HSy
    du.HSy = e*gy*HSy -D*HSy -dHSy*HSy +dHSx*HSx
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


let
    @info "Create figure S9"

    ###################### time span of the simulations
    tspan_start = (0.0, 10_000.0)
    tsave_tstart = 9_500.0:10_000.0
    tspan = (0.0, 7500.0)
    tsave = 0.0:7500.0
    ######################

    ###################### parameters that do not change
    p = (
        S = 4.8,
        D = 0.3,
        Nh = 1.5,
        rmax = 0.7,
        h = 0.53,
        e = 0.33,
        dN = 1.0,
        dA = 0.001,
        aS = 1.2,
        aI = 1.0,
        AcritS = 0.3/(1.2*(0.33-0.53*0.3)),
        AcritI = 0.3/(1.0*(0.33-0.53*0.3)),
        kS = 0.0,
        kI = 2.0,
        dmaxS = 0.5,
        dmaxI = 0.2
    )
    ######################


    ###################### inital conditions
    u0_onlyI = ca.ComponentArray(
        Nx=7.0, Ny=2.0,
        Ax=1.0, Ay=1.0,
        HIx=0.2, HIy=0.1,
    )

    u0_onlyS = ca.ComponentArray(
        Nx=7.0, Ny=2.0,
        Ax=1.0, Ay=1.0,
        HSx=0.2, HSy=0.1,
    )
    ######################

    ###################### Callback
    cb = local_cb(; local_threshold=1e-30)
    ######################

    prob_onlyI = ODEProblem(
        only_inferior!,
        u0_onlyI,
        tspan_start,
        p;)
    sol_onlyI = solve(prob_onlyI, Vern9();
        saveat=tsave_tstart,
        reltol=1e-12,
        abstol=1e-12,
        maxiters=1e20);

    prob_onlyS = ODEProblem(
        only_superior!,
        u0_onlyS,
        tspan_start,
        p;)
    sol_onlyS = solve(prob_onlyS, Vern9();
        saveat=tsave_tstart,
        reltol=1e-12,
        abstol=1e-12,
        maxiters=1e20);

    ############### prepare inital conditions
    u_onlyI = Tables.columntable(sol_onlyI.u);
    u0_S_invader = ca.ComponentArray(;
        Nx=u_onlyI.Nx[end],
        Ny=u_onlyI.Ny[end],
        Ax=u_onlyI.Ax[end],
        Ay=u_onlyI.Ay[end],
        HIx=u_onlyI.HIx[end],
        HIy=u_onlyI.HIy[end],
        HSx=0.001, HSy=0.0001
    )

    u_onlyS = Tables.columntable(sol_onlyS.u);
    u0_I_invader = ca.ComponentArray(;
        Nx=u_onlyS.Nx[end],
        Ny=u_onlyS.Ny[end],
        Ax=u_onlyS.Ax[end],
        Ay=u_onlyS.Ay[end],
        HSx=u_onlyS.HSx[end],
        HSy=u_onlyS.HSy[end],
        HIx=0.001, HIy=0.0001
    )

    prob_S_invader = ODEProblem(
        both_competitors!,
        u0_S_invader,
        tspan,
        p;
        callback=cb)
    sol_S_invader = solve(prob_S_invader, Vern9();
        saveat=tsave,
        reltol=1e-12,
        abstol=1e-12,
        maxiters=1e20);


    prob_I_invader = ODEProblem(
        both_competitors!,
        u0_I_invader,
        tspan,
        p;
        callback=cb)
    sol_I_invader = solve(prob_I_invader, Vern9();
        saveat=tsave,
        reltol=1e-12,
        abstol=1e-12,
        maxiters=1e20);

    ############### store outputs
    u_S_invader = Tables.columntable(sol_S_invader.u);

    fig = Figure(resolution = (500, 400))
    ax1 = Axis(fig[1,1];
        yscale = log10,
        xticklabelsvisible = false,
        topspinevisible = false,
        rightspinevisible = false,
        xminorticksvisible = true,
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(9),
        yticks = (
            10.0 .^ (1:-1:-10),
            ["", "10⁰", "", "10⁻²", "", "10⁻⁴", "", "10⁻⁶", "", "10⁻⁸", "", "10⁻¹⁰"]),
        limits=(-500, 7500, 1e-10, 1e1))


    lines!(sol_S_invader.t, u_S_invader.Ax; color=:green,
        label=L"A_x")
    lines!(sol_S_invader.t, u_S_invader.Ay; color=:green,
        linestyle=:dash, label=L"A_y")

    lines!(sol_S_invader.t, u_S_invader.HSx; color=:red,
        label=L"S_x")
    lines!(sol_S_invader.t, u_S_invader.HSy; color=:red,
        linestyle=:dash, label=L"S_y")

    lines!(sol_S_invader.t, u_S_invader.HIx; color=:orange,
        label=L"I_x")
    lines!(sol_S_invader.t, u_S_invader.HIy; color=:orange,
        linestyle=:dash, label=L"I_y")


    t = length(u_onlyI.Ax)
    ts = -1:-1:-t
    lines!(ts, u_onlyI.Ax; color=:green)
    lines!(ts, u_onlyI.Ay; color=:green, linestyle=:dash)
    lines!(ts, u_onlyI.HIx; color=:orange)
    lines!(ts, u_onlyI.HIy; color=:orange, linestyle=:dash)

    text!("A";
        fontsize = 26,
        font = "Computer Modern Sans Serif 14 Bold",
        position = (0, 1e-8),
        align = (:left, :top),
        color = :black
    )

    Legend(fig[1:2,2], ax1; framevisible=true, framecolor=:grey)



    ax1 = Axis(fig[2,1];
        yscale = log10,
        xlabel ="time",
        topspinevisible = false,
        rightspinevisible = false,
        xminorticksvisible = true,
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(9),
        yticks = (
            10.0 .^ (1:-1:-10),
            ["", "10⁰", "", "10⁻²", "", "10⁻⁴", "", "10⁻⁶", "", "10⁻⁸", "", "10⁻¹⁰"]),
        limits=(-500, 7500, 1e-10, 1e1))

    u_I_invader = Tables.columntable(sol_I_invader.u);
    lines!(sol_I_invader.t, u_I_invader.Ax; color=:green,
        label=L"A_x")
    lines!(sol_I_invader.t, u_I_invader.Ay; color=:green,
        linestyle=:dash, label=L"A_y")

    lines!(sol_I_invader.t, u_I_invader.HSx; color=:red,
        label=L"S_x")
    lines!(sol_I_invader.t, u_I_invader.HSy; color=:red,
        linestyle=:dash, label=L"S_y")

    lines!(sol_I_invader.t, u_I_invader.HIx; color=:orange,
        label=L"I_x")
    lines!(sol_I_invader.t, u_I_invader.HIy; color=:orange,
        linestyle=:dash, label=L"I_y")


    t = length(u_onlyS.Ax)
    ts = -1:-1:-t
    lines!(ts, u_onlyS.Ax; color=:green)
    lines!(ts, u_onlyS.Ay; color=:green, linestyle=:dash)
    lines!(ts, u_onlyS.HSx; color=:red)
    lines!(ts, u_onlyS.HSy; color=:red, linestyle=:dash)


    text!("B";
        fontsize = 26,
        font = "Computer Modern Sans Serif 14 Bold",
        position = (0, 1e-8),
        align = (:left, :top),
        color = :black
    )

    Label(fig[1:2, 0],
        "density", rotation=pi/2)

    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 1, 5)

    save("figures/S9_Sinvader.pdf", fig;)
    display(fig)

    nothing
end
