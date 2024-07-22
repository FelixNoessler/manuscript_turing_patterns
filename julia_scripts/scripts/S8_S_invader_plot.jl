using CairoMakie
using JLD2

let
    @info "Create figure S8"
    ##############################
    ti_attack = [1.0, 1.3]
    ti1 = [0.308134,  0.252933]
    ti2 = [0.048367,  0.080921]
    ##############################

    color_osc = :dodgerblue1
    color_static = :lightblue
    coex_thresh = 1e-10
    cv_tresh = 0.2

    #### S invader
    sim_result = load("simulation_results/S8_S_invader.jld2")
    HS_survived = sim_result["H_density"][:, :, 1] .> coex_thresh
    HI_survived = sim_result["H_density"][:, :, 2] .> coex_thresh
    coexistence = HS_survived .&& HI_survived

    cvs = sim_result["cvs"][:, :]
    cvs[isnan.(cvs)] .= 0.0

    osc_coexistence = ones(size(cvs))
    osc_coexistence[cvs .< cv_tresh .|| .! coexistence] .= NaN
    static_coexistence = ones(size(cvs))
    static_coexistence[cvs .> cv_tresh .|| .! coexistence] .= NaN

    #### I invader
    sim_result_I_invader = load("simulation_results/03_dS_dI_pattern.jld2")

    HS_survived_I_invader = sim_result_I_invader["H_density"][:, :, 1] .> coex_thresh
    HI_survived_I_invader  = sim_result_I_invader["H_density"][:, :, 2] .> coex_thresh
    coexistence_I_invader = HS_survived_I_invader .&& HI_survived_I_invader

    #### diff between results
    coexistence_only_I_invader = fill(NaN, size(cvs))
    coexistence_only_I_invader[.! coexistence .& coexistence_I_invader] .= 1.0

    coexistence_only_S_invader = fill(NaN, size(cvs))
    coexistence_only_S_invader[coexistence .& .! coexistence_I_invader] .= 1.0


    ####
    dS = sim_result["dS_vals"]
    dI = sim_result["dI_vals"]

    fig = Figure(; fontsize = 20)

    axis_size = 650
    xticklabels = ["     10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹   "]
    yticklabels = ["10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹"]

    Label(fig[2, 2], "                  osc. TI"; halign = :left)
    Label(fig[2, 2], "static TI                "; halign = :right)
    Label(fig[2, 2], "          no TI"; halign = :center)

    Label(fig[3, 3], "        static TI"; valign = :top, rotation = -pi/2)
    Label(fig[3, 3], "osc. TI                  "; valign = :bottom, rotation = -pi/2)
    Label(fig[3, 3], "no TI        "; valign = :center, rotation = -pi/2)

    Axis(fig[3, 2];
        width = axis_size, height = axis_size,
        topspinevisible = false,
        rightspinevisible = false,
        xscale = log10,
        yscale = log10,
        xticks = (10.0 .^ (-3.0:1.0), xticklabels),
        yticks = (10.0 .^ (-3.0:1.0), yticklabels),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(9),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(9),
        xgridvisible = false, ygridvisible = false,
        limits = (10 ^ -2.5, 10 ^ 0.5, 10 ^ - 2.5,  10 ^ 0.5))

    heatmap!(
        dS, dI,
        osc_coexistence',
        colormap=[color_osc])
    heatmap!(
        dS, dI,
        static_coexistence',
        colormap=[color_static])

    vlines!(ti1[findfirst(ti_attack .== 1.3)]; color = :black, linestyle = :dash)
    vlines!(ti2[findfirst(ti_attack .== 1.3)]; color = :black, linestyle = :dash)
    hlines!(ti1[findfirst(ti_attack .== 1.0)]; color = :black, linestyle = :dot)
    hlines!(ti2[findfirst(ti_attack .== 1.0)]; color = :black, linestyle = :dot)
    text!([1e-2, 1], [0.45, 1]; text = ["bet-hedging", "maladaptive\ndispersal"],
        align = (:center, :center))

    heatmap!(
        dS, dI,
        coexistence_only_I_invader',
        colormap=[:orange])
    heatmap!(
        dS, dI,
        coexistence_only_S_invader',
        colormap=[:red])

    Legend(fig[3, 4],
           [[MarkerElement(color = :red, marker = :rect, markersize = 30),
             MarkerElement(color = :orange, marker = :rect, markersize = 30)],
            [MarkerElement(color = color_static, marker = :rect, markersize = 30),
             MarkerElement(color = color_osc, marker = :rect, markersize = 30)],
            [LineElement(linestyle = :dash), LineElement(linestyle = :dot)]],
           [["the superior competitor\nis the invader", "the inferior competitor\nis the invader"],
            ["with static dynamics", "with oscillatory dynamics"],
            ["of superior competitor",
            "of inferior competitor"]],
            ["Coexistence only if" ,"Coexistence in both\ninvasion scenarios", "Turing boundaries"],
           framevisible = true,
           gridshalign = :left,
           titlehalign = :left,
           groupgap = 20, rowgap = 5, patchlabelgap = 10, tellheight = true)


    ylabellayout = GridLayout(fig[3, 1])
    xlabellayout = GridLayout(fig[4, 2])

    Label(ylabellayout[2, 1], "dispersal rate of the inferior competitor",
        rotation = pi/2)
    Label(ylabellayout[1, 1], L"d_I",
        rotation = pi/2)
    Label(xlabellayout[1, 1], "dispersal rate of the superior competitor")
    Label(xlabellayout[1, 2], L"d_S")

    rowgap!(fig.layout, 1, 0)
    rowgap!(fig.layout, 2, 0)
    rowgap!(fig.layout, 3, 5)
    colgap!(fig.layout, 1, 5)
    colgap!(fig.layout, 2, 0)
    colgap!(fig.layout, 3, 25)

    rowgap!(ylabellayout, 1, 7)
    colgap!(xlabellayout, 1, 7)


    resize_to_layout!(fig)

    save("figures/S8_S_invader.png", fig; px_per_unit = 10)

    display(fig)

    nothing
end
