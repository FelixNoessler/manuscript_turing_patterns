using CairoMakie
using JLD2

let
    @info "Create figure 3"
    ##############################
    ti_attack = [1.0, 1.3]
    ti1 = [0.308134,  0.252933]
    ti2 = [0.048367,  0.080921]
    ##############################

    color_osc = :dodgerblue1
    color_static = :lightblue

    sim_result = load("simulation_results/03_dS_dI_pattern.jld2")
    env_result = load("simulation_results/03_dS_dI_env_het.jld2")
    dS = sim_result["dS_vals"]
    dI = sim_result["dI_vals"]

    coex_thresh = 1e-10
    cv_tresh = 0.05

    #### pattern formation
    HS_survived = sim_result["H_density"][:, :, 1] .> coex_thresh
    HI_survived = sim_result["H_density"][:, :, 2] .> coex_thresh
    coexistence = HS_survived .&& HI_survived

    cvs = sim_result["cvs"][:, :]
    cvs[isnan.(cvs)] .= 0.0

    osc_coexistence = ones(size(cvs))
    osc_coexistence[cvs .< cv_tresh .|| .! coexistence] .= NaN
    static_coexistence = ones(size(cvs))
    static_coexistence[cvs .> cv_tresh .|| .! coexistence] .= NaN


    #### environmental heterogeneity
    env_HS_survived = env_result["H_density"][:, :, 1] .> coex_thresh
    env_HI_survived = env_result["H_density"][:, :, 2] .> coex_thresh
    env_coexistence = env_HS_survived .&& env_HI_survived

    env_cvs = env_result["cvs"][:, :]
    env_cvs[isnan.(env_cvs)] .= 0.0

    env_osc_coexistence = ones(size(env_cvs))
    env_osc_coexistence[env_cvs .< cv_tresh .|| .! env_coexistence] .= NaN
    env_static_coexistence = ones(size(env_cvs))
    env_static_coexistence[env_cvs .> cv_tresh .|| .! env_coexistence] .= NaN

    env_coexistence_plot = fill(1, size(env_cvs)...)
    env_coexistence_plot[.! env_coexistence .&& coexistence] .= 2

    fig = Figure(; fontsize = 20)


    axis_size = 650
    xticklabels = ["     10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹   "]
    yticklabels = ["10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹"]

    Label(fig[1, 2], "Self-organised pattern formation",
        font = "TeX Gyre Heros Makie Bold", fontsize = 22)
    Label(fig[2, 2], "                       osc. TI"; halign = :left)
    Label(fig[2, 2], "static TI               "; halign = :right)
    Label(fig[2, 2], "          no TI"; halign = :center)

    Label(fig[3, 3], "              static TI"; valign = :top, rotation = -pi/2)
    Label(fig[3, 3], "osc. TI               "; valign = :bottom, rotation = -pi/2)
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
    contour!(
        dS, dI,
        env_coexistence_plot',
        colormap=[:red])
    contourf!(
        dS, dI,
        env_coexistence_plot',
        colormap=[(:black, 0.0), (:black, 0.3)])
    vlines!(ti1[findfirst(ti_attack .== 1.3)]; color = :black, linestyle = :dash)
    vlines!(ti2[findfirst(ti_attack .== 1.3)]; color = :black, linestyle = :dash)
    hlines!(ti1[findfirst(ti_attack .== 1.0)]; color = :black, linestyle = :dot)
    hlines!(ti2[findfirst(ti_attack .== 1.0)]; color = :black, linestyle = :dot)
    text!([1e-2, 1], [0.45, 1]; text = ["bet-hedging", "maladaptive\ndispersal"],
        align = (:center, :center))
    text!([0.0035], [3]; text = "A",
        fontsize = 26,
        font = "Computer Modern Sans Serif 14 Bold",
        align = (:left, :top),
        color = :black)

    ##############################
    #A: d_S=0.005, d_I=0.018, B : d_S=1, d_I=0.47, C: d_S=0.005, d_I=0.1, D: d_S=1, d_I=0.13.
    markersize = 20; color = :black; strokecolor = :white; strokewidth = 1;
    scatter!([0.005], [0.1]; color, strokecolor, strokewidth, markersize, marker = :rect)
    scatter!([0.005], [0.018]; color, strokecolor,strokewidth,  markersize, marker = :pentagon)
    scatter!([1], [0.47]; color, strokecolor, strokewidth, markersize, marker = :diamond)
    scatter!([1], [0.13]; color, strokecolor, strokewidth, markersize, marker = :utriangle)
    ##############################

    legend_layout = GridLayout(fig[3, 4]; valign = :top)
    Legend(legend_layout[1, 1],
        [MarkerElement(color = (:black, 0.3), marker = :rect, markersize = 30,
            markerstrokewidth = 1.5, markerstrokecolor = :red)],
        ["Coexistence only possible\nwith self-organised\npattern formation due to\nheterogeneity modulation"],
        framevisible = false,
        patchlabelgap = 10, tellheight = true)
    Legend(legend_layout[2, 1],
           [[MarkerElement(color = color_static, marker = :rect, markersize = 30),
           MarkerElement(color = color_osc, marker = :rect, markersize = 30)],
           [LineElement(linestyle = :dash),
            LineElement(linestyle = :dot)]],
           [["with static dynamics", "with oscillatory dynamics"],
            ["of the superior competitor",
            "of the inferior competitor"]],
            ["Coexistence", "Turing boundaries"],
           framevisible = false,
           gridshalign = :left,
           titlehalign = :left,
           groupgap = 20, rowgap = 5, patchlabelgap = 10, tellheight = true)
    Box(legend_layout[1:2, 1], color = :transparent, strokecolor = :black)

    Axis(fig[3, 3:4];
        title = "Environmental habitat\nheterogeneity",
        alignmode = Outside(),
        width = axis_size / 3, height = axis_size / 3,
        valign = :bottom,
        topspinevisible = false,
        rightspinevisible = false,
        xscale=log10,
        yscale=log10,
        xticks = (10.0 .^ (-3.0:1.0), xticklabels),
        yticks = (10.0 .^ (-3.0:1.0), yticklabels),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(9),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(9),
        xgridvisible = false, ygridvisible = false,
        limits=(10 ^ -2.5, 10 ^ 0.5, 10 ^ - 2.5,  10 ^ 0.5),
        xlabel = L"d_S",
        ylabel = L"d_I")
    heatmap!(
        dS, dI,
        env_osc_coexistence',
        colormap=[color_osc])
    heatmap!(
        dS, dI,
        env_static_coexistence',
        colormap=[color_static])
    vlines!(ti1[findfirst(ti_attack .== 1.3)]; color = :black, linestyle = :dash)
    vlines!(ti2[findfirst(ti_attack .== 1.3)]; color = :black, linestyle = :dash)
    hlines!(ti1[findfirst(ti_attack .== 1.0)]; color = :black, linestyle = :dot)
    hlines!(ti2[findfirst(ti_attack .== 1.0)]; color = :black, linestyle = :dot)
    text!([0.004], [3]; text = "B",
        fontsize = 24,
        font = "Computer Modern Sans Serif 14 Bold",
        align = (:left, :top),
        color = :black)


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
    rowgap!(legend_layout, 1, 5)

    resize_to_layout!(fig)

    save("figures/03_dS_dI.pdf", fig)

    display(fig)

    nothing
end
