using CairoMakie
using JLD2


let
    @info "Create figure 2"
    ##############################
    ti_attack = [1,        1.1,      1.2,      1.25,     1.3,      1.4,      1.5]
    ti1 = [0.195527, 0.217102, 0.236287, 0.244931, 0.252933, 0.26703,  0.278672]
    ti2 = [0.08636,  0.083894, 0.082123, 0.081458, 0.080921, 0.080188, 0.07984 ]
    ##############################

    sim_result = load("simulation_results/02_dS_dI_pattern.jld2")
    env_result = load("simulation_results/02_dS_dI_env_het.jld2")
    dS = sim_result["dS_vals"]
    dI = sim_result["dI_vals"]

    coex_thresh = 1e-10
    cv_tresh = 0.2

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

    Label(fig[0, 1:2], "Self-organised pattern formation",
        font = "TeX Gyre Heros Makie Bold", fontsize = 22)
    Label(fig[1, 1], "osc. TI"; halign = :left)
    Label(fig[1, 1], "static TI"; halign = :right)
    Label(fig[1, 1], "          no TI"; halign = :center)

    Label(fig[2, 2], "static TI"; valign = :top, rotation = -pi/2)
    Label(fig[2, 2], "osc. TI"; valign = :bottom, rotation = -pi/2)
    Label(fig[2, 2], "no TI        "; valign = :center, rotation = -pi/2)


    Axis(fig[2, 1];
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
        limits = (10 ^ -2.5, 10 ^ 0.5, 10 ^ - 2.5,  10 ^ 0.5),
        xlabel = "dispersal rate of the superior competitor dS",
        ylabel = "dispersal rate of the inferior competitor dI")
    heatmap!(
        dS, dI,
        osc_coexistence',
        colormap=[:blue])
    heatmap!(
        dS, dI,
        static_coexistence',
        colormap=[:lightblue])
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
    text!([1e-2, 1], [0.45, 1]; text = ["Bet hedging", "Maladaptive\ndispersal"],
        align = (:center, :center))
    text!([0.0035], [3]; text = "A",
        fontsize = 26,
        font = "Computer Modern Sans Serif 14 Bold",
        align = (:left, :top),
        color = :black)

    Legend(fig[1:2, 3],
           [[MarkerElement(color = (:black, 0.3), marker = :rect, markersize = 30,
                markerstrokewidth = 1.5, markerstrokecolor = :red)],
           [MarkerElement(color = :lightblue, marker = :rect, markersize = 30),
           MarkerElement(color = :blue, marker = :rect, markersize = 30)],
           [LineElement(linestyle = :dash), LineElement(linestyle = :dot)]],
           [["Coexistence only possible\nwith self-organised\npattern formation"],
           ["with static dynamic", "with oscillatory dynamic"],
            ["of superior competitor",
            "of inferior competitor"]],
            ["", "Coexistence", "Turing boundaries"],
           framevisible = true,
           gridshalign = :left, valign = :top,
           titlehalign = :left,
           groupgap = 20, rowgap = 5, patchlabelgap = 10)

    Axis(fig[2, 2:3];
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
        colormap=[:blue])
    heatmap!(
        dS, dI,
        env_static_coexistence',
        colormap=[:lightblue])
    vlines!(ti1[findfirst(ti_attack .== 1.3)]; color = :black, linestyle = :dash)
    vlines!(ti2[findfirst(ti_attack .== 1.3)]; color = :black, linestyle = :dash)
    hlines!(ti1[findfirst(ti_attack .== 1.0)]; color = :black, linestyle = :dot)
    hlines!(ti2[findfirst(ti_attack .== 1.0)]; color = :black, linestyle = :dot)
    text!([0.004], [3]; text = "B",
        fontsize = 24,
        font = "Computer Modern Sans Serif 14 Bold",
        align = (:left, :top),
        color = :black)

    rowgap!(fig.layout, 1, 0)
    rowgap!(fig.layout, 2, 0)
    colgap!(fig.layout, 1, 0)
    colgap!(fig.layout, 2, 25)
    resize_to_layout!(fig)

    display(fig)

    save("figures/02_dS_dI.png", fig; px_pet_unit = 4)

    nothing
end
