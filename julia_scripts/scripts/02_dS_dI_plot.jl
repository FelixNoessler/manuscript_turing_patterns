using CairoMakie
using JLD2


let
    @info "Create figure 2"
    ##############################
    # turing_bounds = [
    #     # a = 0.8 [k=0], [k=2]
    #     [[0.0122, 1.204],[0.0123]],
    #     # a = 1.0 [k=0], [k=2]
    #     [[0.0364, 0.5391],[0.0364]],
    #     # a = 1.2 [k=0], [k=2]
    #     [[0.0518, 0.4406], [0.0518]],
    #     # a = 1.3 [k=0], [k=2]
    #     [[0.0575, 0.4181],[0.0575]],
    #     # a = 1.5 [k=0], [k=2]
    #     [[0.0519, 0.4295],[0.1099]]
    # ]
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

    # HS_survived_env = env_result["H_density"][:, :, 1] .> coex_thresh
    # HI_survived_env = env_result["H_density"][:, :, 2] .> coex_thresh
    # coexistence_env = HS_survived_env .&& HI_survived_env

    # pattern_coex = coexistence .&& .! coexistence_env
    # only_hetcoex = coexistence_env .&& .! coexistence

    fig = Figure(; fontsize=22)


    axis_size = 600
    xticklabels = ["     10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹   "]
    yticklabels = ["10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹"]
    Axis(fig[1, 1];
        width = axis_size, height = axis_size,
        title = "Self-organised pattern formation",
        backgroundcolor=(:white,1),
        xscale=log10,
        yscale=log10,
        xticks = (10.0 .^ (-3:1), xticklabels),
        yticks = (10.0 .^ (-3:1), yticklabels),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(9),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(9),
        xgridvisible=false, ygridvisible=false,
        limits=(1e-3, 1e1, 1e-3, 1e1),
        xlabel = "dispersal rate of the superior competitor dS",
        ylabel = "dispersal rate of the inferior competitor dI",
        aspect = DataAspect()
    )

    heatmap!(
        dS, dI,
        osc_coexistence',
        colormap=[:blue])
    heatmap!(
        dS, dI,
        static_coexistence',
        colormap=[:lightblue])

    Axis(fig[1, 2];
        width = axis_size, height = axis_size,
        title = "Predefined habitat heterogeneity",
        backgroundcolor=(:white,1),
        xscale=log10,
        yscale=log10,
        xticks = (10.0 .^ (-3:1), xticklabels),
        yticks = (10.0 .^ (-3:1), yticklabels),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(9),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(9),
        xgridvisible=false, ygridvisible=false,
        limits=(1e-3, 1e1, 1e-3, 1e1),
        xlabel = "dispersal rate of the superior competitor dS",
        ylabel = "dispersal rate of the inferior competitor dI",
        aspect = DataAspect()
    )

    heatmap!(
        dS, dI,
        env_osc_coexistence',
        colormap=[:blue])
    heatmap!(
        dS, dI,
        env_static_coexistence',
        colormap=[:lightblue])


    Legend(fig[1, 3],
           [MarkerElement(color = :blue, marker = :rect, markersize = 30),
           MarkerElement(color = :lightblue, marker = :rect, markersize = 30)],
           [" Coexistence\n with oscillatory dynamic",
            " Coexistence\n with static dynamic"],
           framevisible = false)


    resize_to_layout!(fig)
    # heatmap!(
    #     dS,
    #     dI,
    #     .! pattern_coex',
    #     colormap=(:greys, 0.2))
    # contour!(
    #     dS,
    #     dI,
    #     .! pattern_coex',
    #     levels=2,
    #     color=(:white, 1))

    ###### superior competitor
    # S_bounds = turing_bounds[aS_val .== [0.8, 1.0, 1.2, 1.3, 1.5]][1][kS_val .== [0.0, 2.0]][1]

    # for b in S_bounds
    #     lines!([b,b], [1e-3, 1e1];
    #         color=:grey,
    #         linewidth=2)
    # end
    # band!(S_bounds, [1e-3, 1e-3], [1e1, 1e1];
    #     color=(:grey, 0.3),
    #     linewidth=4)

    # text!(facet_letter[u];
    #     fontsize = 26,
    #     font = "Computer Modern Sans Serif 14 Bold",
    #     position = (1.2e-3, 8),
    #     align = (:left, :top),
    #     color = :black
    # )



    ##########
    # Label(fig[1, 2:3], "Sensitivity of the heterotrophs")
    # top_labels = [L"k_S = 0,\; k_I = 0", L"k_S = 0,\;k_I = 2"]
    # for i in 1:2
    #     Box(fig[2,i+1], color = :gray95, strokevisible=true)
    #     Label(fig[2, i+1], top_labels[i],
    #         alignmode=Outside(5),
    #         tellwidth=false,
    #         fontsize=24)
    # end

    # ##########
    # Label(fig[3, 1], L"\text{ }\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad d_{max,I}",
    #     rotation=pi/2,
    #     tellheight=false,
    #     fontsize=24)
    # Label(fig[3, 1], "Max. dispersal rate of the inferior competitor        ",
    #     tellheight=false,
    #     rotation=pi/2)

    # ##########
    # Label(fig[4, 2:3],
    #     L"\text{ }\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\;\;\; d_{max,S}",
    #     fontsize=24)
    # Label(fig[4, 2:3], "Maximal dispersal rate of the superior competitor     ")

    # [rowgap!(fig.layout, i, s) for (i,s) in enumerate([5, 0, 10])]
    # [colgap!(fig.layout, i, s) for (i,s) in enumerate([10, 5])]
    # rowsize!(fig.layout, 3, Aspect(2, 1))

    # resize_to_layout!(fig)
    display(fig)

    save("figures/02_dS_dI.png", fig;)

    nothing
end
