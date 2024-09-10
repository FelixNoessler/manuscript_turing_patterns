using CairoMakie
using JLD2
import CSV
import Tables

let
    @info "Create figure S6"

    ##############################
    scen1 = load("simulation_results/S6_aI_dI_part1.jld2")
    scen2 = load("simulation_results/S6_aI_dI_part2.jld2")
    oTI = CSV.read("simulation_results/oTI_boundary.txt", Tables.matrix)
    sTI = CSV.read("simulation_results/sTI_boundary.txt", Tables.matrix)
    ##############################

    color_osc = :dodgerblue1
    color_static = :lightblue

    axis_size = 500
    yticklabels = ["10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹"]

    coex_thresh = 1e-10
    cv_tresh = 0.05
    dS_vals = [0.005, 1.0]

    fig = Figure(; fontsize = 20)

    for (i,scen) in enumerate([scen1, scen2])
        aI = scen["aI_vals"]
        dI = scen["dI_vals"]

        S_survived = scen["H_density"][:, :, 1] .> coex_thresh
        I_survived = scen["H_density"][:, :, 2] .> coex_thresh
        coexistence = S_survived .&& I_survived

        cvs = scen["cvs"][:, :]
        cvs[isnan.(cvs)] .= 0.0

        osc_coexistence = ones(size(cvs))
        osc_coexistence[cvs .< cv_tresh .|| .! coexistence] .= NaN

        static_coexistence = ones(size(cvs))
        static_coexistence[cvs .> cv_tresh .|| .! coexistence] .= NaN

        only_S_survided = fill(NaN, size(cvs))
        only_S_survided[S_survived .& .! I_survived] .= 1

        Axis(fig[3, i+1];
            width = axis_size, height = axis_size,
            yscale = log10,
            xticks = 0.3:0.2:1.3,
            yticks = (10.0 .^ (-3.0:1.0), yticklabels),
            yticklabelsvisible = i == 1 ? true : false,
            xminorticksvisible = true,
            yminorticksvisible = true,
            yminorticks = IntervalsBetween(9),
            xgridvisible = false, ygridvisible = false,
            limits = (0.2, maximum(aI), minimum(dI), maximum(dI)))


        global box_element = poly!(Point2f[(0.0, 1e-3), (0.0, 3), (1.3, 3), (1.3, 1e-3)],
                                   color = Makie.LinePattern())

        heatmap!(aI, dI, only_S_survided, colormap=[:white])
        heatmap!(aI, dI, osc_coexistence, colormap=[color_osc])
        heatmap!(aI, dI, static_coexistence, colormap=[color_static])

        lines!(oTI[:, 1], oTI[:, 2]; color = :black, linestyle = :dot,
               linewidth = 3.0)
        lines!(sTI[:, 1], sTI[:, 2]; color = :black, linestyle = :dot,
               linewidth = 3.0)


        text!([0.22], [2.5]; text = i == 1 ? "A" : "B",
            fontsize = 26,
            font = "Computer Modern Sans Serif 14 Bold",
            align = (:left, :top),
            color = :black)
    end

    Label(fig[3, 4], "static TI"; valign = :top, rotation = -pi/2)
    Label(fig[3, 4], "osc. TI"; valign = :bottom, rotation = -pi/2)
    Label(fig[3, 4], "no TI        "; valign = :center, rotation = -pi/2)

    Legend(fig[3, 5],
        [[MarkerElement(color = color_static, marker = :rect, markersize = 30),
          MarkerElement(color = color_osc, marker = :rect, markersize = 30)],
          [MarkerElement(color = :white, marker = :rect, markersize = 30,
                         strokecolor = :black, strokewidth = 0.1),
          box_element],
         [LineElement(linestyle = :dot, linewidth = 3.0)]],
        [["with static dynamics", "with oscillatory dynamics"],
         ["the inferior competitor", "the superior competitor"],
         ["of the inferior competitor"]],
        ["Coexistence", "Extinction of", "Turing boundaries"],
        framevisible = true,
        gridshalign = :left,
        titlehalign = :left,
        groupgap = 20, rowgap = 5, patchlabelgap = 10, tellheight = true)


    Label(fig[1, 2:3], "dispersal rate of the superior competitor")
    Box(fig[2, 2])
    Label(fig[2, 2], L"d_S = 0.005", padding = (5,5,5,5))
    Box(fig[2, 3])
    Label(fig[2, 3], L"d_S = 1.0", padding = (5,5,5,5))

    ylabellayout = GridLayout(fig[3, 1])
    Label(ylabellayout[1, 1], L"d_I"; rotation = pi/2)
    Label(ylabellayout[2, 1], "dispersal rate of the inferior competitor"; rotation = pi/2)
    rowgap!(ylabellayout, 1, 7)

    xlabellayout = GridLayout(fig[4, 2:3])
    Label(xlabellayout[1, 1], "attack rate of the inferior competitor")
    Label(xlabellayout[1, 2], L"a_I")
    colgap!(xlabellayout, 1, 7)

    rowgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 2, 0)
    rowgap!(fig.layout, 3, 5)

    colgap!(fig.layout, 1, 5)
    colgap!(fig.layout, 2, 15)
    colgap!(fig.layout, 3, 0)

    resize_to_layout!(fig)

    save("figures/S6_aI_dI.pdf", fig)

    display(fig)

    nothing
end
