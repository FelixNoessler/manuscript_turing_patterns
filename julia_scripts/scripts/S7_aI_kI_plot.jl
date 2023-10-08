using CairoMakie
using JLD2

let
    @info "Create figure S7"

    ##############################
    turing_bounds = [
        # a = 0.8 [k=0], [k=2]
        [[0.0122, 1.204],[0.0123]],
        # a = 1.0 [k=0], [k=2]
        [[0.0364, 0.5391],[0.0364]],
        # a = 1.2 [k=0], [k=2]
        [[0.0518, 0.4406], [0.0518]],
        # a = 1.3 [k=0], [k=2]
        [[0.0575, 0.4181],[0.0575]],
        # a = 1.5 [k=0], [k=2]
        [[0.0519, 0.4295],[0.1099]]
    ]
    ##############################
    sim_result = load("simulation_results/S7_aI.jld2")

    aI = sim_result["aI_vals"]
    kI = sim_result["kI_vals"]

    dmaxS = sim_result["dmaxS_vals"]
    dmaxI = sim_result["dmaxI_vals"]

    facet_letter = ["A", "B", "C"]

    fig = Figure(;
    resolution=(1000, 440),
    fontsize=22)

    xticklabels = ["    0.0", "0.4", "0.8", "1.2     "]
    for sce in eachindex(dmaxS)
        Axis(fig[3, 1+sce];
            yticklabelsvisible = sce == 1 ? true : false,
            xgridvisible=false, ygridvisible=false,
            # xminorticksvisible = true,
            # xminorticks = IntervalsBetween(2),
            xticks = (0.0:0.4:1.2, xticklabels))

        coex_thresh = 1e-12
        HS_survived = sim_result["H_density"][:, :, sce, 1] .> coex_thresh
        HI_survived = sim_result["H_density"][:, :, sce, 2] .> coex_thresh
        coexistence = HS_survived .&& HI_survived

        cvs = sim_result["cvs"][:, :, sce]
        cvs[isnan.(cvs)] .= 0.0
        cv_tresh = 0.3

        osc_coexistence = ones(size(cvs))
        osc_coexistence[cvs .> cv_tresh] .= 0
        osc_coexistence[.! coexistence] .= 1
        osc_heat = ones(size(osc_coexistence))
        osc_heat[iszero.(osc_coexistence)] .= NaN

        static_coexistence = ones(size(cvs))
        static_coexistence[cvs .< cv_tresh] .= 0
        static_coexistence[.! coexistence] .= 1
        static_heat = ones(size(static_coexistence))
        static_heat[iszero.(static_coexistence)] .= NaN

        contourf!(
            aI, kI,
            osc_coexistence,
            levels=2,
            colormap=[Makie.RGB(45/255,175/255,20/255)])
        contourf!(
            aI, kI,
            static_coexistence,
            levels=2,
            colormap=[:blue1])
        text!(facet_letter[sce];
            fontsize = 26,
            font = "Computer Modern Sans Serif 14 Bold",
            position = (0.02, 1.95),
            align = (:left, :top),
            color = :black
        )

        if sce == 2
            scatter!(1.16, 0.1; marker=:star5, color=:black)
        end

    end

    ##########
    Label(fig[1, 2:4], "Maximal dispersal rates of the heterotrophs")
    top_labels = [
    L"d_{max,S} = 0.001,\; d_{max,I} = 1",
    L"d_{max,S} = 1,\; d_{max,I} = 1",
    L"d_{max,S} = 1,\; d_{max,I} = 0.001"]
    for i in 1:3
        Box(fig[2,i+1], color = :gray95, strokevisible=true)
        Label(fig[2, i+1], top_labels[i],
            alignmode=Outside(5),
            tellwidth=false,
            fontsize=24)
    end


    ##########
    Label(fig[3, 1], L"\qquad\qquad\qquad\;\;\;\;\;\; k_I",
        rotation=pi/2,
        tellheight=false,
        fontsize=24)
    Label(fig[3, 1], "Sensitivity of the inferior competitor  ",
        tellheight=false,
        rotation=pi/2)

    ##########
    Label(fig[4, 2:4], L"\qquad\qquad\qquad\;\;\;\;\;\; a_I",
        fontsize=24)
    Label(fig[4, 2:4], "Attack rate of the inferior competitor  ",
        )

    [rowgap!(fig.layout, i, s) for (i,s) in enumerate([5, 0, 10])]
    [colgap!(fig.layout, i, s) for (i,s) in enumerate([10, 5, 5])]
    rowsize!(fig.layout, 3, Aspect(2, 1))

    # resize_to_layout!(fig)
    display(fig)
    save("figures/S7_aI.pdf", fig;)
    nothing
end
