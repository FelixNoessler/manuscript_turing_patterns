using CairoMakie
using JLD2

let
    @info "Create figure S6"

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
    kS_0_result = load("simulation_results/02_dmax_dmax.jld2")
    kS_2_result = load("simulation_results/S6_kS_kI.jld2")

    dmaxS = kS_0_result["dmaxS_vals"]
    dmaxI = kS_0_result["dmaxI_vals"]

    aS_val = 1.2
    kS_val = [0.0, 2.0]

    xticklabels = ["     10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹   "]
    yticklabels = ["10⁻³", "10⁻²", "10⁻¹", "10⁰", "\n10¹"]

    fig = Figure(;
        resolution=(900, 1),
        fontsize=22)

    sce = 0
    facet_letter = ["A", "B", "C", "D"]

    sim_result = nothing
    # y = 1
    for x in 1:2
        for y in 1:2
            sim_result = isone(y) ? kS_0_result : kS_2_result
            sce = x

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

            Axis(fig[2+y, 1+x];
                backgroundcolor=(:white,1),
                xscale=log10,
                yscale=log10,
                xticks = (10.0 .^ (-3:1), xticklabels),
                yticks = (10.0 .^ (-3:1), yticklabels),
                xminorticksvisible = true,
                xminorticks = IntervalsBetween(9),
                yminorticksvisible = true,
                yminorticks = IntervalsBetween(9),
                yticklabelsvisible = x == 1 ? true : false,
                xticklabelsvisible = y == 1 ? false : true,
                xgridvisible=false, ygridvisible=false,
                limits=(1e-3, 1e1, 1e-3, 1e1)
            )

            contourf!(
                dmaxS, dmaxI,
                osc_coexistence',
                levels=2,
                colormap=[Makie.RGB(45/255,175/255,20/255)])
            contourf!(
                dmaxS, dmaxI,
                static_coexistence',
                levels=2,
                colormap=[:blue1])


            # ###### superior competitor
            aS = 1.2
            kS = y == 1 ? 0.0 : 2.0
            S_bounds = turing_bounds[aS .== [0.8, 1.0, 1.2, 1.3, 1.5]][1][kS .== [0.0, 2.0]][1]

            for b in S_bounds
                lines!([b,b], [1e-3, 1e1];
                    color=:grey,
                    linewidth=2)
            end
            if isone(length(S_bounds))
                push!(S_bounds, 1e1)
            end
            band!(S_bounds, [1e-3, 1e-3], [1e1, 1e1];
                color=(:grey, 0.3),
                linewidth=4)

            letter_index = x
            letter_index = y == 2 ? x + 2 : x

            text!(facet_letter[letter_index];
                fontsize = 26,
                font = "Computer Modern Sans Serif 14 Bold",
                position = (1.2e-3, 8),
                align = (:left, :top),
                color = :black
            )
        end
    end

    ##########
    Label(fig[1, 2:3], "Sensitivity of the inferior competitor")
    top_labels = [L"k_I = 0", L"k_I = 2"]
    for i in 1:2
        Box(fig[2,i+1], color = :gray95, strokevisible=true)
        Label(fig[2, i+1], top_labels[i],
            alignmode=Outside(5),
            tellwidth=false,
            fontsize=24)
    end

    ##########
    Label(fig[3:4, 5], "Sensitivity of the superior competitor",
        rotation=-pi/2)
    right_labels = [L"k_S = 0", L"k_S = 2"]
    for i in 1:2
        Box(fig[2+i, 4], color = :gray95, strokevisible=true)
        Label(fig[2+i, 4], right_labels[i],
            alignmode=Outside(5),
            tellwidth=true,
            fontsize=24,
            rotation=-pi/2)
    end


    ##########
    Label(fig[3:4, 1], L"\qquad\qquad\qquad\qquad\;\;\;\;\, d_{max,I}",
        rotation=pi/2,
        tellheight=false,
        fontsize=24)
    Label(fig[3:4, 1], "Max. dispersal rate of the inferior competitor        ",
        tellheight=false,
        rotation=pi/2)

    ##########
    Label(fig[5, 2:3], L"\qquad\qquad\qquad\qquad\qquad\;\; d_{max,S}",
        fontsize=24)
    Label(fig[5, 2:3], "Maximal dispersal rate of the superior competitor     ",
        )

    [rowgap!(fig.layout, i, s) for (i,s) in enumerate([5, 0, 5, 10])]
    [colgap!(fig.layout, i, s) for (i,s) in enumerate([10, 5, 0, 5])]
    rowsize!(fig.layout, 3, Aspect(2, 1))
    rowsize!(fig.layout, 4, Aspect(2, 1))

    resize_to_layout!(fig)
    display(fig)
    save("figures/S6_kS_kI.pdf", fig;)
    nothing

end
