using CairoMakie
using JLD2


let
    @info "Create figure S2"
    sim_result = load("simulation_results/S2_bifurc_a.jld2")

    a = sim_result["a_vals"]
    minresult = sim_result["minresult"]
    maxresult = sim_result["maxresult"]

    fig = Figure(resolution=(600, 350))
    ax1 = Axis(fig[1,2];
        xminorticksvisible = false,
        yminorticksvisible = false,
        xticks=0:0.5:2,
        xgridvisible=false,
        ygridvisible=false,
        topspinevisible=false,
        rightspinevisible=false,
        limits=(-0.05, 2.02, -0.1, 4.6))


    lineprop = Dict([
        :linewidth => 3])

    labels = ["Nutrients", "Autotrophs", "Heterotroph"]
    colors = [:grey, :green, :red]
    for i in 1:3
        lines!(a, minresult[:, i], label = labels[i],
            color=colors[i], linewidth=3)
        lines!(a, maxresult[:, i],
            color=colors[i], linewidth=3)
    end

    ############ Hopf-bifurcation
    lines!([1.33, 1.327], [-0.1, 5], linestyle=:dot, color=:black, linewidth=2)
    text!(["stable", "limit cycle"],
        position=[(1, 3.7), (1.82, 3.7)],
        align=(:center, :center))

    ############ legend
    Legend(fig[1,3], ax1; framevisible=false)

    ############ ylabel
    Label(fig[1,1], "Minima and maxima of densities",
        rotation=pi/2,
        tellheight=false)

    ############ xlabel
    Label(fig[2, 2], L"\qquad\qquad\;\;\;\;\,a_H",
        tellwidth=false,
        fontsize=22)
    Label(fig[2, 2], "Attack rate of the heterotroph  ",
        tellwidth=false)

    colgap!(fig.layout, 1, 5)
    colgap!(fig.layout, 2, 5)
    rowgap!(fig.layout, 1, 5)

    display(fig)

    save("figures/S2_bifurc_a.pdf", fig)

    nothing
end
