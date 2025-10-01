# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

]activate ../

include("../src/networks.jl")
include("../src/utils.jl")

gr(fontfamily="Helvetica",foreground_color_legend = nothing,background_color_legend = nothing, titlefontsize=11, tickfontsize=10, legendfontsize=8,labelfontsize=10,grid=true, tick_direction=:out,size=(400,300))

# illustration
Random.seed!(123)
gsmall = union(binet(500, [50,10], 4, 0, EXPONENTS), binet(500, [50,50], 0, 0.5, EXPONENTS))
@time gp = graphplot(gsmall,curves=false,markercolor = repeat([1,4],inner=nv(gsmall)÷2),markersize=0.1,size=(450,400))
gp|>savefigname("../figs/Sfigs/raw/network1.svg", save=false)

Random.seed!(123)
gsmall = union(binet(500, [50,10], 2, 0,  EXPONENTS), binet(500, [50,50], 0,　0.5, EXPONENTS))
@time gp = graphplot(gsmall,curves=false,markercolor = repeat([1,4],inner=nv(gsmall)÷2),markersize=0.1,size=(450,400))
gp|>savefigname("../figs/Sfigs/raw/network2.svg", save=false)

# ## Network creation and calculate network characteristics

include("../src/networks.jl")

Random.seed!(123)
binets1 = create_networks1(100)
binets2 = create_networks2(100);

ns1, ns2 = analyse_networks_and_summarise_results(binets1, binets2);

include("../src/networks.jl")
plot_simple_adj_ngm_R0(ns1, ns2) |> savefigname("../figs/Sfigs/network1_2_spectrum.svg", save = true)
plot_simple_rel_R0(ns1, ns2) |> savefigname("../figs/Sfigs/network1_2_persentsexual.svg", save = true)

Random.seed!(123)
binets1_assort = create_networks1_assort(100)
binets2_assort = create_networks2_assort(100);

ns3, ns4 = analyse_networks_and_summarise_results(binets1_assort, binets2_assort);

# +
include("../src/networks.jl")
pl_R01 = plot_adj_ngm_eigen_R0(ns1, ns2; lim=(8,18))
pl_rel1 = plot_rel_R0_sex(ns1, ns2; );
pl_R02 = plot_adj_ngm_eigen_R0(ns3, ns4; lim=(8,18), version=2)
pl_rel2 = plot_rel_R0_sex(ns3, ns4; version=2);
pls = [pl_R01, pl_rel1, pl_R02, pl_rel2]

for (i,pl) in enumerate(pls)
     pl |> savefigname("../figs/Sfigs/network1_6_panel$(i).svg", save = true)
end
plot(pls..., layout=(2,2), size=(880, 800), bottom_margin=5Plots.mm, left_margin=5Plots.mm) 
# -

# ## Proportion of each agent

include("../src/networks.jl")
groupedbar_compare_prop_ngm_adj(ns1, ns2) |> 
    savefigname("../figs/Sfigs/networkeigen.svg", save = true)

groupedbar_compare_prop_ngm_adj(ns3, ns4; version=2)

# ## Weighted networks

include("../src/networks.jl")
Random.seed!(123)
@time binets1w = create_weighted_network1(100; weight_factor=1.0, community_weight=1.0/7.0)
@time binets2w = create_weighted_network2(100; weight_factor=1.0, community_weight=1.0/7.0);

ns1w, ns2w = analyse_networks_and_summarise_results(binets1w, binets2w);

include("../src/networks.jl")
pl_R01w = plot_adj_ngm_eigen_R0(ns1w, ns2w; lim=(8,18), version=3)
pl_rel1w = plot_rel_R0_sex(ns1w, ns2w; lim=(0.0, 61), version=3);
for (i,pl) in enumerate([pl_R01w, pl_rel1w])
     pl |> savefigname("../figs/Sfigs/network1_6_panel$(i+4).svg", save = true)
end
plot(pl_R01w, pl_rel1w, size=(800,400), bottom_margin=5Plots.mm, left_margin=5Plots.mm) 











