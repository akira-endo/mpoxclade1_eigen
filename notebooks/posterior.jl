# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.6
#   kernelspec:
#     display_name: Julia Multithreads 1.9.3
#     language: julia
#     name: julia-multithreads-1.9
# ---

using Pkg
Pkg.activate("../")

gr(fontfamily="Helvetica",foreground_color_legend = nothing,background_color_legend = nothing, titlefontsize=11, tickfontsize=10, legendfontsize=8,labelfontsize=10,grid=true, tick_direction=:out,size=(400,300))

# load functions
include("../src/eigen.jl")
include("../src/eigen_setup.jl")
include("../src/eigen_output.jl");

loadfit = load("../outputs/sexualfit_main.jld2")

for key in keys(loadfit)
    println(key)
    println("empirical contact matrix")
    loadfit[key][1]|>collect.|>chainof.|>ess .|>display
    println("synthetic contact matrix")
    loadfit[key][3]|>collect.|>chainof .|>ess.|>display
end

# ## Model averaging

loadwaic = load("../outputs/sexualfit_main_waic.jld2")

waic_kivu = waic.(loadwaic["kivu2024_fit"][1]|>collect)

waic_kamituga = waic.(loadwaic["kamituga2024_fit"][1]|>collect)

waic_otherhz = waic.(loadwaic["otherhz2024_fit"][1]|>collect)

waic_burundi = waic.(loadwaic["burundi2024_fit"][1]|>collect)

endemic_ll = [cm.misc[:opt].nll(mean(chainof(cm)|>Array,dims=1)) for cm in (loadwaic["endemic2015_24_fit"].zmb_fit|>collect.|>first)]

modelweights = endemic_ll .+first.(waic_kivu)./2
normalize(exp.(.-(modelweights.-minimum(modelweights))),1)

modelweights = endemic_ll .+first.(waic_kivu)./2
normalize(exp.(.-(modelweights.-minimum(modelweights))),1)

waic_drc = (first.(waic_kivu).+first.(waic.(loadwaic["kamituga2024_fit"][1]|>collect)).+first.(waic.(loadwaic["otherhz2024_fit"][1]|>collect)))./2

modelweights_drc = endemic_ll .+first.(waic_drc)./2
normalize(exp.(.-(modelweights_drc.-minimum(modelweights_drc))),1)

waic_allIb = (first.(waic_kivu).+first.(waic.(loadwaic["kamituga2024_fit"][1]|>collect)).+first.(waic.(loadwaic["otherhz2024_fit"][1]|>collect)))./2 .+ first.(waic.(loadwaic["burundi2024_fit"][1]|>collect))

modelweights_allIb = endemic_ll .+first.(waic_allIb)./2
normalize(exp.(.-(modelweights_allIb.-minimum(modelweights_allIb))),1)

# ## Plot posterior

# +
Random.seed!(123)
parlabel =["p[".* string.([15,20,30,40]).*"–".*string.([19,29,39,49]).*"]"; "q[".* string.([15,20,30,40]).*"–".*string.([19,29,39,49]).*"]"; "w_F / λ_C"; "w_M / λ_C"]
modelweights=[0.21,0.33,0.29,0.17]
dictkeys = ["kivu2024_fit", "kamituga2024_fit", "otherhz2024_fit", "burundi2024_fit"]
mixposteriors = [MixtureModel(posteriorjointdist.(cms),modelweights) for cms in getindex.(Ref(loadfit),dictkeys).|>first.|>collect]

for (i,key) in enumerate(dictkeys)#, ctype in (:all, :phys, :home, :physhome)
    figname="correlation_".*key.*"_average.svg"
corrplot(rand(mixposteriors[i],1000)|>permutedims,size=(2000,1800),tickfontsize=15,labelfontsize=20,label=parlabel,left_margin=10Plots.mm,xrotation=90,markercolor=cgrad([:gray20, :gray20]),fc=:thermal)|>savefigname("../figs/Sfigs/raw/correlation/".*figname,save=true)
end

# +
parlabel =["p[".* string.([15,20,30,40]).*"–".*string.([19,29,39,49]).*"]"; "q[".* string.([15,20,30,40]).*"–".*string.([19,29,39,49]).*"]"; "w_F"; "w_M"]

for key in keys(loadfit), ctype in (:all, :phys, :home, :physhome)
    figname="correlation_".*key.*"_".*string(ctype).*".svg"
corrplot((chainof(getfield(loadfit[key][1],ctype))|>Array),size=(2000,2000),label=parlabel,left_margin=10Plots.mm,xrotation=90,markercolor=cgrad([:gray20, :gray20]),fc=:thermal)|>savefigname("../figs/Sfigs/raw/correlation/".*figname,save=false)
end
# -


