# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
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

# load functions
include("../src/eigen.jl")
@suppress include("../src/eigen_setup.jl")
include("../src/eigen_output.jl");

# +
drc_endemic = aggregateagegroups(readJSON(:drc_endemic_Aug2024),[0:5:15;20:10:50])
drc_kamituga = aggregateagegroups(readJSON(:drc_kamituga_Aug2024),[0:5:20;30:10:50])
burundi = aggregateagegroups(readJSON(:burundi_Aug2024),[0:5:15;20:10:50])
burundi2 = aggregateagegroups(readJSON(:burundi_midSep2024),[0:5:15;20:10:50])
burundi3= aggregateagegroups(readJSON(:burundi_Oct2024),[0:5:15;20:10:50])
burundi2m3=aggregateagegroups(readJSON(:burundi_Oct2024),[0:5:15;20:10:50])
burundi2m3.cases .-= burundi2.cases
drc_otherhz = aggregateagegroups(readJSON(:drc_otherhz_Aug2024),[0:5:15;20:10:50])
drc_kivu = aggregateagegroups(readJSON(:drc_kivu_Aug2024),[0:5:15;20:10:50])

function sexratio(p::Pyramid, ages=p.ageinterval)
    ag = in.(p.ageinterval, Ref(ages))
    cases = sum.(getindex.(p.cases,Ref(ag)))
    cases, normalize(cases,1)
end
# -

ages = (15,20,30,40)
@show sexratio(drc_endemic,ages)
@show sexratio(drc_kivu,ages)
@show sexratio(drc_kamituga,ages)
@show sexratio(drc_otherhz,ages)
@show sexratio(burundi,ages)
@show sexratio(burundi2,ages)
@show sexratio(burundi3,ages)
@show sexratio(burundi2m3,ages)

ages = (0,5,10,50)
@show sexratio(drc_endemic,ages)
@show sexratio(drc_kivu,ages)
@show sexratio(drc_kamituga,ages)
@show sexratio(drc_otherhz,ages)
@show sexratio(burundi,ages)
@show sexratio(burundi2,ages)
@show sexratio(burundi3,ages)
@show sexratio(burundi2m3,ages)

@show sexratio(drc_endemic)
@show sexratio(drc_kivu)
@show sexratio(drc_kamituga)
@show sexratio(drc_otherhz)
@show sexratio(burundi)
@show sexratio(burundi2)
@show sexratio(burundi3)
@show sexratio(burundi2m3)

gr(fontfamily="Helvetica",foreground_color_legend = nothing,background_color_legend = nothing, titlefontsize=11, tickfontsize=10, legendfontsize=8,labelfontsize=10,grid=true, tick_direction=:out,size=(400,300))

# +
module Baseline
bayesian=true
include("../src/sensitivity_setup.jl")
include("../src/sensitivity_output.jl")
end
module ActivePop_pu
include("../src/sensitivity_setup.jl")
bcond=[0.2,0.01]
include("../src/sensitivity_output.jl")
end
module ActivePop_pl
include("../src/sensitivity_setup.jl")
bcond=[0.05,0.01]
include("../src/sensitivity_output.jl")
end
module ActivePop_qu
include("../src/sensitivity_setup.jl")
bcond=[0.1,0.05]
include("../src/sensitivity_output.jl")
end
module ActivePop_ql
include("../src/sensitivity_setup.jl")
bcond=[0.1,0.002]
include("../src/sensitivity_output.jl")
end
module Assortative
include("../src/sensitivity_setup.jl")

for cm in zmb2024_sexual
    cm.misc[:modifier!]=propmixsensitivity!
end
for cm in zmb2024_sexual_b
    cm.misc[:modifier!]=propmixsensitivity!
end
assort=[0.357,0.393]
include("../src/sensitivity_output.jl")
end
module Separate
bayesian=true
include("../src/sensitivity_setup.jl")

for cm in zmb2024_sexual
    cm.misc[:modifier!]=propmixseparate!
end
for cm in zmb2024_sexual_b
    cm.misc[:modifier!]=propmixseparate!
end
assort=[0.1]
bcond=[0.01,0.001]
include("../src/sensitivity_output.jl")
end
module ByLocationVE
include("../src/sensitivity_setup.jl")
estkeys = [Symbol.(["1_" "2_"],"addmat",4:7)|>vec;:addmat_v;:addmat_w;:s_vax]
include("../src/sensitivity_output.jl")
end
# -

modulenames = ["Baseline", "ActivePop_pu", "ActivePop_pl", "ActivePop_qu", "ActivePop_ql","Assortative","Separate","ByLocationVE"]

#modules = getfield.(Ref(Main), Symbol.(modulenames))
#outputdict=Dict(modulenames.=>getfield.(modules, :frac_R0_samples))
#save("../outputs/sensitivity.jld2",outputdict)
#=for (n,m) in zip(modulenames,modules)
    outputdict = Dict(["kamituga2024_fit","otherhz2024_fit","kivu2024_fit","burundi2024_fit"].=>[m.kamituga2024_fit,m.otherhz2024_fit,m.kivu2024_fit,m.burundi2024_fit])
@suppress(save("../outputs/sexualfit_"*n*".jld2",outputdict))
end=#

loadfrac=load("../outputs/sensitivity.jld2")

medfracR0=[[median.(loc) for loc in loadfrac[key]] for key in modulenames]
xlabels=["baseline", "2x high-activity males", "0.5x high-activity males", "5x high-activity females", "0.2x high-activity females", "age assortativity", "0.1x crossover transmission", "location-specific coverage"]
plot([(reduce(hcat,(getindex.(medfracR0,x))))' for x in 1:4],color=[1 1 2 2],linestyle=[:solid :dash],ylim=(0,1),legend=:topright,label=[["all contacts" "physical only " "at home" "physical & home"].*" + sexual" fill(:none,1,12)],xtickfontsize=9,bottom_margin=4Plots.mm,xtick=(1:8,xlabels),xrotation=45,ylabel="proportion attributable to sexual contacts",size=(500,400),linealpha=0.9)
sensplot=annotate!(1.5, first.(plot_data).+0.07.+[0,0,-0.13,-0.02], text.(["N & S Kivu", "Kamituga", "Other HZ", "Burundi"],8,"Helvetica")) 

# Baseline weights
Ib_weights = pweights(normalize([0.21, 0.33,0.29,0.17],1))
base_frac=[begin avgmodel = Baseline.MixtureModel(posteriordist.(frac),Ib_weights)
        (median(avgmodel),quantile.(avgmodel,Ref([0.025,0.975])))
    end
for frac in Baseline.frac_R0_samples]

# DRC weights
Ib_weights = pweights(normalize([0.23, 0.43,0.24,0.09],1))
drc_frac=[begin avgmodel = Baseline.MixtureModel(posteriordist.(frac),Ib_weights)
        (median(avgmodel),quantile.(avgmodel,Ref([0.025,0.975])))
    end
for frac in Baseline.frac_R0_samples]

# AllIb weights
Ib_weights = pweights(normalize([5e-5, 0.01,0.07,0.92],1))
all_frac=[begin avgmodel = Baseline.MixtureModel(posteriordist.(frac),Ib_weights)
        (median(avgmodel),quantile.(avgmodel,Ref([0.025,0.975])))
    end
for frac in Baseline.frac_R0_samples]

overlay_xticks!(p::Plots.Plot, new_ticks)=xticks!(vcat.(xticks(p)[1],new_ticks[1]))
sensplot2=deepcopy(sensplot)
scatter!(sensplot2,[8.66],first.(base_frac)',series=:scatter,marker=:x,color=:black,label=:none)
scatter!(sensplot2,[9],first.(drc_frac)',series=:scatter,marker=:x,color=3,label=:none)
scatter!(sensplot2,[9.33],first.(all_frac)',series=:scatter,marker=:x,color=4,label=:none, labelfontsize=10,bottom_margin=-4Plots.mm)
overlay_xticks!(sensplot2,[Tuple(([9-1/3:1/3:9+1/3;9],["(baseline)","(weights 1)", "(weights 2)", "model weights:                                    "]))])|>savefigname("../figs/Sfigs/raw/sensitivity2.svg",save=false)

# ByLocatoinVE: estimates
Ib_weights = pweights([0.21, 0.33,0.29,0.17])
[begin
    sus_vax=ByLocationVE.MixtureModel(ByLocationVE.posteriordist.(fit[1]|>collect,11),Ib_weights)
        (median(sus_vax),quantile(sus_vax,[0.025,0.975]))|>display
    end
for fit in getfield.(Ref(ByLocationVE),[:kivu2024_fit,:kamituga2024_fit,:otherhz2024_fit,:burundi2024_fit])]

