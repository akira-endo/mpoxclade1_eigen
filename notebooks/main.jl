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
include("../src/eigen_setup.jl")
include("../src/eigen_output.jl")

const serial=(Int[],Int[])

# +
# write data
# include("../src/data.jl");

# load data
#datasets = readJSON(dir = "../data/pyramid/")

# load MCMC results
loadfit = load("../outputs/sexualfit_main.jld2")
# -

gr(fontfamily="Helvetica",foreground_color_legend = nothing,background_color_legend = nothing, titlefontsize=11, tickfontsize=10, legendfontsize=8,labelfontsize=10,grid=true, tick_direction=:out,size=(400,300))

# Fit to 2011-15 data in Tshuapa
tshuapaplot = plot(tshuapa_h2hag,color=:black)
tshuapa2015_fit = output_fit(
    tshuapa_h2hag,
    zmb_skeleton = zmb2015,
    drc_skeleton = drc2015,
    dataplots = tshuapaplot,
    bayesian=true);
zmb2015_fit = tshuapa2015_fit.zmb_fit;

# Fit to 2024 data in DRC (only required for validation in the next cell)
endemicplot = plot(plot(drc_endemic_ag),ylim=(0,0.35),xtickfontsize=9);
endemic2024_fit = output_fit(
    drc_endemic_ag,
    zmb_skeleton = zmb2024,
    drc_skeleton = drc2024,
    dataplots = endemicplot,
    bayesian=true
    );

endemicplot = plot(plot(drc_endemic_ag),ylim=(0,0.35),xtickfontsize=9);
endemic2024_validate = output_validate(
    drc_endemic_ag,
    zmb_skeleton = endemic2024_fit.zmb_fit,
    zmb_ref = tshuapa2015_fit.zmb_fit,
    drc_skeleton = endemic2024_fit.drc_fit,
    drc_ref = tshuapa2015_fit.drc_fit,
    dataplots = endemicplot,
    bayesian=true
    );


# Fit to both 2015 and 2024 clade Ia data
Random.seed!(1)
endemic2015_24_fit = output_fit(
    [tshuapa_h2hag,drc_endemic_ag];
    zmb_skeleton = [zmb2015,zmb2024],
    drc_skeleton = [drc2015,drc2024],
    dataplots = [tshuapaplot,endemicplot],
    bayesian=true);
zmb2015_24_fit = endemic2015_24_fit.zmb_fit
drc2015_24_fit = endemic2015_24_fit.drc_fit;

# +
drc_kamituga = aggregateagegroups(readJSON(:drc_kamituga_Aug2024),[0:5:20;30:10:50])
burundi = aggregateagegroups(readJSON(:burundi_Aug2024),[0:5:15;20:10:50])
drc_otherhz = aggregateagegroups(readJSON(:drc_otherhz_Aug2024),[0:5:15;20:10:50])
otherhzplot = plot(plot(drc_otherhz,color = [1 2]),ylim=(0,0.5),xtickfontsize=9)
burundiplot = plot(plot(drc_otherhz,color = [1 2]),ylim=(0,0.5),xtickfontsize=9)
kamitugaplot = plot(plot(drc_kamituga,color=[1 2]),ylim=(0,0.5),size=(600,300))

drc_kivu = aggregateagegroups(readJSON(:drc_kivu_Aug2024),[0:5:15;20:10:50])
kivuplot = plot(plot(drc_kivu),ylim=(0,0.5),xtickfontsize=9);
# -


zmb2024_sexual = addsexualcontact!(deepcopy(zmb2024),[15;20:10:40]; modifier! = propmix!, countrycode = "COD",year = 2024);
drc2024_sexual = addsexualcontact!(deepcopy(drc2024),[15;20:10:40]; modifier! = propmix!, countrycode = "COD",year = 2024);

kamituga2024_fit = output_sexual_fit(
    drc_kamituga,
    zmb_sexual_skeleton = zmb2024_sexual,
    zmb_ref=last.(zmb2015_24_fit|>collect),
    drc_sexual_skeleton = drc2024_sexual,
    drc_ref=last.(drc2015_24_fit|>collect),
    dataplots = collapseplot(drc_kamituga),
    estkeys = [Symbol.(["1_" "2_"],"addmat",4:7)|>vec;:addmat_v;:addmat_w],
    bayesian=true,
    load = loadfit["kamituga2024_fit"]
  );

kivu2024_fit = output_sexual_fit(
    drc_kivu,
    zmb_sexual_skeleton = zmb2024_sexual,
    zmb_ref=last.(zmb2015_24_fit|>collect),
    drc_sexual_skeleton = drc2024_sexual,
    drc_ref=last.(drc2015_24_fit|>collect),
    dataplots = collapseplot(drc_kivu),
    estkeys = [Symbol.(["1_" "2_"],"addmat",4:7)|>vec;:addmat_v;:addmat_w],
    bayesian=true,
    load=loadfit["kivu2024_fit"]
    );

otherhz2024_fit = output_sexual_fit(
    drc_otherhz,
    zmb_sexual_skeleton = zmb2024_sexual,
    zmb_ref=last.(zmb2015_24_fit|>collect),
    drc_sexual_skeleton = drc2024_sexual,
    drc_ref=last.(drc2015_24_fit|>collect),
    dataplots = collapseplot(drc_otherhz),
    estkeys = [Symbol.(["1_" "2_"],"addmat",4:7)|>vec;:addmat_v;:addmat_w],
    bayesian=true,
    load=loadfit["otherhz2024_fit"]
    );

# +
zmb2024_sexual_b = addsexualcontact!(deepcopy(bdi2024),[15;20:10:40]; modifier! = propmix!, countrycode = "BDIC",year = 2024);
drc2024_sexual_b = addsexualcontact!(deepcopy(bdi_s2024),[15;20:10:40]; modifier! = propmix!, countrycode = "BDIC",year = 2024);

burundi = aggregateagegroups(readJSON(:burundi_Oct2024),[0:5:15;20:10:50])
burundi0 = aggregateagegroups(readJSON(:burundi_midSep2024),[0:5:15;20:10:50])
burundi.cases .-= broadcast.(min,burundi0.cases,burundi.cases)

zmb_ref = deepcopy(last.(zmb2015_24_fit|>collect))
drc_ref = deepcopy(last.(drc2015_24_fit|>collect))

for cmt in zmb_ref
    cmt.parameters[:s_vax].=1-(1-cmt.parameters[:s_vax][])*1.0
end
for cmt in drc_ref
    cmt.parameters[:s_vax].=1-(1-cmt.parameters[:s_vax][])*1.0
end

burundi2024_fit = output_sexual_fit(
    burundi,
    zmb_sexual_skeleton = zmb2024_sexual_b,
    zmb_ref=zmb_ref,
    drc_sexual_skeleton = drc2024_sexual_b,
    drc_ref=drc_ref,
    dataplots = collapseplot(burundi),
    estkeys = [Symbol.(["1_" "2_"],"addmat",4:7)|>vec;:addmat_v;:addmat_w],
    bayesian=true,
    load=loadfit["burundi2024_fit"]
    );
# -
kamituga2024_skeleton=deepcopy(kamituga2024_fit)
for cm in kamituga2024_skeleton.zmb_fit for el in cm.addmat el.*=0. end end
for cm in kamituga2024_skeleton.drc_fit for el in cm.addmat el.*=0. end end
kamituga2024_validate = output_validate(
    drc_kamituga,
    zmb_skeleton = kamituga2024_skeleton.zmb_fit,
    zmb_ref = zmb2015_24_fit|>collect.|>first,
    drc_skeleton = kamituga2024_skeleton.drc_fit,
    drc_ref = drc2015_24_fit|>collect.|>first,
    dataplots = drc_kamituga|>collapseplot,
    bayesian=true
    );


otherhz2024_skeleton=deepcopy(otherhz2024_fit)
for cm in otherhz2024_skeleton.zmb_fit for el in cm.addmat el.*=0. end end
for cm in otherhz2024_skeleton.drc_fit for el in cm.addmat el.*=0. end end
otherhz2024_validate = output_validate(
    drc_otherhz,
    zmb_skeleton = otherhz2024_skeleton.zmb_fit,
    zmb_ref = zmb2015_24_fit|>collect.|>first,
    drc_skeleton = otherhz2024_skeleton.drc_fit,
    drc_ref = drc2015_24_fit|>collect.|>first,
    dataplots = drc_otherhz|>collapseplot,
    bayesian=true
    );

kivu2024_skeleton=deepcopy(kivu2024_fit)
for cm in kivu2024_skeleton.zmb_fit for el in cm.addmat el.*=0. end end
for cm in kivu2024_skeleton.drc_fit for el in cm.addmat el.*=0. end end
kivu2024_validate = output_validate(
    drc_kivu,
    zmb_skeleton = kivu2024_skeleton.zmb_fit,
    zmb_ref = zmb2015_24_fit|>collect.|>first,
    drc_skeleton = kivu2024_skeleton.drc_fit,
    drc_ref = drc2015_24_fit|>collect.|>first,
    dataplots = drc_kivu|>collapseplot,
    bayesian=true
    );

burundi2024_skeleton=deepcopy(burundi2024_fit)
for cm in burundi2024_skeleton.zmb_fit for el in cm.addmat el.*=0. end end
for cm in burundi2024_skeleton.drc_fit for el in cm.addmat el.*=0. end end
burundi2024_validate = output_validate(
    burundi,
    zmb_skeleton = burundi2024_skeleton.zmb_fit,
    zmb_ref = zmb_ref,
    drc_skeleton = burundi2024_skeleton.drc_fit,
    drc_ref = drc_ref,
    dataplots = burundi|>collapseplot,
    bayesian=true
    );


# save results
outputdict = Dict(["endemic2015_24_fit","kamituga2024_fit","otherhz2024_fit","kivu2024_fit","burundi2024_fit"].=>[endemic2015_24_fit, kamituga2024_fit,otherhz2024_fit,kivu2024_fit,burundi2024_fit])
#@suppress(save("../outputs/sexualfit_main2.jld2",outputdict))

# ## Main figures

plot(endemic2015_24_fit.zmb_plot[1], title="Tshuapa, 2011-2015") |>savefigname("../figs/fig1/raw/tshuapa.svg", save = false)

plot(endemic2015_24_fit.zmb_plot[2], title="Endemic provinces, 2024",legend=:none) |>savefigname("../figs/fig1/raw/endemicprovs.svg", save = false)

plot(kivu2024_fit.zmb_plot, kivu2024_validate.zmb_fit|>collect, title="North and South Kivu, 2024",size=(780,300),xrotation=0,bottom_margin=4Plots.PlotMeasures.mm,ylim=(0,0.301),
color =[1 1 2 2],linestyle=[:solid :dash],linealpha=[fill(0.2,7);0],label=:none,left_margin=3Plots.PlotMeasures.mm,xlabel="age (male, female)") |>savefigname("../figs/fig1/raw/kivu.svg", save = false)

plot(kamituga2024_fit.zmb_plot, kamituga2024_validate.zmb_fit|>collect, title="Kamituga health zone",size=(400,300),xrotation=45,bottom_margin=3Plots.PlotMeasures.mm,legend=:none,ylim=(0,0.301),
    color =[1 1 2 2],linestyle=[:solid :dash],linealpha=[fill(0.2,7);0],label=:none,xlabel="age (male, female)") |>savefigname("../figs/fig1/raw/kamituga.svg", save = false)

plot(otherhz2024_fit.zmb_plot, otherhz2024_validate.zmb_fit|>collect, title="Other health zones",size=(400,300),xrotation=45,bottom_margin=3Plots.PlotMeasures.mm,legend=:none,ylim=(0,0.301),
        color =[1 1 2 2],linestyle=[:solid :dash],linealpha=[fill(0.2,7);0],label=:none, xlabel="age (male, female)") |>savefigname("../figs/fig1/raw/otherhz.svg", save = false)

plot(burundi2024_fit.zmb_plot, burundi2024_validate.zmb_fit|>collect, title="Burundi, 2024",size=(800,300),xrotation=0,bottom_margin=4Plots.PlotMeasures.mm,legend=:none,ylim=(0,0.201),
        color =[1 1 2 2],linestyle=[:solid :dash],linealpha=[fill(0.2,7);0],label=:none,xlabel="age (male, female)",left_margin=3Plots.PlotMeasures.mm)|>savefigname("../figs/fig1/raw/burundi.svg",save=false)

# ## Results summary

# +
est_kivu = CrI(kivu2024_fit.zmb_fit)
est_kamituga = CrI(kamituga2024_fit.zmb_fit)
est_otherhz = CrI(otherhz2024_fit.zmb_fit)
est_burundi = CrI(burundi2024_fit.zmb_fit);

# model weights
Ib_weights = pweights([0.21, 0.33,0.29,0.17])
# -

eig_kamituga=b_eigenanalysis(kamituga2024_fit.zmb_fit)
eig_kivu=b_eigenanalysis(kivu2024_fit.zmb_fit)
eig_otherhz=b_eigenanalysis(otherhz2024_fit.zmb_fit)
eig_burundi=b_eigenanalysis(burundi2024_fit.zmb_fit);

# +
# model averaged susceptibility estimates
sus_0_5=MixtureModel(posteriordist.(zmb2015_24_fit|>collect.|>first,1),Ib_weights)
(median(sus_0_5),quantile(sus_0_5,[0.025,0.975]))|>display

sus_vax=MixtureModel(posteriordist.(zmb2015_24_fit|>collect.|>first,2),Ib_weights)
(median(sus_vax),quantile(sus_vax,[0.025,0.975]))|>display

# -

[quantile.(sus_0_5.components,q) for q in (0.025,0.5,0.975)]|>display
[quantile.(sus_vax.components,q) for q in (0.025,0.5,0.975)]|>display

# fraction of R0 attributable to sexual transmisssion
frac_R0_samples =[ (broadcast.(-,1, broadcast.(/,eig.eigval0, eig.eigval))) for eig in [eig_kivu, eig_kamituga, eig_otherhz, eig_burundi]]
frac_R0 = [MixtureModel([MixtureModel(Normal.(post,0)) for post in loc],Ib_weights) for loc in frac_R0_samples]
quantile.(frac_R0,Ref([0.025,0.5,0.975]))|>display
[quantile.(loc.components,Ref([0.025,0.5,0.975])) for loc in frac_R0]

# fraction of R0 attributable to sexual transmisssion: Synthetic matrix
eig_kamituga=b_eigenanalysis(kamituga2024_fit.drc_fit)
eig_kivu=b_eigenanalysis(kivu2024_fit.drc_fit)
eig_otherhz=b_eigenanalysis(otherhz2024_fit.drc_fit)
eig_burundi=b_eigenanalysis(burundi2024_fit.drc_fit);
frac_R0_samples =[ (broadcast.(-,1, broadcast.(/,eig.eigval0, eig.eigval))) for eig in [eig_kivu, eig_kamituga, eig_otherhz, eig_burundi]]
frac_R0 = [MixtureModel([MixtureModel(Normal.(post,0)) for post in loc],pweights([0,1])) for loc in frac_R0_samples]
quantile.(frac_R0,Ref([0.025,0.5,0.975]))|>display

# +
# transmission frequency heatmap
eig_endemic = (ngm=endemic2015_24_fit.zmb_fit|>collect.|>last.|>ngm, eigcases=endemic2015_24_fit.zmb_fit|>collect.|>last.|>dominanteigvec.|>(x->[x]), eigval0 = endemic2015_24_fit.zmb_fit|>collect.|>first.|>dominanteigval,eigval = endemic2015_24_fit.zmb_fit|>collect.|>last.|>dominanteigval)

heats=[begin mixingmatrices = (broadcast.(.*,eig.ngm,(eig.eigcases')|>vec).|>sum)
    R0scale=0.82./(eig_endemic.eigval0)
    mixingmatrices.*=R0scale
    meanmixmat=exp.(mean(broadcast.(log,mixingmatrices),Ib_weights))
    heatmap(meanmixmat,color=cgrad(:inferno,scale=(x->exp(4x))),clim=(0,0.3)) end for eig in [eig_endemic, eig_kivu, eig_kamituga, eig_otherhz, eig_burundi]]

plot(heats[[1,2,5]]...,ticks=(1:8,makeagegrouplabel([0:5:20;30:10:50])),size=(800,650),xrotation=45, 
    xlabel = "age (infector)", title = ["Endemic provinces" "North and South Kivu" "Burundi"], ylabel = ["age (infectee)" "" "age (infectee)"],
    bottom_margin = 2Plots.PlotMeasures.mm) |>savefigname("../figs/fig2/raw/heatmaps.svg", save = false)
# +
# vaccine impact heatmaps
popvaxrange, fswvaxrange = 0:0.5:100, 0:0.5:100
kivu2024_project=deepcopy(kivu2024_fit.zmb_fit)
kivuRmap=vaccineRmap.(kivu2024_project|>collect,eig_endemic.eigval0,ve=0.86,poprange = popvaxrange./100, fswrange = fswvaxrange./100)
burundi2024_project=deepcopy(burundi2024_fit.zmb_fit)
burundiRmap=vaccineRmap.(burundi2024_project|>collect,eig_endemic.eigval0,ve=0.86,poprange = popvaxrange./100, fswrange = fswvaxrange./100);


# -

hm_kivu = heatmap(fswvaxrange,popvaxrange,mean(kivuRmap.|>first,Ib_weights),clim=(0.5,1.5), xticks = (0:20:100), color=cgrad(:coolwarm, [0.0, 0.5, 0.501, 1],categorical=false))
contour!(hm_kivu, fswvaxrange,popvaxrange,mean(kivuRmap.|>first,Ib_weights),contour_labels=true, levels = 0.5:0.1:1.5,color=:black)
hm_burundi = heatmap(fswvaxrange,popvaxrange,mean(burundiRmap.|>first,Ib_weights),clim=(0.5,1.5), xticks = (0:20:100),color=cgrad(:coolwarm, [0.0, 0.5, 0.501, 1],categorical=false))
contour!(hm_burundi,fswvaxrange,popvaxrange,mean(burundiRmap.|>first,Ib_weights),contour_labels=true, levels = 0.5:0.05:1.5,color=:black)
hm2_kivu = heatmap(fswvaxrange,popvaxrange,mean(100 .*(kivuRmap.|>last),Ib_weights),clim=(0.0,60), xticks = (0:20:100), color=cgrad(:Blues_3,rev=true))
contour!(hm2_kivu, fswvaxrange,popvaxrange,mean(100 .*(kivuRmap.|>last),Ib_weights),contour_labels=true, levels = 0.0:5:100,color=:black)
hm2_burundi = heatmap(fswvaxrange,popvaxrange,mean(100 .*(burundiRmap.|>last),Ib_weights),clim=(0.0,60), xticks = (0:20:100),color=cgrad(:Blues_3,rev=true))
contour!(hm2_burundi,fswvaxrange,popvaxrange,mean(100 .*(burundiRmap.|>last),Ib_weights),contour_labels=true, levels = 0.0:5:100,color=:black);

plot(hm_kivu,title = "reproduction number" ,xlabel = "vaccine uptake % (FSW)", ylabel = "vaccine uptake % (adult aged 20–39)",size=(400,350)) |>savefigname("../figs/fig2/raw/kivu_vax.svg",save=false)
plot(hm_burundi,title = "reproduction number",xlabel = "vaccine uptake % (FSW)", ylabel = "vaccine uptake % (adult aged 20–39)",size=(400,350)) |>savefigname("../figs/fig2/raw/burundi_vax.svg",save=false)
plot(hm2_kivu,title = "% reduction in reproduction number" ,xlabel = "vaccine uptake % (FSW)", ylabel = "vaccine uptake % (adult aged 20–39)",size=(400,350)) |>savefigname("../figs/fig2/raw/kivu_red.svg",save=false)
plot(hm2_burundi,title = "% reduction in reproduction number",xlabel = "vaccine uptake % (FSW)", ylabel = "vaccine uptake % (adult aged 20–39)",size=(400,350)) |>savefigname("../figs/fig2/raw/burundi_red.svg",save=false)

# ## Supple figures

# Eigenvector without accounting for susceptibility
plot(tshuapaplot, zmb2015_original|>collect; label = ["all contacts" "physical only" "at home" "physical & home"],legend=(0.1,0.96),ylim=(0,0.4),color = [1 1 2 2],linestyle = [:solid :dash],title="Tshuapa, 2011-2015 (null-model)")|>savefigname("../figs/Sfigs/raw/null1.svg", save=false)
plot(tshuapaplot, drc2015_original|>collect; label = ["all contacts" "physical only" "at home" "physical & home"],legend=(0.1,0.96),ylim=(0,0.4),color=[:royalblue :firebrick],title="Tshuapa, 2011-2015 (null-model)") |>savefigname("../figs/Sfigs/raw/null2.svg", save=false)

# training vs validation
plot(tshuapa2015_fit.zmb_plot, title="Tshuapa, 2011-2015 (training)",legend=:none)|>savefigname("../figs/Sfigs/raw/train1.svg", save = false)
plot(endemic2024_validate.zmb_plot, title="Endemic provinces, 2024 (validation)",legend=:none) |>savefigname("../figs/Sfigs/raw/train2.svg", save=false)
plot(tshuapa2015_fit.drc_plot, title="Tshuapa, 2011-2015 (training)",legend=:none)|>savefigname("../figs/Sfigs/raw/validation1.svg", save=false)
plot(endemic2024_validate.drc_plot, title="Endemic provinces, 2024 (validation)",legend=:none) |>savefigname("../figs/Sfigs/raw/validation2.svg", save=false)

#DRC synthetic matrix
plot(endemic2015_24_fit.drc_plot[1], title="Tshuapa, 2011-2015") |>savefigname("../figs/Sfigs/raw/syn_tshuapa.svg", save = false)
plot(endemic2015_24_fit.drc_plot[2], title="Endemic provinces, 2024",legend=:none) |>savefigname("../figs/Sfigs/raw/syn_endemicprovs.svg", save = false)

plot(kivu2024_fit.drc_plot, kivu2024_validate.drc_fit|>collect, title="North and South Kivu, 2024",size=(800,300),xrotation=0,bottom_margin=4Plots.PlotMeasures.mm,left_margin=3Plots.PlotMeasures.mm,ylim=(0,0.301),xlabel="age (male, female)",
        color =[:royalblue :firebrick],linealpha=[fill(0.2,7);0],label=:none) |>savefigname("../figs/Sfigs/raw/syn_kivu.svg", save = false)

plot(kamituga2024_fit.drc_plot, kamituga2024_validate.drc_fit|>collect, title="Kamituga health zone",size=(400,300),xrotation=45,bottom_margin=3Plots.PlotMeasures.mm,legend=:none,ylim=(0,0.301),
    color =[:royalblue :firebrick],linealpha=[fill(0.2,7);0],label=:none,xlabel="age (male, female)") |>savefigname("../figs/Sfigs/raw/syn_kamituga.svg", save = false)

plot(otherhz2024_fit.drc_plot, otherhz2024_validate.drc_fit|>collect, title="Other health zones",size=(400,300),xrotation=45,bottom_margin=3Plots.PlotMeasures.mm,legend=:none,ylim=(0,0.301),
        color =[:royalblue :firebrick],linealpha=[fill(0.2,7);0],label=:none,xlabel="age (male, female)") |>savefigname("../figs/Sfigs/raw/syn_otherhz.svg", save =false)

plot(burundi2024_fit.drc_plot, burundi2024_validate.drc_fit|>collect, title="Burundi, 2024",size=(800,300),xrotation=0,bottom_margin=4Plots.PlotMeasures.mm,legend=:none,ylim=(0,0.201),
        color =[:royalblue :firebrick],linealpha=[fill(0.2,7);0],label=:none,xlabel="age (male, female)") |>savefigname("../figs/Sfigs/raw/syn_burundi.svg", save =false)

# NGM convergence check
ngm_tshuapa2015 = ngm.(endemic2015_24_fit.zmb_fit|>collect.|>first).|>first
tshuapa_zoonotic = aggregatecategories(readJSON(:tshuapa2015_zoonotic))
casegens = normalize.((ngm_tshuapa2015).^[0 1 2 100].*tshuapa_zoonotic.cases,1)
plot([plot(casegens[x,:], xticks=(1:6,makeagegrouplabel(tshuapa_zoonotic)),ylim=(0,0.4),xlabel="age",ylabel="proportion",label=["animal-exposed" "H-to-H gen ".*string.((1:2)') "eigenvector"],
    color = [:black :lightblue :steelblue :darkred],linestyle = [:dash :solid :solid :solid]) for x in 1:4]...,size = (600,400),legendfontsize=7, legend=[(0.67,0.95) :none :none :none], title = ["all contacts" "physical contacts" "all contacts at home" "physical contacts at home"]) |>savefigname("../figs/Sfigs/raw/convergence1.svg", save =false)

ngm_kivu2024 = collapseblockmat.(ngm.(kivu2024_fit.zmb_fit|>collect))
initcases=Int.(in.((1:32), Ref([4,20])))
casegens = normalize.(reduce.(vcat,splitsum.(ngm_kivu2024.^[(1:4)' 100].*Ref(initcases),4,Ref([1,3]))),1)
plot([plot(casegens[x,:], xticks=(1:16,repeat(makeagegrouplabel(drc_kivu),2)),xrotation=45,ylim=(0,0.4),linealpha=Float64.((1:16).!=8),xlabel="age (male, female)",ylabel="proportion",label=["gen ".*string.((1:4)') "eigenvector"],tickfontsize = 8,labelfontsize=9, bottom_margin=5Plots.PlotMeasures.mm,
    color = Gray.((range(0,0.9,length=5))')) for x in 1:4]...,size = (600,400),legendfontsize=6, legend=[(0.8,1) :none :none :none], title = ["all contacts" "physical contacts" "all contacts at home" "physical contacts at home"]) |>savefigname("../figs/Sfigs/raw/convergence2.svg", save =false)

# +
# Estimated proportion of sexually active individuals
dictkeys = ["kivu2024_fit", "kamituga2024_fit", "otherhz2024_fit", "burundi2024_fit"]
mixposteriors = [MixtureModel(posteriorjointdist.(cms),Ib_weights) for cms in getindex.(Ref(loadfit),dictkeys).|>first.|>collect]

drc_na = loadfit["kivu2024_fit"][1][1].misc[:pop]
bdi_na = loadfit["burundi2024_fit"][1][1].misc[:pop]
Random.seed!(123)
postsamples = rand.(mixposteriors,1000) # posterior samples from averaged model

activepops = [begin
    m_active = cumsum(na[1:4].*samples[1:4,:],dims=1)
    m_active./=m_active[end:end,:]
    f_active = cumsum(na[5:8].*samples[5:8,:],dims=1)
    f_active./=f_active[end:end,:]
    medians=median([m_active;f_active],dims=2)
    errors = abs.(quantile.(eachrow([m_active;f_active]),[0.025 0.975]).-medians)
    [diff([0.;medians[1:4];0.;medians[5:8]],dims=1)[[1:4;6:end]] [zeros(1,2); errors[1:end-1,:]]]
    end for (samples, na) in zip(postsamples, [fill(drc_na,3);[bdi_na]])]
m_medians=reduce(hcat,getindex.(activepops,Ref(1:4),1))'
f_medians=reduce(hcat,getindex.(activepops,Ref(5:8),1))'
m_error=(reduce(hcat,getindex.(activepops,Ref(1:4),3))',reduce(hcat,getindex.(activepops,Ref(1:4),2))')
f_error=(reduce(hcat,getindex.(activepops,Ref(5:8),3))',reduce(hcat,getindex.(activepops,Ref(5:8),2))');

# +
colors=repeat([collect(palette(:devon,8))[7:-1:4] collect(palette(:lipari,9))[8:-1:5]]|>permutedims,inner=(4,1))
groupedbar([m_medians; f_medians] ,bar_position=:stack,color=colors,legend=:none,label=false,ms = 0)
plot!((1:8).+0.1.*[0,0,1,1,1,1,1,1],1 .-cumsum([m_medians; f_medians],dims=2)[:,1:end-1], markershape=:none,color=:black, seriestype=:scatter, yerror=([m_error[1];f_error[1]][:,2:end].*[0 0 1],[m_error[2];f_error[2]][:,2:end].*[0 0 1]),ms=0,label=false)
plot!((1:8).-0.1.*[0,0,1,0,0,0,1,0],1 .-cumsum([m_medians; f_medians],dims=2)[:,1:end-1], markershape=:none,color=:black, seriestype=:scatter, yerror=([m_error[1];f_error[1]][:,2:end].*[1 0 0],[m_error[2];f_error[2]][:,2:end].*[1 0 0]),ms=0,label=false)
plot!((1:8),1 .-cumsum([m_medians; f_medians],dims=2)[:,1:end-1], markershape=:none,color=:black, seriestype=:scatter, yerror=([m_error[1];f_error[1]][:,2:end].*[0 1 0],[m_error[2];f_error[2]][:,2:end].*[0 1 0]),ms=0,label=false)

labels=string.([15,20,30,40]).*"–".*string.([19,29,39,49]).*[" M" " F"]|>vec|>permutedims
scatter!(fill(-1,8,8),markershape=:square,ylim=(0,1),color=colors[[1,5],:]|>permutedims|>vec|>permutedims,legend=:outertopright,label=labels,size=(500,300),
ylabel="proportion",xticks=(1:8,repeat(["N&S Kivu","Kamituga","Other HZ", "Burundi"],2)),xrotation=45,bottom_margin=4Plots.mm,xtickfontsize=9)|>savefigname("../figs/Sfigs/raw/activepop.svg",save=true)
# -
[quantile.(eachrow(postsamples[x][9:10,:]),Ref([0.025,0.5,0.975])) for x in 1:4]

[quantile(sqrt.(prod(postsamples[x][9:10,:],dims=1)|>vec),[0.025,0.5,0.975]) for x in 1:4]

[quantile(sqrt.(diff(log.(postsamples[x][9:10,:]),dims=1).|>exp|>vec),[0.025,0.5,0.975]) for x in 1:4]


