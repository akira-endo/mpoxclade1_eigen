# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

using Pkg
Pkg.activate("../")

# +
#Pkg.instantiate()
# -

# load functions
include("../src/eigen.jl");
include("../src/eigen_setup.jl");
include("../src/eigen_output.jl")

gr(fontfamily="Helvetica",foreground_color_legend = nothing,background_color_legend = nothing, titlefontsize=11, tickfontsize=10, legendfontsize=8,labelfontsize=10,grid=true, tick_direction=:out,size=(400,300))

# +
endemicplot = plot(plot(drc_endemic_ag),ylim=(0,0.35),xtickfontsize=9)
tshuapaplot = plot(tshuapa_h2hag,color=:black)
endemic2015_24_fit = output_fit(
    [tshuapa_h2hag,drc_endemic_ag];
    zmb_skeleton = [zmb2015,zmb2024],
    drc_skeleton = [drc2015,drc2024],
    dataplots = [tshuapaplot,endemicplot], bayesian=true
    );

zmb2015_24_fit = endemic2015_24_fit.zmb_fit
drc2015_24_fit = endemic2015_24_fit.drc_fit;
# -

# functions to create contact matrices for a given year
filters = [nothing, (phys_contact=1,),(cnt_home=1,),(phys_contact=1,cnt_home=1)]
drc_matrices(year,filters=filters) = @suppress (;zip([:all, :phys, :home, :physhome] ,getfield.(contactmatrix.(:zimbabwe_survey, Ref(drc_endemic_ag), "COD", filters,year=year,refyear=2012,refcountrycode="MAN"),:matrix).|>first)...)
bdi_matrices(year,filters=filters) = @suppress (;zip([:all, :phys, :home, :physhome] ,getfield.(contactmatrix.(:zimbabwe_survey, Ref(drc_endemic_ag), "BDI", filters,year=year,refyear=2012,refcountrycode="MAN"),:matrix).|>first)...)
# note: For 2024 census data use "BDIC" instead of "BDI". year has to be 2024, otherwise it will throw an error

# +
# update community contact matrix
kivuCMs = loadfit["kivu2024_fit"][1] # vector of ContactMatrix [all, phys, home, physhome]
for (cm, mat) in zip(kivuCMs, drc_matrices(2025))
    # replace community contact elements of cm with mat.
    # replacing [1,2] only will work ([1,4] is a shallow copy)
    cm.matrix[1,2].=mat./2 # /2 because sex-specific
end

burundiCMs = loadfit["burundi2024_fit"][1] # vector of ContactMatrix [all, phys, home, physhome]
for (cm, mat) in zip(burundiCMs, bdi_matrices(2025))
    # replace community contact elements of cm with mat.
    # replacing [1,2] only will work ([1,4] is a shallow copy)
    cm.matrix[1,2].=mat./2 # /2 because sex-specific
end
# -

kivu_eigvals = MCMCiterate.(Ref(dominanteigval), kivuCMs|>collect)
burundi_eigvals = MCMCiterate.(Ref(dominanteigval), burundiCMs|>collect)

kivuCMs_nosexual = deepcopy(zmb2015_24_fit)
for (cm, mat) in zip(kivuCMs_nosexual, drc_matrices(2025))
    # replace community contact elements of cm with mat.
    cm[2].matrix[1].=mat
end
burundiCMs_nosexual = deepcopy(zmb2015_24_fit)
for (cm, mat) in zip(burundiCMs_nosexual, bdi_matrices(2025))
    # replace community contact elements of cm with mat.
    cm[2].matrix[1].=mat
end
function dominanteigval_nosexual(cm)
    cm.addmat.=zero(cm.addmat)
    dominanteigval(cm)
end

# ### Analysis

loadfit = load("../outputs/sexualfit_main.jld2")
kivu_chain_sex = collect(values(loadfit["kivu2024_fit"][1]))
burundi_chain_sex = collect(values(loadfit["burundi2024_fit"][1]))
zmb_fit = zmb2015_24_fit
kivu_chain_nosex = [zmb_fit.all[2], zmb_fit.phys[2], zmb_fit.home[2], zmb_fit.physhome[2]]
burundi_chain_nosex = deepcopy(kivu_chain_nosex)

# +
include("../src/Reff_projection_utils.jl");

# compute eigenvalue
nonsex_kivu = compute_reff(kivu_chain_nosex, drc_matrices; sexual=false)
total_kivu = compute_reff(kivu_chain_sex, drc_matrices; sexual=true)
nonsex_buri = compute_reff(burundi_chain_nosex, bdi_matrices; sexual=false)
total_buri = compute_reff(burundi_chain_sex, bdi_matrices; sexual=true)

# beta
med2015 = quantile(nonsex_kivu[2015], 0.5)
β = 0.82 / med2015
@printf("\nβ  = %.4f (2015 non-sexual Kivu Reff = 0.82)\n", β)

# compute reff 
scaled_nonsex_kivu = Dict(yr => nonsex_kivu[yr] .* β for yr in keys(nonsex_kivu))
scaled_total_kivu  = Dict(yr => total_kivu[yr]  .* β for yr in keys(total_kivu))
scaled_nonsex_buri = Dict(yr => nonsex_buri[yr] .* β for yr in keys(nonsex_buri))
scaled_total_buri  = Dict(yr => total_buri[yr]  .* β for yr in keys(total_buri))

# output
summarize(scaled_nonsex_kivu; name="Kivu non-sexual Reff (Ia)")
summarize(scaled_total_kivu;  name="Kivu total Reff (Ib)")
summarize(scaled_nonsex_buri; name="Burundi non-sexual Reff (Ia)")
summarize(scaled_total_buri;  name="Burundi total Reff (Ib)")

# +
# plot
years_vec = [2010, 2015, 2020, 2024, 2030]

R0_drc = [quantile(scaled_nonsex_kivu[y], 0.5)   for y in years_vec]
R0_drc_low = [quantile(scaled_nonsex_kivu[y], 0.025) for y in years_vec]
R0_drc_up = [quantile(scaled_nonsex_kivu[y], 0.975) for y in years_vec]

R0_bdi = [quantile(scaled_nonsex_buri[y], 0.5)   for y in years_vec]
R0_bdi_low = [quantile(scaled_nonsex_buri[y], 0.025) for y in years_vec]
R0_bdi_up = [quantile(scaled_nonsex_buri[y], 0.975) for y in years_vec]

R0_drc_Ib = [quantile(scaled_total_kivu[y], 0.5)    for y in years_vec]
R0_drc_Ib_low = [quantile(scaled_total_kivu[y], 0.025)  for y in years_vec]
R0_drc_Ib_up = [quantile(scaled_total_kivu[y], 0.975)  for y in years_vec]

R0_bdi_Ib = [quantile(scaled_total_buri[y], 0.5)    for y in years_vec]
R0_bdi_Ib_low = [quantile(scaled_total_buri[y], 0.025)  for y in years_vec]
R0_bdi_Ib_up = [quantile(scaled_total_buri[y], 0.975)  for y in years_vec]

drc_halfwidth = (R0_drc .- R0_drc_low, R0_drc_up .- R0_drc)
bdi_halfwidth = (R0_bdi .- R0_bdi_low, R0_bdi_up .- R0_bdi)
drc_Ib_halfwidth = (R0_drc_Ib .- R0_drc_Ib_low, R0_drc_Ib_up .- R0_drc_Ib)
bdi_Ib_halfwidth = (R0_bdi_Ib .- R0_bdi_Ib_low, R0_bdi_Ib_up .- R0_bdi_Ib)

default_color = palette(:auto)[1]
colours = [:royalblue, 1, 2, :firebrick,
           RGBA(default_color*0.8,1.0),
           RGBA(65/255*0.8,105/255*0.8,225/255*0.8,1.0)]
bottom_margin = 10Plots.PlotMeasures.mm

sel = findall(x->x in (2024,2030), years_vec)

plt = plot([2015,2015], [0.79,0.85],
           color=:black, xlabel="year", ylabel="reproduction number",
           label="", ylim=(0,1.7), size=(600,300),
           bottom_margin=bottom_margin)

# clade Ia (DRC)
plot!(plt, years_vec, R0_drc,
      ribbon=drc_halfwidth, fillalpha=0.15,
      color=colours[2], label="clade Ia (DRC)")

# clade Ia (Burundi)
plot!(plt, years_vec, R0_bdi,
      ribbon=bdi_halfwidth, fillalpha=0.15,
      color=colours[3], label="clade Ia (Burundi)")

# clade Ib (DRC)
plot!(plt, years_vec[sel], R0_drc_Ib[sel],
      ribbon=(drc_Ib_halfwidth[1][sel], drc_Ib_halfwidth[2][sel]),
      fillalpha=0.15, color=colours[1], label="clade Ib (DRC)")

# clade Ib (Burundi)
plot!(plt, years_vec[sel], R0_bdi_Ib[sel],
      ribbon=(bdi_Ib_halfwidth[1][sel], bdi_Ib_halfwidth[2][sel]),
      fillalpha=0.15, color=colours[4], label="clade Ib (Burundi)")

hline!(plt, [1], color=RGBA(0.5,0.5,0.5,0.5),
       linestyle=:dash, label="")

scatter!(plt, [2015],[0.82], color=:black,
         marker=:circle, markersize=4, markerstrokewidth=0,
         label=nothing)
plot!(plt, [2013,2017],[0.82,0.82], color=:black, label=nothing)
annotate!(plt, 2015,0.7, text("Tshuapa","Helvetica",:black,7))

plot!(plt, [2023.1,2023.1], [0.76,0.92], color=colours[5], label=nothing)
scatter!(plt, [2023.1],[0.84], marker=:circle, color=colours[5],
         markersize=4, markerstrokewidth=0, label=nothing)
plot!(plt, [2023.3,2023.3], [0.77,0.93], color=colours[5], label=nothing)
scatter!(plt, [2023.3],[0.85],
         marker=:circle, markercolor=:white,
         markerstrokecolor=colours[5],
         markersize=4, markerstrokewidth=1, label=nothing)
annotate!(2022.1, 0.84, text("pre-Ib 1", "Helvetica", :black, 7), fontfamily="Helvetica")
annotate!(2024.3, 0.85, text("pre-Ib 2", "Helvetica", :black, 7), fontfamily="Helvetica")

plot!(plt, [2024.1,2024.1], [1.46,1.55], color=colours[6], label=nothing)
scatter!(plt, [2024.1],[1.50], marker=:circle, color=colours[6],
         markersize=4, markerstrokewidth=0, label=nothing)
plot!(plt, [2024.3,2024.3], [1.33,1.41], color=colours[6], label=nothing)
scatter!(plt, [2024.3],[1.37],
         marker=:circle, markercolor=:white,
         markerstrokecolor=colours[6],
         markersize=4, markerstrokewidth=1, label=nothing)
annotate!(plt, 2023.3,1.50, text("SK 1","Helvetica",:black,7))
annotate!(plt, 2023.5,1.37, text("SK 2","Helvetica",:black,7))

#savefig(plt, "fig_2b.svg")
display(plt)
