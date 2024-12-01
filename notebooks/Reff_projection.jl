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
#     display_name: Julia 1.8.3
#     language: julia
#     name: julia-1.8
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

endemic2015_24_fit = output_fit(
    [tshuapa_h2hag,drc_endemic_ag];
    zmb_skeleton = [zmb2015,zmb2024],
    drc_skeleton = [drc2015,drc2024],
    dataplots = endemicplot
    );
zmb2015_24_fit = endemic2015_24_fit.zmb_fit
drc2015_24_fit = endemic2015_24_fit.drc_fit;
# -

keys(zmb2015_24_fit)

dominanteigval(zmb2015_24_fit.all[2]) # [2] is ContactMatrix for 2024 data 

zmb2015_24_fit.all[2].susceptibility

zmb2015_24_fit.all[2].parameters

# add new parameters for the group just under partially vaccinated
for cm_vec in zmb2015_24_fit
    for cm in cm_vec
        cm.parameters[:s_postvax]=fill(0.5)
        cm.susceptibility[1][end-2]=cm.parameters[:s_postvax]
    end
end 

zmb2015_24_fit.all[2].parameters

zmb2015_24_fit.all[2].susceptibility

zmb2015_24_fit.all[2].parameters[:s_postvax].=1

zmb2015_24_fit.all[2].susceptibility

zmb2015_24_fit.all[2].parameters

zmb2015_24_fit.all[2].susceptibility[1]

# ### DRC

include("../src/Reff_projection_utils.jl");

# +
years = [2010, 2015, 2020, 2024, 2030]
s_infant_values = [2.2, 1.9, 1.5, 1.3]
s_vax_values = [0.26, 0.28, 0.22, 0.23]

all_susceptibilities = Dict{Tuple{Float64, Float64}, Dict}()
# get the correct combination of susceptibility by age group
for (s_infant, s_vax) in zip(s_infant_values, s_vax_values)
    s_partvax = (1 + s_vax) / 2  
    susceptibilities = create_susceptibilities_for_years(years, s_infant, s_partvax, s_vax)
    all_susceptibilities[(s_infant, s_vax)] = susceptibilities
end

for (key, susceptibilities) in all_susceptibilities
    println("For s_infant = $(key[1]), s_vax = $(key[2]):")
    println(susceptibilities)
end
# -

dominant_eigvals_all = []
dominant_eigvals_phys = []
dominant_eigvals_home = []
dominant_eigvals_physhome = []
include("../src/Reff_projection_setup.jl")

# all contact
println("Dominant Eigenvalues: ", dominant_eigvals_all)
println("Used combination: s_infant = $(first_key[1]), s_vax = $(first_key[2])")

# physical contact
println("Dominant Eigenvalues: ", dominant_eigvals_phys)
println("Used combination: s_infant = $(second_key[1]), s_vax = $(second_key[2])")

# home contact
println("Dominant Eigenvalues: ", dominant_eigvals_home)
println("Used combination: s_infant = $(third_key[1]), s_vax = $(third_key[2])")

# physical home contact
println("Dominant Eigenvalues: ", dominant_eigvals_physhome)
println("Used combination: s_infant = $(forth_key[1]), s_vax = $(forth_key[2])")

# +
# R0 projection with community cotacts only
# model weights
w1 = 0.22
w2 = 0.38
w3 = 0.21
w4 = 0.17

# model averaging
weighted_avg = []
for i in 1:length(dominant_eigvals_all)
    weighted_value = w1 * dominant_eigvals_all[i] + w2 * dominant_eigvals_phys[i] + w3 * dominant_eigvals_home[i] + w4 * dominant_eigvals_physhome[i]
    push!(weighted_avg, weighted_value)
end
println("Weighted averages: ", weighted_avg)
# -

β=0.82/weighted_avg[2]
R₀ = [x * β for x in weighted_avg]
println("R₀: ", R₀)
println("β: ", β)

# Nouth and South Kivu (community+sexual contact)
s1 = 1 - 0.39
s2 = 1 - 0.4
s3 = 1 - 0.42
s4 = 1 - 0.42
weighted_avg_sexual_nskivu = model_averaging_sexual(s1, s2, s3, s4, dominant_eigvals_all, dominant_eigvals_phys, dominant_eigvals_home, dominant_eigvals_physhome)
println("Weighted averages (sexual): ", weighted_avg_sexual_nskivu)

sR₀_kivu = [x * β for x in weighted_avg_sexual_nskivu]
println("sR₀: ", sR₀_kivu)
println("sβ: ", β)

nskivu_plot=plot(years, R₀, xlabel="Year", ylabel="R₀", label="community contacts", ylim=(0,:auto), title="Nouth and South Kivu")
plot!([2024, 2030], [sR₀_kivu[4], sR₀_kivu[5]], label="community & sexual contacts")
plot!([2015, 2015], [0.79, 0.85], color=:black, label=nothing)
scatter!([2015], [0.82], color=:black, marker=:circle, markersize=3, markerstrokewidth=0, label=nothing)
hline!([1], color=:black, linestyle=:dashdot, label="")

# ### Burundi

# +
years = [2010, 2015, 2020, 2024, 2030]
s_infant_values = [2.2, 1.9, 1.5, 1.3]
s_vax_values = 1　.- 0.7 .* (1 .- [0.26, 0.28, 0.22, 0.23])

all_susceptibilities = Dict{Tuple{Float64, Float64}, Dict}()
# get the correct combination of susceptibility by age group
for (s_infant, s_vax) in zip(s_infant_values, s_vax_values)
    s_partvax = (1 + s_vax)/2   
    susceptibilities = create_susceptibilities_for_years(years, s_infant, s_partvax, s_vax)
    all_susceptibilities[(s_infant, s_vax)] = susceptibilities
end

for (key, susceptibilities) in all_susceptibilities
    println("For s_infant = $(key[1]), s_vax = $(key[2]):")
    println(susceptibilities)
end
# -

dominant_eigvals_all = []
dominant_eigvals_phys = []
dominant_eigvals_home = []
dominant_eigvals_physhome = []
include("../src/Reff_projection_setup.jl")

println("Dominant Eigenvalues: ", dominant_eigvals_all)
println("Used combination: s_infant = $(first_key[1]), s_vax = $(first_key[2])")

println("Dominant Eigenvalues: ", dominant_eigvals_phys)
println("Used combination: s_infant = $(second_key[1]), s_vax = $(second_key[2])")

println("Dominant Eigenvalues: ", dominant_eigvals_home)
println("Used combination: s_infant = $(third_key[1]), s_vax = $(third_key[2])")

println("Dominant Eigenvalues: ", dominant_eigvals_physhome)
println("Used combination: s_infant = $(forth_key[1]), s_vax = $(forth_key[2])")

# +
# R0 projection with community cotacts only
# model weights
w1 = 0.22
w2 = 0.38
w3 = 0.21
w4 = 0.17

# model averaging
weighted_avg = []
for i in 1:length(dominant_eigvals_all)
    weighted_value = w1 * dominant_eigvals_all[i] + w2 * dominant_eigvals_phys[i] + w3 * dominant_eigvals_home[i] + w4 * dominant_eigvals_physhome[i]
    push!(weighted_avg, weighted_value)
end
println("Weighted averages: ", weighted_avg)
# -

R₀_burundi = [x * β for x in weighted_avg]
println("R₀ (Burundi): ", R₀_burundi)
println("β : ", β)

# +
# Burundi (community+sexual contact)
s1 = 1 - 0.05
s2 = 1 - 0.06
s3 = 1 - 0.05
s4 = 1 - 0.04

weighted_avg_sexual_burundi = model_averaging_sexual(s1, s2, s3, s4, dominant_eigvals_all, dominant_eigvals_phys, dominant_eigvals_home, dominant_eigvals_physhome)

println("Weighted averages (sexual_burundi): ", weighted_avg_sexual_burundi)
# -

sR₀_burundi = [x * β for x in weighted_avg_sexual_burundi]
println("sR₀: ", sR₀_burundi)

# ### Plot

# +
using Colors
default_color = palette(:auto)[1]
colours = [:royalblue, 1, 2, :firebrick, RGBA(red(default_color) * 0.8, green(default_color) * 0.8, blue(default_color) * 0.8, 1.0), RGBA(65/255 * 0.8, 105/255 * 0.8, 225/255 * 0.8, 1.0)]  
bottom_margin = 10 * Plots.PlotMeasures.mm

plot([2015, 2015], [0.79, 0.85], color=:black, xlabel="year", ylabel="reproduction number", label="", ylim=(0,1.7), title="", size=(600,300), bottom_margin=bottom_margin)
scatter!([2015], [0.82], color=:black, marker=:black, markersize=4, markerstrokewidth=0, label=nothing)
plot!([2013, 2017], [0.82, 0.82], color=:black, label=nothing)
annotate!(2015, 0.7, text("Tshuapa", "Helvetica", :black, 7))

plot!(years, R₀, color=colours[2], label="clade Ia (DRC)")
plot!([2024, 2030], [sR₀_kivu[4], sR₀_kivu[5]], label="clade Ib (DRC)", color=colours[1])
plot!(years, R₀_burundi, label="clade Ia (Burundi)", color=colours[3])
plot!([2024, 2030], [sR₀_burundi[4], sR₀_burundi[5]], label="clade Ib (Burundi)", color=colours[4])

plot!([2023.1, 2023.1], [0.78, 0.93], color=colours[5], label=nothing)
scatter!([2023.1], [0.85], marker=:circle, color=colours[5], markersize=4, markerstrokewidth=0, label=nothing)
plot!([2023.3, 2023.3], [0.76, 0.91], color=colours[5], label=nothing)
scatter!([2023.3], [0.83], marker=:circle, markercolor=:white, markerstrokecolor= colours[5], markersize=4, markerstrokewidth=1, label=nothing)
annotate!(2022.1, 0.83, text("pre-Ib 1", "Helvetica", :black, 7), fontfamily="Helvetica")
annotate!(2024.3, 0.81, text("pre-Ib 2", "Helvetica", :black, 7), fontfamily="Helvetica")

plot!([2024.1, 2024.1], [1.46, 1.55], color=colours[6], label=nothing)
scatter!([2024.1], [1.50], color=colours[6], marker=:circle, markersize=4, markerstrokewidth=0, label=nothing)
plot!([2024.3, 2024.3], [1.33, 1.41], color=colours[6], label=nothing)
scatter!([2024.3], [1.37], marker=:circle, markersize=4, markercolor=:white, markerstrokewidth=1, markerstrokecolor=colours[6], label=nothing)
annotate!(2023.3, 1.50, text("SK 1", "Helvetica", :black, 7))
annotate!(2023.5, 1.37, text("SK 2", "Helvetica",:black, 7))

hline!([1], color=RGBA(0.5, 0.5, 0.5, 0.5), linestyle=:dash, label="")
