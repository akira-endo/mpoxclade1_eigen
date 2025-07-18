# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %%
using Pkg
#Pkg.activate("../")

# %%
# Set the plotting backend and styles
using Plots
gr(fontfamily="Helvetica", foreground_color_legend=nothing, background_color_legend=nothing,
   titlefontsize=11, tickfontsize=10, legendfontsize=8, labelfontsize=10,
   grid=true, tick_direction=:out, size=(400, 300))


# %%
# Load the data
using CSVFiles, DataFrames
merged_data = DataFrame(load("../data/intermediate/merged_data.csv")) #|> first

# %%
colourmap = cgrad(:inferno)
# margin
left_margin = 10 * Plots.PlotMeasures.mm
right_margin = 10 * Plots.PlotMeasures.mm
top_margin = 10 * Plots.PlotMeasures.mm
bottom_margin = 10 * Plots.PlotMeasures.mm
# create and filter df for the plot (cladeII R0-criticalR0 cladeI)
plot_data1 = DataFrame(
    cR0 = merged_data.clade1_cR0_14,
    R0clade2 = merged_data.clade2_R0_14,
    closest_SAR = merged_data.closest_SAR_14,
    cR0_10 = merged_data.clade1_cR0_14 .* (1/(1-0.86*0.1)),
    cR0_30 = merged_data.clade1_cR0_14 .* (1/(1-0.86*0.3)),
    cR0_50 = merged_data.clade1_cR0_14 .* (1/(1-0.86*0.5)),
    cR0_70 = merged_data.clade1_cR0_14 .* (1/(1-0.86*0.7))
)
plot_data1 = filter(row -> row.R0clade2 <= 3, plot_data1)

# Plot (cladeII R0-criticalR0 cladeI)
plot_1 = plot(plot_data1.R0clade2, plot_data1.cR0, linecolor = colourmap[1], linewidth = 1.5, alpha = 0.8,
    xlabel = "clade II R₀ in the previous epidemic", ylabel = "clade I R₀ for epidemic takeoff",
    legend = (0.75,0.5), label = "0%", legendtitle = "\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0vaccine uptake", legendtitlefonthalign=:right, legendtitlefontsize=8,
    left_margin = left_margin, right_margin = right_margin, top_margin = top_margin, bottom_margin = bottom_margin, xlim = (1,3),
    yticks = ([0, 1, 2, 4, 6, 8], string.([0, 1, 2, 4, 6, 8])), ylim = (0, maximum(plot_data1.cR0))
)
plot!(plot_data1.R0clade2, plot_data1.cR0_10, linecolor = colourmap[50], linewidth = 1.5, alpha = 0.8, label = "10%")
plot!(plot_data1.R0clade2, plot_data1.cR0_30, linecolor = colourmap[100], linewidth = 1.5, alpha = 0.8, label = "30%")
plot!(plot_data1.R0clade2, plot_data1.cR0_50, linecolor = colourmap[150], linewidth = 1.5, alpha = 0.8, label = "50%")
plot!(plot_data1.R0clade2, plot_data1.cR0_70, linecolor = colourmap[200], linewidth = 1.5, alpha = 0.8, label = "70%")

# %%
colourmap = cgrad(:inferno)
# Creat and filter df
plot_data2 = DataFrame(
    cR0 = merged_data.clade1_cR0_14,
    effectiveS = 1 ./ merged_data.clade1_cR0_14,
    effectiveS_10 = (1-0.86*0.1) ./ (merged_data.clade1_cR0_14),
    effectiveS_30 = (1-0.86*0.3) ./ (merged_data.clade1_cR0_14),
    effectiveS_50 = (1-0.86*0.5) ./ (merged_data.clade1_cR0_14),
    effectiveS_70 = (1-0.86*0.7) ./ (merged_data.clade1_cR0_14),
    R0clade2 = merged_data.clade2_R0_14,
    closest_SAR = merged_data.closest_SAR_14
)
plot_data2 = filter(row -> row.R0clade2 <= 3, plot_data2)

# Plot (cladeIIR0-effective susceptible prop)
plot_2 = plot(plot_data2.R0clade2, plot_data2.effectiveS .* 100, linecolor = colourmap[1], linewidth = 1.5, alpha = 0.8,
    xlabel = "clade II R₀ in the previous epidemic", ylabel = "effective susceptible proportion %", legend = (0.75,0.9),label = "0%", legendtitle = "\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0vaccine uptake", legendtitlefonthalign=:right,legendtitlefontsize=8,
    left_margin=left_margin, right_margin=right_margin, top_margin=top_margin, bottom_margin=bottom_margin, xlim=(1, 3),
     yticks=0:20:100, ylim=(0, 100)
)
plot!(plot_data2.R0clade2, plot_data2.effectiveS_10 .* 100, linecolor = colourmap[50], linewidth = 1.5, alpha = 0.8, label = "10%")
plot!(plot_data2.R0clade2, plot_data2.effectiveS_30 .* 100, linecolor = colourmap[100], linewidth = 1.5, alpha = 0.8, label = "30%")
plot!(plot_data2.R0clade2, plot_data2.effectiveS_50 .* 100, linecolor = colourmap[150], linewidth = 1.5, alpha = 0.8, label = "50%")
plot!(plot_data2.R0clade2, plot_data2.effectiveS_70 .* 100, linecolor = colourmap[200], linewidth = 1.5, alpha = 0.8, label = "70%")

# %%
# Combine the two plots
combined_plot = plot(plot_2, plot_1, layout = (1, 2), size = (800, 300)) |> display
#savefig(combined_plot, "fig_s7.svg")
