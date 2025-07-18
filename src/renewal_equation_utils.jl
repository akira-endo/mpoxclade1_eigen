using CSVFiles
using Dates
using DataFrames
using DataFramesMeta
using Distributions
using MCMCChains
using Pipe
using Plots
using ProgressMeter
using Random
using StatsBase
using StatsPlots
using Turing

PATH_SUSPECTED = "../data/suspected_south_kivu_DRC.csv"
PATH_WHO_CONFIRMED = "../data/weekly_AFR_cases_by_country_as_of_15_September_2024.csv"

function default_plot_setting()
	gr(fontfamily = "Helvetica",
		foreground_color_legend = nothing,
		background_color_legend = nothing,
		titlefontsize = 11, tickfontsize = 10,
		legendfontsize = 8, labelfontsize = 10,
		grid = true, tick_direction = :out,
		size = (600, 450))
end

function end_of_week_given_week(year::Int, week_number::Int)
	# Calculate the first day of the year
	start_date = Date(year, 1, 1)
	# Calculate the first Sunday of the year
	days_to_first_sunday = 7 - dayofweek(start_date) + 1
	first_sunday = start_date + Day(days_to_first_sunday - 1)
	# Calculate the Sunday of the given week number
	# Subtract 1 because the first week already includes the first Sunday
	end_of_week = first_sunday + Week(week_number - 1)
	return end_of_week
end

function read_suspected_DRC_SK()
	df = DataFrame(load(PATH_SUSPECTED))
	df.date = end_of_week_given_week.(df.year, df.week)
	df.date = @. ifelse(df.year == 2024, df.date - Dates.Day(7), df.date)
	df.xvalue .= 1:nrow(df)

	df_SK = @subset df :year .== 2024
	df_SK.new_suspect_SK = parse.(Int64, df_SK.new_suspect_SK)
	df_SK.xvalue .= 1:nrow(df_SK)
	return (df, df_SK)
end

function read_WHO_confirmed()::DataFrame
	df_all = @pipe DataFrame(load(PATH_WHO_CONFIRMED)) |>
				   sort(_, :week_end_date)
	df_DRC = @pipe filter(x -> x.country == "Democratic Republic of the Congo", df_all) |>
				   sort(_, :week_end_date)

	latest_date = Date("2024-08-11")
	df_conf = @chain df_DRC begin
		@subset :iso3 .== "COD"
		@orderby :week_end_date
		@subset (:week_end_date .<= latest_date) .& (:week_end_date .>= Date("2023-01-01"))
	end
	df_conf.xvalue .= 1:nrow(df_conf)

	#N_SPILLOVER = 2.53 # per week
	#df_conf.spillover .= N_SPILLOVER
	N_SPILLOVER = 6 # per week
	@transform!(df_conf,
		:spillover = ifelse.(:new_confirmed_cases .< N_SPILLOVER,
			:new_confirmed_cases, N_SPILLOVER))
	return df_conf
end

abstract type UnivariateOffset end

Distributions.mean(d::UnivariateOffset) = mean(d.d) + d.offset
Base.length(x::UnivariateOffset) = 1
Base.iterate(x::UnivariateOffset) = (x, nothing)
Base.iterate(x::UnivariateOffset, nothing) = nothing

function Distributions.pdf(d::UnivariateOffset, x::Real)
	if x < d.offset
		return 0.0
	else
		return pdf(d.d, x - d.offset)
	end
end

function Distributions.cdf(d::UnivariateOffset, x::Real)
	if x < d.offset
		return 0.0
	else
		return cdf(d.d, x - d.offset)
	end
end

struct GammaOffset <: UnivariateOffset
	shape::Float64
	scale::Float64
	offset::Float64
	d::UnivariateDistribution
end
GammaOffset(shape, scale, offset) = GammaOffset(shape, scale, offset, Gamma(shape, scale))

struct WeibullOffset <: UnivariateOffset
	shape::Float64
	scale::Float64
	offset::Float64
	d::UnivariateDistribution
end
WeibullOffset(shape, scale, offset) = WeibullOffset(shape, scale, offset, Weibull(shape, scale))
