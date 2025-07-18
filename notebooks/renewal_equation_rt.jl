using Pkg
Pkg.activate("../")

include("../src/renewal_equation.jl")
default_plot_setting()

df_DRC, df_SK = read_suspected_DRC_SK()
# For week display
df_conf = read_WHO_confirmed()
df_conf2023 = @subset df_conf (:week_end_date .< Date("2024-01-01") ) .& (:week_end_date .>= Date("2023-01-07"))
df_conf2023Aug = @subset df_conf2023 (:week_end_date .< Date("2023-09-01"))
nothing

# +
# Number of confirmed cases
@show df_conf2023Aug.new_confirmed_cases |> sum
@show Ispill = sum(df_conf2023Aug.new_confirmed_cases) * 0.177 / (52 * 2/3)
gens = [Valentina_wb_ser1, Valentina_wb_ser2]

res_DRC = estimate_DRC_RenewalRt(df_conf2023Aug, gens, Ispill)
res_SK = estimate_SK_RenewalRt(df_SK)
res_Rt = vcat(res_DRC, res_SK)
display(res_Rt)
nothing

# +
include("../src/renewal_equation.jl")
R_means = res_Rt[:, :Rt_m]
R_confs = [[res_Rt[i, :Rt_l], res_Rt[i, :Rt_h]] for i in 1:nrow(res_Rt)]

Random.seed!(12)
Ct = df_conf2023Aug.new_confirmed_cases
df_conf2023Aug.spillover = truncated.(Poisson(Ispill), 0, Ct) .|> rand

plot_epicurve_SK_DRC(df_DRC, df_SK, df_conf, df_conf2023, df_conf2023Aug, R_means, R_confs)
# -











