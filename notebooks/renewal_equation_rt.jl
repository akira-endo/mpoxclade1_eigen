include("../src/renewal_equation.jl")
default_plot_setting()

df_DRC, df_SK = read_suspected_DRC_SK()
# For week display
df_conf = read_WHO_confirmed()
df_conf2023 = @subset df_conf (:week_end_date .< Date("2024-01-01") ) .& (:week_end_date .>= Date("2023-01-07")) 
df_conf2023Sep = @subset df_conf2023 (:week_end_date .< Date("2023-09-01"))
nothing

RtHols = estimate_multiple_RenewalRt(df_conf2023Sep, df_SK)
R_means, R_confs = merge_means_confs(RtHols)
plot_epicurve_SK_DRC(df_DRC, df_SK, df_conf, df_conf2023, df_conf2023Sep, R_means, R_confs)


