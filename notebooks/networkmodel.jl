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
#     display_name: Julia 14 Threads 1.9.3
#     language: julia
#     name: julia-14-threads-1.9
# ---

]activate ../

include("../src/networks.jl")
include("../src/utils.jl")

gr(fontfamily="Helvetica",foreground_color_legend = nothing,background_color_legend = nothing, titlefontsize=11, tickfontsize=10, legendfontsize=8,labelfontsize=10,grid=true, tick_direction=:out,size=(400,300))

# illustration
Random.seed!(123)
gsmall1 = union(binet(500, [50,10], 100,0,3),binet(500, [50,50], 0,0.5,3))
@time gp = graphplot(gsmall,curves=false,markercolor = repeat([1,4],inner=nv(gsmall)÷2),markersize=0.1,size=(450,400)) 
gp|>savefigname("../figs/Sfigs/raw/network1.svg", save=false)

Random.seed!(123)
gsmall = union(bipropnet(500, [50,10], 50,0,3),bipropnet(500, [50,50], 0,0.5,3))
@time gp = graphplot(gsmall,curves=false,markercolor = repeat([1,4],inner=nv(gsmall)÷2),markersize=0.1,size=(450,400))
gp|>savefigname("../figs/Sfigs/raw/network2.svg", save=false)

# +
Random.seed!(123)
gh = binet.(fill(5000,100), Ref([500,100]),1000, 0., 3)
gl = binet.(fill(5000,100), Ref([500,500]), 0,0.5,3)
gc = mixednet.(fill(10000,100),0,0,7.5)
binets1=(gh,gl,gc,union.(gh,gl))

gh = binet.(fill(5000,100), Ref([500,100]),500, 0., 3)
gl = binet.(fill(5000,100), Ref([500,500]), 0,0.5,3)
gc = mixednet.(fill(10000,100),0,0,7.5)
binets2=(gh,gl,gc,union.(gh,gl));

# +
# Size of 1st and 2nd giant components
ccs = connected_components.(binets1[4])
componentsizes = sort.(broadcast.(length,ccs),rev=true)
@show getindex.(componentsizes,Ref(1:2)).|>Base.Fix2(getindex,1)|>quantile
@show getindex.(componentsizes,Ref(1:2)).|>Base.Fix2(getindex,2)|>quantile;

ccs = connected_components.(binets2[4])
componentsizes = sort.(broadcast.(length,ccs),rev=true)
@show getindex.(componentsizes,Ref(1:2)).|>Base.Fix2(getindex,1)|>quantile
@show getindex.(componentsizes,Ref(1:2)).|>Base.Fix2(getindex,2)|>quantile;
# -

# degree distribution
degdist = mean(get.(Ref(degree_histogram(binets1[4][x])),1:100, 0) for x in 1:length(binets1[4]))
degdist=degdist[degdist.>0]
plot([1;(1 .-min.(1,cumsum(degdist./sum(degdist))))[1:end-1]],xscale=:log10,yscale=:log10,legend = :none,xlabel = "degree", ylabel = "CCDF")
degdist = mean(get.(Ref(degree_histogram(binets2[4][x])),1:100, 0) for x in 1:length(binets2[4]))
degdist=degdist[degdist.>0]
plot!([1;(1 .-min.(1,cumsum(degdist./sum(degdist))))[1:end-1]],xscale=:log10,yscale=:log10,legend = :none,xlabel = "degree", ylabel = "CCDF")|>display

assortativity.(binets1[4])|>Base.Fix2(quantile,(0.,0.025,0.5,0.975,1))

@time R0 = adjacency_spectrum.(gs).|>maximum

# compute adjacency spectrum (~1h)
BLAS.set_num_threads(1)
adjR0s1 = zeros(length(gh))
adjR0s2 = zeros(length(gh))
totalgraphs=union.(binets1[3],binets1[4])
@time Threads.@threads for i in 1:length(totalgraphs)
    adjR0s1[i]=totalgraphs[i]|>adjacency_spectrum|>last
end
totalgraphs=union.(binets2[3],binets2[4])
@time Threads.@threads for i in 1:length(totalgraphs)
    adjR0s2[i]=totalgraphs[i]|>adjacency_spectrum|>last
end

BLAS.set_num_threads(7)
cR0s1 = zeros(length(gh))
cR0s2 = zeros(length(gh))
totalgraphs=union.(binets1[3],binets1[2])
@time for i in 1:length(totalgraphs)
    cR0s1[i]=totalgraphs[i]|>adjacency_spectrum|>last
end
totalgraphs=union.(binets2[3],binets2[2])
@time for i in 1:length(totalgraphs)
    cR0s2[i]=totalgraphs[i]|>adjacency_spectrum|>last
end

approxR0s1 = approxR0.(binets1[1:3]...)
approxR0s2 = approxR0.(binets2[1:3]...)

# save files
 CSV.write("../data/intermediate/networkmodel_spectrum1.csv",DataFrame((adj = adjR0s1, ngm = approxR0s1, com = cR0s1)))
 CSV.write("../data/intermediate/networkmodel_spectrum2.csv",DataFrame((adj = adjR0s2, ngm = approxR0s2, com = cR0s2)))

# read files
ns1 = CSV.read("../data/intermediate/networkmodel_spectrum1.csv", DataFrame)
ns2 = CSV.read("../data/intermediate/networkmodel_spectrum2.csv", DataFrame);

scatter(ns1[:,1],ns1[:,2],xlimit=(0,20),ylimit=(0,20),markersize=1.5,markerstrokewidth=0, label="network 1")
scatter!(ns2[:,1],ns2[:,2],xlimit=(5,20),ylimit=(5,20),markersize=1.5,markerstrokewidth=0, label="network 2")
plot!([0,20],[0,20],color=:black,linewidth=0.75,label="0% error")
plot!([0,20],[0,20*0.95],color=:grey,linestyle=:dash,linewidth=0.75,label="5% error")
plot!([0,20],[0,20*0.9],color=:grey,linestyle=:dashdot,linewidth=0.75,label="10% error")
plot!(xlabel="dominant eigenvalue of adjacency matrix",ylabel="dominant eigenvalue of modelled NGM",size=(400,400)) |> savefigname("../figs/Sfigs/raw/networkspectrum.svg", save = false)

@show quantile(ns1[:,1],(0.025,0.5,0.975))
@show quantile(ns2[:,1],(0.025,0.5,0.975));

@show quantile(1 .-ns1[:,3]./ns1[:,1],(0.025,0.5,0.975))
@show quantile(1 .-ns2[:,3]./ns2[:,1],(0.025,0.5,0.975));

@time am=adjacency_matrix(gs[1])

@time [am*normalize((am^n*fill(1e-4,10000)),1)|>sum for n in 1:10]|>plot



highriskind=((adjacency_matrix(gs[1]).+diagm(fill(1,10000)))*in.(1:10000,Ref(5001:5500)).!=0)

g1=deepcopy(gs[1])
rem_vertices!(g1,(1:10000)[.!highriskind])|>sort

degree(g1)|>histogram

@time am=adjacency_matrix(g1)

@time [am*normalize((am^n*fill(1e-4,1123)),1)|>sum for n in 1:20]|>plot


