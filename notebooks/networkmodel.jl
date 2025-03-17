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

]activate .

# +
using Graphs
using StatsBase
using Plots
using GraphRecipes
using Random
function mixednet(size, size_high, meandeg_high , meandeg_random, exponent = 3)
    sf = static_scale_free(size_high, meandeg_high*size_high÷2|>Int, exponent)#barabasi_albert(size_high, size_high÷10, 3)
    er =erdos_renyi(size, meandeg_random/(size-1))
    union(sf, er)
end
function mixednet(size, size_high, meandeg_random)
    sf = static_scale_free(size_high, 3size_high, 2.7)#barabasi_albert(size_high, size_high÷10, 3)
    er =erdos_renyi(size, meandeg_random/(size-1))
    union(sf, er)
end

function bimixednet(size, size_high, meandeg_random) # generates 2x size
    nets = [mixednet(size, size_high, meandeg_random) for _ in 1:2]
    binet = blockdiag(nets[1], nets[2]) # join nets
    for e in collect(edges(binet))
        sd = [src(e),dst(e)]
        rewire = @view sd[sample(1:2)]
        rewire .+= rewire[]>size ? -size : +size
        rem_edge!(binet, e)
        add_edge!(binet, sd[1],sd[2])
    end
    binet
end
function bipropnet(sizes, sizes_high::AbstractVector{<:Integer}, nedges_high::Integer, meandeg_random::Real,exponent=3)
    nets = mixednet.(sizes, sizes_high, nedges_high./sizes_high, meandeg_random,exponent)
    edgepairs = zip(shuffle(edges(nets[1])|>collect),shuffle(edges(nets[2])|>collect))
    binet = SimpleGraph(sum(nv.(nets)),0)
    for ep in edgepairs
        sd = ([src(ep[2]),dst(ep[2])].+sizes[1])|>shuffle
        add_edge!(binet, src(ep[1]),sd[1])
        add_edge!(binet, dst(ep[1]),sd[2])
    end
    binet
end
# -

gr(fontfamily="Helvetica",foreground_color_legend = nothing,background_color_legend = nothing, titlefontsize=11, tickfontsize=10, legendfontsize=8,labelfontsize=10,grid=true, tick_direction=:out,size=(400,300))

# illustration
Random.seed!(123)
gsmall = bimixednet(500, 50, 0.5)
@time gp = graphplot(gsmall,curves=false,markercolor = repeat([1,4],inner=nv(gsmall)÷2),markersize=0.1,size=(450,400))



# +
#savefig(gp,"../figs/Sfigs/raw/network.svg")
# -

Random.seed!(123)
@time gs = bimixednet.(fill(5000,1), 500, 0.5)

# Size of 1st and 2nd giant components
ccs = connected_components.(gs)
componentsizes = sort.(broadcast.(length,ccs),rev=true)
@show getindex.(componentsizes,Ref(1:2)).|>Base.Fix2(getindex,1)|>quantile
@show getindex.(componentsizes,Ref(1:2)).|>Base.Fix2(getindex,2)|>quantile;

# degree distribution
degdist = mean(get.(Ref(degree_histogram(gs[x])),1:100, 0) for x in 1:length(gs))
degdist=degdist[degdist.>0]
plot([1;(1 .-min.(1,cumsum(degdist./sum(degdist))))[1:end-1]],xscale=:log10,yscale=:log10,legend = :none,xlabel = "degree", ylabel = "CCDF")

assortativity.(gs)|>Base.Fix2(quantile,(0.,0.025,0.5,0.975,1))

@time R0 = adjacency_spectrum.(gs).|>maximum

degs = degree(gs[1])
@time ordered=sort(1:nv(gs[1]),by=x->degs[x],rev=true)
[mean(degree(gs[x])[ordered[1:1000]],StatsBase.weights(degree(gs[x])[ordered[1:1000]])) for x in 1:1]

[mean(degree(gs[x])[[1:500;5001:5500]],StatsBase.weights(degree(gs[x])[[1:500;5001:5500]])) for x in 1:1]

gcopy=deepcopy(gs)

rem_vertices!.(gcopy,Ref([501:5000;5501:10000]))

R0_h = adjacency_spectrum.(gcopy).|>maximum

scatter(R0,R0_h, xlimit=(0,7),ylimit=(0,7))

@time am=adjacency_matrix(gs[1])

@time [am*normalize((am^n*fill(1e-4,10000)),1)|>sum for n in 1:10]|>plot

degree(gs[1])[1:1000]|>histogram

highriskind=((adjacency_matrix(gs[1]).+diagm(fill(1,10000)))*in.(1:10000,Ref(5001:5500)).!=0)

g1=deepcopy(gs[1])
rem_vertices!(g1,(1:10000)[.!highriskind])|>sort

degree(g1)|>histogram

@time am=adjacency_matrix(g1)

@time [am*normalize((am^n*fill(1e-4,1123)),1)|>sum for n in 1:20]|>plot

Random.seed!(123)
gsmall = union(bipropnet(500, [50,10], 90,0,[2.7,3.1]),bipropnet(500, [50,10], 0,0.5,3))
@time gp = graphplot(gsmall,curves=false,markercolor = repeat([1,4],inner=nv(gsmall)÷2),markersize=0.1,size=(450,400))

Random.seed!(123)
gsmall = union(bipropnet(500, [50,10], 50,0,[2.7,3.1]),bipropnet(500, [50,10], 0,0.5,3))
@time gp = graphplot(gsmall,curves=false,markercolor = repeat([1,4],inner=nv(gsmall)÷2),markersize=0.1,size=(450,400))

function approxR0(gh, gl, gc, sizes_high = (500, 100))
    C = degree(union(gl,gc))|>mean
    gs = union(gh,gl)
    d_M = degree(gs)[1:sizes_high[1]]
    nd_M = mean(d_M,StatsBase.weights(d_M)) # neighbour degree for male
    d_F = degree(gs)[5000 .+ (1:sizes_high[2])]
    nd_F = mean(d_F,StatsBase.weights(d_F)) # neighbour degree for female
    ngm = [0 0 0 nd_F mean(d_F) 0
        fill(C*0.1,1,6)
        fill(C*0.9,1,6)
nd_M mean(d_M) 0 0 0 0
 fill(C*0.02,1,6)
 fill(C*0.98,1,6)]
    eigen(ngm).values[end]|>Real
end
 

scatter(trueR0s,aR0,xlimit=(5,8),ylimit=(5,8))
plot!([0,10],[0,10])

@time aR0 = approxR0.(gh,gl,gc)

using LinearAlgebra
trueR0s = zeros(length(gh))
totalgraphs=union(union(gh[i],gl[i]),gc[i])
@time Threads.@threads for (i, g) in enumerate(totalgraphs)
    trueR0s[i]=g|>adjacency_spectrum|>last
end
trueR0s

Random.seed!(123)
@time gh = bipropnet.(fill(5000,100), Ref([500,100]),900, 0., Ref([2.7,3.1]))
gl = bipropnet.(fill(5000,100), Ref([500,100]), 0,0.5,3)
gc = mixednet.(fill(10000,100),0,0,2.5)

Random.seed!(123)
@time gh = bipropnet.(fill(5000,10), Ref([500,100]),500, 0., Ref([2.7,3.1]))
gl = bipropnet.(fill(5000,10), Ref([500,100]), 0,0.5,3)
gc = mixednet.(fill(10000,10),0,0,2.5)
trueR0s_500 = zeros(length(gh))
totalgraphs=union.(union.(gh,gl),gc)
@time Threads.@threads for i in 1:length(totalgraphs)
    trueR0s_500[i]=totalgraphs[i]|>adjacency_spectrum|>last
end
trueR0s_500

Random.seed!(123)
@time gh = bipropnet.(fill(5000,10), Ref([500,100]),900, 0., Ref([2.7,3.1]))
gl = bipropnet.(fill(5000,10), Ref([500,100]), 0,0.5,3)
gc = mixednet.(fill(10000,10),0,0,5)
trueR0s_500 = zeros(length(gh))
totalgraphs=union.(union.(gh,gl),gc)
@time Threads.@threads for i in 1:length(totalgraphs)
    trueR0s_500[i]=totalgraphs[i]|>adjacency_spectrum|>last
end
trueR0s_500

@time aR0_500 = approxR0.(gh,gl,gc)
scatter(trueR0s_500,ylimit=(0,8))
scatter!(aR0_500)

@time R0 = adjacency_spectrum(union(gh,gl))

adjacency_spectrum(gh)

adjacency_spectrum(union(gl,gc))

@show degree(union(gl,gc))|>mean
@show degree(union(gl,gh))|>mean
@show mean(degree(union(gl,gh))[1:500],StatsBase.weights(degree(union(gl,gh))[1:500]))
@show mean(degree(union(gl,gh))[5001:5100],StatsBase.weights(degree(union(gl,gh))[5001:5100]))
mean(degree(gh),StatsBase.weights(degree(gh)))

[0 0 0 11.28 8.92 0
 fill(3*0.1,1,6)
fill(3*0.9,1,6)
2.68 1.784 0 0 0 0
 fill(3*0.02,1,6)
 fill(3*0.98,1,6)]|>eigen


