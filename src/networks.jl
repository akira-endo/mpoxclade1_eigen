using CSV
using Graphs
using StatsBase
using Plots
using GraphRecipes
using Random
using LinearAlgebra
using DataFrames

function mixednet(size, size_high, meandeg_high , meandeg_random, exponent = 3)
    sf = static_scale_free(size_high, meandeg_high*size_high÷2|>Int, exponent)
    er =erdos_renyi(size, meandeg_random/(size-1))
    union(sf, er)
end
function binet(sizes, sizes_high::AbstractVector{<:Integer}, nedges_high::Integer, meandeg_random::Real ,exponent=3)
    if sizes_high[1]<sizes_high[2] throw("Error: 1st element of sizes_high is smaller than the 2nd") end
    hsizeratio=Int(sizes_high[1]/sizes_high[2])
    nets = mixednet.(sizes, maximum(sizes_high), fill(2nedges_high./maximum(sizes_high),2), meandeg_random,exponent)

    edgepairs = zip(shuffle(edges(nets[1])|>collect),shuffle(edges(nets[2])|>collect))
    binet = SimpleGraph(sum(nv.(nets)),0)
    for ep in edgepairs
        sd = (sizes[1].+ (([src(ep[2]),dst(ep[2])].-1).÷hsizeratio.+1))|>shuffle
        add_edge!(binet, src(ep[1]),sd[1]) # divide by hsizeratio to collapse pseudo vertices
        add_edge!(binet, dst(ep[1]),sd[2])
    end
    binet
end

function approxR0(gh, gl, gc, sizes_high = (500, 100))
    if !(nv(gh)==nv(gl)==nv(gc)) throw("sizes do not match") end
    size = nv(gh)÷2
    C = (degree(union(gh,union(gl,gc)))|>mean)/2
    gs = union(gh,gl)
    d_M = degree(gs)[1:sizes_high[1]]
    nd_M = mean(d_M,StatsBase.weights(d_M)) # neighbour degree for male
    d_F = degree(gs)[size .+ (1:sizes_high[2])]
    nd_F = mean(d_F,StatsBase.weights(d_F)) # neighbour degree for female
    ngm = [0 0 nd_F mean(d_F)*sizes_high[2]/size
        fill(C,1,4)
nd_M mean(d_M)*sizes_high[1]/size 0 0
 fill(C,1,4)]
    eigen(ngm).values[end]|>Real
end