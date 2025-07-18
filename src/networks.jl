using CSVFiles
using Graphs
using StatsBase
using Plots
using GraphRecipes
using Random
using LinearAlgebra
using DataFrames
using MCMCChains

function mixednet(size, size_high, meandeg_high , meandeg_random, exponent = 3)
    sf = static_scale_free(size_high, meandeg_high*size_high÷2|>Int, exponent)
    er =erdos_renyi(size, meandeg_random/(size-1))
    union(sf, er)
end

function binet(sizes, sizes_high::AbstractVector{<:Integer}, nedges_high::Integer, meandeg_random::Real ,exponent=3, assort = 0)
    if sizes_high[1]<sizes_high[2] throw("Error: 1st element of sizes_high is smaller than the 2nd") end
    hsizeratio=Int(sizes_high[1]/sizes_high[2])
    nets = mixednet.(sizes, maximum(sizes_high), fill(2nedges_high./maximum(sizes_high),2), meandeg_random,exponent)
    d1=degree(nets[1]); d2=degree(nets[2])

    edgepairs = zip(shuffle(edges(nets[1])|>collect),shuffle(edges(nets[2])|>collect))
    binet = SimpleGraph(sum(nv.(nets)),0)
    for ep in edgepairs
        sd = [src(ep[2]),dst(ep[2])]|>shuffle # edge in female sf network
        # degree assortativity
        if rand()<assort && (d1[src(ep[1])]-d1[dst(ep[1])])*(d2[sd[1]]-d2[sd[2]])<0
            reverse!(sd)
        end

        sd .= (sizes[1].+ ((sd.-1).÷hsizeratio.+1)) # divide by hsizeratio to collapse female vertices
        add_edge!(binet, src(ep[1]),sd[1])
        add_edge!(binet, dst(ep[1]),sd[2])
    end
    binet
end

function binet2(sizes, sizes_high::AbstractVector{<:Integer}, nedges_high::Integer, meandeg_random::Real ,exponent=3, assort = 0)
    if sizes_high[1]<sizes_high[2] throw("Error: 1st element of sizes_high is smaller than the 2nd") end
    hsizeratio=Int(sizes_high[1]/sizes_high[2])
    nets = mixednet.(sizes, maximum(sizes_high), fill(2nedges_high./maximum(sizes_high),2), meandeg_random,exponent)
    d1=degree(nets[1]); d2=groupsplit(degree(nets[2]),sizes[1]÷hsizeratio).|>sum

    edgepairs = zip(shuffle(edges(nets[1])|>collect),shuffle(edges(nets[2])|>collect))
    binet = SimpleGraph(sum(nv.(nets)),0)
    for ep in edgepairs
        sd = [src(ep[2]),dst(ep[2])]|>shuffle # edge in female sf network
        sd .= (sizes[1].+ ((sd.-1).÷hsizeratio.+1)) # divide by hsizeratio to collapse female vertices
        # degree assortativity
        if rand()<assort && (d1[src(ep[1])]-d1[dst(ep[1])])*(d2[sd[1]-sizes[1]]-d2[sd[2]-sizes[1]])<0
            reverse!(sd)
        end
        add_edge!(binet, src(ep[1]),sd[1])
        add_edge!(binet, dst(ep[1]),sd[2])
    end
    binet
end


function approxR0(gh, gl, gc, sizes_high = (500, 100), nd_M=nothing, nd_F=nothing)
    if !(nv(gh)==nv(gl)==nv(gc)) throw("sizes do not match") end
    size = nv(gh)÷2
    C = (degree(union(gh,union(gl,gc)))|>mean)/2
    gs = union(gh,gl)
    d_M = degree(gs)[1:sizes_high[1]]
    if isnothing(nd_M) nd_M = mean(d_M,StatsBase.weights(d_M)) end # neighbour degree for male
    d_F = degree(gs)[size .+ (1:sizes_high[2])]
    if isnothing(nd_F) nd_F = mean(d_F,StatsBase.weights(d_F)) end # neighbour degree for female
    ngm = [0 0 nd_F mean(d_F)*sizes_high[2]/size
        fill(C,1,4)
nd_M mean(d_M)*sizes_high[1]/size 0 0
 fill(C,1,4)]
    eigen(ngm).values[end]|>Real
end

adjacency_dominant(g)=mean(degree(g),pweights(eigenvector_centrality(g)))
function bipartiteassortativity(g)
    newg = SimpleDiGraph(nv(g),0)
    for e in edges(g) add_edge!(newg,src(e),dst(e)) end
    assortativity(newg)
end