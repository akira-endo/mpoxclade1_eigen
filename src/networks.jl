#using CSV
using CategoricalArrays
using CSVFiles
using Colors
using DataFrames
using DataFramesMeta
using FStrings
using Graphs
using GraphRecipes
using LinearAlgebra
using MCMCChains
using Plots
using Pipe
using Random
using SimpleWeightedGraphs
using StatsBase
using StatsPlots

const BinetTuple = NTuple{4, Vector{SimpleGraph{Int64}}}
const BinetWeightedTuple = NTuple{4, Vector{SimpleWeightedGraph{Int64, Float64}}}
const UnionBinetTuple = Union{BinetTuple, BinetWeightedTuple}
const UnionGraph = Union{SimpleGraph, SimpleWeightedGraph}

"""Create a mixture of scale-free and erdos_renyi graphs.

Args:
- size: the size of network for erdos_rency part.
- size_high: the size of network for scale-free part.
- meandeg_high: the mean degree of scale-free part.
- meandeg_random: the mean degree of erdos_renyi part.
- exponent: the exponent of scale-free part.
"""
function mixednet(size::Int64, size_high::Int64,
	meandeg_high::Real, meandeg_random::Real, exponent::Real = 3)
	sf = static_scale_free(size_high, meandeg_high*size_high÷2|>Int, exponent; finite_size_correction = false)
	er = erdos_renyi(size, meandeg_random/(size-1))
	union(sf, er)
end


"""Create a bipartite network from two mixed networks.

Args:
- size: the size of network.
- sizes_high: a vector of two integers, the latter elements is used to determine the ratio.
- meandeg_high: the number of edges in scale-free part.
- meandeg_random: the mean degree of erdos_renyi part.
- exponent: the exponent of scale-free part.
- assort: the degree assortativity, between 0 and 1. 0 for no assortativity.
"""
function binet(size::Int64, sizes_high::AbstractVector{<:Integer},
	meandeg_high::Real, meandeg_random::Real,
	exponent::AbstractVector{<:Real} = [3], assort::Real = 0)
	if sizes_high[1]<sizes_high[2]
		throw("Error: 1st element of sizes_high is smaller than the 2nd")
	end
	hsizeratio=Int(sizes_high[1]/sizes_high[2])
	nets = mixednet.(size, maximum(sizes_high), fill(meandeg_high, 2), meandeg_random, exponent)
	d1=degree(nets[1]);
	d2=degree(nets[2])

	edgepairs = zip(shuffle(edges(nets[1])|>collect), shuffle(edges(nets[2])|>collect))
	g_binet = SimpleGraph(sum(nv.(nets)), 0)
	for ep in edgepairs
		sd = [src(ep[2]), dst(ep[2])]|>shuffle # edge in female sf network
		# degree assortativity
		if rand()<assort && (d1[src(ep[1])]-d1[dst(ep[1])])*(d2[sd[1]]-d2[sd[2]])<0
			reverse!(sd)
		end

		sd .= (size .+ ((sd .- 1) .÷ hsizeratio .+ 1)) # divide by hsizeratio to collapse female vertices
		add_edge!(g_binet, src(ep[1]), sd[1])
		add_edge!(g_binet, dst(ep[1]), sd[2])
	end
	g_binet
end

function binet2(sizes, sizes_high::AbstractVector{<:Integer},
	meandeg_high::Real, meandeg_random::Real,
	exponent::AbstractVector{<:Real} = [3], assort = 0)
	if sizes_high[1]<sizes_high[2]
		throw("Error: 1st element of sizes_high is smaller than the 2nd")
	end
	hsizeratio = Int(sizes_high[1] / sizes_high[2])
	nets = mixednet.(sizes, maximum(sizes_high), fill(meandeg_high, 2), meandeg_random, exponent)
	d1=degree(nets[1]);
	d2=groupsplit(degree(nets[2]), sizes[1]÷hsizeratio) .|> sum

	edgepairs = zip(shuffle(edges(nets[1])|>collect), shuffle(edges(nets[2])|>collect))
	binet = SimpleGraph(sum(nv.(nets)), 0)
	for ep in edgepairs
		sd = [src(ep[2]), dst(ep[2])]|>shuffle # edge in female sf network
		sd .= (sizes[1] .+ ((sd .- 1) .÷ hsizeratio .+ 1)) # divide by hsizeratio to collapse female vertices
		# degree assortativity
		if rand()<assort && (d1[src(ep[1])]-d1[dst(ep[1])])*(d2[sd[1]-sizes[1]]-d2[sd[2]-sizes[1]])<0
			reverse!(sd)
		end
		add_edge!(binet, src(ep[1]), sd[1])
		add_edge!(binet, dst(ep[1]), sd[2])
	end
	binet
end

function remove_edges_randomly(g::SimpleGraph, rem_p::Float64)::SimpleGraph
	n_vert = nv(g)
	edge_pair = g |> edges |> collect
	g_new = SimpleGraph(nv(g), 0)
	for ep in edge_pair
		if rand() < rem_p
			continue
		end
		add_edge!(g_new, src(ep), dst(ep))
	end
	g_new
end

function get_dominant_eigenvalue_vector(ngm::Matrix)::Tuple{<:Real, Vector{Float64}}
	eg = eigen(ngm)
	dom_evec = @pipe eg.vectors[:, end] .|> abs |> normalize(_, 1)
	(eg.values[end], dom_evec)
end

vec_vec_to_df(v) = @pipe mapreduce(permutedims, vcat, v) |> DataFrame(_, :auto)
function eigens_to_vec_and_df(eg)
	e_vals, e_vecs = zip(eg...)
	(e_vals |> collect, vec_vec_to_df(e_vecs))
end

function adj_props_groups(binets::UnionBinetTuple, adj_eigen::Vector)::DataFrame
	ind_m = vcat(1:5000)
	ind_f = vcat(5001:10000)

    adj_props = []
    for (gc, gh, ghl, evec) in zip(binets[3], binets[1], binets[4], adj_eigen)
        gchl = union(gc, ghl)
		adj_whole = adjacency_matrix(gchl)
		adj_sex = adjacency_matrix(gh)
		adj_com = adj_whole .- adj_sex

        # Calculate proportions of infectees (not infectors) through sexual and community
        prop_sex = adj_sex * evec
        prop_com = adj_com * evec
        prop_mH = sum(prop_sex[ind_m])
        prop_mL = sum(prop_com[ind_m])
        prop_fH = sum(prop_sex[ind_f])
        prop_fL = sum(prop_com[ind_f])
        props = normalize([prop_mH, prop_mL, prop_fH, prop_fL], 1)
        push!(adj_props, props)
    end
    adj_df = vec_vec_to_df(adj_props)
    @rename!(adj_df, :adj_mH = :x1, :adj_mL = :x2, :adj_fH = :x3, :adj_fL = :x4)
    return adj_df
end

"""
Usage:
```
approxR0s1 = approxR0.(binets1[1:3]...)
R0 = approxR0(binets1[1][1], binets1[1][2], binets1[1][3])
````
"""
function approxR0(gh::SimpleGraph, gl, gc, ghl, sizes_high = (500, 100),
	nd_M = nothing, nd_F = nothing;
	verbose = false,
)::Tuple{<:Real, Vector{Float64}}
	if !(nv(gh)==nv(gl)==nv(gc))
		throw("sizes do not match")
	end
	size = nv(gh)÷2
	C = (get_degree(union(gh, union(gl, gc)))|>mean)/2
	gs = union(gh, gl)
	d_M = get_degree(gs)[1:sizes_high[1]]
	if isnothing(nd_M)
		nd_M = mean(d_M, StatsBase.weights(d_M))
	end # neighbour degree for male
	d_F = get_degree(gs)[size .+ (1:sizes_high[2])]
	if isnothing(nd_F)
		nd_F = mean(d_F, StatsBase.weights(d_F))
	end # neighbour degree for female
	ngm = [
		0 0 nd_F mean(d_F)*sizes_high[2]/size
		fill(C, 1, 4)
		nd_M mean(d_M)*sizes_high[1]/size 0 0
		fill(C, 1, 4)]
	if verbose == true
		println(f"""Degrees: nd_M = {nd_M:.2f}, nd_F = {nd_F:.2f},
		mean_d_F = {mean(d_F):.2f}, mean_d_M = {mean(d_M):.2f},
		C = {C:.2f}
		""")
	end
	get_dominant_eigenvalue_vector(ngm)
end

function approxR0_com(gh::SimpleGraph, gl, gc, ghl)::Tuple{<:Real, Vector{Float64}}
	if !(nv(gl)==nv(gc))
		throw("sizes do not match")
	end
	size = nv(gl)÷2
	C = (get_degree(union(gh, union(gl, gc)))|>mean)/2
	ngm = [
		0 0 0 0
		fill(C, 1, 4)
		0 0 0 0
		fill(C, 1, 4)]
	get_dominant_eigenvalue_vector(ngm)
end

function approxR0_eigen(binets::UnionBinetTuple)
	n_sim = length(binets[1])
	nd_M = zeros(n_sim)
	nd_F = zeros(n_sim)
	for i in 1:n_sim
		ev = normalize(binets[1][i] |> eigenvector_centrality, 1)
		nd_M[i] = mean(get_degree(binets[1][i])[1:500], pweights(ev[1:500]))
		nd_F[i] = mean(get_degree(binets[1][i])[5001:5100], pweights(ev[5001:5100]))
	end
	approxR0.(binets[1:4]..., Ref((500, 100)), nd_M, nd_F)
end

function approxR0(gh::SimpleWeightedGraph, gl::SimpleWeightedGraph,
	gc::SimpleWeightedGraph, ghl::SimpleWeightedGraph,
	sizes_high = (500, 100),
	nd_M = nothing, nd_F = nothing,
	verbose = false,
)::Tuple{<:Real, Vector{Float64}}
	if !(nv(gh)==nv(gl)==nv(gc)==nv(ghl))
		throw("sizes do not match")
	end
	size = nv(gh)÷2
	# Since community contacts do not separate contact types.
	gc_w = edges(gc) |> collect |> (x -> x[1]) |> weight
	C = (get_degree2(union(gc, ghl))|>mean)/2 * gc_w

	d_M = get_degree(ghl)[1:sizes_high[1]]
	if isnothing(nd_M)
		nd_M = mean(d_M, StatsBase.weights(d_M))
	end # neighbour degree for male
	d_F = get_degree(ghl)[size .+ (1:sizes_high[2])]
	if isnothing(nd_F)
		nd_F = mean(d_F, StatsBase.weights(d_F))
	end # neighbour degree for female
	ngm = [
		0 0 nd_F mean(d_F)*sizes_high[2]/size
		fill(C, 1, 4)
		nd_M mean(d_M)*sizes_high[1]/size 0 0
		fill(C, 1, 4)]
	if verbose == true
		println(f"""Degrees: nd_M = {nd_M:.2f}, nd_F = {nd_F:.2f},
		mean_d_F = {mean(d_F):.2f}, mean_d_M = {mean(d_M):.2f},
		C = {C:.2f}
		""")
	end
	get_dominant_eigenvalue_vector(ngm)
end

function approxR0_com(gh::SimpleWeightedGraph, gl::SimpleWeightedGraph,
	gc::SimpleWeightedGraph, ghl::SimpleWeightedGraph,
)::Tuple{<:Real, Vector{Float64}}
	if !(nv(gl)==nv(gc)==nv(ghl))
		throw("sizes do not match")
	end
	size = nv(gl)÷2
	# Since community contacts do not separate contact types.
	gc_w = edges(gc) |> collect |> (x -> x[1]) |> weight
	C = (get_degree2(union(gc, ghl))|>mean)/2 * gc_w
	ngm = [
		0 0 0 0
		fill(C, 1, 4)
		0 0 0 0
		fill(C, 1, 4)]
	get_dominant_eigenvalue_vector(ngm)
end

get_degree(g::SimpleGraph)::Vector = degree(g)
get_degree(g::SimpleWeightedGraph)::Vector = degree_matrix(g; dir = :both) |> diag |> Vector
get_degree2(g::SimpleGraph)::Vector = degree(g)
get_degree2(g::SimpleWeightedGraph)::Vector = degree(g)
function adjacency_dominant(g::UnionGraph; eig_vec = nothing)
	if isnothing(eig_vec)
		eig_vec = eigenvector_centrality(g)
	end
	mean(get_degree(g), pweights(eig_vec))
end

function bipartiteassortativity(g)
	newg = SimpleDiGraph(nv(g), 0)
	for e in edges(g)
		add_edge!(newg, src(e), dst(e))
	end
	assortativity(newg)
end

EXPONENTS = [3.00, 3.00] # [3.00, 3.00] # [2.70, 2.90]
function create_networks(n_sim::Int64,
	meandeg_high, meandeg_low, meandeg_com)
	exponents = EXPONENTS
	gh = binet.(fill(5000, n_sim), Ref([500, 100]), meandeg_high, 0.0, Ref(exponents))
	gl = binet.(fill(5000, n_sim), Ref([500, 500]), 0.0, meandeg_low, Ref(exponents))
	gc = mixednet.(fill(10000, n_sim), 0, 0, meandeg_com)
	(gh, gl, gc, union.(gh, gl))
end

create_networks1(n_sim::Int64 = 100) = create_networks(n_sim, 4.0, 0.5, 7.5)
create_networks2(n_sim::Int64 = 100) = create_networks(n_sim, 2.0, 0.5, 7.5)

function create_networks1_assort(n_sim::Int64 = 100)
	exponents = EXPONENTS
	gh = binet2.(fill(5000, n_sim), Ref([500, 100]), 4.0, 0.0, Ref(exponents), 1)
	gl = binet.(fill(5000, n_sim), Ref([500, 500]), 0, 0.5, Ref(exponents))
	gc = mixednet.(fill(10000, n_sim), 0, 0, 7.5)
	binets = (gh, gl, gc, union.(gh, gl))
	binets
end

function create_networks2_assort(n_sim::Int64 = 100)
	exponents = EXPONENTS
	gh = binet2.(fill(5000, n_sim), Ref([500, 100]), 2.0, 0.0, Ref(exponents), 1)
	gl = binet.(fill(5000, n_sim), Ref([500, 500]), 0, 0.5, Ref(exponents))
	gc = mixednet.(fill(10000, n_sim), 0, 0, 7.5)
	binets = (gh, gl, gc, union.(gh, gl))
	binets
end

function create_weighted_network1(n_sim::Int64 = 100;
	weight_factor = 1.0, community_weight = 1.0)
	binets = create_networks(n_sim, 4.0/weight_factor, 0.5/weight_factor, 7.5/community_weight)
	binets_to_weighted(binets, weight_factor; community_weight = community_weight);
end

function create_weighted_network2(n_sim::Int64 = 100;
	weight_factor = 1.0, community_weight = 1.0)
	binets = create_networks(n_sim, 2.0/weight_factor, 0.25/weight_factor, 7.5/community_weight)
	binets_to_weighted(binets, weight_factor; community_weight = community_weight);
end

function check_first_second_giant_components(binets::UnionBinetTuple)
	println("GCC for sexual contact network")
	ccs = connected_components.(binets[4])
	componentsizes = sort.(broadcast.(length, ccs), rev = true)
	@show getindex.(componentsizes, Ref(1:2)) .|> Base.Fix2(getindex, 1)|>quantile
	@show getindex.(componentsizes, Ref(1:2)) .|> Base.Fix2(getindex, 2)|>quantile;

	println("GCC for the whole contact network")
	ccs = connected_components.(union.(binets[3], binets[4]))
	componentsizes = sort.(broadcast.(length, ccs), rev = true)
	@show getindex.(componentsizes, Ref(1)) .|> Base.Fix2(getindex, 1)|>quantile
	#@show getindex.(componentsizes, Ref(1:2)) .|> Base.Fix2(getindex, 2)|>quantile;
end

function plot_degree(binets1::UnionBinetTuple, binets2::UnionBinetTuple; title = "", x_max = 150)
	g = union.(binets1[4], binets1[3])
	degdist = mean(get.(Ref(degree_histogram(g[x])), 1:x_max, 0) for x in 1:length(g))
	degdist=degdist[degdist .> 0]
	yticks = ([1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5])
	xlim = [1, x_max]
	pl = plot(
		scale = :log10, yscale = :log10,
		xlabel = "degree", ylabel = "CCDF",
		yticks = yticks, xlim = xlim, title = title,
	)
	plot!(pl, [1; (1 .- min.(1, cumsum(degdist ./ sum(degdist))))[1:(end-1)]],
		label = "Network 1")

	g = union.(binets2[4], binets2[3])
	degdist = mean(get.(Ref(degree_histogram(g[x])), 1:x_max, 0) for x in 1:length(g))
	degdist=degdist[degdist .> 0]
	plot!(pl, [1; (1 .- min.(1, cumsum(degdist ./ sum(degdist))))[1:(end-1)]],
		label = "Network 2")
	display(pl)
end

function get_adjacency_R0(binets::UnionBinetTuple)::Tuple{Vector, Vector}
	n_sim = length(binets[1])
	adjR0s = zeros(n_sim)
	adj_eigen = []
	totalgraphs = union.(binets[3], binets[4])
	for i in 1:n_sim
		eig_vec = @pipe totalgraphs[i] |> eigenvector_centrality |> normalize(_, 1)
		push!(adj_eigen, eig_vec)
		adjR0s[i] = adjacency_dominant(totalgraphs[i]; eig_vec=eig_vec)
	end
	(adjR0s, adj_eigen)
end

function get_adjacency_R0_com(binets::UnionBinetTuple)::Tuple{Vector, Vector}
	n_sim = length(binets[1])
	adjR0s_com = zeros(n_sim)
	adj_eigen = []
	#totalgraphs = union.(binets[2], binets[3])
	totalgraphs = binets[3]
	for i in 1:n_sim
		eig_vec = @pipe totalgraphs[i] |> eigenvector_centrality |> normalize(_, 1)
		push!(adj_eigen, eig_vec)
		adjR0s_com[i] = adjacency_dominant(totalgraphs[i]; eig_vec=eig_vec)
	end
	(adjR0s_com, adj_eigen)
end

function plot_simple_adj_ngm_R0(ns1::DataFrame, ns2::DataFrame; lim = (8, 18), kwds...)
	lim_max = lim[2]
	pl = plot(;xlim = lim, ylim = lim, kwds...)
	scatter!(pl, ns1[:, :adj], ns1[:, :ngm],
		markersize = 1.5, markerstrokewidth = 0, label = "network 1")
	scatter!(ns2[:, :adj], ns2[:, :ngm],
		markersize = 1.5, markerstrokewidth = 0, label = "network 2")
	plot!([0, lim_max], [0, lim_max], color = :black, linewidth = 0.75, label = "0% error")
	plot!([0, lim_max], [0, lim_max*0.95], color = :grey, linestyle = :dash, linewidth = 0.75, label = "5% error")
	plot!([0, lim_max], [0, lim_max*0.9], color = :grey, linestyle = :dashdot, linewidth = 0.75, label = "10% error")
	plot!(xlabel = "dominant eigenvalue of adjacency matrix", ylabel = "dominant eigenvalue of modelled NGM",
		size = (400, 400))
	pl
end

function plot_simple_rel_R0(ns1::DataFrame, ns2::DataFrame; lim = (0, 61), kwds...)
	pl = plot(;xlim = lim, ylim = lim, kwds...)
	scatter!(pl, ns1[:, :adj_relRsex], ns1[:, :ngm_relRsex],
		markersize = 1.5, markerstrokewidth = 0, label = "network 1")
	scatter!(ns2[:, :adj_relRsex], ns2[:, :ngm_relRsex],
		markersize = 1.5, markerstrokewidth = 0, label = "network 2")
	add_pp_error!(pl, lim[2])
	plot!(pl, xlabel = "% sexual contribution in adjacency matrix", ylabel = "% sexual contribution in modelled NGM",
		size = (400, 400))
end

function add_pp_error!(pl::Plots.Plot, lim_max)
	plot!(pl, [0, lim_max], [0, lim_max], color = :black, linewidth = 0.75, label = "0 pp error")
	#plot!(pl, [0, lim_max], [0, lim_max*0.80], color = :grey, linestyle = :dash, linewidth = 0.75, label = "")
	plot!(pl, [0, lim_max], [0 - 5, lim_max - 5], color = :grey, linestyle = :dash, linewidth = 0.75, label = "5 pp error")
	plot!(pl, [0, lim_max], [0 - 10, lim_max - 10], color = :grey, linestyle = :dash, linewidth = 0.75, label = "10 pp error")
	#plot!(pl, [0, lim_max], [0 + 5, lim_max + 5], color = :grey, linestyle = :dash, linewidth = 0.75, label = "")
	#plot!(pl, [0, lim_max], [0 + 10, lim_max + 10], color = :grey, linestyle = :dashdot, linewidth = 0.75, label = "")
end

function plot_rel_R0_sex(ns1::DataFrame, ns2::DataFrame; version = 1, lim = (0, 61))
	lim_max = lim[2]
	nw1 = (version==1) ? "network 1" : (version == 2) ? "network 3" : "network 5"
	nw2 = (version==1) ? "network 2" : (version == 2) ? "network 4" : "network 6"
	pl = plot(xlimit = lim, ylimit = lim)
	markersize = 1.9
	markersize_adj = 0.8

	scatter!(pl, ns1[:, :adj_relRsex], ns1[:, :ngm_relRsex],
		markersize = markersize, markerstrokewidth = 0, color = 1,
		label = "$nw1 (neighbour degree)", alpha = 1.0)
	scatter!(pl, ns1[:, :adj_relRsex], ns1[:, :ngm_e_relRsex],
		markersize = markersize + markersize_adj, markerstrokewidth = 0.5, color = colorant"navy",
		marker = :+, label = "$nw1 (eigenvalue)", alpha = 1.0)
	scatter!(pl, ns2[:, :adj_relRsex], ns2[:, :ngm_relRsex],
		markersize = markersize, markerstrokewidth = 0, color = 2,
		label = "$nw2 (neighbour degree)")
	scatter!(pl, ns2[:, :adj_relRsex], ns2[:, :ngm_e_relRsex],
		markersize = markersize + markersize_adj, markerstrokewidth = 0.5, color = colorant"firebrick4",
		marker = :+, label = "$nw2 (eigenvalue)")
	add_pp_error!(pl, lim_max)
	plot!(pl, xlabel = "% sexual contribution in adjacency matrix", ylabel = "% sexual contribution in modelled NGM",
		size = (400, 400))
	pl
end

function plot_adj_ngm_eigen_R0(ns1::DataFrame, ns2::DataFrame; version = 1, lim = (0, 20), kwds...)
	lim_max = lim[2]
	nw1 = (version==1) ? "network 1" : (version == 2) ? "network 3" : "network 5"
	nw2 = (version==1) ? "network 2" : (version == 2) ? "network 4" : "network 6"
	markersize = 1.9
	markersize_adj = 0.8

	pl = plot(; kwds...)
	scatter!(pl, ns1[:, :adj], ns1[:, :ngm], xlimit = lim, ylimit = lim,
		markersize = markersize, markerstrokewidth = 0, color = 1,
		label = "$nw1 (neighbour degree)")
	scatter!(ns1[:, :adj], ns1[:, :ngm_e], xlimit = lim, ylimit = lim,
		markersize = markersize + markersize_adj, markerstrokewidth = 0.5, color = :navy,
		marker = :+, label = "$nw1 (eigenvalue)", alpha = 1.0)
	scatter!(ns2[:, :adj], ns2[:, :ngm], xlimit = lim, ylimit = lim,
		markersize = markersize, markerstrokewidth = 0, color = 2,
		label = "$nw2 (neighbour degree)")
	scatter!(ns2[:, :adj], ns2[:, :ngm_e], xlimit = lim, ylimit = lim,
		markersize = markersize + markersize_adj, markerstrokewidth = 0.5, color = :firebrick4,
		marker = :+, label = "$nw2 (eigenvalue)", alpha = 1.0)

	plot!([0, lim_max], [0, lim_max], color = :black, linewidth = 0.75, label = "0% error")
	plot!([0, lim_max], [0, lim_max*0.95], color = :grey, linestyle = :dash, linewidth = 0.75, label = "5% error")
	plot!([0, lim_max], [0, lim_max*0.9], color = :grey, linestyle = :dashdot, linewidth = 0.75, label = "10% error")
	plot!(xlabel = "dominant eigenvalue of adjacency matrix", ylabel = "dominant eigenvalue of modelled NGM",
		size = (400, 400))
	pl
end

function modify_results(
	adjR0s, adjR0s_com,
	approxR0s, approxR0s_com, approxeig,
	ngm_prop, ngm_e_prop, adj_prop,
)
	@rename!(ngm_prop, :ngm_mH = :x1, :ngm_mL = :x2, :ngm_fH = :x3, :ngm_fL = :x4)
	@rename!(ngm_e_prop, :ngm_e_mH = :x1, :ngm_e_mL = :x2, :ngm_e_fH = :x3, :ngm_e_fL = :x4)
	ns = DataFrame(
		(adj = adjR0s,
		adj_com = adjR0s_com,
		ngm = approxR0s,
		ngm_com = approxR0s_com,
		ngm_e = approxeig)
	)
	ns = hcat(ns, ngm_prop, ngm_e_prop, adj_prop)
	ns = @chain ns begin
		@transform :adj_relRsex = (:adj .- :adj_com) ./ :adj .* 100
		@transform :ngm_relRsex = (:ngm .- :ngm_com) ./ :ngm .* 100
		@transform :ngm_e_relRsex = (:ngm_e .- :ngm_com) ./ :ngm_e .* 100
	end
	ns
end

function get_summarised_R0_values(binets1::UnionBinetTuple, binets2::UnionBinetTuple)
	adjR0s1, adj_eigen1 = get_adjacency_R0(binets1);
	adjR0s2, adj_eigen2 = get_adjacency_R0(binets2);
	adjR0s1_com, _ = get_adjacency_R0_com(binets1);
	adjR0s2_com, _ = get_adjacency_R0_com(binets2);

	approxR0s1, ngm_prop1 = approxR0.(binets1[1:4]...) |> eigens_to_vec_and_df
	approxR0s2, ngm_prop2 = approxR0.(binets2[1:4]...) |> eigens_to_vec_and_df
	approxR0s1_com, _ = approxR0_com.(binets1[1:4]...) |> eigens_to_vec_and_df
	approxR0s2_com, _ = approxR0_com.(binets1[1:4]...) |> eigens_to_vec_and_df
	approxeig1, ngm_e_prop1 = approxR0_eigen(binets1) |> eigens_to_vec_and_df
	approxeig2, ngm_e_prop2 = approxR0_eigen(binets2) |> eigens_to_vec_and_df

	adj_prop1 = adj_props_groups(binets1, adj_eigen1)
	adj_prop2 = adj_props_groups(binets2, adj_eigen2)
	ns1 = modify_results(
		adjR0s1, adjR0s1_com,
		approxR0s1, approxR0s1_com, approxeig1,
		ngm_prop1, ngm_e_prop1, adj_prop1,
	)
	ns2 = modify_results(
		adjR0s2, adjR0s2_com,
		approxR0s2, approxR0s2_com, approxeig2,
		ngm_prop2, ngm_e_prop2, adj_prop2,
	)
	cols = [:adj_relRsex, :ngm_relRsex, :ngm_e_relRsex]
	@show ns1[:, cols] |> describe
	@show ns2[:, cols] |> describe;
	(ns1, ns2)
end

function analyse_networks_and_summarise_results(binets1::BinetTuple, binets2::BinetTuple)
	# Size of 1st and 2nd giant components
	check_first_second_giant_components(binets1)
	check_first_second_giant_components(binets2);
	# degree distribution
	plot_degree(binets1, binets2)

	# Check bipartite assortativity
	println("Assortativity for bipartite network")
	bipartiteassortativity.(binets1[1])|>Base.Fix2(quantile, (0.0, 0.025, 0.5, 0.975, 1)) |> display
	bipartiteassortativity.(binets2[1])|>Base.Fix2(quantile, (0.0, 0.025, 0.5, 0.975, 1)) |> display
	# Check assortativity
	println("Assortatativity for the whole network")
	assortativity.(binets1[1])|>Base.Fix2(quantile, (0.0, 0.025, 0.5, 0.975, 1))|>display
	assortativity.(binets2[1])|>Base.Fix2(quantile, (0.0, 0.025, 0.5, 0.975, 1)) |> display

	println("Mean degree for the high sexual activity network")
	@pipe degree.(binets1[1]) .|> mean |> quantile(_, [0.1, 0.5, 0.9]) |> (x -> x .* 10) |> println
	@pipe degree.(binets2[1]) .|> mean |> quantile(_, [0.1, 0.5, 0.9]) |> (x -> x .* 10) |> println

	get_summarised_R0_values(binets1, binets2)
end

function analyse_networks_and_summarise_results(binets1::BinetWeightedTuple, binets2::BinetWeightedTuple)
	# Size of 1st and 2nd giant components
	check_first_second_giant_components(binets1)
	check_first_second_giant_components(binets2);
	# degree distribution
	plot_degree(binets1, binets2)

	println("Mean degree for the high sexual activity network")
	@pipe get_degree.(binets1[1]) .|> mean |> quantile(_, [0.1, 0.5, 0.9]) |> (x -> x .* 10) |> println
	@pipe get_degree.(binets2[1]) .|> mean |> quantile(_, [0.1, 0.5, 0.9]) |> (x -> x .* 10) |> println

	get_summarised_R0_values(binets1, binets2)
end

function get_med_conf(ns::DataFrame)
	ns_qs = describe(ns,
		:median,
		(x -> quantile(x, 0.025)) => :q025,
		(x -> quantile(x, 0.975)) => :q975,
	)
	@transform!(ns_qs,
		:error_h = :q975 .- :median,
		:error_l = :median .- :q025
	);
	d_med = Dict(r[:variable] => r[:median] for r in eachrow(ns_qs))
	d_eH = Dict(r[:variable] => r[:error_h] for r in eachrow(ns_qs))
	d_eL = Dict(r[:variable] => r[:error_l] for r in eachrow(ns_qs))
	(d_med, d_eH, d_eL)
end

function layout_matrix(d1::Dict, d2::Dict)
	[
		d1[:adj_mH] d1[:adj_fH] d1[:adj_mL] d1[:adj_fL];
		d1[:ngm_mH] d1[:ngm_fH] d1[:ngm_mL] d1[:ngm_fL];
		d2[:adj_mH] d2[:adj_fH] d2[:adj_mL] d2[:adj_fL];
		d2[:ngm_mH] d2[:ngm_fH] d2[:ngm_mL] d2[:ngm_fL];
	]
end

function groupedbar_compare_prop_ngm_adj(ns1::DataFrame, ns2::DataFrame; version = 1)
	nw1 = (version==1) ? "Network 1" : (version == 2) ? "Network 3" : "Network 5"
	nw2 = (version==1) ? "Network 2" : (version == 2) ? "Network 4" : "Network 6"

	ns_qs1 = get_med_conf(ns1);
	ns_qs2 = get_med_conf(ns2);
	mat_qs = layout_matrix.(ns_qs1, ns_qs2)
	devon = palette(:devon, 8)
	lipari = palette(:lipari, 9)
	colors = repeat([devon[4], lipari[6], devon[6], lipari[8]][:, :] |> permutedims, inner = (4, 1))

	pl = plot(; ylabel = "proportion in eigenvector")
	groupedbar!(pl, mat_qs[1], color = colors, bar_position = :stack, legend = :outertopright,
		labels = ["Male (sexual acquired)" "Female (sexual acquired)" "Male (community acquired)" "Female (community acquired)"])

	plot!(pl, 1.0 .- cumsum(mat_qs[1], dims = 2)[:, 1:(end-1)], seriestype = :scatter, label = "", markershape = :none, color = :black, ms = 0,
		yerror = (mat_qs[2][:, 1:(end-1)], mat_qs[2][:, 1:(end-1)]),
	)
	xticks_labels = ["adjacency matrix", "NGM", "adjacency matrix", "NGM"]
	yticks_ = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
	plot!(pl, xticks = ([1, 2, 3, 4], xticks_labels), yticks = (yticks_, yticks_),
		xrotation = 0, size = (700, 400), bottom_margin = 7Plots.mm, left_margin = 5Plots.mm)
	annotate!(pl, 1.5, -0.12, (nw1, 10, :center))
	annotate!(pl, 3.5, -0.12, (nw2, 10, :center))
	plot!(pl)
	pl
end

function Graphs.union(g::SimpleWeightedGraph, h::SimpleWeightedGraph)::SimpleWeightedGraph
	adj = adjacency_matrix(g) .+ adjacency_matrix(h)
	SimpleWeightedGraph(adj)
end

function add_weights(g::SimpleGraph, w::Float64)::SimpleWeightedGraph
	adj = Graphs.LinAlg.adjacency_matrix(g) .* w
	SimpleWeightedGraph(adj)
end

"""
Args:
- sexual_weight: relative risk of sexual contacts to community contacts.
"""
function binets_to_weighted(binets::BinetTuple,
	sexual_weight::Float64; community_weight::Float64 = 1.0)::BinetWeightedTuple
	gh = add_weights.(binets[1], sexual_weight)
	gl = add_weights.(binets[2], sexual_weight)
	gc = add_weights.(binets[3], community_weight)
	# Duplicated should be removed.
	ghl = @pipe union.(binets[1], binets[2]) |> add_weights.(_, sexual_weight)
	binets = (gh, gl, gc, ghl)
end

function plot_for_compare(ns1, ns2, ns1w, ns2w; title = "")
	plot_adj_ngm_R0(ns1, ns2)
	pl_R01 = plot_adj_ngm_eigen_R0(ns1, ns2; lim = (5, 30));
	pl_rel1 = plot_rel_R0_sex(ns1, ns2; lim = (0, 90));

	plot_adj_ngm_R0(ns1w, ns2w; lim = (5, 30))
	pl_R01w = plot_adj_ngm_eigen_R0(ns1w, ns2w; lim = (5, 30))
	pl_rel1w = plot_rel_R0_sex(ns1w, ns2w; lim = (0, 90));

	plot(pl_R01, pl_rel1, pl_R01w, pl_rel1w, layout = (2, 2), size = (880, 800),
		bottom_margin = 5Plots.mm, left_margin = 5Plots.mm,
		title = title) |> display
end

function add_annotate!(pl::Plots.Plot, s)
    annotate!(pl,(-0.1,1.1),
            text(s, :black, :left, :center, 14) )
end
