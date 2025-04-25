ENV["TMPDIR"] = "/tmp"
using LinearAlgebra
using Countries
using CSV
using DataFrames
using Distributions
using Plots
using RCall
using Suppressor
using StructTypes
using JSON3
using Optim
using FiniteDiff
using Printf
using Turing
using AdvancedHMC
using Random
using LogDensityProblemsAD
using LogDensityProblems
using FiniteDifferences
using Retry

@rimport socialmixr as smr
@rimport base as r
@rimport hhh4contacts as h4c

Scalar=Array{Float64,0}
Base.adjoint(x::Scalar)=x
Base.zero(x::Scalar)=fill(0.)
Base.zero(m::AbstractArray{<:AbstractArray})=zero.(m)
Base.:+(f::Float64, s::Scalar)=f.+s
Base.:+(s::Scalar,f::Float64)=f.+s
Base.:*(s::Scalar,f::Float64)=f.*s
Base.:*(f::Float64,s::Scalar)=f.*s

Base.convert(a::Type{Scalar},v::Float64)=fill(v)


PDict = Dict{Symbol,Scalar}

## structs
mutable struct Pyramid{A<:AbstractArray,S1, S2}
    cases::A
    casecategories::S1
    ageinterval::Vector{Int}
    graphlabel::S2
    misc::Dict{Symbol,Any}
end
Base.:+(p1::Pyramid,p2::Pyramid)=Pyramid(p1.cases.+p2.cases, p1.casecategories,p1.ageinterval,p1.graphlabel,p1.misc)
Base.length(p::Pyramid)=1
Base.iterate(p::Pyramid) = (p, nothing)
Base.iterate(p::Pyramid,nothing) = nothing
Pyramid(cases,casecategories,ageinterval::AbstractVector{<:Integer},graphlabel,n::Nothing)=Pyramid(cases,casecategories,convert(Vector{Int},ageinterval),graphlabel,Dict{Symbol,Any}())
Pyramid(cases,casecategories,ageinterval::AbstractVector{<:Integer},graphlabel)=Pyramid(cases,casecategories,convert(Vector{Int},ageinterval),graphlabel,Dict{Symbol,Any}())

mutable struct ContactMatrix
    matrix::Matrix{Matrix{Float64}}
    ageinterval::Vector{Int}#V1
    susceptibility::Vector{Vector{Scalar}}
    parameters::PDict
    addmat::Matrix{Matrix{Float64}}
    misc::Dict{Symbol,Any}
end
ContactMatrix(m::Vector,a,s,p,am,mi)=ContactMatrix(m[:,:],a,s,p,am,mi)
ContactMatrix(m,a,s,p,am::Vector,mi)=ContactMatrix(m,a,s,p,am[:,:],mi)
ContactMatrix(m::Vector,a,s,p,am::Vector,mi)=ContactMatrix(m[:,:],a,s,p,am[:,:],mi)
Base.length(cm::ContactMatrix)=1
Base.iterate(cm::ContactMatrix,p) = nothing
Base.iterate(cm::ContactMatrix) = (cm, nothing)
#Base.iterate(cm::ContactMatrix,p) = (p, nothing)
#Base.iterate(cm::ContactMatrix) = nothing

## functions
include("../src/utils.jl")
# pyramid
function aggregatecategories(p::Pyramid, categories = p.casecategories)
    ind=findall(Base.Fix2(in,categories),p.casecategories)
    Pyramid([reduce((.+),getindex(p.cases,ind))], Symbol(categories...), p.ageinterval, p.graphlabel)
end

function aggregateagegroups(p::Pyramid, newageinterval)
    if !(newageinterval ⊆ p.ageinterval) error("new age intervals are not subset of existing intervals") end
    ind = findfirst.((.==)(newageinterval),Ref(p.ageinterval))
    Pyramid(groupsum.(p.cases,Ref(ind)),p.casecategories,newageinterval,p.graphlabel)
end

function combine(p1::Pyramid,p2::Pyramid)
    if p1.ageinterval!=p2.ageinterval error("age intervals do not match") end
    Pyramid(vcat(p1.cases,p2.cases), vcat(p1.casecategories,p2.casecategories), p1.ageinterval, vcat(p1.graphlabel,p2.graphlabel))
end

function makeagegrouplabel(ageinterval)
    string.(ageinterval).*"–".* [string.(ageinterval[begin+1:end].-1);""]
end
makeagegrouplabel(p::Pyramid)=makeagegrouplabel(p.ageinterval)

# contact matrices
function partvax_cov(countrycode, year, cessation = 1980)
    unvax_age = year-cessation
    vaxageint = [0,unvax_age÷10*10,unvax_age+1,(unvax_age÷10+1)*10]
    popsize(vaxageint,countrycode=string(countrycode),year=year)[2:3]|>Base.Fix2(normalize,1)|>last
end
function contactmatrix(surveyname::Symbol, ageinterval::AbstractVector, countries = nothing, filter = nothing, sus_func = x->one(0.);year, refyear = 2013, refcountrycode="ZWE")
    partcov=partvax_cov(countries, year)
    demogchange = popsize(ageinterval, countrycode=countries, year = year)./popsize(ageinterval,countrycode=refcountrycode,year=refyear)
    if refcountrycode=="MAN" refcountrycode = "ZWE" end
    cmt = socialmixr_eig(surveyname, ageinterval, get_country(refcountrycode).name,filter;susceptibility=demogchange)
    ContactMatrix([cmt.matrix], ageinterval, convert(Vector{Vector{Union{Scalar,Float64}}},[map.(sus_func,ageinterval)]),PDict(),[zero(cmt.matrix)],Dict{Symbol,Any}(:issynthetic=>false,:partcov=>partcov))
end
function contactmatrix(syntheticdata::AbstractDict, ageinterval::AbstractVector, countrycode_s, sus_func = x->one(0.);year=2024,refyear = 2020)
    partcov=partvax_cov(countrycode_s, year)
    refcountry = countrycode_s
    if countrycode_s==:BDIC refcountry = :BDI end
    demogchange = popsize(ageinterval,countrycode=string(countrycode_s),year=year)./popsize(ageinterval,countrycode=string(refcountry),year=refyear)
    if countrycode_s==:BDIC countrycode_s = :BDI end
    cmt = synthetic_eig(syntheticdata, ageinterval, countrycode_s;year=year,susceptibility=demogchange)
    ContactMatrix([cmt.matrix], ageinterval, convert(Vector{Vector{Union{Scalar,Float64}}},[map.(sus_func,ageinterval)]),PDict(),[zero(cmt.matrix)],Dict{Symbol,Any}(:issynthetic=>true,:partcov=>partcov))
end

contactmatrix(surveyname::Symbol, p::Pyramid, countries = nothing, filter = nothing, sus_func = x->one(0.);year,refyear = 2013, refcountrycode="ZWE")= contactmatrix(surveyname, p.ageinterval, countries, filter, sus_func;year=year,refyear = refyear, refcountrycode=refcountrycode)

contactmatrix(syntheticdata::AbstractDict, p::Pyramid, countrycode_s, sus_func = x->one(0.);year=2024,refyear=2020)= contactmatrix(syntheticdata, p.ageinterval, countrycode_s, sus_func;year=year,refyear=refyear)

function Base.:+(c1::ContactMatrix, m::AbstractArray)
    if iszero(c1.addmat) ContactMatrix(c1.matrix,c1.ageinterval,c1.susceptibility,c1.parameters,m,c1.misc)
    else ContactMatrix(c1.matrix+m,c1.ageinterval,c1.susceptibility,c1.parameters,c1.addmat,c1.misc) end
end

setsusceptibility!(cm::ContactMatrix, sus_func::Union{Function,AbstractArray{Function}})=cm.susceptibility.= broadcast.(sus_func, [cm.ageinterval],Ref(Ref(cm)))#[sus_func.(cm.ageinterval,Ref(cm))]
setsusceptibility!(cm::ContactMatrix, susceptibility::AbstractArray)=cm.susceptibility.=[susceptibility]

Base.merge(d::Dict,n::Nothing)=d
Base.merge(n::Nothing,d::Dict)=d
Base.merge(n::Nothing,m::Nothing)=nothing
function Base.vcat(cm1::ContactMatrix, cm2::ContactMatrix)
    if cm1.ageinterval!=cm2.ageinterval error("age intervals do not match") end
    susceptibility = ifelse(cm1.susceptibility==cm2.susceptibility, cm1.susceptibility,[cm1.susceptibility; cm2.susceptibility])
    ContactMatrix([cm1.matrix;cm2.matrix],cm1.ageinterval,susceptibility,merge(cm1.parameters,cm2.parameters),[cm1.addmat;cm2.addmat],merge(cm1.misc,cm2.misc))
end
function Base.hcat(cm1::ContactMatrix, cm2::ContactMatrix)
    if cm1.ageinterval!=cm2.ageinterval error("age intervals do not match") end
    susceptibility = ifelse(cm1.susceptibility==cm2.susceptibility, cm1.susceptibility,[cm1.susceptibility  cm2.susceptibility])
    ContactMatrix([cm1.matrix cm2.matrix],cm1.ageinterval,susceptibility,merge(cm1.parameters,cm2.parameters),[cm1.addmat cm2.addmat],merge(cm1.misc,cm2.misc))
end

# contact matrix data
wpp2024=CSV.read("../data/wpp2024/populationdata.csv",DataFrame)[:,[:Location,:Iso3,:AgeStart,:Time,:Value]]
rename!(wpp2024,:Value=>:population,:AgeStart=>Symbol("lower.age.limit"))
function popsize(ageinterval; countrycode, year,data=wpp2024) 
fildata = filter(:Time=>==(year),filter(:Iso3=>==(countrycode),data))
rcopy(smr.pop_age(fildata, ageinterval.|>Int)).population|>Base.Fix2(normalize,1)
end


function socialmixr_eig(surveyname, ageinterval, countries = nothing, filter=nothing; susceptibility = 1, addmat = zeros(fill(size(ageinterval)[1],2)...))
    smixr = @suppress( rcopy(r.suppressWarnings(smr.contact_matrix(surveyname, countries = countries, var"age.limits" = ageinterval, filter = filter,symmetric=true,var"return.demography"=true,var"estimated.contact.age"="sample"))) )
    cmt = smixr[:matrix]
    if countries=="Zimbabwe" cmt./=2 end # as Zimbabwe contact matrix contain two days of contacts per participant
    cmt .+= addmat
    cmt = susceptibility'.* cmt
    ρ = eigvals(cmt')[end]|>abs
    ev = abs.(normalize(eigvecs(cmt')[:,end],1))
    (eigval = ρ, eigvec = ev, matrix = cmt, demography = smixr[:demography], ageinterval = ageinterval)
end


function synthetic_eig(contactdata, ageinterval, countrycodes::AbstractArray;susceptibility = 1,year=2024)
    out = synthetic_eig.(Ref(contactdata),Ref(ageinterval),countrycodes,Ref(susceptibility);year=year)
    (countrycode = getfield.(out,:countrycode), eigval = getfield.(out,:eigval), eigvec = getfield.(out,:eigvec), matrix = getfield.(out,:matrix))
end
function synthetic_eig(contactdata, ageinterval, countrycode::Symbol; susceptibility = 1,year=2024)
    cmt = contactdata[countrycode]
    # merge contact matrix using {hhh4contacts}
    @rput cmt
    R"library(hhh4contacts); rownames(cmt)<-0:15*5"
    pop = popsize((0:15).*5,countrycode=string(countrycode),year=year)
    @rput pop
    groupmap = (;(Symbol.([@sprintf("%02d",a) for a in ageinterval]).=>[string.(x) for x in collect.(range.(ageinterval,[ageinterval[2:end];76].-1,step=5))])...)
    @rput groupmap
    cmt = rcopy(R"aggregateC(cmt,groupmap,pop)") #should rearrange rows
    cmt = susceptibility'.* cmt
    ρ = eigvals(cmt')[end]|>abs
    ev = abs.(normalize(eigvecs(cmt')[:,end],1))
    (countrycode=countrycode,eigval = ρ, eigvec = ev, matrix = cmt)
end
function synthetic_eig(contactdata, ageinterval, countrycode::Nothing; susceptibility = 1,year=2024)
    synthetic_eig(contactdata, ageinterval, keys(contactdata)|>collect, susceptibility;year=year)
end

# next generation matrix
function ngm(cm::ContactMatrix, R0 = nothing)
    if haskey(cm.parameters,:s_partvax) cm.parameters[:s_partvax] .= 1-(1-cm.parameters[:s_vax][])*cm.misc[:partcov] end # set partvax
    cmt = broadcast.(*,cm.susceptibility', broadcast.(+,cm.matrix, cm.addmat))
    if !isnothing(R0) cmt./=dominanteigval(cm) end
    cmt|>transpose
end
dominanteigval(cm::ContactMatrix) = cm|>ngm|>collapseblockmat|>eigvals.|>abs|>maximum
function dominanteigvec(cm::ContactMatrix) # now assumes it's synthetic matrix if ageinterval and matrix doesn't match: to be improved in future
    colngm = ngm(cm)|>collapseblockmat
    dev = normalize(eigvecs(colngm)[:,end],1).|>abs # dominant eigenvector, normalised
    if length(cm.ageinterval) > size(first(cm.matrix),1) error("finer age intervals are specified than the contact data") end
    if length(cm.ageinterval) < size(ngm(cm)|>collapseblockmat,1)
        if !iszero(size(colngm,1) % length(cm.ageinterval)) error("contact matrix size not multiple of age interval") end
        folds = size(first(colngm),1) % length(cm.ageinterval)
        dev = splitsum(dev,4,[1,3])|>collapseblockmat|>vec # ad hoc for now
    end
    dev
end


# inference
function likelihood(p::Pyramid, cm::ContactMatrix)
    if p.ageinterval!=cm.ageinterval error("age intervals do not match") end
    if haskey(p.misc,:weight) w = p.misc[:weight] else w = 1 end
    logpdf(Multinomial(sum(sum.(p.cases)),dominanteigvec(cm)),vcat(p.cases...))*w
end

struct NLL
    p::Pyramid
    cm::ContactMatrix
    mutateparameters::Vector{Scalar}
    modifier! ::Function
    penalty::Float64
end

struct NLLs{P<:Pyramid,C<:ContactMatrix}
    p::Vector{P}
    cm::Vector{C}
    mutateparameters::Vector{Scalar}
    modifier! ::Function
    penalty::Float64
    sync::Vector{Symbol}
end
voidmodifier(x...)=0.
function (nll::NLL)(x)
        for (parameter, el) in zip(nll.mutateparameters,x) parameter.=el end
        modifier_ll = nll.modifier!(nll.p, nll.cm) * nll.penalty
        -likelihood(nll.p,nll.cm)-modifier_ll
end
function (nlls::NLLs)(x)
        for (parameter, el) in zip(nlls.mutateparameters,x) parameter.=el end
        modifier_ll = mean(nlls.modifier!.(nlls.p, nlls.cm)) * nlls.penalty
        if !(isempty(nlls.sync)) for i in 1:length(nlls.cm)-1 overwriteparameters!(nlls.cm[i+1], nlls.cm[i],nlls.sync) end end # assume cm is in growing order in terms of # parameters
        -sum(likelihood.(nlls.p,nlls.cm))-modifier_ll
end

function optimise!(mutateparameters::Vector{Scalar}, p::Pyramid, cm::ContactMatrix; sync=false,modifier!::Function = voidmodifier) # sync is only placeholder in this method
    nll=NLL(p,cm,mutateparameters,modifier!, 10000.)
    @time opt = optimize(nll, zeros(length(mutateparameters)).+1e-6,fill(500.,length(mutateparameters)),getindex.(mutateparameters),Fminbox(LBFGS()),Optim.Options(g_tol=1e-5, x_tol=1e-8))
    hess = FiniteDiff.finite_difference_hessian(nll,opt.minimizer)
    nll(opt.minimizer)
    (minimizer = opt.minimizer, minimum = opt.minimum,  hessian = hess ,result=opt,nll=nll)
end

function optimise!(mutateparameters::Vector{Scalar}, p::AbstractArray{<:Pyramid}, cm::AbstractArray{<:ContactMatrix}; sync = Symbol[], modifier!::Function = voidmodifier)
    nlls=NLLs(p,cm,mutateparameters,modifier!,10000.,sync)
    @time opt = optimize(nlls, zeros(length(mutateparameters)).+1e-6,fill(100.,length(mutateparameters)),getindex.(mutateparameters),Fminbox(LBFGS()), Optim.Options(g_tol=1e-5, x_tol=1e-8, time_limit=1800.))
    hess = FiniteDiff.finite_difference_hessian(nlls,opt.minimizer)
    nlls(opt.minimizer)#for (parameter, el) in zip(mutateparameters,opt.minimizer) parameter.=el end
    (minimizer = opt.minimizer, minimum = opt.minimum,  hessian = hess,result=opt,nll=nlls)
end

## optimise! using MCMC
function b_optimise!(mutateparameters::Vector{Scalar}, p::Pyramid, cm::ContactMatrix; sync=false,modifier!::Function = voidmodifier) # sync is only placeholder in this method
    opt = optimise!(mutateparameters, p, cm; sync=sync,modifier! = modifier!)
    init = opt.minimizer
    hess = opt.hessian
    nll=NLL(p,cm,mutateparameters,modifier!,10000.)

    if applicable(modifier!,nothing) # i.e. if modifier was not set
        println("method: importance sampling resampling")
        return runISR(nll, init, hess, 2000)
    else
        println("method: No-U-turn sapler")
        return runmcmc(nll,init, 2000, 500)
    end
end
function b_optimise!(mutateparameters::Vector{Scalar}, p::AbstractArray{<:Pyramid}, cm::AbstractArray{<:ContactMatrix}; sync = Symbol[], modifier!::Function = voidmodifier)
    opt = optimise!(mutateparameters, p, cm; sync = sync, modifier! = modifier!)
    init = opt.minimizer
    hess = opt.hessian
    nlls=NLLs(p,cm,mutateparameters,modifier!,10000.,sync)
    if  applicable(modifier!,nothing) # i.e. if modifier was not set
        println("method: importance sampling resampling")
        return runISR(nlls, init, hess, 2000)
    else
        println("method: No-U-turn sapler")
        return runmcmc(nlls,init, 2000, 500)
    end
end
function runISR(nll, init, hess, n_samples = 1000)
    props = rand(MvNormal(init, inv(hess)),n_samples*10)|>eachcol
    lpprops = logpdf.(Ref(MvNormal(init, inv(hess))), props)
    ld = .-nll.(props)
    lw = ld.-lpprops
    w = pweights(exp.(lw.-maximum(lw)))
    res = sample(props, w, n_samples, replace = true)
    chain = setinfo(Chains(res), (logdensity=ld,))
    med = quantile(chain,q=[0.5]).nt[2]
    nll(med) # update contact matrix via nll
    (med = med, ld = ld, chain=chain, nll=nll, ess = sum(w)^2/sum(w.^2))
end
function runmcmc(nll, init, n_samples = 1000, n_adapts = 500)
    @model function poistrick(x=0, nll_input = nll, len = length(nll.mutateparameters))
        lpar ~ filldist(Normal(0,5),8)
        lvw ~ filldist(Normal(0,2),len-8)
        transpar = sqrt.(2exp.(lpar)).+1
        x~Poisson(nll_input([transpar;exp.(lvw)]))
        Turing.@addlogprob! sum(lpar)+sum(lvw)
        return exp.([lpar;lvw])
    end
    lp = LogDensityFunction(poistrick())
    model = AdvancedHMC.LogDensityModel(
            LogDensityProblemsAD.ADgradient(Val(:FiniteDifferences), lp,;fdm=FiniteDifferences.forward_fdm(2,1))
    )
    # run MCMC. Retry if stopped due to hitting NaN
    seedint = isempty(serial[1]) ? 1 : pop!(serial[1])
    AHMCchain=nothing
    @repeat 20 try
        Random.seed!(seedint)
        AHMCchain = AbstractMCMC.sample(
            model,
            AdvancedHMC.NUTS(0.8),
            n_samples+n_adapts;
            n_adapts = n_adapts,
            initial_params = [log.((init[1:8].-1).^2 ./2);log.(init[9:end])].+0.01,
        )
    catch e
        @retry if typeof(e) == ArgumentError seedint+=1 end
    end
    push!(serial[2], seedint)
    ## post processing
    transθ = getfield.(getfield.(AHMCchain,:z),:θ)[n_adapts.+(1:n_samples)]
    len = transθ|>first|>length
    par, vw = broadcast.(exp,getindex.(transθ, Ref(1:8))), broadcast.(exp,getindex.(transθ, Ref(len-1:len)))

    # rescale par to meet boundary conditions
    na = convert(Vector{Float64},nll.cm.misc[:pop])./2
    denomweights=Ref(na[4:7])./([sum(na[4:7]),1].*nll.cm.misc[:bcond])
    pa = ((@view el[1:(8)÷2]) for el in par)
    qa = ((@view el[(8)÷2+1:end]) for el in par)
    for (p,q) in zip(pa,qa)
        p./= sum(denomweights[1].* p)
        q./= sum(denomweights[2].* q)
    end

    ld = getfield.(getfield.(getfield.(AHMCchain,:z),:ℓπ),:value)
    chain = setinfo(Chains(vcat.(par,vw)), (logdensity=ld,))
    med = quantile(chain,q=[0.5]).nt[2]
    LogDensityProblems.logdensity(lp, log.(med)) # update contact matrix via nll
    (med = med, ld = ld, chain=chain,nll=nll)
end




function estimateparameters!(cms, p::Union{Pyramid,AbstractArray{<:Pyramid}}, parameters;sync=false,modifier!::Function = voidmodifier, bayesian=false)
    res = Vector{Any}(undef, length(cms))
    if length(cms)>2
    Threads.@threads for i in 1:length(cms)
        (cmt, parms)= (zip(cms, parameters)|>collect)[i]
        firstel = typeof(cmt) <: ContactMatrix ? cmt : first(cmt) # if cmt is an array of ContactMatrix take the first
    opt = !bayesian ? optimise!(parms, p, cmt,sync=sync,modifier! = modifier!) : b_optimise!(parms, p, cmt,sync=sync,modifier! = modifier!)
    for el in cmt el.misc[:opt] = opt end
        res[i]=opt
        end
    else
        for i in 1:length(cms)
        (cmt, parms)= (zip(cms, parameters)|>collect)[i]
        firstel = typeof(cmt) <: ContactMatrix ? cmt : first(cmt) # if cmt is an array of ContactMatrix take the first
    opt = !bayesian ? optimise!(parms, p, cmt,sync=sync,modifier! = modifier!) : b_optimise!(parms, p, cmt,sync=sync,modifier! = modifier!)
    for el in cmt el.misc[:opt] = opt end
        res[i]=opt
        end
    end
    res
end

function overwriteparameters!(c1::ContactMatrix,c2::ContactMatrix, parnames = true) # all if true, vector of parameter names to only ovewrwrite specific parameters
    for key in intersect(keys(c1.parameters),keys(c2.parameters), parnames == true ? keys(c1.parameters) : parnames)
        c1.parameters[key].=c2.parameters[key]
    end
end

## Handling MCMC chain
chainof(cm::ContactMatrix)=cm.misc[:opt].chain
chainof(cm::AbstractVector{<:ContactMatrix})=cm[1].misc[:opt].chain
function MCMCiterate(f::Function, cm::ContactMatrix, chain::Chains=chainof(cm))
    if haskey(cm.misc,:opt) med=copy(cm.misc[:opt].med) end
    len=size(chain)[2]
    out = [begin
            parv=collect(parvec)
            if len>=10 parv[1:8].=sqrt.(2parv[1:8]).+1 end
            cm.misc[:opt].nll(parv) # update cm with MCMC slice
            f(cm)
    end for parvec in (chain|>Array|>eachrow)]
    
    if len>=10 med[1:8].=sqrt.(2med[1:8]).+1 end
    cm.misc[:opt].nll(med) # revert cm to median estimates
    out
end

function pwlikelihood(cm)
    p = cm.misc[:opt].nll.p
    ps = (Pyramid(cases,p.casecategories,p.ageinterval,p.graphlabel,p.misc) for cases in onehot_vov(p.cases))
    likelihood.(ps,Ref(cm))
end
    
function waic(cm::ContactMatrix)
    nll = cm.misc[:opt].nll
    casevec=vcat(nll.p.cases...)
    pwl = MCMCiterate(pwlikelihood,cm)
    mcmclen=length(pwl)
    pwlmat = reduce(hcat,pwl)'
    waic_k = [-2*(Turing.logsumexp(col) - log(mcmclen) - var(col)) for col in eachcol(pwlmat)]
    waic = sum(waic_k,Weights(casevec)) - 2*(Turing.logfactorial(sum(casevec))-sum(Turing.logfactorial.(casevec))) # 2nd term: factorials in mutlnomial pdf
    se_k = [-2√(var(LogNormal(0,std(col)))/mcmclen  + 2var(col)^2/(mcmclen-1)) for col in eachcol(pwlmat)]
    se_waic = √sum(se_k.^2,Weights(casevec)) #2sqrt(sum([var(col,corrected=false) for col in eachcol(pwlmat)],Weights(casevec)))
    (waic, se_waic)
end

## Data
#load contact survey data
R"zimbabwe_survey <- readRDS('../data/zimbabwe.rds')"

# load synthetic contact data
R"base::load('../data/synthetic/contact_all.rdata')"
R"base::load('../data/synthetic/contact_home.rdata')"
contact_home = rcopy(R"contact_home")
contact_all = rcopy(R"contact_all");

# load case data
StructTypes.StringType(::Type{Pyramid}) = StructTypes.Mutable()
function readJSON(tags::AbstractArray{<:Union{String,Symbol}}=String[]; dir =  "../data/pyramid/")
    if isempty(tags) tags = filter(Base.Fix2(endswith,".json"),readdir(dir)) end
    tags = first.(splitext.(string.(tags)))
    Dict(Symbol.(tags).=>[JSON3.read(read(dir*string(tag)*".json"), Pyramid) for tag in tags])
end
function readJSON(tag::Union{String,Symbol}; dir = "../data/pyramid/")
    JSON3.read(read(dir*string(tag)*".json"), Pyramid)
end
function readJSON!(p::Pyramid,tag::Union{String,Symbol}; dir = "../data/pyramid/")
    JSON3.read!(read(dir*string(tag)*".json"), p)
end