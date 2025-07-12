include("../src/eigen.jl")
@suppress include("../src/eigen_setup.jl")
include("../src/eigen_output.jl")
tshuapaplot = plot(tshuapa_h2hag,color=:black)
endemicplot = plot(plot(drc_endemic_ag),ylim=(0,0.35),xtickfontsize=9);
Random.seed!(1)
endemic2015_24_fit = output_fit(
    [tshuapa_h2hag,drc_endemic_ag];
    zmb_skeleton = [zmb2015,zmb2024],
    drc_skeleton = [drc2015,drc2024],
    dataplots = [tshuapaplot,endemicplot],
    bayesian=false);
zmb2015_24_fit = endemic2015_24_fit.zmb_fit
drc2015_24_fit = endemic2015_24_fit.drc_fit;

drc_kamituga = aggregateagegroups(readJSON(:drc_kamituga_Aug2024),[0:5:20;30:10:50])
drc_otherhz = aggregateagegroups(readJSON(:drc_otherhz_Aug2024),[0:5:15;20:10:50])

burundi = aggregateagegroups(readJSON(:burundi_Oct2024),[0:5:15;20:10:50])
burundi0 = aggregateagegroups(readJSON(:burundi_midSep2024),[0:5:15;20:10:50])
burundi.cases .-= broadcast.(min,burundi0.cases,burundi.cases)

otherhzplot = plot(plot(drc_otherhz,color = [1 2]),ylim=(0,0.5),xtickfontsize=9)
burundiplot = plot(plot(drc_otherhz,color = [1 2]),ylim=(0,0.5),xtickfontsize=9)
kamitugaplot = plot(plot(drc_kamituga,color=[1 2]),ylim=(0,0.5),size=(600,300))

drc_kivu = aggregateagegroups(readJSON(:drc_kivu_Aug2024),[0:5:15;20:10:50])
kivuplot = plot(plot(drc_kivu),ylim=(0,0.5),xtickfontsize=9);

zmb2024_sexual = addsexualcontact!(deepcopy(zmb2024),[15;20:10:40]; modifier! = propmix!, countrycode = "COD",year = 2024);
drc2024_sexual = addsexualcontact!(deepcopy(drc2024),[15;20:10:40]; modifier! = propmix!, countrycode = "COD",year = 2024);
zmb2024_sexual_b = addsexualcontact!(deepcopy(bdi2024),[15;20:10:40]; modifier! = propmix!, countrycode = "BDIC",year = 2024);
drc2024_sexual_b = addsexualcontact!(deepcopy(bdi_s2024),[15;20:10:40]; modifier! = propmix!, countrycode = "BDIC",year = 2024);

zmb_ref = deepcopy(last.(zmb2015_24_fit|>collect))
drc_ref = deepcopy(last.(drc2015_24_fit|>collect))

function propmixsensitivity!(p::Pyramid,cm::ContactMatrix,assort::AbstractVector = assort)
    @assert all(0 .≤assort.≤1) "assort must be 0-1" # assort = [MF, FM]
    ll=propmix!(p,cm)
    for rowid in 1:size(cm.addmat,1)
        for mat in cm.addmat[rowid,:]
            aid=rowid<3 ? 1 : 2
            rowsum = sum(mat,dims=2)
            mat.*=(1-assort[aid])
            mat.+=diagm((rowsum.*assort[aid])|>vec)
        end
    end
    ll
end

function propmixseparate!(p::Pyramid,cm::ContactMatrix,assort::AbstractVector = assort)
    @assert all(0 .≤assort.≤1) "assort must be 0-1" # assort = [MF, FM]
    ll=propmix!(p,cm)
    for rowid in [1,3]
        for mat in cm.addmat[rowid,[2,4]]
            cmat = cm.matrix[1,2]
            mat.=.-(1-assort[1]).*cmat
        end
    end
    ll
end

function b_optimise!(mutateparameters::Vector{Scalar}, p::Pyramid, cm::ContactMatrix; sync=false,modifier!::Function = voidmodifier) # sync is only placeholder in this method
    opt = optimise!(mutateparameters, p, cm; sync=sync,modifier! = modifier!)
    init = opt.minimizer
    hess = opt.hessian
    nll=NLL(p,cm,mutateparameters,modifier!,10000.)

    if applicable(modifier!,nothing) # i.e. if modifier was not set
        println("method: importance sampling resampling")
        return runISR(nll, init, inv(hess), 2000)
    else
        println("method: No-U-turn sapler")
        return runmcmc(nll,init, 2000, 500)
    end
end
function estimateparameters!(cms, p::Union{Pyramid,AbstractArray{<:Pyramid}}, parameters;sync=false,modifier!::Function = voidmodifier, bayesian=false)
    res = Vector{Any}(undef, length(cms))
    if length(cms)>2
    for i in 1:length(cms)
        (cmt, parms)= (zip(cms, parameters)|>collect)[i]
        firstel = typeof(cmt) <: ContactMatrix ? cmt : first(cmt) # if cmt is an array of ContactMatrix take the first
    opt = !bayesian ? optimise!(parms, p, cmt,sync=sync,modifier! = modifier!) : b_optimise!(parms, p, cmt,sync=sync,modifier! = modifier!)
    for el in cmt el.misc[:opt] = opt end
        res[i]=opt
        end
    else
        bayesian=false
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
        
    