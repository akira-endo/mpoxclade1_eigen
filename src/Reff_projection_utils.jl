scalar(x) = (isa(x,AbstractArray) && ndims(x)==0) ? x[] : x
function build_susceptibility_vec(cm::ContactMatrix, yr::Int)
    ai = cm.ageinterval; n=length(ai)
    pinf  = scalar(cm.parameters[:s_infant])
    ppart = scalar(cm.parameters[:s_partvax])
    pvax  = scalar(cm.parameters[:s_vax])
    v = ones(n); v[1]=pinf
    last3=(n-2):n
    if   yr==2010       v[last3].=pvax
    elseif yr==2015     v[last3[1]]=ppart; v[last3[2:3]].=pvax
    elseif yr==2020     v[last3[2:3]].=pvax
    elseif yr==2024     v[last3[2]]=ppart; v[last3[3]]=pvax
    elseif yr==2030     v[last3[3]]=pvax
    else error("") end
    return v
end

function replace_community!(cm::ContactMatrix, mat::Matrix{Float64}; sexual::Bool=true)
    B=(cm.matrix); (r,c)=size(B)
    if !sexual
        @assert (r,c)==(1,1)
        B[1]=mat
    else
        if     (r,c)==(1,1)
            B[1]=mat
        elseif (r,c)==(1,4)
            B[1,2].=mat./2
        elseif (r,c)==(2,2)
            B[1,2].=mat./2; B[2,1].=mat./2
        else error("") end
    end
end
strip_sexual!(cm::ContactMatrix) = cm.addmat .= zero(cm.addmat)

# collapseblockmat from original dominant, for non-sex we'll use this:
function dominanteigval_safe(cm::ContactMatrix)
    full  = broadcast(+, cm.matrix, cm.addmat)
    small = collapseblockmat(full)
    sus   = scalar.(first(cm.susceptibility))
    @assert length(sus)==size(small,1)
    M     = sus' .* small
    maximum(abs.(eigvals(M')))
end

const FILTERS = [nothing, (phys_contact=1,), (cnt_home=1,), (phys_contact=1,cnt_home=1)]
function drc_matrices(ageint::Vector{Int}, yr::Int; filters=FILTERS)
    cms = @suppress [
      contactmatrix(:zimbabwe_survey, ageint, "COD", flt, (_)->one(0.0);
                    year=yr, refyear=2012, refcountrycode="MAN")
      for flt in filters
    ]
    mats = getfield.(cms, :matrix) .|> first
    (; all=mats[1], phys=mats[2], home=mats[3], physhome=mats[4])
end
function bdi_matrices(ageint::Vector{Int}, yr::Int; filters=FILTERS)
    cms = @suppress [
      contactmatrix(:zimbabwe_survey, ageint, "BDI", flt, (_)->one(0.0);
                    year=yr, refyear=2012, refcountrycode="MAN")
      for flt in filters
    ]
    mats = getfield.(cms, :matrix) .|> first
    (; all=mats[1], phys=mats[2], home=mats[3], physhome=mats[4])
end

# reff sampler
function compute_reff(chain::Vector{ContactMatrix}, mat_fn::Function;
                      sexual::Bool, weights=[0.21, 0.33, 0.29, 0.17]) # model weights
    W      = weights ./ sum(weights)
    years  = [2010,2015,2020,2024,2030]
    ageint = chain[1].ageinterval
    N      = length(chain[1].misc[:opt].chain)
    out    = Dict{Int,Vector{Float64}}()

    # choose eigen function
    eigfn = sexual ? dominanteigval : dominanteigval_safe

    for yr in years
        cms  = [deepcopy(chain[i]) for i in 1:4]
        mats = values(mat_fn(ageint, yr)) |> collect
        for i in 1:4
            replace_community!(cms[i], mats[i]; sexual=sexual)
            setsusceptibility!(cms[i], ones(length(cms[i].ageinterval)))
            sexual || strip_sexual!(cms[i])
        end

        samples = Vector{Float64}(undef, N)
        for s in 1:N
            acc = 0.0
            for i in 1:4
                cm = cms[i]
                θ  = vec(Array(cm.misc[:opt].chain)[s, :, 1])
                cm.misc[:opt].nll(θ)
                setsusceptibility!(cm, build_susceptibility_vec(cm, yr))
                acc += W[i] * eigfn(cm)
            end
            samples[s] = acc
        end
        out[yr] = samples
    end

    return out
end

function summarize(dict::Dict{Int,Vector{Float64}}; name="")
    println("\n--- $name ---")
    for yr in sort(collect(keys(dict)))
        q = quantile(dict[yr], [0.5,0.025,0.975])
        @printf("%4d : %.3f (95%% CrI %.3f–%.3f)\n",
                yr, q[1], q[2], q[3])
    end
end