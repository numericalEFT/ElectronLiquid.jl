using JLD2
using Measurements
# using Gaston
using PyCall
using PyPlot
using ElectronLiquid
using ElectronGas

const filenameZ = "data_Z.jld2"
const filenameK = "data_K.jld2"
const parafilename = "para_wn_1minus0.csv"
const Zrenorm = false     # turn off Z renormalization 
para = ParaMC(rs=1.0, beta=40.0, Fs=-0.0, order=3, mass2=1.0, isDynamic=false)

function zfactor(data, β)
    return @. (imag(data[2, 1]) - imag(data[1, 1])) / (2π / β)
end

function mu(data)
    return real(data[1, 1])
end

function derive_zmu(datatuple, isSave)
    df = CounterTerm.fromFile(parafilename)
    para, ngrid, kgrid, data = datatuple
    printstyled(UEG.short(para), color=:yellow)
    println()

    for p in sort([k for k in keys(data)])
        println("$p: μ = $(mu(data[p]))   z = $(zfactor(data[p], para.β))")
    end

    _mu = Dict()
    for (p, val) in data
        _mu[p] = mu(val)
    end
    _z = Dict()
    for (p, val) in data
        _z[p] = zfactor(val, para.β)
    end

    dzi, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _z)
    println("dz: ", dzi)

    z = []
    for i in eachindex(dzi)
        push!(z, 1.0 / (1.0 + sum(dzi[1:i])))
    end

    ############# save to csv  #################
    # println(df)
    for o in keys(data)
        # println(o)
        # global df
        paraid = UEG.paraid(para)
        df = CounterTerm.appendDict(df, paraid, Dict("order" => o, "μ" => _mu[o].val, "μ.err" => _mu[o].err, "Σw" => _z[o].val, "Σw.err" => _z[o].err); replace=true)
    end

    # println("new dataframe\n$df")
    isSave && CounterTerm.toFile(df, parafilename)

    return z
end

function loaddata(para, FileName=filename)
    key = UEG.short(para)
    f = jldopen(FileName, "r")
    # println(key)
    # println(keys(f))
    p, ngrid, kgrid, sigma = f[key]
    # println(sigma)
    order = p.order
    _partition = UEG.partition(para.order)
    rdata, idata = Dict(), Dict()
    for (ip, p) in enumerate(_partition)
        rdata[p] = real(sigma[p][:, :])
        idata[p] = imag(sigma[p][:, :])
    end
    return ngrid, kgrid, rdata, idata
end

function renormalize(para, sigma, Zrenorm)
    mu, sw = CounterTerm.getSigma(para, parafile=parafilename)
    ############ z renormalized  ##########################
    dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)
    # println(para.order)
    # println(dmu)
    sigma = CounterTerm.chemicalpotential_renormalization(para.order, sigma, dmu)
    return sigma
end

function process(para, Zrenorm, filename)
    ngrid, kgrid, rdata, idata = loaddata(para, filename)

    return renormalize(para, rdata, Zrenorm), renormalize(para, idata, Zrenorm), kgrid, ngrid
end


if abspath(PROGRAM_FILE) == @__FILE__
    # isSave = false
    # if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
    #     # the second parameter may be set to save the derived parameters
    #     isSave = true
    # end
    isSave = true

    fz = jldopen(filenameZ, "r")
    key = keys(fz)[1]
    if UEG.paraid(fz[key][1]) == UEG.paraid(para)
        zmu = derive_zmu(fz[key], isSave)
    end

    println("zfactor: ", zmu)

    rSw_k, iSw_k, kgrid, ngrid = process(para, Zrenorm, filenameK)

    # println(rSw_k)
    # println(iSw_k)
    # println(kgrid)

    kF_label = searchsortedfirst(kgrid, para.kF)

    rS_dk = []
    for val in rSw_k
        push!(rS_dk, @. (val[1, :] - val[1, kF_label]) / (kgrid - kgrid[kF_label]) * para.me / para.kF)
        # push!(rS_dk, @. (val[1, :] - val[1, kF_label]) / (kgrid.grid - kgrid.grid[kF_label]) * para.me / para.kF)
    end

    i_δk = 4
    dm = []
    for o in 1:para.order
        y = sum([v.val for v in rS_dk[j]] for j in 1:o)
        e = [v.err for v in rS_dk[o]]
        push!(dm, measurement((y[kF_label+i_δk] + y[kF_label-i_δk]) / 2, (e[kF_label+i_δk] + e[kF_label-i_δk]) / 2 + abs(y[kF_label+i_δk] - y[kF_label-i_δk]) / 2))
    end

    meff = @. 1.0 / (zmu * (1.0 - dm))
    println(meff)
end