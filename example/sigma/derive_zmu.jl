using ElectronLiquid
using ElectronGas
using Measurements
using Printf
using JLD2

# rs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0]
rs = [4.0,]
mass2 = [1.0,]
Fs = [-0.0,]
beta = [25.0,]
order = [4,]
isDynamic = false

const filename = "data_Z.jld2"

# const parafilename = "para.csv"
const parafilename = "para_wn_1minus0.csv"

function zfactor(data, β)
    return @. (imag(data[2, 1]) - imag(data[1, 1])) / (2π / β)
end

function mu(data)
    return real(data[1, 1])
end

function process(para, datatuple, isSave)
    df = CounterTerm.fromFile(parafilename)
    ngrid, kgrid, data = datatuple
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
    for i in eachindex(dzi)
        println("z[$i]: ", 1.0 / (1.0 + sum(dzi[1:i])))
    end

    for i in eachindex(dmu)
        println("dmu[$i]: ", dmu[i])
    end

    ############# save to csv  #################
    # println(df)
    for P in keys(data)
        # println(P)
        # global df
        paraid = UEG.paraid(para)
        df = CounterTerm.appendDict(df, paraid, Dict("partition" => P, "μ" => _mu[P].val, "μ.err" => _mu[P].err, "Σw" => _z[P].val, "Σw.err" => _z[P].err); replace=true)
    end

    # println("new dataframe\n$df")
    isSave && CounterTerm.toFile(df, parafilename)
end

if abspath(PROGRAM_FILE) == @__FILE__

    # @assert length(ARGS) >= 1 "One argument for the data file name is required!"
    # filename = ARGS[1]
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    f = jldopen(filename, "r")

    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic)
        kF = para.kF
        for key in keys(f)
            loadpara = UEG.ParaMC(key)
            if UEG.paraid(loadpara) == UEG.paraid(para)
                process(para, f[key], isSave)
            end
        end
    end
end
