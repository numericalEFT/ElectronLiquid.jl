using ElectronLiquid
using Measurements
using JLD2

rs = [5.0,]
mass2 = [0.01, ]
Fs = [-0.0, ]
beta = [100, ]
order = [2,]

const filename = "data_Z.jld2"

function zfactor(data, β)
    return @. (imag(data[2, 1]) - imag(data[1, 1])) / (2π / β)
end

function mu(data)
    return real(data[1, 1])
end

function process(datatuple, isSave)
    df = CounterTerm.fromFile()
    para, ngrid, kgrid, data = datatuple

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
    println("zfactor: ", dzi)

    ############# save to csv  #################
    # println(df)
    for o in keys(data)
        # global df
        paraid = UEG.paraid(para)
        df = CounterTerm.appendDict(df, paraid, Dict("order" => o, "μ" => _mu[o].val, "μ.err" => _mu[o].err, "Σw" => _z[o].val, "Σw.err" => _z[o].err); replace=true)
    end

    # println("new dataframe\n$df")
    isSave && CounterTerm.toFile(df)
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
        para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=true)
        kF = para.kF
        for key in keys(f)
            if UEG.paraid(f[key][1]) == UEG.paraid(para)
                process(f[key], isSave)
            end
        end
    end
end