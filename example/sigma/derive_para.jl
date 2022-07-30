using ElectronLiquid
using Measurements
using JLD2

const filename = "data_Z.jld2"

function zfactor(data, β)
    return @. (imag(data[2, 1]) - imag(data[1, 1])) / (2π / β)
end

function mu(data)
    return real(data[1, 1])
end

function process(FileName, isSave)
    f = jldopen(FileName, "r")
    df = CounterTerm.fromFile()
    for key in keys(f)
        println(key)
        value = f[key]
        para, ngrid, kgrid, data = value[1], value[2], value[3], value[4]

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

        ############# save to csv  #################
        # println(df)
        for o in keys(data)
            # global df
            paraid = UEG.paraid(para)
            df = CounterTerm.appendDict(df, paraid, Dict("order" => o, "μ" => _mu[o].val, "μ.err" => _mu[o].err, "Σw" => _z[o].val, "Σw.err" => _z[o].err); replace=true)
        end
    end
    println("new dataframe\n$df")
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
    process(filename, isSave)
end