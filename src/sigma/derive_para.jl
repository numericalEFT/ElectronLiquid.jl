"""
By definition, the sigma renormalization is defined as
Σ1 = Σ11
Σ2 = Σ20+Σ11*δμ1
Σ3 = Σ30+Σ11*δμ2+Σ12*δμ1^2+Σ21*δμ1
Σ4 = Σ40+Σ11*δμ3+Σ12*(2*δμ1*δμ2)+Σ13*δμ1^3+Σ21*δμ2+Σ22*δμ1^2+Σ31*δμ1
"""

include("../common/para_builder.jl")
using .UEG
include("../common/counterterm.jl")
using .CounterTerm

using Measurements
# using DataFrames
using JLD2

const filename = "dataCT.jld2"

function zfactor(idata, para)
    return (idata[2] - idata[1]) / (2π / para.β)
end

function mu(rdata)
    return rdata[1]
end

function loaddata(para, avg, std)
    order = para.order
    _partition = UEG.partition(para.order)
    rdata, idata = Dict(), Dict()
    for (ip, p) in enumerate(_partition)
        rdata[p] = [measurement(real(avg[ip, wi]), real(std[ip, wi])) for wi in 1:length(avg[ip, :])]
        idata[p] = [measurement(imag(avg[ip, wi]), imag(std[ip, wi])) for wi in 1:length(avg[ip, :])]
    end
    return order, _partition, rdata, idata
end

function process(FileName, isSave)
    f = jldopen(FileName, "r")
    df = fromFile()
    for key in keys(f)
        println(key)
        value = f[key]
        para, avg, std = value[1], value[2], value[3]
        _order, _partition, rdata, idata = loaddata(para, avg, std)

        for p in sort([k for k in keys(rdata)])
            println("$p: μ = $(mu(rdata[p]))   z = $(zfactor(idata[p], para))")
        end

        _mu = Dict()
        for (p, val) in rdata
            _mu[p] = mu(val)
        end
        _z = Dict()
        for (p, val) in idata
            _z[p] = zfactor(val, para)
        end

        ############# save to csv  #################
        # println(df)
        for o in keys(rdata)
            # global df
            paraid = UEG.paraid(para)
            df = appendDict(df, paraid, Dict("order" => o, "μ" => _mu[o].val, "μ.err" => _mu[o].err, "Σw" => _z[o].val, "Σw.err" => _z[o].err); replace=true)
        end
    end
    println("new dataframe\n$df")
    isSave && toFile(df)
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

    # para =
    # p = ParaMC(rs=5.0, beta=100.0, Fs=-0.0, order=Order, mass2=1e-5)
    # _order, _partition, rdata, idata = loaddata(filename)
    # println("original data up to the Order = $_order")
    # for p in sort([k for k in keys(rdata)])
    #     println("$p: μ = $(mu(rdata[p]))   z = $(zfactor(idata[p]))")
    # end
    # # rdata = mergeInteraction(rdata)
    # # idata = mergeInteraction(idata)
    # # println("merged data: ")
    # # for p in sort([k for k in keys(rdata)])
    # #     println("$p: μ = $(mu(rdata[p]))   z = $(zfactor(idata[p]))")
    # # end

    # _mu = Dict()
    # for (p, val) in rdata
    #     _mu[p] = mu(val)
    # end
    # _z = Dict()
    # for (p, val) in idata
    #     _z[p] = zfactor(val)
    # end

    # ############# save to csv  #################
    # df = fromFile()
    # println(df)
    # for o in keys(rdata)
    #     global df
    #     df = appendDict(df, paraid, Dict("order" => o, "μ" => _mu[o].val, "μ.err" => _mu[o].err, "Σw" => _z[o].val, "Σw.err" => _z[o].err))
    # end
    # println(df)
    # isSave && toFile(df)

    # _μ, δμ = chemicalpotential(_order, _mu, isFock)
    # _z = chemicalpotential_renormalization(_order, _z, δμ)
    # println("z factor and counterterm")
    # for o in 1:_order
    #     println("order $o: z = $(1 / (1 + sum(_z[1:o])))  δz = $(_z[o])")
    # end

end