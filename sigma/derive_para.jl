"""
By definition, the sigma renormalization is defined as
Σ1 = Σ11
Σ2 = Σ20+Σ11*δμ1
Σ3 = Σ30+Σ11*δμ2+Σ12*δμ1^2+Σ21*δμ1
Σ4 = Σ40+Σ11*δμ3+Σ12*(2*δμ1*δμ2)+Σ13*δμ1^3+Σ21*δμ2+Σ22*δμ1^2+Σ31*δμ1
"""

include("../common/parameter.jl")
include("../common/counterterm.jl")

using .CounterTerm
using Measurements
# using DataFrames
using JLD2

function zfactor(idata)
    return (idata[2] - idata[1]) / (2π / β)
end

function mu(rdata)
    return rdata[1]
end

function loaddata(FileName)
    f = jldopen(FileName, "r")
    avg, std = f["avg"], f["std"]
    order = f["order"]
    _partition = f["partition"]
    rdata, idata = Dict(), Dict()
    for (ip, p) in enumerate(_partition)
        rdata[p] = [measurement(real(avg[ip, wi]), real(std[ip, wi])) for wi in 1:length(avg[ip, :])]
        idata[p] = [measurement(imag(avg[ip, wi]), imag(std[ip, wi])) for wi in 1:length(avg[ip, :])]
    end
    return order, _partition, rdata, idata
end

if abspath(PROGRAM_FILE) == @__FILE__

    @assert length(ARGS) >= 1 "One argument for the data file name is required!"
    filename = ARGS[1]
    _order, _partition, rdata, idata = loaddata(filename)
    println("original data up to the Order = $_order")
    for p in sort([k for k in keys(rdata)])
        println("$p: μ = $(mu(rdata[p]))   z = $(zfactor(idata[p]))")
    end
    rdata = mergeInteraction(rdata)
    idata = mergeInteraction(idata)
    println("merged data: ")
    for p in sort([k for k in keys(rdata)])
        println("$p: μ = $(mu(rdata[p]))   z = $(zfactor(idata[p]))")
    end

    _mu = Dict()
    for (p, val) in rdata
        _mu[p] = mu(val)
    end
    _μ, δμ = chemicalpotential(_order, _mu, isFock)

    _z = Dict()
    for (p, val) in idata
        _z[p] = zfactor(val)
    end
    _z = chemicalpotential_renormalization(_order, _z, δμ)
    println("z factor and counterterm")
    for o in 1:_order
        println("order $o: z = $(1 / (1 + sum(_z[1:o])))  δz = $(_z[o])")
    end

    ############# save to csv  #################
    df = fromFile()
    for o in 1:_order
        global df
        e = deepcopy(paraid)
        e["order"] = o
        e["δμ"] = δμ[o].val
        e["δμ.err"] = δμ[o].err
        e["δz"] = _z[o].val
        e["δz.err"] = _z[o].err

        appendDict!(df, e)
    end
    println(df)
    toFile(df)

    # println(1 / (1 + _z[2]), ", ", 1 / (1 + _z[2] + _z[3]), ", ", 1 / (1 + _z[2] + _z[3] + _z[4]))
    # println(1 / (1 + _z[2]), ", ", 1 / (1 + _z[2] + _z[3]))
end