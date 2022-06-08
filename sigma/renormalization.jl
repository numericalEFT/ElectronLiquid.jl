"""
By definition, the sigma renormalization is defined as
Σ1 = Σ11
Σ2 = Σ20+Σ11*δμ1
Σ3 = Σ30+Σ11*δμ2+Σ12*δμ1^2+Σ21*δμ1
Σ4 = Σ40+Σ11*δμ3+Σ12*(2*δμ1*δμ2)+Σ21*δμ2+Σ22*δμ1^2+Σ31*δμ1
"""

include("../common/parameter.jl")

using Measurements
using JLD2

const Order = 3
const FileName = "data.jld2"

partition = [(1, 0, 0),  # order 1
    (2, 0, 0), (1, 1, 0), (1, 0, 1),  #order 2
    (3, 0, 0), (2, 1, 0), (2, 0, 1), (1, 1, 1), (1, 2, 0), (1, 0, 2), #order 3
]

partition = [p for p in sort(partition) if p[1] + p[2] + p[3] <= Order]

println("Diagram set: ", partition)

function zfactor(idata)
    return (idata[2] - idata[1]) / (2π / β)
    # return imag(avg[2] - avg[1]) / (2π / β), (abs(imag(std[2])) + abs(imag(std[1]))) / (2π / β)
    # return imag(avg[1]) / (π / β), imag(std[1]) / (π / β)
end

function mu(rdata)
    return rdata[1]
end

"""
Merge interaction order and the main order
(tot, go, io) --> (tot+io, go)
"""
function mergeInteraction(data)
    res = Dict()
    for (p, val) in data
        println(p)
        mp = (p[1] + p[3], p[2])
        if haskey(res, mp)
            res[mp] += val
        else
            res[mp] = val
        end
    end
    return res
end

function load()
    f = jldopen(FileName, "r")
    avg, std = f["avg"], f["std"]
    rdata, idata = Dict(), Dict()
    for (ip, p) in enumerate(partition)
        rdata[p] = [measurement(real(avg[ip, wi]), real(std[ip, wi])) for wi in 1:length(avg[ip, :])]
        idata[p] = [measurement(imag(avg[ip, wi]), imag(std[ip, wi])) for wi in 1:length(avg[ip, :])]
    end
    return rdata, idata
end

function chemicalpotential(rdata)
    _partition = sort([k for k in keys(rdata)])
    # println(_partition)
    _mu = Dict()
    for (p, val) in rdata
        _mu[p] = mu(val)
    end
    δμ = Vector{Any}(undef, Order)
    μ = Vector{Any}(undef, Order)
    if Order >= 1
        μ[1] = _mu[(1, 0)]
        # δμ[1] = -μ[1] #for the Fock-renormalized G scheme only
        δμ[1] = 0.0 #for the Fock-renormalized G scheme only
    end
    if Order >= 2
        μ[2] = _mu[(2, 0)] + δμ[1] * _mu[(1, 1)]
        δμ[2] = -μ[2]
    end
    if Order >= 3
        # Σ3 = Σ30+Σ11*δμ2+Σ12*δμ1^2+Σ21*δμ1
        μ[3] = _mu[(3, 0)] + δμ[1] * _mu[(2, 1)] + δμ[1]^2 * _mu[(1, 2)] + δμ[2] * _mu[(1, 1)]
        δμ[3] = -μ[3]
    end
    if Order >= 4
        # Σ4 = Σ40+Σ11*δμ3+Σ12*(2*δμ1*δμ2)+Σ21*δμ2+Σ22*δμ1^2+Σ31*δμ1
        μ[4] = _mu[(4, 0)] + δμ[1] * _mu[(3, 1)] + δμ[1]^2 * _mu[(2, 2)] + δμ[2] * _mu[(2, 1)] + 2 * δμ[1] * δμ[2] * _mu[(1, 2)] + δμ[3] * _mu[(1, 1)]
        δμ[4] = -μ[4]
    end

    return μ, δμ
end

if abspath(PROGRAM_FILE) == @__FILE__
    rdata, idata = load()
    println("original data: ")
    for (p, val) in rdata
        println("$p: μ = $(mu(rdata[p]))   z = $(zfactor(idata[p]))")
    end
    rdata = mergeInteraction(rdata)
    idata = mergeInteraction(idata)
    println("merged data: ")
    for (p, val) in rdata
        println("$p: μ = $(mu(rdata[p]))   z = $(zfactor(idata[p]))")
    end
    _μ, δμ = chemicalpotential(rdata)
    println(_μ)
    println(δμ)
end