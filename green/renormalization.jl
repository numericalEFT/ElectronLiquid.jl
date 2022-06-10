"""
By definition, the sigma renormalization is defined as
Σ1 = Σ11
Σ2 = Σ20+Σ11*δμ1
Σ3 = Σ30+Σ11*δμ2+Σ12*δμ1^2+Σ21*δμ1
Σ4 = Σ40+Σ11*δμ3+Σ12*(2*δμ1*δμ2)+Σ13*δμ1^3+Σ21*δμ2+Σ22*δμ1^2+Σ31*δμ1
"""

include("../common/parameter.jl")

using Measurements
using JLD2

const Order = 4
const FileName = "data.jld2"

partition = [(0, 0, 0), (0, 1, 0), (0, 2, 0),
    (1, 0, 0),  # order 1
    (2, 0, 0), (1, 1, 0), (1, 0, 1),  #order 2
    (3, 0, 0), (2, 1, 0), (2, 0, 1), (1, 1, 1), (1, 2, 0), (1, 0, 2), #order 3
    (4, 0, 0), (3, 1, 0), (3, 0, 1), (2, 1, 1), (2, 2, 0), (2, 0, 2), (1, 3, 0), (1, 0, 3), (1, 2, 1), (1, 1, 2) #order 4
]

partition = [p for p in sort(partition) if p[1] + p[2] + p[3] <= Order]

println("Diagram set: ", partition)

"""
Merge interaction order and the main order
(tot, go, io) --> (tot+io, go)
"""
function mergeInteraction(data)
    res = Dict()
    for (p, val) in data
        # println(p)
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
    data = Dict()
    for (ip, p) in enumerate(partition)
        data[p] = measurement(avg[ip], std[ip])
    end
    return data
end

function chemicalpotential(rdata)
    _partition = sort([k for k in keys(rdata)])
    # println(_partition)
    _mu = Dict()
    for (p, val) in rdata
        _mu[p] = val
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
        δμ[2] = -μ[2]/_mu[(0,1)]
    end
    if Order >= 3
        # Σ3 = Σ30+Σ11*δμ2+Σ12*δμ1^2+Σ21*δμ1
        μ[3] = _mu[(3, 0)] + δμ[1] * _mu[(2, 1)] + δμ[1]^2 * _mu[(1, 2)] + δμ[2] * _mu[(1, 1)]
        println(_mu[(3, 0)], ", ", _mu[(1, 1)], ", ", δμ[2])
        δμ[3] = -μ[3]/_mu[(0,1)]
    end
    if Order >= 4
        # Σ4 = Σ40+Σ11*δμ3+Σ12*(2*δμ1*δμ2)+Σ13*δμ1^3+Σ21*δμ2+Σ22*δμ1^2+Σ31*δμ1
        # μ[4] = _mu[(4, 0)] + δμ[1] * _mu[(3, 1)] + δμ[1]^2 * _mu[(2, 2)] + δμ[2] * _mu[(2, 1)]+ (δμ[1])^3 * _mu[(1, 3)] + 2 * δμ[1] * δμ[2] * _mu[(1, 2)] + δμ[3] * _mu[(1, 1)]+ δμ[2]^2 * _mu[(0, 2)]
        # μ[4] = _mu[(4, 0)]  + δμ[2] * _mu[(2, 1)] + δμ[3] * _mu[(1, 1)]+ δμ[2]^2 * _mu[(0, 2)]
        μ[4] = _mu[(4, 0)]  + δμ[2] * _mu[(2, 1)] + δμ[3] * _mu[(1, 1)]- δμ[2]^2 * _mu[(0, 2)]
        println((_mu[(4, 0)]  - δμ[2]^2 * _mu[(0, 2)])/_mu[(0, 1)])
        δμ[4] = -μ[4]/_mu[(0,1)]
    end

    return μ, δμ
end

if abspath(PROGRAM_FILE) == @__FILE__
    data = load()
    println("original data: ")
    for p in sort([k for k in keys(data)])
        println("$p: density = $(data[p])")
    end
    data = mergeInteraction(data)
    println("merged data: ")
    for p in sort([k for k in keys(data)])
        println("$p: density = $(data[p])")
    end
    _μ, δμ = chemicalpotential(data)
    println(_μ)
    println(δμ)
end