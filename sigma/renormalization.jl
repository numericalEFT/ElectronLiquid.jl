include("io.jl")

function chemicalpotential(Order, rdata)
    # _partition = sort([k for k in keys(rdata)])
    # println(_partition)
    _mu = Dict()
    for (p, val) in rdata
        _mu[p] = mu(val)
    end
    δμ = Vector{Any}(undef, Order)
    μ = Vector{Any}(undef, Order)
    if Order >= 1
        μ[1] = _mu[(1, 0)]
        if isFock
            δμ[1] = 0.0 #for the Fock-renormalized G scheme only
        else
            δμ[1] = -μ[1] #for the Fock-renormalized G scheme only
        end
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
        # Σ4 = Σ40+Σ11*δμ3+Σ12*(2*δμ1*δμ2)+Σ13*δμ1^3+Σ21*δμ2+Σ22*δμ1^2+Σ31*δμ1
        μ[4] = _mu[(4, 0)] + δμ[1] * _mu[(3, 1)] + δμ[1]^2 * _mu[(2, 2)] + δμ[2] * _mu[(2, 1)] + (δμ[1])^3 * _mu[(1, 3)] + 2 * δμ[1] * δμ[2] * _mu[(1, 2)] + δμ[3] * _mu[(1, 1)]
        # μ[4] = _mu[(4, 0)]  + δμ[2] * _mu[(2, 1)] + δμ[3] * _mu[(1, 1)]
        δμ[4] = -μ[4]
    end

    return μ, δμ
end

if abspath(PROGRAM_FILE) == @__FILE__

    @assert length(ARGS) >= 1 "One argument for the data file name is required!"
    filename = ARGS[1]
    order, _partition, rdata, idata = load(filename)
    println("original data up to the Order = $order")
    for p in sort([k for k in keys(rdata)])
        println("$p: μ = $(mu(rdata[p]))   z = $(zfactor(idata[p]))")
    end
    rdata = mergeInteraction(rdata)
    idata = mergeInteraction(idata)
    println("merged data: ")
    for p in sort([k for k in keys(rdata)])
        println("$p: μ = $(mu(rdata[p]))   z = $(zfactor(idata[p]))")
    end
    _μ, δμ = chemicalpotential(order, rdata)
    println(_μ)
    println(δμ)

    _z = Dict()
    for (p, val) in idata
        _z[p] = zfactor(val)
    end
    _z = chemicalpotential_renormalization(order, _z, δμ)
    println(_z)
    for o in 1:order
        println("order $o: z = ", 1 / (1 + sum(_z[1:o])))
    end
    # println(1 / (1 + _z[2]), ", ", 1 / (1 + _z[2] + _z[3]), ", ", 1 / (1 + _z[2] + _z[3] + _z[4]))
    # println(1 / (1 + _z[2]), ", ", 1 / (1 + _z[2] + _z[3]))
end