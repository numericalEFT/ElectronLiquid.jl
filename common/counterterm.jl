
"""
Merge interaction order and the main order
(normal_order, G_order, W_order) --> (normal+W_order, G_order)
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

"""
    function chemicalpotential_renormalization(order, data, δμ)
    
    merge different diagrammatic orders with proper chemical potential renormalization

    By definition, the chemical potential renormalization is defined as
    Σ1 = Σ11
    Σ2 = Σ20+Σ11*δμ1
    Σ3 = Σ30+Σ11*δμ2+Σ12*δμ1^2+Σ21*δμ1
    Σ4 = Σ40+Σ11*δμ3+Σ12*(2*δμ1*δμ2)+Σ13*δμ1^3+Σ21*δμ2+Σ22*δμ1^2+Σ31*δμ1

# Arguments
order : total order
data  : Dict{Order_Tuple, Actual_Data}, where Order_Tuple is a tuple of two integer Tuple{Normal_Order+W_Order, G_Order}
δμ    : chemical potential renormalization for each order
"""
function chemicalpotential_renormalization(order, data, δμ)
    # _partition = sort([k for k in keys(rdata)])
    # println(_partition)
    @assert order <= 4 "Order $order hasn't been implemented!"
    @assert length(δμ) <= order
    d = data
    z = Vector{eltype(values(d))}(undef, order)
    if order >= 1
        z[1] = d[(1, 0)]
    end
    if order >= 2
        z[2] = d[(2, 0)] + δμ[1] * d[(1, 1)]
    end
    if order >= 3
        # Σ3 = Σ30+Σ11*δμ2+Σ12*δμ1^2+Σ21*δμ1
        z[3] = d[(3, 0)] + δμ[1] * d[(2, 1)] + δμ[1]^2 * d[(1, 2)] + δμ[2] * d[(1, 1)]
    end
    if order >= 4
        # Σ4 = Σ40+Σ11*δμ3+Σ12*(2*δμ1*δμ2)+Σ13*δμ1^3+Σ21*δμ2+Σ22*δμ1^2+Σ31*δμ1
        # z[4] = _z[(4, 0)] + δμ[1] * _z[(3, 1)] + δμ[1]^2 * _z[(2, 2)] + δμ[2] * _z[(2, 1)]+ (δμ[1])^3 * _z[(1, 3)] + 2 * δμ[1] * δμ[2] * _z[(1, 2)] + δμ[3] * _z[(1, 1)]
        z[4] = d[(4, 0)] + δμ[2] * d[(2, 1)] + δμ[3] * d[(1, 1)]
    end
    return z
end

function readDeltaMu(_paraid=paraid, parafile=parafileName)
    println("read: ")
    df = DataFrame(CSV.File(parafileName))
    for (key, value) in _paraid
        # println("select: ", key, ", ", value, " type, ", typeof(value))
        # println(df)
        if value isa Real
            df = filter(row -> (row[key] ≈ value), df)
        else
            df = filter(row -> row[key] == value, df)
        end
    end
    println("Selected para: ", df)
    println(df["δμ"])
    println(df["δμ.err"])
end
