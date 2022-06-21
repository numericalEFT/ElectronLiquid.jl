module CounterTerm
# using CSV
using DataFrames
using DelimitedFiles
# using PyCall
using Measurements
export mergeInteraction, chemicalpotential_renormalization, fromFile, toFile, appendDict, chemicalpotential
export muCT

const parafileName = joinpath(@__DIR__, "para.csv") # ROOT/common/para.csv

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
    # println("size: ", size(d[(1, 0)]))
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

"""
    function chemicalpotential(order, Σ, isfock)

    Derive the chemicalpotential shift and the chemical potential counterterm for each order from the self-energy

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
function chemicalpotential(order, Σ, isfock)
    @assert order <= 4 "Order $order hasn't been implemented!"

    δμ = Vector{Any}(undef, order)
    μ = Vector{Any}(undef, order)
    if order >= 1
        μ[1] = Σ[(1, 0)]
        if isfock
            δμ[1] = 0.0 #for the Fock-renormalized G scheme only
        else
            δμ[1] = -μ[1] #for the Fock-renormalized G scheme only
        end
    end
    if order >= 2
        μ[2] = Σ[(2, 0)] + δμ[1] * Σ[(1, 1)]
        δμ[2] = -μ[2]
    end
    if order >= 3
        # Σ3 = Σ30+Σ11*δμ2+Σ12*δμ1^2+Σ21*δμ1
        μ[3] = Σ[(3, 0)] + δμ[1] * Σ[(2, 1)] + δμ[1]^2 * Σ[(1, 2)] + δμ[2] * Σ[(1, 1)]
        δμ[3] = -μ[3]
    end
    if order >= 4
        # Σ4 = Σ40+Σ11*δμ3+Σ12*(2*δμ1*δμ2)+Σ13*δμ1^3+Σ21*δμ2+Σ22*δμ1^2+Σ31*δμ1
        μ[4] = Σ[(4, 0)] + δμ[1] * Σ[(3, 1)] + δμ[1]^2 * Σ[(2, 2)] + δμ[2] * Σ[(2, 1)] + (δμ[1])^3 * Σ[(1, 3)] + 2 * δμ[1] * δμ[2] * Σ[(1, 2)] + δμ[3] * Σ[(1, 1)]
        # μ[4] = _mu[(4, 0)]  + δμ[2] * _mu[(2, 1)] + δμ[3] * _mu[(1, 1)]
        δμ[4] = -μ[4]
    end

    println("Chemical Potential shift and counterterm:")
    for o in 1:order
        println("order $o:  μ = $(μ[o])  δμ = $(δμ[o])")
    end
    return μ, δμ
end

function fromFile(parafile=parafileName)
    try
        data, header = readdlm(parafile, ',', header=true)
        return DataFrame(data, vec(header))
    catch e
        println(e)
        println("Failed to load from $parafile. We will initialize the file instead")
        return nothing
    end
end

function toFile(df, parafile=parafileName)
    if isnothing(df)
        @warn "Empty dataframe $df, nothing to save"
        return
    end
    println("Save the parameters to the file $parafile")
    writedlm(parafile, Iterators.flatten(([names(df)], eachrow(df))), ',')
end

compareRow(row, _dict) = all(value isa AbstractFloat ? row[key] ≈ value : row[key] == value for (key, value) in _dict)

function appendDict(df::Union{Nothing,DataFrame}, paraid::Dict, data::Dict)
    # if isempty(filter(row -> compareRow(row, paraid), df)) == false
    #     #duplicated paraid, but may 
    #     if replace == false
    #         return df
    #     end
    #     @warn "Add new row with the existing paraid! \n $(paraid)"
    # end
    d = merge(paraid, data)
    if isnothing(df) || isempty(df) || nrow(df) == 0
        return DataFrame(d)
    else
        if isempty(filter(row -> compareRow(row, d), df))
            #duplicated paraid and data, then simply skip
            df = deepcopy(df)
            append!(df, d)
        end
        return df
    end
end

"""
    function muCT(df, paraid, order)

    Extract the chemical potential counterterm for the given parameters.
    If there are multiple counterterms of the same order, then the counterterm with the smallest errorbar will be returned

# Arguments
df     : DataFrame
paraid : Dictionary of parameter names and values
order  : the truncation order
"""
function muCT(df::DataFrame, paraid::Dict, order::Int)
    # println(df)
    df = filter(row -> compareRow(row, paraid), df)
    sort!(df, "δμ.err") #sort error from small to large
    @assert isempty(df) == false "no data exist for $paraid"
    # println(df)

    δμ = [filter(r -> r["order"] == o, df)[1, "δμ"] for o in 1:order-1]
    err = [filter(r -> r["order"] == o, df)[1, "δμ.err"] for o in 1:order-1]
    # err = [mu for mu in df[!, "δμ.err"]]
    return measurement.(δμ, err)
end

end