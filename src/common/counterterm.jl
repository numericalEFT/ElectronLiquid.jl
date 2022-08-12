module CounterTerm
# using CSV
using DataFrames
using DelimitedFiles
using TaylorSeries
using Printf

using ..UEG

# using PyCall
using ..Measurements
export mergeInteraction, fromFile, toFile, appendDict
# export muCT, zCT
export z_renormalization, chemicalpotential_renormalization, renormalization
export sigmaCT
export getSigma

const parafileName = joinpath(@__DIR__, "para.csv") # ROOT/common/para.csv

function partition(order::Int)
    # normal order, G order, W order
    par = [(1, 0, 0),  # order 1
        (2, 0, 0), (1, 1, 0), (1, 0, 1),  #order 2
        (3, 0, 0), (2, 1, 0), (2, 0, 1), (1, 1, 1), (1, 2, 0), (1, 0, 2), #order 3
        (4, 0, 0), (3, 1, 0), (3, 0, 1), (2, 1, 1), (2, 2, 0), (2, 0, 2), (1, 3, 0), (1, 0, 3), (1, 2, 1), (1, 1, 2) #order 4
    ]
    return sort([p for p in par if p[1] + p[2] + p[3] <= order])
end

"""
Merge interaction order and the main order
(normal_order, G_order, W_order) --> (normal+W_order, G_order)
"""
function mergeInteraction(data)
    if data isa Dict && all(x -> length(x) == 3, keys(data))
        res = Dict()
        for (p, val) in data
            @assert length(p) == 3
            # println(p)
            mp = (p[1] + p[3], p[2])
            if haskey(res, mp)
                res[mp] += val
            else
                res[mp] = val
            end
        end
        return res
    else # nothing to merge
        return data
    end
end

"""
    function renormalization(order, data, δμ, δz=nothing; nbody=1, zrenorm=true)
    
    First perform the chemical potential renormalization, then perform the z-factor renormalization

# Arguments
order : total order
data  : Dict{Order_Tuple, Actual_Data}, where Order_Tuple is a tuple of two integer Tuple{Normal_Order+W_Order, G_Order}
δμ    : chemical potential counterterm
δz    : z-factor counterterm

zrenorm : turn on or off the z-factor renormalization
nbody : nbody=1 for the one-body vertex function (self-energy, or Γ3) and nbody=2 for the two-body vertex function
"""
function renormalization(order, data, δμ, δz=nothing; nbody=1, zrenorm=true)
    data = chemicalpotential_renormalization(order, data, δμ)
    println(data)
    if zrenorm && (isnothing(δz) == false)
        data = z_renormalization(order, data, δz, nbody)
    end
    return data
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
    @assert length(δμ) + 1 >= order
    data = mergeInteraction(data)
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
        z[4] = d[(4, 0)] + δμ[1] * d[(3, 1)] + δμ[1]^2 * d[(2, 2)] + δμ[2] * d[(2, 1)] + (δμ[1])^3 * d[(1, 3)] + 2 * δμ[1] * δμ[2] * d[(1, 2)] + δμ[3] * d[(1, 1)]
        # z[4] = d[(4, 0)] + δμ[2] * d[(2, 1)] + δμ[3] * d[(1, 1)]
    end
    return z
end

"""
    function z_renormalization(order, data, δz, nbody::Int)

    By defintion, z=1+δz1+δz2+...

    Then the z-factor renormalization resuffles the power series in the following way:
    1. For the one-body vertex function
        DR = z*D ==> (D1+D2+D3+D4+...)*(1+δz1+δz2+δz3+...) = D1 + (D2+D1*δz1) + (D3+D2*δz1+D1*δz2) + (D4+D3*δz1+D2*δz2+D1*δz3) + ...
    2. For the two-body vertex function
        DR = z^2*D ==> (D1+D2+D3+D4+...)*(1+δz1+δz2+δz3+...)^2 = D1 + (D2+2*D1*δz1) + (D3+2*D2*δz1+D1*(δz1^2+2*δz2) + (D4+2*D3*δz1+D2*(δz1^2+δz2)+2*D1*(δz1*δz2+δz3)) + ...
    
    Note that the z-factor renormalization doesn't alter the W-order and the G-order, therefore, D_{m,n,k} renormalizes just as D_{m}.

# Arguments
order : total order
data  : Vector of data, or Dict{Order_Tuple, Actual_Data}, where Order_Tuple is a tuple of two integer Tuple{Normal_Order+W_Order, G_Order} or three integer Tuple{Normal_Order, G_Order, W_Order}
δz    : z-factor renormalization for each order
nbody : nbody=1 for the one-body vertex function (self-energy, or Γ3) and nbody=2 for the two-body vertex function
"""
function z_renormalization(order, data, δz, nbody::Int)
    # @assert order <= 2 "Order $order hasn't been implemented!"
    @assert order <= length(δz) + 1
    data = mergeInteraction(data)
    function addOneZ(order, data::AbstractDict, δz)
        nd = Dict()
        for orders in keys(data)
            if sum(orders) > order
                continue
            end
            nd[orders] = deepcopy(data[orders])
            maxO = orders[1] #the highest order
            others = orders[2:end]
            for o in 1:maxO-1
                lower = (maxO - o, others...)
                nd[orders] += data[lower] * δz[o]
            end
        end
        return nd
    end
    function addOneZ(order, data::AbstractVector, δz)
        # println(data)
        # println(length(data))
        # println(order)
        @assert order <= length(data)
        nd = deepcopy(data[1:order])
        for _order in 1:order
            maxO = _order #the highest order
            for o in 1:maxO-1
                nd[_order] += data[maxO-o] * δz[o]
            end
        end
        return nd
    end
    #nobody vertex function requires z^n factor
    for n in 1:nbody
        data = addOneZ(order, data, δz)
    end
    return data
end

"""
    function sigmaCT(order, μ, sw=Dict(key => 0.0 for key in keys(μ)); isfock=false)

    Derive the chemicalpotential and z-factor counterterm for each order from the self-energy.

# Arguments
order : total order
μ     : ReΣ(kF, w=0),     Dict{Order_Tuple, Actual_Data}, where Order_Tuple is a tuple of two integer Tuple{Normal_Order+W_Order, G_Order}, or three integer Tuple{Normal_Order, W_Order, G_Order}
sw    : dImΣ(kF, w=0)/dw, Dict{Order_Tuple, Actual_Data}, where Order_Tuple is a tuple of two integer Tuple{Normal_Order+W_Order, G_Order}, or three integer Tuple{Normal_Order, W_Order, G_Order}

# Return (δzi, δμ, δz)
The convention is the following:
δzi : 1/z = 1+δzi_1+δzi_2+... 
δμ  : chemical_potential_shift_without_z_renormalization = δμ_1+δμ_2+...
δz  : z = 1+δz_1+δz_2+...

# Remark:
The chemical potential shift is the chemical potential shift without z-renormalization.
"""
function sigmaCT(order, μ, sw=Dict(key => 0.0 for key in keys(μ)); isfock=false, verbose=0)
    swtype = typeof(collect(values(sw))[1])
    mutype = typeof(collect(values(μ))[1])
    δzi = zeros(swtype, order)
    δz = zeros(swtype, order)
    δμ = zeros(mutype, order)
    for o in 1:order
        # println("zR: ", zR)
        μR = mergeInteraction(μ)
        swR = mergeInteraction(sw)

        # println(swR)
        swR = chemicalpotential_renormalization(o, swR, δμ)
        μR = chemicalpotential_renormalization(o, μR, δμ)

        if isfock && o == 1
            δμ[o] = zero(swtype)
            δz[o] = zero(mutype)
            δzi[o] = zero(mutype)
        else
            δμ[o] = -μR[o]
            # δz[o] = swR[o] # somehow, the current scheme misses a factor of -1 in the z-factor counterterm
            δzi[o] = swR[o]
        end
    end

    zi = Taylor1([1.0, δzi...], order)
    z = 1 / zi
    δz = [getcoeff(z, o) for o in 1:order]


    if verbose > 0
        printstyled(@sprintf("%8s  %24s  %24s  %24s\n", "order", "δzi", "δμ", "δz"), color=:green)
        for o in 1:order
            @printf("%8d  %24s  %24s  %24s\n", o, "$(δzi[o])", "$(δμ[o])", "$(δz[o])")
        end
    end
    return δzi, δμ, δz
end


function fromFile(parafile=parafileName)
    try
        data, header = readdlm(parafile, ',', header=true)
        df = DataFrame(data, vec(header))
        sortdata!(df)
        return df
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
    sortdata!(df)
    println("Save the parameters to the file $parafile")
    writedlm(parafile, Iterators.flatten(([names(df)], eachrow(df))), ',')
end

compareRow(row, _dict) = all(value isa AbstractFloat ? row[key] ≈ value : row[key] == value for (key, value) in _dict)

sortdata!(df::DataFrame) = sort!(df, ["dim", "spin", "isDynamic", "isFock", "rs", "beta", "Fs", "Fa", "mass2"])

function compactOrder(orders)
    @assert length(orders) == 3
    o1, o2, o3 = orders
    @assert o1 < 10 && o2 < 10 && o3 < 10
    return o1 * 100 + o2 * 10 + o3
end

function appendDict(df::Union{Nothing,DataFrame}, paraid::Dict, data::Dict; replace=true)
    # if isempty(existing) == false
    #     #     #duplicated paraid, but may 
    #     #     if replace == false
    #     #         return df
    #     #     end
    #     #     @warn "Add new row with the existing paraid! \n $(paraid)"
    #     # sort!(existing, "μ.err") #sort error from small to large
    #     # sort!(existing, "Σw.err") #sort error from small to large
    #     # if (data["Σw.err"] < existing["Σw.err"][end]) && (data["μ.err"] < existing["μ.err"][end])
    #     #     delete!(existing, length(existing))
    #     # end
    #     for row in eachrow(existing)
    #         if (data["Σw.err"] < existing["Σw.err"][end]) && (data["μ.err"] < existing["μ.err"][end])
    #             delete!(existing, row)
    #         end
    #     end
    # end
    d = merge(paraid, data)
    if haskey(d, "order")
        if length(d["order"]) == 3
            d["order"] = compactOrder(d["order"])
        end
    end
    if isnothing(df) || isempty(df) || nrow(df) == 0
        return DataFrame(d)
    else
        # println("data\n$d")
        # remove the entries with larger error
        # df = filter(row -> (!(compareRow(row, paraid) && (row["order"] == d["order"]) && row["Σw.err"] > d["Σw.err"] && row["μ.err"] > d["μ.err"])), df)
        # if replace
        olddf = filter(row -> ((compareRow(row, paraid) && (row["order"] == d["order"]))), df)
        if isempty(olddf)
            # println("to save")
            append!(df, d)
            sortdata!(df)
        else
            bigerrdf = filter(row -> ((compareRow(row, paraid) && (row["order"] == d["order"]) && row["Σw.err"] > d["Σw.err"] && row["μ.err"] > d["μ.err"])), df)
            if isempty(bigerrdf) == false
                # replace only if the the new data has better quality for all quantitites
                # println("to replace with the new one")
                println(bigerrdf)
                df = filter(row -> (!(compareRow(row, paraid) && (row["order"] == d["order"]) && row["Σw.err"] > d["Σw.err"] && row["μ.err"] > d["μ.err"])), df)
                append!(df, d)
                sortdata!(df)
            end
        end
        return df
    end
end

"""
    function getSigma(para::ParaMC; order=para.order, parafile=parafileName)

    read self-energy parameters from the file

# Arguments
df     : DataFrame
paraid : Dictionary of parameter names and values
order  : the truncation order
"""
function getSigma(para::ParaMC; order=para.order, parafile=parafileName)
    df = fromFile(parafile)
    return getSigma(df, UEG.paraid(para), order)
end
function getSigma(df::DataFrame, paraid::Dict, order::Int)
    if order == 0
        return [], []
    end
    _partition = partition(order) # to construct the counterterms, we only need to calculate order-1 sigma
    # println(df)
    df = filter(row -> compareRow(row, paraid), df)
    sort!(df, "μ.err") #sort error from small to large
    @assert isempty(df) == false "no data exist for $paraid"
    # println(df)

    mu = Dict()
    for o in _partition
        v = filter(r -> r["order"] == compactOrder(o), df)[1, "μ"]
        err = filter(r -> r["order"] == compactOrder(o), df)[1, "μ.err"]
        mu[o] = measurement(v, err)
    end

    sort!(df, "Σw.err") #sort error from small to large
    @assert isempty(df) == false "no data exist for $paraid"

    sw = Dict()
    for o in _partition
        v = filter(r -> r["order"] == compactOrder(o), df)[1, "Σw"]
        err = filter(r -> r["order"] == compactOrder(o), df)[1, "Σw.err"]
        sw[o] = measurement(v, err)
    end

    return mu, sw
end

# """
#     function muCT(df, paraid, order)

#     Extract the chemical potential counterterm for the given parameters.
#     If there are multiple counterterms of the same order, then the counterterm with the smallest errorbar will be returned

# # Arguments
# df     : DataFrame
# paraid : Dictionary of parameter names and values
# order  : the truncation order
# """
# function muCT(df::DataFrame, paraid::Dict, order::Int)
#     if order == 1
#         return [], []
#     end
#     _partition = partition(order)
#     # println(df)
#     df = filter(row -> compareRow(row, paraid), df)
#     sort!(df, "δμ.err") #sort error from small to large
#     @assert isempty(df) == false "no data exist for $paraid"
#     # println(df)

#     δμ = [filter(r -> r["order"] == o, df)[1, "δμ"] for o in 1:order-1]
#     err = [filter(r -> r["order"] == o, df)[1, "δμ.err"] for o in 1:order-1]
#     # δz = [filter(r -> r["order"] == o, df)[1, ""] for o in 1:order-1]
#     # err = [filter(r -> r["order"] == o, df)[1, "δμ.err"] for o in 1:order-1]
#     # err = [mu for mu in df[!, "δμ.err"]]
#     return measurement.(δμ, err)
# end

# """
#     function zCT(df, paraid, order)

#     Extract the z-factor counterterm for the given parameters.
#     If there are multiple counterterms of the same order, then the counterterm with the smallest errorbar will be returned

# # Arguments
# df     : DataFrame
# paraid : Dictionary of parameter names and values
# order  : the truncation order
# """
# function zCT(df::DataFrame, paraid::Dict, order::Int)
#     if order == 1
#         return []
#     end
#     # println(df)
#     df = filter(row -> compareRow(row, paraid), df)
#     sort!(df, "δz.err") #sort error from small to large
#     @assert isempty(df) == false "no data exist for $paraid"
#     # println(df)

#     δz = [filter(r -> r["order"] == o, df)[1, "δz"] for o in 1:order-1]
#     err = [filter(r -> r["order"] == o, df)[1, "δz.err"] for o in 1:order-1]
#     return measurement.(δz, err)
# end

end