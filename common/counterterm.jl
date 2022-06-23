module CounterTerm
# using CSV
using DataFrames
using DelimitedFiles
# using PyCall
using Measurements
export mergeInteraction, chemicalpotential_renormalization, fromFile, toFile, appendDict, chemicalpotential
export muCT, zCT

const parafileName = joinpath(@__DIR__, "para.csv") # ROOT/common/para.csv

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
    @assert order <= 2 "Order $order hasn't been implemented!"
    @assert order <= length(δz) + 1
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
    function chemicalpotential(order, Σ, isfock)

    Derive the chemicalpotential shift and the chemical potential counterterm for each order from the self-energy

    By definition, the chemical potential renormalization is defined as
    Σ1 = Σ10
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
    Σ = mergeInteraction(Σ)

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

"""
    function chemicalpotential(order, Σ, isfock)

    Derive the chemicalpotential shift and the chemical potential counterterm for each order from the self-energy

    By definition, the chemical potential renormalization is defined as
    Σ1 = Σ10
    Σ2 = Σ20+Σ11*δμ1
    Σ3 = Σ30+Σ11*δμ2+Σ12*δμ1^2+Σ21*δμ1
    Σ4 = Σ40+Σ11*δμ3+Σ12*(2*δμ1*δμ2)+Σ13*δμ1^3+Σ21*δμ2+Σ22*δμ1^2+Σ31*δμ1

# Arguments
order : total order
data  : Dict{Order_Tuple, Actual_Data}, where Order_Tuple is a tuple of two integer Tuple{Normal_Order+W_Order, G_Order}
δμ    : chemical potential renormalization for each order
"""
function derive_onebody_parameter_from_sigma(order, μ, z=zeros(order); isfock=false)
    δz = zeros(order)
    δμ = zeros(order)
    μR = mergeInteraction(μ)
    zR = mergeInteraction(z)
    for o in 1:order
        # println(zR)
        zR = z_renormalization(o, zR, δz, 1)
        # println(zR)
        zR = chemicalpotential_renormalization(o, zR, δμ)

        μR = z_renormalization(o, μR, δz, 1)
        μR = chemicalpotential_renormalization(o, μR, δμ)
        if isfock && o == 1
            δμ[o] = 0.0
            δz[o] = 0.0
        else
            δμ[o] = -μR[o]
            δz[o] = zR[o]
        end
    end

    println("Onebody counterterm:")
    for o in 1:order
        println("order $o:  δμ = $(δμ[o]), δz = $(δz[o])")
    end
    return δμ, δz
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
    if order == 1
        return []
    end
    # println(df)
    df = filter(row -> compareRow(row, paraid), df)
    sort!(df, "δμ.err") #sort error from small to large
    @assert isempty(df) == false "no data exist for $paraid"
    # println(df)

    δμ = [filter(r -> r["order"] == o, df)[1, "δμ"] for o in 1:order-1]
    err = [filter(r -> r["order"] == o, df)[1, "δμ.err"] for o in 1:order-1]
    # δz = [filter(r -> r["order"] == o, df)[1, ""] for o in 1:order-1]
    # err = [filter(r -> r["order"] == o, df)[1, "δμ.err"] for o in 1:order-1]
    # err = [mu for mu in df[!, "δμ.err"]]
    return measurement.(δμ, err)
end

"""
    function zCT(df, paraid, order)

    Extract the z-factor counterterm for the given parameters.
    If there are multiple counterterms of the same order, then the counterterm with the smallest errorbar will be returned

# Arguments
df     : DataFrame
paraid : Dictionary of parameter names and values
order  : the truncation order
"""
function zCT(df::DataFrame, paraid::Dict, order::Int)
    if order == 1
        return []
    end
    # println(df)
    df = filter(row -> compareRow(row, paraid), df)
    sort!(df, "δz.err") #sort error from small to large
    @assert isempty(df) == false "no data exist for $paraid"
    # println(df)

    δz = [filter(r -> r["order"] == o, df)[1, "δz"] for o in 1:order-1]
    err = [filter(r -> r["order"] == o, df)[1, "δz.err"] for o in 1:order-1]
    return measurement.(δz, err)
end

end