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

# const parafileName = joinpath(@__DIR__, "para.csv") # ROOT/common/para.csv
const parafileName = "para.csv" # ROOT/common/para.csv

"""
Hard-coded counterterm partitions for the self-energy in the form (n_loop, n_μ, n_λ).
"""
function partition(order::Int)
    # normal order, G order, W order
    # NOTE: partitions of the form (0, nμ, nλ) vanish for Σ diagrams,
    #       since there is no interaction line at zeroth loop order
    par = [
        # order 1
        (1, 0, 0),
        # order 2
        (2, 0, 0), (1, 1, 0), (1, 0, 1),
        # order 3
        (3, 0, 0), (2, 1, 0), (2, 0, 1),
        (1, 1, 1), (1, 2, 0), (1, 0, 2),
        # order 4
        (4, 0, 0), (3, 1, 0), (3, 0, 1), (2, 1, 1), (2, 2, 0),
        (2, 0, 2), (1, 3, 0), (1, 0, 3), (1, 2, 1), (1, 1, 2),
        #order 5
        (5, 0, 0), (4, 1, 0), (4, 0, 1), (3, 2, 0), (3, 1, 1), (3, 0, 2), (2, 3, 0), (2, 2, 1),
        (2, 1, 2), (2, 0, 3), (1, 4, 0), (1, 3, 1), (1, 2, 2), (1, 1, 3), (1, 0, 4),
        #order 6
        (6, 0, 0), (5, 1, 0), (5, 0, 1), (4, 2, 0), (4, 1, 1), (4, 0, 2), (3, 3, 0), (3, 2, 1),
        (3, 1, 2), (3, 0, 3), (2, 4, 0), (2, 3, 1), (2, 2, 2), (2, 1, 3), (2, 0, 4), (1, 5, 0),
        (1, 4, 1), (1, 3, 2), (1, 2, 3), (1, 1, 4), (1, 0, 5),
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
- `order` : total order
- `data`  : Dict{Order_Tuple, Actual_Data}, where Order_Tuple is a tuple of two integer Tuple{Normal_Order+W_Order, G_Order}
- `δμ`    : chemical potential counterterm
- `δz`    : z-factor counterterm

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
- `order` : total order
- `data`  : Dict{Order_Tuple, Actual_Data}, where Order_Tuple is a tuple of two integer Tuple{Normal_Order+W_Order, G_Order}
- `δμ`    : chemical potential renormalization for each order
- `offset` (Int, optional): the first order (=normal+W_order) offset (defaults to 0).
"""
function chemicalpotential_renormalization(order, data, δμ; offset::Int=0)
    # _partition = sort([k for k in keys(rdata)])
    # println(_partition)
    @assert order <= 5 "Order $order hasn't been implemented!"
    @assert length(δμ) + 1 >= order
    data = mergeInteraction(data)
    d = data
    # println("size: ", size(d[(1, 0)]))
    # z = Vector{eltype(values(d))}(undef, order)
    sample = collect(values(d))[1]
    z = [zero(sample) for i in 1:order]
    # z = Vector{eltype(values(d))}[]
    # z = []
    # println(typeof(z))
    if order >= 1
        z[1] = d[(1 + offset, 0)]
    end
    if order >= 2
        z[2] = d[(2 + offset, 0)] + δμ[1] .* d[(1 + offset, 1)]
    end
    if order >= 3
        # Σ3 = Σ30+Σ11*δμ2+Σ12*δμ1^2+Σ21*δμ1
        z[3] = d[(3 + offset, 0)] + δμ[1] .* d[(2 + offset, 1)] + δμ[1] .^ 2 .* d[(1 + offset, 2)] + δμ[2] .* d[(1 + offset, 1)]
    end
    if order >= 4
        # Σ4 = Σ40+Σ11*δμ3+Σ12*(2*δμ1*δμ2)+Σ13*δμ1^3+Σ21*δμ2+Σ22*δμ1^2+Σ31*δμ1
        z[4] = d[(4 + offset, 0)] +
               d[(3 + offset, 1)] .* δμ[1] +
               d[(2 + offset, 2)] .* δμ[1] .^ 2 +
               d[(2 + offset, 1)] .* δμ[2] +
               d[(1 + offset, 3)] .* (δμ[1]) .^ 3 +
               d[(1 + offset, 2)] .* 2 .* δμ[1] .* δμ[2] +
               d[(1 + offset, 1)] .* δμ[3]
    end
    if order >= 5
        # Σ5 = Σ50 + Σ41*δμ1 + ...
        z[5] =
            d[(5 + offset, 0)] +
            d[(4 + offset, 1)] .* δμ[1] +
            d[(3 + offset, 2)] .* δμ[1] .^ 2 +
            d[(2 + offset, 3)] .* δμ[1] .^ 3 +
            d[(1 + offset, 4)] .* δμ[1] .^ 4 +
            d[(3 + offset, 1)] .* δμ[2] +
            d[(2 + offset, 2)] .* 2 .* δμ[1] .* δμ[2] +
            d[(1 + offset, 3)] .* 3 .* δμ[1] .^ 2 .* δμ[2] +
            d[(1 + offset, 2)] .* (δμ[2] .^ 2 + 2 * δμ[1] .* δμ[3]) +
            d[(2 + offset, 1)] .* δμ[3] +
            d[(1 + offset, 1)] .* δμ[4]
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
- `order` : total order
- `data`  : Vector of data, or Dict{Order_Tuple, Actual_Data}, where Order_Tuple is a tuple of two integer Tuple{Normal_Order+W_Order, G_Order} or three integer Tuple{Normal_Order, G_Order, W_Order}
- `δz`    : z-factor renormalization for each order
- `nbody` : nbody=1 for the one-body vertex function (self-energy, or Γ3) and nbody=2 for the two-body vertex function
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
                nd[orders] += data[lower] .* δz[o]
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
                nd[_order] += data[maxO-o] .* δz[o]
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
- `order` : total order
- `μ`     : ReΣ(kF, w=0),     Dict{Order_Tuple, Actual_Data}, where Order_Tuple is a tuple of two integer Tuple{Normal_Order+W_Order, G_Order}, or three integer Tuple{Normal_Order, W_Order, G_Order}
- `sw`    : dImΣ(kF, w=0)/dw, Dict{Order_Tuple, Actual_Data}, where Order_Tuple is a tuple of two integer Tuple{Normal_Order+W_Order, G_Order}, or three integer Tuple{Normal_Order, W_Order, G_Order}
- `isfock`: if true (false) Fock renormalization is turned on (off)
- `verbose`: verbosity level (0 (default): no output to stdout, 1: print to stdout)

# Return (δzi, δμ, δz)
The convention is the following:
- `δzi` : 1/z = 1+δzi_1+δzi_2+... 
- `δμ`  : chemical_potential_shift_without_z_renormalization = δμ_1+δμ_2+...
- `δz`  : z = 1+δz_1+δz_2+...

# Remark:
The chemical potential shift is the chemical potential shift without z-renormalization.
"""
function sigmaCT(order, μ, sw=Dict(key => 0.0 for key in keys(μ)); isfock=false, verbose=0)
    println(sw)
    # swtype = typeof(collect(values(sw))[1])
    # mutype = typeof(collect(values(μ))[1])
    sw1 = collect(values(sw))[1]
    mu1 = collect(values(μ))[1]

    δzi = [zero(sw1) for i in 1:order]
    δz = [zero(sw1) for i in 1:order]
    δμ = [zero(mu1) for i in 1:order]
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

    # _δzi = [i == 0 ? one(sw1) : δzi[i] for i in 0:order]
    # zi = Taylor1(_δzi, order)
    # z = 1 / zi
    # δz = [getcoeff(z, o) for o in 1:order]
    δz = _inverse(δzi)

    if verbose > 0
        printstyled(@sprintf("%8s  %24s  %24s  %24s\n", "order", "δzi", "δμ", "δz"), color=:green)
        for o in 1:order
            @printf("%8d  %24s  %24s  %24s\n", o, "$(δzi[o])", "$(δμ[o])", "$(δz[o])")
        end
    end
    return δzi, δμ, δz
end

function _inverse(z::AbstractVector{T}) where {T}
    order = length(z)
    zi = [zero(z[1]) for i in 1:length(z)]
    # zi = zeros(T, order)
    if order >= 1
        zi[1] = -z[1]
    end
    if order >= 2
        zi[2] = z[1] .^ 2 - z[2]
    end
    if order >= 3
        zi[3] = -z[1] .^ 3 + 2z[1] .* z[2] - z[3]
    end
    if order >= 4
        zi[4] = z[1] .^ 4 - 3z[1] .^ 2 .* z[2] + z[2] .^ 2 + 2z[1] .* z[3] - z[4]
    end
    if order >= 5
        zi[5] = -z[1] .^ 5 + 4z[1] .^ 3 .* z[2] - 3z[1] .* z[2] .^ 2 - 3z[1] .^ 2 .* z[3] + 2z[2] .* z[3] + 2z[1] .* z[4] - z[5]
    end
    if order >= 6
        error("order must be <= 5")
    end
    return zi
end


"""
    function densityCT(order, nw; isfock=false)
Derive the chemicalpotential shift from the density.
# Arguments
- `order` : total order
- `nw`     : ∫G(k,0⁻)dk,   Dict{Order_Tuple, Actual_Data}, where Order_Tuple is a tuple of two integer Tuple{Normal_Order+W_Order, G_Order}, or three integer Tuple{Normal_Order, W_Order, G_Order}
- `isfock`: if true (false) Fock renormalization is turned on (off)
- `verbose`: verbosity level (0 (default): no output to stdout, 1: print to stdout)
# Return δμ
The convention is the following:
- `δμ`  : chemical_potential_shift_with_renormalization_condition = δμ_1+δμ_2+...
"""
function densityCT(order, nw; isfock=false, verbose=0)
    nwtype = typeof(collect(values(nw))[1])
    δμ = zeros(nwtype, order)
    nR = mergeInteraction(nw)

    if isfock
        δμ[1] = zero(nwtype)
        # if order >= 3
        #     δμ[2] = -nR[(3,0)]/nR[(1,1)]
        # end
        # if order >= 4
        #     δμ[3] = -(nR[(4,0)] + nR[(2,1)] * δμ[2]) / nR[(1,1)]
        # end
        # if order >= 5
        #     δμ[4] = -(nR[(5,0)] + nR[(3,1)] * δμ[2] + nR[(2,1)] * δμ[3] + nR[(1,2)] * δμ[2]^2) / nR[(1,1)]
        # end
    else
        if order >= 1
            δμ[1] = -nR[(1, 0)] / nR[(0, 1)]
        end
    end
    if order >= 2
        δμ[2] = -(nR[(2, 0)] + nR[(1, 1)] * δμ[1] + nR[(0, 2)] * δμ[1]^2) / nR[(0, 1)]
    end
    if order >= 3
        δμ[3] = -(nR[(3, 0)] + nR[(2, 1)] * δμ[1] + nR[(1, 1)] * δμ[2] +
                  nR[(1, 2)] * δμ[1]^2 + nR[(0, 3)] * δμ[1]^3 + nR[(0, 2)] * 2 * δμ[2] * δμ[1]) / nR[(0, 1)]
    end
    if order >= 4
        δμ[4] = -(nR[(4, 0)] + nR[(3, 1)] * δμ[1] + nR[(2, 1)] * δμ[2] + +nR[(1, 1)] * δμ[3] +
                  nR[(2, 2)] * δμ[1]^2 + nR[(1, 3)] * δμ[1]^3 + nR[(1, 2)] * 2 * δμ[2] * δμ[1] +
                  nR[(0, 4)] * δμ[1]^4 + nR[(0, 3)] * 3 * δμ[1]^2 * δμ[2] + nR[(0, 2)] * (2 * δμ[1] * δμ[3] + δμ[2]^2)) / nR[(0, 1)]
    end
    if verbose > 0
        printstyled(@sprintf("%8s  %24s \n", "order", "δμ"), color=:green)
        for o in 1:order
            @printf("%8d  %24s \n", o, "$(δμ[o])")
        end
    end
    return δμ
end


"""
    function fromFile(parafile=parafileName; root_dir=@__DIR__)

Loads self-energy counterterm data from a CSV file `parafile` and returns a DataFrame.

# Arguments
- `parafile` : name of the CSV file to load from
- `root_dir` : the root directory of `parafile` (default: `<ElectronLiquid_Root>/common/`)
- `verbose`  : verbosity level (0: no output to stdout, 1 (default): print to stdout)
"""
function fromFile(parafile=parafileName; root_dir=@__DIR__, verbose=1)
    parafile = joinpath(root_dir, parafile) # ROOT/common/para.csv
    verbose > 0 && println("Reading para from $parafile")
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

"""
    function fromFile(parafile=parafileName; root_dir=@__DIR__)

Saves self-energy counterterm data specified by a DataFrame `df` to a CSV file `parafile`.

# Arguments
- `df`       : DataFrame containing the self-energy counterterm data
- `parafile` : name of the CSV file to save to
- `root_dir` : the root directory of `parafile` (default: `<ElectronLiquid_Root>/common/`)
- `verbose`  : verbosity level (0: no output to stdout, 1 (default): print to stdout)
"""
function toFile(df, parafile=parafileName; root_dir=@__DIR__, verbose=1)
    parafile = joinpath(root_dir, parafile) # ROOT/common/para.csv
    if isnothing(df)
        @warn "Empty dataframe $df, nothing to save"
        return
    end
    sortdata!(df)
    verbose > 0 && println("Save the parameters to the file $parafile")
    writedlm(parafile, Iterators.flatten(([names(df)], eachrow(df))), ',')
end

compareRow(row, _dict) = all(value isa AbstractFloat ? row[key] ≈ value : row[key] == value for (key, value) in _dict)

sortdata!(df::DataFrame) = sort!(df, ["order", "dim", "spin", "isDynamic", "isFock", "rs", "beta", "Fs", "Fa", "mass2"])

function compactPartition(P)
    @assert length(P) == 3
    P1, P2, P3 = P
    @assert P1 < 10 && P2 < 10 && P3 < 10
    return P1 * 100 + P2 * 10 + P3
end

function appendDict(df::Union{Nothing,DataFrame}, paraid::Dict, data::Dict; replace=true, verbose=1)
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
    if haskey(d, "partition")
        if length(d["partition"]) == 3
            d["partition"] = compactPartition(d["partition"])
        end
    end
    if isnothing(df) || isempty(df) || nrow(df) == 0
        return DataFrame(d)
    else
        # println("data\n$d")
        # remove the entries with larger error
        # df = filter(row -> (!(compareRow(row, paraid) && (row["partition"] == d["partition"]) && row["Σw.err"] > d["Σw.err"] && row["μ.err"] > d["μ.err"])), df)
        # if replace
        olddf = filter(row -> ((compareRow(row, paraid) && (row["partition"] == d["partition"]))), df)
        P = d["partition"]
        if isempty(olddf)
            verbose > 0 && println("will save $P for $paraid")
            append!(df, d)
            sortdata!(df)
        else
            bigerrdf = filter(row -> ((compareRow(row, paraid) && (row["partition"] == d["partition"]) && row["Σw.err"] >= d["Σw.err"] && row["μ.err"] >= d["μ.err"])), df)
            if isempty(bigerrdf) == false
                # replace only if the the new data has better quality for all quantitites
                verbose > 0 && println("will replace $P for $paraid")
                # println(bigerrdf)
                df = filter(row -> (!(compareRow(row, paraid) && (row["partition"] == d["partition"]) && row["Σw.err"] >= d["Σw.err"] && row["μ.err"] >= d["μ.err"])), df)
                append!(df, d)
                sortdata!(df)
            end
        end
        return df
    end
end

"""
    function getSigma(para::ParaMC; order=para.order, parafile=parafileName)

Derives the counterterms `mu` and `sw` from self-energy data stored in the CSV file `parafile`.

# Arguments
- `para`   : the physical parameter set to load data for
- `order`  : the maximum simulation order
- `parafile` : name of the CSV file to load from
- `root_dir` : the root directory of `parafile` (default: `<ElectronLiquid_Root>/common/`)
"""
function getSigma(para::ParaMC; order=para.order, parafile=parafileName, root_dir=@__DIR__)
    # println(parafile)
    df = fromFile(parafile; root_dir=root_dir)
    @assert isnothing(df) == false "file $parafile failed to load"
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
    for P in _partition
        v = filter(r -> r["partition"] == compactPartition(P), df)[1, "μ"]
        err = filter(r -> r["partition"] == compactPartition(P), df)[1, "μ.err"]
        mu[P] = measurement(v, err)
    end

    sort!(df, "Σw.err") #sort error from small to large
    @assert isempty(df) == false "no data exist for $paraid"

    sw = Dict()
    for P in _partition
        v = filter(r -> r["partition"] == compactPartition(P), df)[1, "Σw"]
        err = filter(r -> r["partition"] == compactPartition(P), df)[1, "Σw.err"]
        sw[P] = measurement(v, err)
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

#     δμ = [filter(r -> r["partition"] == o, df)[1, "δμ"] for o in 1:order-1]
#     err = [filter(r -> r["partition"] == o, df)[1, "δμ.err"] for o in 1:order-1]
#     # δz = [filter(r -> r["partition"] == o, df)[1, ""] for o in 1:order-1]
#     # err = [filter(r -> r["partition"] == o, df)[1, "δμ.err"] for o in 1:order-1]
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

#     δz = [filter(r -> r["partition"] == o, df)[1, "δz"] for o in 1:order-1]
#     err = [filter(r -> r["partition"] == o, df)[1, "δz.err"] for o in 1:order-1]
#     return measurement.(δz, err)
# end

end