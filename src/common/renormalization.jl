module Renorm
using DataFrames
using DelimitedFiles
using TaylorSeries

using ..Measurements
export mergeInteraction, fromFile, toFile, appendDict, chemicalpotential
export muCT, zCT
export z_renormalization, chemicalpotential_renormalization
export derive_onebody_parameter_from_sigma
export getSigma

function set_params(; order::Int, numct::Int)
    set_variables(:ϵ; order=order, numvars=numct + 1)
end

term(idx::Int) = get_variables()[idx]

function renormalize(data::TaylorN, idx::Int, series::TaylorN)
    newterm = [k == idx ? series : term(k) for k in 1:get_numvars()]
    println("newterm:", newterm)
    return evaluate(data, newterm)
end

function renormalize(data::AbstractArray, idx::Int, series::TaylorN)
    @assert eltype(data) <: TaylorN

    result = similar(data)
    newterm = [k == idx ? series : term(k) for k in 1:get_numvars()]
    for i in eachindex(data)
        result[i] = evaluate(data[i], newterm)
    end
    return result
end

function renormalize(data::AbstractArray, idx::Int, series::AbstractArray)
    @assert eltype(data) <: TaylorN
    @assert eltype(series) <: TaylorN

    @assert size(data) == size(series)
    result = similar(data)
    for i in eachindex(data)
        newterm = [k == idx ? series[i] : term(idx) for k in 1:get_numvars()]
        result[i] = evaluate(data[i], newterm)
    end
    return result
end

# function evaluate(series::TaylorN, x0=ones(get_numvars()))
#     return TaylorSeries.evaluate(series, x0)
# end

# function evaluate(series::AbstractArray, x0=ones(get_numvars()))
#     @assert eltype(series) <: TaylorN
#     result = similar(series)
#     for i in eachindex(series)
#         result[i] = TaylorSeries.evaluate(series[i], x0)
#     end
#     return result
# end

function merge(series, idx::Int; target::Int=1)
    return renormalize(series, idx, term(target))
end

function dict2series(datadict::AbstractDict)
    ks = collect(keys(datadict))
    d1 = datadict[ks[1]]
    @assert all(x -> size(x) == size(d1), values(datadict))
    @assert all(x -> eltype(x) == eltype(d1), values(datadict))


    if d1 isa Number
        type = typeof(d1)
        data = zero(TaylorN{type})
    elseif d1 isa AbstractArray
        type = eltype(d1)
        data = similar(d1, TaylorN{type})
        data .= 0.0
    else
        error("data with $(typeof(d1)) has unsupported type")
    end

    order = get_order()
    numvars = get_numvars()

    for (key, val) in datadict
        @assert length(key) == numvars "a tuple of $numvars integers are expected for key of the data"
        @assert sum(key) <= order "sum of orders should be smaller than the highest order $order"
        power = reduce(*, [term(i)^k for (i, k) in enumerate(key)])

        if val isa Number
            data += val * power
        else
            data .+= val * power
        end
    end
    return data
end

end