module Green

using JLD2
using Printf, LinearAlgebra
using ..CompositeGrids
using ..ElectronGas
using ..MCIntegration
using ..Lehmann

using ..FeynmanDiagram
import ..FeynmanDiagram.FrontEnds: Filter, NoHartree
import ..FeynmanDiagram.Parquet: DiagPara, GreenDiag
using ..Measurements

using ..UEG
using ..Propagator
using ..Diagram

function diagram(paramc::ParaMC, _partition::Vector{T};
    filter=[NoHartree], extK=nothing, optimize_level=1
) where {T}
    return Diagram.diagram_parquet_noresponse(:green, paramc, _partition;
        filter=filter, extK=extK, optimize_level=optimize_level)
end

@inline function phase(varT, extT, l, β)
    tin, tout = varT[extT[1]], varT[extT[2]]
    return exp(1im * π * (2l + 1) / β * (tout - tin))
end

function default_partition(para::ParaMC)
    return [p for p in UEG.partition(para.order, offset=0) if !(p[1] == 0 && p[3] > 0)]
end

function _prepare_parquetad(para::ParaMC, diagram)
    partition, diagpara, FeynGraphs, extT_labels = diagram
    maxMomNum = maximum(p[1] for p in partition) + 1

    funcGraphs! = Dict{Int,Function}()
    leaf_maps = Vector{Dict{Int,Graph}}()
    for (i, key) in enumerate(partition)
        funcGraphs![i], leafmap = Compilers.compile(FeynGraphs[key])
        push!(leaf_maps, leafmap)
    end

    leafStat, loopbasis = FeynmanDiagram.leafstates(leaf_maps, maxMomNum)
    momLoopPool = FrontEnds.LoopPool(:K, para.dim, loopbasis)
    root = zeros(Float64, maximum(length.(extT_labels)))

    return maxMomNum, funcGraphs!, leafStat, leaf_maps, momLoopPool, root
end

function _eval_leafvalues!(idx, varK, varT, para::ParaMC, leafStat, leaf_maps, momLoopPool, isLayered2D::Bool)
    leafval, leafType, leafOrders, leafτ_i, leafτ_o, leafMomIdx = leafStat
    dim, β, me, μ = para.dim, para.β, para.me, para.μ
    tau_num = para.isDynamic ? 2 : 1

    for (i, lftype) in enumerate(leafType[idx])
        if lftype == 0
            continue
        elseif lftype == 1
            τ = varT[leafτ_o[idx][i]] - varT[leafτ_i[idx][i]]
            kq = FrontEnds.loop(momLoopPool, leafMomIdx[idx][i])
            ϵ = dot(kq, kq) / (2me) - μ
            leafval[idx][i] = Propagator.green_derive(τ, ϵ, β, leafOrders[idx][i][1])
        elseif lftype == 2
            diagid = leaf_maps[idx][i].properties
            τ1, τ2 = varT[leafτ_i[idx][i]], varT[leafτ_o[idx][i]]
            kq = FrontEnds.loop(momLoopPool, leafMomIdx[idx][i])
            leafval[idx][i] = Propagator.interaction_derive(
                τ1, τ2, kq, para, leafOrders[idx][i];
                idtype=diagid.type, tau_num=tau_num, isLayered=isLayered2D,
            )
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    return leafval[idx]
end

function _result_dict(partition, result, transform)
    datadict = Dict{eltype(partition),Any}()
    for (o, key) in enumerate(partition)
        avg, std = result.mean[o], result.stdev[o]
        datadict[key] = transform(avg, std)
    end
    return datadict
end

function _save_data(filename::Union{String,Nothing}, para::ParaMC, grid, kgrid, data)
    isnothing(filename) && return
    jldopen(filename, "a+") do f
        key = "$(UEG.short(para))"
        if haskey(f, key)
            @warn("replacing existing data for $key")
            delete!(f, key)
        end
        f[key] = (grid, kgrid, data)
    end
end

include("greenKT.jl")
include("greenKW.jl")

function MC(para::ParaMC;
    kgrid=[para.kF],
    tgrid=[para.β - 1e-8],
    neval=1e6,
    filename::Union{String,Nothing}=nothing,
    partition=default_partition(para),
    isLayered2D=false,
    filter=[NoHartree],
    extK=nothing,
    optimize_level=1,
    verbose=-1,
    kwargs...
)
    diagram_info = diagram(para, partition; filter=filter, extK=extK, optimize_level=optimize_level)
    data, result = KT(para, diagram_info;
        kgrid=kgrid, tgrid=tgrid, neval=neval, print=verbose,
        isLayered2D=isLayered2D, kwargs...)
    _save_data(filename, para, tgrid, kgrid, data)
    return data, result
end

function MC_KW(para::ParaMC;
    kgrid=[para.kF],
    ngrid=[0],
    neval=1e6,
    filename::Union{String,Nothing}=nothing,
    partition=default_partition(para),
    isLayered2D=false,
    filter=[NoHartree],
    extK=nothing,
    optimize_level=1,
    verbose=-1,
    kwargs...
)
    diagram_info = diagram(para, partition; filter=filter, extK=extK, optimize_level=optimize_level)
    data, result = KW(para, diagram_info;
        kgrid=kgrid, ngrid=ngrid, neval=neval, print=verbose,
        isLayered2D=isLayered2D, kwargs...)
    _save_data(filename, para, ngrid, kgrid, data)
    return data, result
end

end
