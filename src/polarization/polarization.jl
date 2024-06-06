module Polarization
using JLD2, CSV

using Printf, LinearAlgebra, DataFrames
using ..StaticArrays
using ..Parameters
using ..CompositeGrids
using ..ElectronGas
using ..MCIntegration
using ..Lehmann

using ..FeynmanDiagram
import ..FeynmanDiagram.FrontEnds: TwoBodyChannel, Alli, PHr, PHEr, PPr, AnyChan
import ..FeynmanDiagram.FrontEnds: Filter, NoHartree, NoFock, DirectOnly, Wirreducible, Girreducible, NoBubble, Proper
import ..FeynmanDiagram.FrontEnds: Response, Composite, ChargeCharge, SpinSpin, UpUp, UpDown
import ..FeynmanDiagram.FrontEnds: AnalyticProperty, Instant, Dynamic
import ..FeynmanDiagram.Parquet: DiagPara, Ver4Diag, PolarDiag
using ..Measurements

using ..UEG
using ..Propagator
import ..Propagator: LeafStateAD
using ..Diagram

import ..Weight

@inline function phase(varT, extT, l, β)
    tin, tout = varT[extT[1]], varT[extT[2]]
    return cos(π * 2l / β * (tout - tin))
end

include("polarizationKT.jl")
include("polarizationKW.jl")

function MC_Clib(para; kgrid=[para.kF,], ngrid=[0], neval=1e6, reweight_goal=nothing,
    # spinPolarPara::Float64=0.0, # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
    filename::Union{String,Nothing}=nothing, partition=UEG.partition(para.order),
    isLayered2D=false, # whether to use the screened Coulomb interaction in 2D or not 
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/")
)
    @assert para.spin == 2 "Only spin-unpolarized case is supported for compiled C library"
    kF = para.kF

    if isLayered2D
        @assert (para.dim == 2) "Only 2D systems supports the tanh screened Coulomb interaction"
    end

    _partition = Vector{eltype(partition)}()
    for p in partition
        p[1] == 1 && p[3] > 0 && continue
        push!(_partition, p)
    end

    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for p in partition
            push!(reweight_goal, 8.0^(p[1] - 1))
        end
        push!(reweight_goal, 1.0)
    end
    neighbor = UEG.neighbor(_partition)
    println(reweight_goal)
    # println(neighbor)

    diaginfo = Polarization.diagram_loadinfo(para, _partition, root_dir=root_dir)
    polar, result = Polarization.KW_Clib(para, diaginfo;
        root_dir=root_dir, isLayered2D=isLayered2D,
        neighbor=neighbor, reweight_goal=reweight_goal,
        kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:nothread)

    if isnothing(polar) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (ngrid, kgrid, polar)
            end
        end
        for (ip, key) in enumerate(_partition)
            println("Group ", key)
            @printf("%10s  %10s   %10s   %10s   %10s \n", "q/kF", "real(avg)", "err", "imag(avg)", "err")
            r, i = real(polar[key]), imag(polar[key])
            for (in, n) in enumerate(ngrid)
                println("n = $n")
                for (iq, q) in enumerate(kgrid)
                    @printf("%10.6f  %10.6f ± %10.6f   %10.6f ± %10.6f\n", q[1] / kF, r[in, iq].val, r[in, iq].err, i[in, iq].val, i[in, iq].err)
                end
            end
        end
    end
    return polar, result
end

function MC(para; kgrid=[para.kF,], ngrid=[0], neval=1e6, reweight_goal=nothing,
    # spinPolarPara::Float64=0.0, # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
    filename::Union{String,Nothing}=nothing, partition=UEG.partition(para.order),
    isLayered2D=false, # whether to use the screened Coulomb interaction in 2D or not 
    filter=[NoHartree], extK=nothing, optimize_level=1
)
    kF = para.kF

    if isLayered2D
        @assert (para.dim == 2) "Only 2D systems supports the tanh screened Coulomb interaction"
    end

    diagram = Diagram.diagram_parquet_noresponse(:chargePolar, para, partition, filter=filter, extK=extK, optimize_level=optimize_level)

    partition = diagram[1]
    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder, vOrder) in partition
            push!(reweight_goal, 8.0^(order - 1))
        end
        push!(reweight_goal, 1.0)
    end

    println(reweight_goal)

    polar, result = Polarization.KW(para, diagram;
        isLayered2D=isLayered2D,
        neighbor=neighbor, reweight_goal=reweight_goal,
        kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:nothread)

    if isnothing(polar) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (ngrid, kgrid, polar)
            end
        end
        for (ip, key) in enumerate(partition)
            println("Group ", key)
            @printf("%10s  %10s   %10s   %10s   %10s \n", "q/kF", "real(avg)", "err", "imag(avg)", "err")
            r, i = real(polar[key]), imag(polar[key])
            for (in, n) in enumerate(ngrid)
                println("n = $n")
                for (iq, q) in enumerate(kgrid)
                    @printf("%10.6f  %10.6f ± %10.6f   %10.6f ± %10.6f\n", q[1] / kF, r[in, iq].val, r[in, iq].err, i[in, iq].val, i[in, iq].err)
                end
            end
        end
    end
    return polar, result
end

function diagram_loadinfo(paramc::ParaMC, _partition::Vector{T};
    filter=[NoHartree], transferLoop=nothing,
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/"), filename="extvars_chargePolar.jld2"
) where {T}
    diagpara = Vector{DiagPara}()
    extT_labels = Vector{Vector{Int}}[]

    fname = joinpath(root_dir, filename)
    jldopen(fname, "r") do f
        for p in _partition
            key_str = join(string.(p))
            if key_str in keys(f)
                extT = f[key_str][1]
                push!(diagpara, Diagram.diagPara(PolarDiag, paramc.isDynamic, p[1], paramc.spin, filter, transferLoop))
                push!(extT_labels, extT)
            else
                error("$(key_str) not found in $(fname)")
            end
        end
    end
    return (_partition, diagpara, extT_labels)
end

include("source_codeParquetAD/Cwrapper_chargePolar_ParquetAD.jl")

const evalfuncParquetAD_chargePolar_map = Dict(
    (1, 0, 0) => eval_chargePolar_ParquetAD100!,
    (1, 1, 0) => eval_chargePolar_ParquetAD110!,
    (1, 2, 0) => eval_chargePolar_ParquetAD120!,
    (1, 3, 0) => eval_chargePolar_ParquetAD130!,
    (1, 4, 0) => eval_chargePolar_ParquetAD140!,
    (1, 5, 0) => eval_chargePolar_ParquetAD150!,
    (2, 0, 0) => eval_chargePolar_ParquetAD200!,
    (2, 0, 1) => eval_chargePolar_ParquetAD201!,
    (2, 0, 2) => eval_chargePolar_ParquetAD202!,
    (2, 0, 3) => eval_chargePolar_ParquetAD203!,
    (2, 0, 4) => eval_chargePolar_ParquetAD204!,
    (2, 1, 0) => eval_chargePolar_ParquetAD210!,
    (2, 1, 1) => eval_chargePolar_ParquetAD211!,
    (2, 1, 2) => eval_chargePolar_ParquetAD212!,
    (2, 1, 3) => eval_chargePolar_ParquetAD213!,
    (2, 2, 0) => eval_chargePolar_ParquetAD220!,
    (2, 2, 1) => eval_chargePolar_ParquetAD221!,
    (2, 2, 2) => eval_chargePolar_ParquetAD222!,
    (2, 3, 0) => eval_chargePolar_ParquetAD230!,
    (2, 3, 1) => eval_chargePolar_ParquetAD231!,
    (2, 4, 0) => eval_chargePolar_ParquetAD240!,
    (3, 0, 0) => eval_chargePolar_ParquetAD300!,
    (3, 0, 1) => eval_chargePolar_ParquetAD301!,
    (3, 0, 2) => eval_chargePolar_ParquetAD302!,
    (3, 0, 3) => eval_chargePolar_ParquetAD303!,
    (3, 1, 0) => eval_chargePolar_ParquetAD310!,
    (3, 1, 1) => eval_chargePolar_ParquetAD311!,
    (3, 1, 2) => eval_chargePolar_ParquetAD312!,
    (3, 2, 0) => eval_chargePolar_ParquetAD320!,
    (3, 2, 1) => eval_chargePolar_ParquetAD321!,
    (3, 3, 0) => eval_chargePolar_ParquetAD330!,
    (4, 0, 0) => eval_chargePolar_ParquetAD400!,
    (4, 0, 1) => eval_chargePolar_ParquetAD401!,
    (4, 0, 2) => eval_chargePolar_ParquetAD402!,
    (4, 1, 0) => eval_chargePolar_ParquetAD410!,
    (4, 1, 1) => eval_chargePolar_ParquetAD411!,
    (4, 2, 0) => eval_chargePolar_ParquetAD420!,
    (5, 0, 0) => eval_chargePolar_ParquetAD500!,
    (5, 0, 1) => eval_chargePolar_ParquetAD501!,
    (5, 1, 0) => eval_chargePolar_ParquetAD510!,
    (6, 0, 0) => eval_chargePolar_ParquetAD600!
)

end