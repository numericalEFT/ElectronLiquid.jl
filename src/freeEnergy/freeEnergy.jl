module FreeEnergy

using JLD2, CSV
using Printf, LinearAlgebra, DataFrames
using ..CompositeGrids
using ..ElectronGas
using ..MCIntegration
using ..Lehmann

using ..FeynmanDiagram
import ..FeynmanDiagram.FrontEnds: Filter, NoHartree, NoFock, DirectOnly, Wirreducible, Girreducible, NoBubble, Proper
import ..FeynmanDiagram.FrontEnds: AnalyticProperty, Instant, Dynamic
import ..FeynmanDiagram.Parquet: DiagPara, VacuumDiag
using ..Measurements

using ..UEG
using ..Propagator
import ..Propagator: LeafStateAD
using ..Diagram

include("parquetAD.jl")
include("parquetAD_Clib.jl")

function MC(para; neval=1e6, reweight_goal=nothing,
    # spinPolarPara::Float64=0.0, # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
    filename::Union{String,Nothing}=nothing, partition=UEG.partition(para.order, offset=0),
    isLayered2D=false, # whether to use the screened Coulomb interaction in 2D or not 
    filter=[NoHartree], optimize_level=1, verbose=-1
)
    diagram = Diagram.diagram_GV_freeE(para, partition, filter=filter, optimize_level=optimize_level)

    partition = diagram[1]
    neighbor = UEG.neighbor(partition)
    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder, vOrder) in partition
            reweight_factor = 2.0^(order - 1)
            push!(reweight_goal, reweight_factor)
        end
        push!(reweight_goal, 1.0)
    end

    freeE, result = FreeEnergy.ParquetAD(para, diagram; isLayered2D=isLayered2D, neighbor=neighbor, print=verbose,
        reweight_goal=reweight_goal, neval=neval, parallel=:nothread)

    if isnothing(freeE) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (freeE,)
            end
        end
        for (ip, key) in enumerate(partition)
            println("Group ", key)
            @printf("%10s   %10s \n", "avg", "err")
            r = freeE[key]
            @printf("%10.6f ± %10.6f\n", r.val, r.err)
        end
    end
    return freeE, result
end

function MC_Clib(para; neval=1e6, reweight_goal=nothing,
    # spinPolarPara::Float64=0.0, # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
    filename::Union{String,Nothing}=nothing, partition=UEG.partition(para.order, offset=0),
    isLayered2D=false, # whether to use the screened Coulomb interaction in 2D or not 
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/"), verbose=-1
)
    @assert para.spin == 2 "Only spin-unpolarized case is supported for compiled C library"

    if isLayered2D
        @assert (para.dim == 2) "Only 2D systems supports the tanh screened Coulomb interaction"
    end
    neighbor = UEG.neighbor(partition)
    _partition = Vector{eltype(partition)}()
    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder, vOrder) in partition
            order == 0 && vOrder > 0 && continue
            reweight_factor = 2.0^(order - 1)
            push!(reweight_goal, reweight_factor)
            push!(_partition, (order, sOrder, vOrder))
        end
        push!(reweight_goal, 1.0)
    end

    diagpara = Vector{DiagPara}()
    for p in _partition
        push!(diagpara, Diagram.diagPara(VacuumDiag, para.isDynamic, p[1], para.spin))
    end
    freeE, result = FreeEnergy.ParquetAD_Clib(para, (_partition, diagpara);
        root_dir=root_dir, isLayered2D=isLayered2D,
        neighbor=neighbor, reweight_goal=reweight_goal,
        neval=neval, parallel=:nothread, print=verbose)

    if isnothing(freeE) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (freeE,)
            end
        end
        for key in _partition
            println("Group ", key)
            @printf("%10s   %10s \n", "avg", "err")
            r = freeE[key]
            @printf("%10.6f ± %10.6f\n", r.val, r.err)
        end
    end
    return freeE, result
end

include("source_codeParquetAD/Cwrapper_freeEnergy_ParquetAD.jl")
const evalfunc_freeEnergy_map = Dict(
    (0, 0, 0) => eval_freeEnergy_ParquetAD000!,
    (0, 1, 0) => eval_freeEnergy_ParquetAD010!,
    (0, 2, 0) => eval_freeEnergy_ParquetAD020!,
    (0, 3, 0) => eval_freeEnergy_ParquetAD030!,
    (0, 4, 0) => eval_freeEnergy_ParquetAD040!,
    (0, 5, 0) => eval_freeEnergy_ParquetAD050!,
    (1, 0, 0) => eval_freeEnergy_ParquetAD100!,
    (1, 0, 1) => eval_freeEnergy_ParquetAD101!,
    (1, 0, 2) => eval_freeEnergy_ParquetAD102!,
    (1, 0, 3) => eval_freeEnergy_ParquetAD103!,
    (1, 0, 4) => eval_freeEnergy_ParquetAD104!,
    (1, 1, 0) => eval_freeEnergy_ParquetAD110!,
    (1, 1, 1) => eval_freeEnergy_ParquetAD111!,
    (1, 1, 2) => eval_freeEnergy_ParquetAD112!,
    (1, 1, 3) => eval_freeEnergy_ParquetAD113!,
    (1, 2, 0) => eval_freeEnergy_ParquetAD120!,
    (1, 2, 1) => eval_freeEnergy_ParquetAD121!,
    (1, 2, 2) => eval_freeEnergy_ParquetAD122!,
    (1, 3, 0) => eval_freeEnergy_ParquetAD130!,
    (1, 3, 1) => eval_freeEnergy_ParquetAD131!,
    (1, 4, 0) => eval_freeEnergy_ParquetAD140!,
    (2, 0, 0) => eval_freeEnergy_ParquetAD200!,
    (2, 0, 1) => eval_freeEnergy_ParquetAD201!,
    (2, 0, 2) => eval_freeEnergy_ParquetAD202!,
    (2, 0, 3) => eval_freeEnergy_ParquetAD203!,
    (2, 1, 0) => eval_freeEnergy_ParquetAD210!,
    (2, 1, 1) => eval_freeEnergy_ParquetAD211!,
    (2, 1, 2) => eval_freeEnergy_ParquetAD212!,
    (2, 2, 0) => eval_freeEnergy_ParquetAD220!,
    (2, 2, 1) => eval_freeEnergy_ParquetAD221!,
    (2, 3, 0) => eval_freeEnergy_ParquetAD230!,
    (3, 0, 0) => eval_freeEnergy_ParquetAD300!,
    (3, 0, 1) => eval_freeEnergy_ParquetAD301!,
    (3, 0, 2) => eval_freeEnergy_ParquetAD302!,
    (3, 1, 0) => eval_freeEnergy_ParquetAD310!,
    (3, 1, 1) => eval_freeEnergy_ParquetAD311!,
    (3, 2, 0) => eval_freeEnergy_ParquetAD320!,
    (4, 0, 0) => eval_freeEnergy_ParquetAD400!,
    (4, 0, 1) => eval_freeEnergy_ParquetAD401!,
    (4, 1, 0) => eval_freeEnergy_ParquetAD410!,
    (5, 0, 0) => eval_freeEnergy_ParquetAD500!,
)

end
