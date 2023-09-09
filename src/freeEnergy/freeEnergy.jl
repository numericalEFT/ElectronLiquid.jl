module FreeEnergy

using JLD2
using Printf, LinearAlgebra
using ..CompositeGrids
using ..ElectronGas
using ..MCIntegration
using ..Lehmann

using ..FeynmanDiagram
using ..Measurements

using ..UEG
using ..Propagator

function diagPara(para::ParaMC, order::Int, filter)
    inter = [FeynmanDiagram.Interaction(ChargeCharge, para.isDynamic ? [Instant, Dynamic] : [Instant,]),]  #instant charge-charge interaction
    DiagParaF64(
        type=VacuumDiag,
        innerLoopNum=order + 1,
        hasTau=true,
        loopDim=para.dim,
        spin=para.spin,
        firstLoopIdx=1,
        interaction=inter,
        filter=filter
    )
end

function diagramGV(paramc::ParaMC, _partition::Vector{T}; filter=[FeynmanDiagram.NoHartree]) where {T}
    diagpara = Vector{DiagParaF64}()
    for p in _partition
        push!(diagpara, diagPara(paramc, p[1], filter))
    end
    FeynGraphs, FermiLabel, BoseLabel, mappings = FeynmanDiagram.diagdictGV(:freeEnergy, _partition, paramc.dim)
    return (_partition, diagpara, FeynGraphs, FermiLabel, BoseLabel, mappings)
end

include("freeEnergyGV.jl")


function MC(para; neval=1e6, filename::Union{String,Nothing}=nothing, diagtype=:GV,
    filter=[NoHartree,], partition=UEG.partition(para.order), reweight_goal=nothing, verbose=0)

    diagram = FreeEnergy.diagramGV(para, partition, filter=filter)
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

    freeE, result = FreeEnergy.GV(para, diagram; neighbor=neighbor, print=verbose,
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

end