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
    if order == 0
        return DiagParaF64(type=VacuumDiag, innerLoopNum=order + 1, hasTau=true, loopDim=para.dim, spin=para.spin, firstLoopIdx=1,
            totalTauNum=1,
            interaction=inter, filter=filter
        )
    else
        return DiagParaF64(
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
end

function diagramGV(paramc::ParaMC, _partition::Vector{T};
    filter=[FeynmanDiagram.NoHartree], spinPolarPara::Float64=0.0) where {T}
    diagpara = Vector{DiagParaF64}()
    gkeys = Vector{T}()
    for p in _partition
        p[1] == 0 && p[3] > 0 && continue
        push!(diagpara, diagPara(paramc, p[1], filter))
        push!(gkeys, p)
    end
    FeynGraphs, FermiLabel, BoseLabel, mappings = FeynmanDiagram.diagdictGV(:freeEnergy, gkeys, paramc.dim, spinPolarPara=spinPolarPara)
    return (gkeys, diagpara, FeynGraphs, FermiLabel, BoseLabel, mappings)
end

include("freeEnergyGV.jl")
include("freeEnergyGV_2dscreened.jl")


function MC(para; neval=1e6, filename::Union{String,Nothing}=nothing, reweight_goal=nothing,
    spinPolarPara::Float64=0.0, # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
    isScreened::Bool=false, # whether to use the screened Coulomb interaction
    filter=[NoHartree,], partition=UEG.partition(para.order, offset=0), verbose=-1)

    diagram = FreeEnergy.diagramGV(para, partition, filter=filter, spinPolarPara=spinPolarPara)
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

    if isScreened == false
        freeE, result = FreeEnergy.GV(para, diagram; neighbor=neighbor, print=verbose,
            reweight_goal=reweight_goal, neval=neval, parallel=:nothread)
    else
        # screened 2D Coulomb interaction
        @assert para.dim == 2 "Only 2D is supported for the tanh screened Coulomb interaction"
        freeE, result = FreeEnergy.GVscreened(para, diagram; neighbor=neighbor, print=verbose,
            reweight_goal=reweight_goal, neval=neval, parallel=:nothread)
    end

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