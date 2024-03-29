module Polarization
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
import ..ExprTreeF64

function diagPara(para::ParaMC, order::Int, filter, response::Response)
    inter = [FeynmanDiagram.Interaction(response, para.isDynamic ? [Instant, Dynamic] : [Instant,]),]  #instant charge-charge interaction
    return DiagParaF64(
        type=PolarDiag,
        hasTau=true,
        innerLoopNum=order,
        loopDim=para.dim,
        spin=para.spin,
        firstLoopIdx=2,
        interaction=inter,
        filter=filter
    )
end

function diagram(paramc::ParaMC, _partition::Vector{T};
    filter=[
        FeynmanDiagram.NoHartree,
        FeynmanDiagram.Proper,
    ],
    response::Response=ChargeCharge
) where {T}
    println("Build the polarization diagrams into an expression tree ...")
    println("Diagram set: ", _partition)
    if response != ChargeCharge
        error("$response response not yet implemented!")
    end

    diag = Vector{ExprTreeF64}()
    diagpara = Vector{DiagParaF64}()
    partition = Vector{T}()
    for p in _partition
        # NOTE: Parquet.polarization returns a dataframe; merge the spin channels
        #       and sum over the external spin: Π_chch = 4 Π_s = 2(Π↑↑ + Π↑↓)
        para = diagPara(paramc, p[1], filter, response)
        extK = DiagTree.getK(para.totalLoopNum, 1)
        df = Parquet.polarization(para, extK; name=:Π)
        @assert allequal(df.extT) "All external times must be the same for the polarization diagrams!"
        extT = df.extT[1]

        # NOTE: We explicitly construct the diagram d so that extT is retained in root DiagramIds
        polar_id = PolarId(para, response; k=extK, t=extT)
        d = Diagram{Float64}(
            polar_id,
            Sum(),
            df.diagram;
            name=:Π,
            factor=para.spin
        )
        dp = DiagTree.derivative([d,], BareGreenId, p[2], index=1)
        dpp = DiagTree.derivative(dp, BareInteractionId, p[3], index=2)

        # the Taylor expansion should be d^n f(x) / dx^n / n!, so there is a factor of 1/n! for each derivative
        for d in dpp
            d.factor *= 1 / factorial(p[2]) / factorial(p[3])
        end

        if isempty(dpp) == false
            if paramc.isFock && (p != (1, 0, 0)) # the Fock diagram itself should not be removed
                DiagTree.removeHartreeFock!(dpp)
            end
            push!(diagpara, para)
            push!(partition, p)
            push!(diag, ExprTree.build(dpp))
        else
            @warn("partition $p doesn't have any diagram. It will be ignored.")
        end
    end

    root = [d.root for d in diag] #get the list of root nodes
    #assign the external Tau to the corresponding diagrams
    extT = [[diag[ri].node.object[idx].para.extT::Tuple{Int,Int} for idx in r] for (ri, r) in enumerate(root)]
    #diag: vector of ExprTreeF64
    result = (partition, diagpara, diag, root, extT)
    return result
end

function diagramGV(paramc::ParaMC, _partition::Vector{T}; filter=[FeynmanDiagram.NoHartree],
    response::Response=ChargeCharge, spinPolarPara::Float64=0.0) where {T}
    gkeys = Vector{T}()
    for p in _partition
        if p[1] == 1 && p[3] > 0
            @warn("partition $p doesn't have any diagram. It will be ignored.")
        else
            push!(gkeys, p)
        end
    end

    if response == ChargeCharge
        FeynGraphs, labelProd, mappings = FeynmanDiagram.diagdictGV(:chargePolar, gkeys, spinPolarPara=spinPolarPara)
    elseif response == SpinSpin
        FeynGraphs, labelProd, mappings = FeynmanDiagram.diagdictGV(:spinPolar, gkeys, spinPolarPara=spinPolarPara)
    else
        error("$response response not yet implemented!")
    end
    diagpara = Vector{DiagParaF64}()
    for p in gkeys
        push!(diagpara, diagPara(paramc, p[1], filter, response))
    end
    return (gkeys, diagpara, FeynGraphs, labelProd, mappings)
end

@inline function phase(varT, extT, l, β)
    tin, tout = varT[extT[1]], varT[extT[2]]
    return cos(π * 2l / β * (tout - tin))
end

include("polarizationKT.jl")
include("polarizationKW.jl")
include("polarizationGV.jl")

function MC(para; response::Response=ChargeCharge, kgrid=[0.0, para.kF], ngrid=[0,],
    filter=[NoHartree,], partition=UEG.partition(para.order), diagtype=:GV,
    spinPolarPara::Float64=0.0, # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
    neval=1e6, reweight_goal=nothing, filename::Union{String,Nothing}=nothing)

    if diagtype == :GV
        diagram = Polarization.diagramGV(para, partition, response=response, filter=filter, spinPolarPara=spinPolarPara)
    elseif diagtype == :Parquet
        diagram = Polarization.diagram(para, partition, response=response, filter=filter)
    else
        error("unknown diagrams' generated type")
    end

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

    if diagtype == :GV
        polar, result = Polarization.GV(para, diagram;
            neighbor=neighbor, reweight_goal=reweight_goal,
            kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:nothread)
    elseif diagtype == :Parquet
        diagram = Polarization.diagram(para, partition)
        polar, result = Polarization.KW(para, diagram;
            neighbor=neighbor, reweight_goal=reweight_goal,
            kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:nothread)
    end

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
                    @printf("%10.6f  %10.6f ± %10.6f   %10.6f ± %10.6f\n", q[1] / para.kF, r[in, iq].val, r[in, iq].err, i[in, iq].val, i[in, iq].err)
                end
            end
        end
    end
    return polar, result
end

end