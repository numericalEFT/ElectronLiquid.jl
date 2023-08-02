module Sigma
using Cuba
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

function diagPara(para::ParaMC, order::Int, filter)
    inter = [FeynmanDiagram.Interaction(ChargeCharge, para.isDynamic ? [Instant, Dynamic] : [Instant,]),]  #instant charge-charge interaction
    DiagParaF64(
        type=SigmaDiag,
        innerLoopNum=order,
        hasTau=true,
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
        # Girreducible,
        # Proper,   #one interaction irreduble diagrams or not
        # NoBubble, #allow the bubble diagram or not
    ],
    dR=false # whether to use the derivative of the renormalized interaction, for the RG purpose
) where {T}
    println("Build the sigma diagrams into an expression tree ...")
    println("Diagram set: ", _partition)

    diag = Vector{ExprTreeF64}()
    diagpara = Vector{DiagParaF64}()
    partition = Vector{T}()
    # diagrams = Vector{Diagram{Float64}}()
    for p in _partition
        para = diagPara(paramc, p[1], filter)
        sd::Vector{Diagram{Float64}} = Parquet.sigma(para).diagram
        sdp = DiagTree.derivative(sd, BareGreenId, p[2], index=1)
        sdpp = DiagTree.derivative(sdp, BareInteractionId, p[3], index=2)
        if dR
            sdpp = DiagTree.derivative(sdpp, BareInteractionId, 1, index=3)
        end
        # the Taylor expansion should be d^n f(x) / dx^n / n!, so there is a factor of 1/n! for each derivative
        for d in sdpp
            d.factor *= 1 / factorial(p[2]) / factorial(p[3])
        end
        if isempty(sdpp) == false
            if paramc.isFock && (p != (1, 0, 0)) # the Fock diagram itself should not be removed
                DiagTree.removeHartreeFock!(sdpp)
            end
            push!(diagpara, para)
            push!(partition, p)
            push!(diag, ExprTree.build(sdpp))
            # append!(diagrams, sdpp)
        else
            @warn("partition $p doesn't have any diagram. It will be ignored.")
        end
    end

    # @time diag = [ExprTree.build(diags) for diags in sigma] # DiagTree to ExprTree
    root = [d.root for d in diag] #get the list of root nodes
    #assign the external Tau to the corresponding diagrams
    extT = [[diag[ri].node.object[idx].para.extT::Tuple{Int,Int} for idx in r] for (ri, r) in enumerate(root)]
    result = (partition, diagpara, diag, root, extT)
    # result = (partition, diagpara, diag, root, extT, diagrams)
    return result
end

function diagramGV(paramc::ParaMC, _partition::Vector{T}; filter=[FeynmanDiagram.NoHartree]) where {T}
    diagpara = Vector{DiagParaF64}()
    for p in _partition
        push!(diagpara, diagPara(paramc, p[1], filter))
    end
    FeynGraphs, FermiLabel, BoseLabel, mappings = FeynmanDiagram.diagdictGV(:sigma, _partition, paramc.dim)
    return (_partition, diagpara, FeynGraphs, FermiLabel, BoseLabel, mappings)
end

@inline function phase(varT, extT, l, β)
    tin, tout = varT[extT[1]], varT[extT[2]]
    return exp(1im * π * (2l + 1) / β * (tout - tin))
end

#include("sigma_generic.jl")
include("sigmaKW.jl")
include("sigmaCuba.jl")
include("sigmaVegas.jl")
include("sigmaGV.jl")

function MC(para; kgrid=[para.kF,], ngrid=[-1, 0, 1], neval=1e6, filename::Union{String,Nothing}=nothing, diagtype=:Parquet)
    kF = para.kF
    _order = para.order

    # partition = [(1, 0, 0), (2, 0, 1), (2, 1, 0), (3, 0, 0)]
    partition = UEG.partition(_order)
    neighbor = UEG.neighbor(partition)
    reweight_goal = Float64[]
    for (order, sOrder, vOrder) in partition
        push!(reweight_goal, 4.0^(order + vOrder - 1))
    end
    push!(reweight_goal, 2.0)

    if diagtype == :GV
        diagram = Sigma.diagramGV(para, partition)
        sigma, result = Sigma.GV(para, diagram;
            neighbor=neighbor, reweight_goal=reweight_goal,
            kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:thread)
    elseif diagtype == :Parquet
        diagram = Sigma.diagram(para, partition)
        sigma, result = Sigma.KW(para, diagram;
            neighbor=neighbor, reweight_goal=reweight_goal,
            kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:thread)
    else
        error("unknown diagrams' generated type")
    end

    if isnothing(sigma) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (ngrid, kgrid, sigma)
            end
        end
        for (ip, key) in enumerate(partition)
            println("Group ", key)
            @printf("%10s  %10s   %10s   %10s   %10s \n", "q/kF", "real(avg)", "err", "imag(avg)", "err")
            r, i = real(sigma[key]), imag(sigma[key])
            for (in, n) in enumerate(ngrid)
                println("n = $n")
                for (iq, q) in enumerate(kgrid)
                    @printf("%10.6f  %10.6f ± %10.6f   %10.6f ± %10.6f\n", q[1] / kF, r[in, iq].val, r[in, iq].err, i[in, iq].val, i[in, iq].err)
                end
            end
        end
    end
    return sigma, result
end


end