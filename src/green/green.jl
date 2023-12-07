module Green

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
        type=GreenDiag,
        innerLoopNum=order,
        hasTau=true,
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
    ]
) where {T}
    println("Build the sigma diagrams into an expression tree ...")
    println("Diagram set: ", _partition)

    dim = paramc.dim
    diag = Vector{ExprTreeF64}()
    diagpara = Vector{DiagParaF64}()
    partition = Vector{T}()
    for p in _partition
        para = diagPara(paramc, p[1], filter)
        d = Parquet.green(para)
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
            push!(diag, ExprTree.build(dpp, dim))
        else
            @warn("partition $p doesn't have any diagram. It will be ignored.")
        end
    end

    root = [d.root for d in diag] #get the list of root nodes
    #assign the external Tau to the corresponding diagrams
    #diag: vector of ExprTreeF64
    result = (partition, diagpara, diag, root)
    return result
end

@inline function phase(varT, extT, l, β)
    tin, tout = varT[extT[1]], varT[extT[2]]
    return exp(1im * π * (2l + 1) / β * (tout - tin))
end

include("greenKT.jl")
include("densityKT.jl")
# include("sigmaCuba.jl")
# include("sigmaVegas.jl")

end
