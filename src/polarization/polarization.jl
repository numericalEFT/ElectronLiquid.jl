module Polarization

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
        df = Parquet.polarization(para, name=:Π)
        d = mergeby(df.diagram; operator=Sum(), name=Symbol("Π $response"), factor=para.spin)[1]
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
    #diag: vector of ExprTreeF64
    result = (partition, diagpara, diag, root)
    return result
end

@inline function phase(varT, extT, l, β)
    tin, tout = varT[extT[1]], varT[extT[2]]
    return exp(1im * π * 2l / β * (tout - tin))
end

include("polarizationKT.jl")
# include("polarizationKW.jl")

end