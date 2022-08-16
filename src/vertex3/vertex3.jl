module Ver3

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
import ..Weight

function diagPara(para::ParaMC, order, filter, transferLoop)
    inter = [FeynmanDiagram.Interaction(ChargeCharge, para.isDynamic ? [Instant, Dynamic] : [Instant,]),]  #instant charge-charge interaction
    return DiagParaF64(
        type=Ver3Diag,
        innerLoopNum=order,
        hasTau=true,
        loopDim=para.dim,
        spin=para.spin,
        firstLoopIdx=3,
        interaction=inter,
        filter=filter,
        transferLoop=transferLoop
    )
end

function diagram(paramc::ParaMC, _partition::Vector{T};
    filter=[
        NoHatree,
        # Girreducible,
        Proper,   #one interaction irreduble diagrams or not
        # NoBubble, #allow the bubble diagram or not
    ]
) where {T}
    # println("Build the vertex4 diagrams into an experssion tree ...")
    # _partition = UEG.partition(order)
    # println("Diagram set: ", _partition)

    Kin, Qout = zeros(16), zeros(16)
    Qout[1], Kin[2]=1.0, 1.0
    legK = [Qout, Kin]

    diag = Vector{ExprTreeF64}()
    diagpara = Vector{DiagParaF64}()
    partition = Vector{T}()
    for p in _partition
        para = diagPara(paramc, p[1], filter, Qout)
        d::Vector{Diagram{Float64}} = Parquet.vertex3(para, legK).diagram
        d = DiagTree.derivative(d, BareGreenId, p[2], index=1)
        d = DiagTree.derivative(d, BareInteractionId, p[3], index=2)
        if isempty(d) == false
            if paramc.isFock # remove the Fock subdiagrams
                DiagTree.removeHatreeFock!(d)
            end
            push!(diagpara, para)
            push!(diag, ExprTree.build(d))
            push!(partition, p)
        else
            @warn("partition $p doesn't have any diagram. It will be ignored.")
        end
    end

    # diag = [ExprTree.build(d) for d in ver4]    #experssion tree representation of diagrams 
    rootuu = [[idx for idx in d.root if d.node.object[idx].para.response == UpUp] for d in diag] #select the diagram with upup
    rootud = [[idx for idx in d.root if d.node.object[idx].para.response == UpDown] for d in diag] #select the diagram with updown
    #assign the external Tau to the corresponding diagrams
    extTuu = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootuu)]
    extTud = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootud)]
    return (partition, diagpara, diag, [rootuu, rootud], [extTuu, extTud])
end

@inline function phaseF(varT, extT, nin, nout, β)
    # println(extT)
    tb, tfin, tfout = varT[extT[1]], varT[extT[2]], varT[extT[3]]
    win, wout = π * (2nin + 1) / β, π * (2nout + 1) / β
    wq = win-wout
    return cos(tfin*win-tfout*wout-tb*wq)

    # if (idx == 1)
    #     return cos(π / β * ((2tb) - (tfin + tfout)))
    #     # return cos(π / β * ((2tb) - 3 * tfin + tfout))
    #     # return cos(π / β * ((2tb) - 3 * tfin + tfout))
    # else
    #     return cos(π / β * (tfin - tfout))
    # end
end

@inline function phaseC(varT, extT, nin, nout, β)
    # println(extT)
    tb, tfin, tfout = varT[extT[1]], varT[extT[2]], varT[extT[3]]
    win, wout = π * (2nin + 1) / β, π * (2nout + 1) / β
    wq = win-wout
    return exp(-1im * (tfin * win - tfout * wout - tb*wq))
end

# @inline function phase(varT, extT, n, β)
#     # println(extT)
#     return phase(varT, extT, n[1], n[2], n[3], β)
# end

include("ver3KW.jl")

end