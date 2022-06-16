using Printf, LinearAlgebra
using CompositeGrids
using ElectronGas
using Parameters, Random, DataFrames
using MCIntegration
using Lehmann

using FeynmanDiagram
using StaticArrays
using JLD2

const steps = 1e7
const Order = 4

include("../common/interaction.jl")

# println(dW0)
# exit(0)
const lgrid = [0, 1]
const Nl = length(lgrid)

println("Build the diagrams into an experssion tree ...")

partition = [(0, 0, 0), (0, 1, 0), (0, 2, 0),
    (1, 0, 0),  # order 1
    (2, 0, 0), (1, 1, 0), (1, 0, 1),  #order 2
    (3, 0, 0), (2, 1, 0), (2, 0, 1), (1, 1, 1), (1, 2, 0), (1, 0, 2), #order 3
    (4, 0, 0), (3, 1, 0), (3, 0, 1), (2, 1, 1), (2, 2, 0), (2, 0, 2), (1, 3, 0), (1, 0, 3), (1, 2, 1), (1, 1, 2) #order 4
]

partition = [p for p in sort(partition) if p[1] + p[2] + p[3] <= Order]

println("Diagram set: ", partition)

diagPara(order) = GenericPara(diagType=GreenDiag, innerLoopNum=order, hasTau=true, loopDim=dim, spin=spin, firstLoopIdx=2,
    interaction=[FeynmanDiagram.Interaction(ChargeCharge, [
        Instant,
        # Dynamic
    ]),],  #instant charge-charge interaction
    filter=[
    # Girreducible,
    # Proper,   #one interaction irreduble diagrams or not
    # NoBubble, #allow the bubble diagram or not
    # NoFock,
    ]
)

# dpara = GenericPara(diagType=SigmaDiag, innerLoopNum=4, hasTau=true, loopDim=dim, spin=spin, firstLoopIdx=2,
#     interaction=[FeynmanDiagram.Interaction(ChargeCharge, [
#         Instant,
#         # Dynamic
#     ]),],  #instant charge-charge interaction
#     filter=[
#     # Girreducible,
#     # Proper,   #one interaction irreduble diagrams or not
#     # NoBubble, #allow the bubble diagram or not
#     NoFock,
#     ]
# )
# sigma4 = Parquet.sigma(dpara)
# sigma4 = ExprTree.build(sigma4.diagram)

gdiag = Dict()
for p in partition
    d = Parquet.green(diagPara(p[1]))
    println(p, ",", d)
    d = DiagTree.derivative(d, BareGreenId, p[2])
    if p[3]>0
        d = DiagTree.derivative(d, BareInteractionId, p[3])
    end
    if p == (1, 0, 0)
        gdiag[p] = d
    else
        gdiag[p] = DiagTree.removeHatreeFock!(d)
    end
end
# DiagTree.removeHatreeFock!(sigma[2, 0, 0])
# println(sigma[2, 0, 0])
# plot_tree(sigma[(2, 0, 0)], maxdepth=8)
# exit(0)

gdiag = [gdiag[p] for p in partition]
const diagpara = [diags.id.para for diags in gdiag]
const diag = [ExprTree.build(diags) for diags in gdiag]
const root = [d.root for d in diag] #select the diagram with upup
#assign the external Tau to the corresponding diagrams
# const extT = [[diag[ri].node.object[idx].para.extT for idx in r] for (ri, r) in enumerate(root)]
# println(root)
# println(extT)
# exit(0)

function integrand(config)
    order = config.curr
    varK, varT = config.var[1], config.var[2]

    ExprTree.evalNaive!(diag[order], varK.data, varT.data, eval)
    w = sum(diag[order].node.current[r] for (ri, r) in enumerate(root[order]))
    # w = diag[order].node.current[1]
    # println(l)
    # println(w)
    # println(diag[order].node.current)
    # println(diag[order].node.current[1])
    # exit(0)
    # println(wuu, ",  ", wud)
    # w = 0.5 / β
    # if order == length(diag)
    #     ExprTree.evalNaive!(sigma4, varK.data, varT.data, eval)
    #     wp = sum(sigma4.node.current[r] * phase(varT, extT[order][ri], l) for (ri, r) in enumerate(sigma4.root))
    #     @assert abs(wp-w)<1e-6 "$wp vs $w" 
    #     println(partition[order])
    #     println("$wp vs $w")
    #     exit(0)
    # end
    # println(order, ": ", w)
    # exit(0)
    return w*2/(2*π)^3*β
     #the current implementation of sigma has an additional minus sign compared to the standard defintion
end

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    # println(config.observable[1][1])
    o = config.curr
    weight = integrand(config)
    config.observable[o] += weight / abs(weight) * factor
end

function MC()
    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF)
    # K.data[:, 1] .= [kF, 0.0, 0.0]
    T = MCIntegration.Tau(β, β / 2.0, offset=2)
    T.data[1] = 0.0
    T.data[2] = -1e-8
    # T = MCIntegration.Tau(1.0, 1.0 / 2.0)
    # X = MCIntegration.Discrete(lgrid[1], lgrid[end])

    dof = [[p.innerLoopNum+1, p.totalTauNum - 2] for p in diagpara] # K, T, ExtKidx
    # dof = [[1, 1, 1], [1, 1, 1], [1, 1, 1], [2, 3, 1]]
    # dof = [[1, 0, 1], [1, 0, 1], [1, 0, 1], [2, 1, 1]]
    # println(dof)
    # exit(0)
    obs = zeros(Float64, length(dof)) # observable for the Fock diagram 

    config = MCIntegration.Configuration(steps, (K, T), dof, obs)
    avg, std = MCIntegration.sample(config, integrand, measure; print=0, Nblock=16, reweight=10000)

    if isnothing(avg) == false

        jldsave("data.jld2", avg=avg, std=std)

        open("data.dat", "w") do f
            @printf(f, "#%7s %16s %16s %16s %16s\n", "freq", "real", "real err", "imag", "imag err")
            for o in 1:length(dof)
                write(f, "# $(partition[o])\n")
                @printf(f, "%16.8f %16.8f\n", avg[o], std[o])
            end

        end

    end

end

MC()