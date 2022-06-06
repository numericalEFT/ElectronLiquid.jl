using Printf, LinearAlgebra
using CompositeGrids
using ElectronGas
using Parameters, Random, DataFrames
using MCIntegration
using Lehmann

using FeynmanDiagram
using StaticArrays

const steps = 1e7

include("../common/interaction.jl")

# println(dW0)
# exit(0)
const lgrid = [0, 1, 2, 3, 4, 5]
const Nl = length(lgrid)

println("Build the diagrams into an experssion tree ...")

const Order = 2

diagPara(order) = GenericPara(diagType=GreenDiag, innerLoopNum=order - 1, hasTau=true, loopDim=dim, spin=spin, firstLoopIdx=2,
    interaction=[FeynmanDiagram.Interaction(ChargeCharge, [
        Instant,
        Dynamic
    ]),],  #instant charge-charge interaction
    filter=[
        # Girreducible,
        Proper,   #one interaction irreduble diagrams or not
        NoBubble, #allow the bubble diagram or not
    ]
)

const diagpara = [diagPara(o) for o in 1:Order]
green = [Parquet.green(diagpara[i]) for i in 1:Order]   #diagram of different orders
# println(ver3[2])
# println(diagpara[2].totalTauNum)
# plot_tree(green[2])
# exit(0)
# dver3_w = DiagTree.derivative(ver3[1].diagram, BareInteractionId)
# plot_tree(dver3_w)

green_g = DiagTree.derivative(green[1], BareGreenId)
# println(dver3_g)
# println(typeof(dver3_g))
# plot_tree(green_g)
# dver3_w = DiagTree.derivative(ver3[1].diagram, BareInteractionId)
# plot_tree(dver3_w)
# plot_tree(ver3[1].diagram)
# exit(0)
# plot_tree(ver4uu[1][1])
# plot_tree(ver4[1].diagram, maxdepth = 9)
const diag = [ExprTree.build(green[1]), ExprTree.build(green_g), ExprTree.build(green[2])]    #experssion tree representation of diagrams 
# const diag = [ExprTree.build(dver3_g), ExprTree.build(ver3[2].diagram)]    #experssion tree representation of diagrams 
# println(diag[1].root)
# println(diag[2].root)
# println(length(diag[1].node.current))
const root = [d.root[1] for d in diag] #select the diagram with upup
#assign the external Tau to the corresponding diagrams
const extT = [diag[ri].node.object[r].para.extT for (ri, r) in enumerate(root)]
# println(extT)
# exit(0)


@inline function phase(varT, extT, l)
    # println(extT)
    tin, tout = varT[extT[1]], varT[extT[2]]
    return exp(1im * π * (2l + 1) / β * (tout - tin))
end

function integrand(config)
    order = config.curr
    l = config.var[3][1]
    varK, varT = config.var[1], config.var[2]

    ExprTree.evalNaive!(diag[order], varK.data, varT.data, eval)
    w = sum(diag[order].node.current[r] * phase(varT, extT[order], l) for (ri, r) in enumerate(root[order]))
    # w = diag[order].node.current[1]
    # println(l)
    # println(w)
    # println(diag[order].node.current)
    # println(diag[order].node.current[1])
    # exit(0)
    # println(wuu, ",  ", wud)
    # w = 0.5 / β
    if order == 1
        return w / (β / π)
    elseif order == 2 || order == 3
        return w / (β / π)^2
    else
        error("not implemented!")
    end
end

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    l = config.var[3][1]
    # println(config.observable[1][1])
    o = config.curr
    weight = integrand(config)
    config.observable[o, l+1] += weight / abs(weight) * factor
end

function MC()
    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= [kF, 0.0, 0.0]
    T = MCIntegration.Tau(β, β / 2.0, offset=1)
    T.data[1] = 0.0
    # T = MCIntegration.Tau(1.0, 1.0 / 2.0)
    X = MCIntegration.Discrete(lgrid[1], lgrid[end])

    # dof = [[diagpara[o].innerLoopNum, diagpara[o].totalTauNum - 1, 1] for o in 1:Order] # K, T, ExtKidx
    dof = [[0, 1, 1], [0, 1, 1], [1, 3, 1]]
    # println(dof)
    obs = zeros(ComplexF64, length(dof), Nl) # observable for the Fock diagram 

    config = MCIntegration.Configuration(steps, (K, T, X), dof, obs)
    avg, std = MCIntegration.sample(config, integrand, measure; print=0, Nblock=16, reweight=10000)

    if isnothing(avg) == false
        for o in 1:length(dof)
            println("diagram ", o)
            for li in 1:Nl
                # @printf("%8.4f   %8.4f ±%8.4f\n", lgrid[li], avg[o, li], std[o, li])
                println(lgrid[li], "   ", avg[o, li], "  +-  ", std[o, li])
            end
        end
    end

end

MC()