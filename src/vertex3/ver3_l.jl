using Printf, LinearAlgebra
using CompositeGrids
using ElectronGas
using Parameters, Random, DataFrames
using MCIntegration
using Lehmann

using FeynmanDiagram
using StaticArrays

const steps = 1e6
const isF = true

include("../common/interaction.jl")

# println(dW0)
# exit(0)
const lgrid = [1, 2]
const Nl = length(lgrid)

println("Build the diagrams into an experssion tree ...")

const Order = 2

Qin = [1.0, 0, 0, 0, 0]
Kin = [0, 1.0, 0, 0, 0]
legK = [Qin, Kin]

diagPara(order) = GenericPara(diagType=Ver3Diag, innerLoopNum=order, hasTau=true, loopDim=dim, spin=spin, firstLoopIdx=3,
    interaction=[FeynmanDiagram.Interaction(ChargeCharge, [
        Instant,
        Dynamic
    ]),],  #instant charge-charge interaction
    filter=[
        # Girreducible,
        Proper,   #one interaction irreduble diagrams or not
        # NoBubble, #allow the bubble diagram or not
    ],
    transferLoop=Qin[1:2+order]
)

const diagpara = [diagPara(o) for o in 1:Order]
ver3 = [Parquet.vertex3(diagpara[i], legK) for i in 1:Order]   #diagram of different orders
# println(ver3[2])

# plot_tree(ver3[2].diagram)
# exit(0)
# dver3_w = DiagTree.derivative(ver3[1].diagram, BareInteractionId)
# plot_tree(dver3_w)

# dver3_g = DiagTree.derivative(ver3[1].diagram, BareGreenId)
# println(dver3_g)
# println(typeof(dver3_g))
# plot_tree(dver3_g)
# dver3_w = DiagTree.derivative(ver3[1].diagram, BareInteractionId)
# plot_tree(dver3_w)
# plot_tree(ver3[1].diagram)
# exit(0)
# plot_tree(ver4uu[1][1])
# plot_tree(ver4[1].diagram, maxdepth = 9)
const diag = [ExprTree.build(ver3[o].diagram) for o in 1:Order]    #experssion tree representation of diagrams 
# const diag = [ExprTree.build(dver3_w), ExprTree.build(ver3[2].diagram)]    #experssion tree representation of diagrams 
# println(diag[1].root)
# println(length(diag[1].node.current))
const rootuu = [[idx for idx in d.root if d.node.object[idx].para.response == UpUp] for d in diag] #select the diagram with upup
const rootud = [[idx for idx in d.root if d.node.object[idx].para.response == UpDown] for d in diag] #select the diagram with updown
#assign the external Tau to the corresponding diagrams
const extTuu = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootuu)]
const extTud = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootud)]


@inline function phase(varT, extT)
    # println(extT)
    tb, tfin, tfout = varT[extT[1]], varT[extT[2]], varT[extT[3]]
    if (isF)
        return cos(π / β * ((2tb) - (tfin + tfout)))
        # return cos(π / β * ((2tb) - 3 * tfin + tfout))
    else
        return cos(π / β * (tfin - tfout))
    end
end

function integrand(config)
    order = config.curr
    x = config.var[3][1]
    varK, varT = config.var[1], config.var[2]

    varK.data[:, 2] = [kF * x, kF * sqrt(1 - x^2), 0.0]

    ExprTree.evalNaive!(diag[order], varK.data, varT.data, eval)
    if !isempty(rootuu[order])
        wuu = sum(diag[order].node.current[root] * phase(varT, extTuu[order][ri]) for (ri, root) in enumerate(rootuu[order]))
    else
        wuu = 0.0
    end
    if !isempty(rootud[order])
        wud = sum(diag[order].node.current[root] * phase(varT, extTud[order][ri]) for (ri, root) in enumerate(rootud[order]))
    else
        wud = 0.0
    end
    # println(wuu, ",  ", wud)
    return Weight(wuu / β, wud / β)
end

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    x = config.var[3][1]
    # println(config.observable[1][1])
    o = config.curr
    weight = integrand(config)
    config.observable[o, 1, 1] += weight.d / 2 / abs(weight) * factor
    config.observable[o, 1, 2] += weight.e / 2 / abs(weight) * factor
    config.observable[o, 2, 1] += weight.d * x / 2 / abs(weight) * factor
    config.observable[o, 2, 2] += weight.e * x / 2 / abs(weight) * factor
end

function MC()
    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=2)
    K.data[:, 1] .= [0.0, 0.0, 0.0]
    T = MCIntegration.Tau(β, β / 2.0)
    X = MCIntegration.Continuous([-1.0, 1.0], 0.2) #x=cos(θ)

    dof = [[diagpara[o].innerLoopNum, diagpara[o].totalTauNum, 1] for o in 1:Order] # K, T, ExtKidx
    # println(dof)
    obs = zeros(Order, Nl, 2) # observable for the Fock diagram 

    config = MCIntegration.Configuration(steps, (K, T, X), dof, obs)
    avg, std = MCIntegration.sample(config, integrand, measure; print=0, Nblock=16)

    function info(idx, di)
        return @sprintf("   %8.4f ±%8.4f", avg[idx, di], std[idx, di])
    end

    if isnothing(avg) == false
        # avg *= NF
        # std *= NF
        N = size(lgrid)[1]
        grid = lgrid

        for o in 1:Order
            println("Order ", o)
            println("UpUp ver3: ")
            for li in 1:N
                @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], avg[o, li, 1], std[o, li, 1])
            end
            println("UpDown ver3: ")
            for li in 1:N
                @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], avg[o, li, 2], std[o, li, 2])
            end
            println("Total ver3: ")
            for li in 1:N
                @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], sum(avg[o, li, :]), sum(std[o, li, :]))
            end
        end
    end

end

MC()