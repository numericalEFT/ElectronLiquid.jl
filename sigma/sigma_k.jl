using Printf, LinearAlgebra
using CompositeGrids
using ElectronGas
using Parameters, Random, DataFrames
using MCIntegration
using Lehmann

using FeynmanDiagram
using StaticArrays
using JLD2

include("../common/interaction.jl")

const steps = 1e6
const Order = 1
# const maxK = 3kF
const Nk = 4
const minK = 0.01kF
const order = 4


# println(dW0)
# exit(0)
const lgrid = [0, 1]
const kgrid = CompositeGrid.LogDensedGrid(:cheb, [0.0, maxK], [0.0, kF], Nk, minK, order)
# println(length(kgrid))
# println(kgrid)
# exit(0)
const Nl = length(lgrid)

println("Build the diagrams into an experssion tree ...")

# partition = [(1, 0, 0),  # order 1
#     (2, 0, 0), (1, 1, 0), (1, 0, 1),  #order 2
#     (3, 0, 0), (2, 1, 0), (2, 0, 1), (1, 1, 1), (1, 2, 0), (1, 0, 2), #order 3
#     (4, 0, 0), (3, 1, 0), (3, 0, 1), (2, 1, 1), (2, 2, 0), (2, 0, 2), (1, 3, 0), (1, 0, 3), (1, 2, 1), (1, 1, 2) #order 4
# ]

# partition = [p for p in sort(partition) if p[1] + p[2] + p[3] <= Order]
const _partition = partition(Order)

println("Diagram set: ", _partition)

diagPara(order) = GenericPara(diagType=SigmaDiag, innerLoopNum=order, hasTau=true, loopDim=dim, spin=spin, firstLoopIdx=2,
    interaction=[FeynmanDiagram.Interaction(ChargeCharge, [
        Instant,
        Dynamic
    ]),],  #instant charge-charge interaction
    filter=[
    # Girreducible,
    # Proper,   #one interaction irreduble diagrams or not
    # NoBubble, #allow the bubble diagram or not
    # NoFock,
    ]
)

sigma = Dict()
for p in _partition
    d = Parquet.sigma(diagPara(p[1])).diagram
    d = DiagTree.derivative(d, BareGreenId, p[2])
    d = DiagTree.derivative(d, BareInteractionId, p[3])
    sigma[p] = d
    # if p == (1, 0, 0)
    #     sigma[p] = d
    # else
    #     sigma[p] = DiagTree.removeHatreeFock!(d)
    # end
end
# DiagTree.removeHatreeFock!(sigma[2, 0, 0])
# println(sigma[2, 0, 0])
# plot_tree(sigma[(2, 0, 0)], maxdepth=8)
# exit(0)

sigma = [sigma[p] for p in _partition]
const diagpara = [diags[1].id.para for diags in sigma]
const diag = [ExprTree.build(diags) for diags in sigma]
const root = [d.root for d in diag] #select the diagram with upup
#assign the external Tau to the corresponding diagrams
const extT = [[diag[ri].node.object[idx].para.extT for idx in r] for (ri, r) in enumerate(root)]
# println(root)
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
    w = sum(diag[order].node.current[r] * phase(varT, extT[order][ri], l) for (ri, r) in enumerate(root[order]))
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
    return w #the current implementation of sigma has an additional minus sign compared to the standard defintion
end

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    l = config.var[3][1]
    # println(config.observable[1][1])
    o = config.curr
    weight = integrand(config)
    config.observable[o, l+1] += weight / abs(weight) * factor
end

function zfactor(avg, std)
    return imag(avg[2] - avg[1]) / (2π / β), (abs(imag(std[2])) + abs(imag(std[1]))) / (2π / β)
    # return imag(avg[1]) / (π / β), imag(std[1]) / (π / β)
end

function chemicalpotential(avg, std)
    return real(avg[1]), abs(std[1])
end

function MC()
    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kF
    # K.data[:, 1] .= [kF, 0.0, 0.0]
    T = MCIntegration.Tau(β, β / 2.0, offset=1)
    T.data[1] = 0.0
    # T = MCIntegration.Tau(1.0, 1.0 / 2.0)
    X = MCIntegration.Discrete(lgrid[1], lgrid[end])

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1] for p in diagpara] # K, T, ExtKidx
    # dof = [[1, 1, 1], [1, 1, 1], [1, 1, 1], [2, 3, 1]]
    # dof = [[1, 0, 1], [1, 0, 1], [1, 0, 1], [2, 1, 1]]
    # println(dof)
    # exit(0)
    obs = zeros(ComplexF64, length(dof), Nl) # observable for the Fock diagram 

    config = MCIntegration.Configuration(steps, (K, T, X), dof, obs)
    avg, std = MCIntegration.sample(config, integrand, measure; print=0, Nblock=64, reweight=10000)

    if isnothing(avg) == false

        jldsave("data.jld2", avg=avg, std=std)

        open("data.dat", "w") do f
            @printf(f, "#%7s %16s %16s %16s %16s\n", "freq", "real", "real err", "imag", "imag err")
            for o in 1:length(dof)
                write(f, "# $(_partition[o])\n")
                for li in 1:Nl
                    @printf(f, "%8.4f %16.8f %16.8f %16.8f %16.8f\n", lgrid[li], real(avg[o, li]), real(std[o, li]), imag(avg[o, li]), imag(std[o, li]))
                end
            end

            # write(f, "\n")

            # nz1, ez1 = zfactor(avg[1, :], std[1, :])
            # nzg1, ezg1 = zfactor(avg[2, :], std[2, :])
            # nzw1, ezw1 = zfactor(avg[3, :], std[3, :])
            # nz2, ez2 = zfactor(avg[4, :], std[4, :])

            # mu1, emu1 = chemicalpotential(avg[1, :], std[1, :])
            # # dmu1 = -mu1
            # dmu1 = 0.0
            # dz1 = -nz1

            # write(f, "# chemical potential: \n")
            # write(f, "# mu1 = $mu1 +- $emu1\n")

            # nmug1, emug1 = chemicalpotential(avg[2, :], std[2, :])
            # nmuw1, emuw1 = chemicalpotential(avg[3, :], std[3, :])
            # nmu2, emu2 = chemicalpotential(avg[4, :], std[4, :])
            # mu2 = nmu2 + nmuw1 + dmu1 * nmug1
            # emu2 = emu2 + emuw1 + emug1 * abs(dmu1) + abs(nmug1) * emu1
            # write(f, "# mu2 = $mu2 +- $emu2\n")


            # write(f, "\n")

            # write(f, "# zfactor: \n")
            # write(f, "# nzg1 = $nzg1 +- $ezg1\n")
            # write(f, "# nzw1 = $nzw1 +- $ezw1\n")
            # write(f, "# nz2 = $nz2 +- $ez2\n")
            # write(f, "# dmu1 = $dmu1 +- $emu1    dz1 = $dz1 +- $ez1\n")
            # write(f, "\n")
            # write(f, "# order 1\n")
            # write(f, "# z1 = $z1  +-  $ez1\n")
            # z2 = nz2 + nzw1 + dmu1 * nzg1
            # ez2 = ez2 + ezw1 + abs(dmu1) * ezg1 + emu1 * abs(nzg1)
            # # z2 = nz2 + nzw1 + dmu1 * nzg1 + dz1 * nz1
            # # dz2 = -z2
            # # ez2 = ez2 + ezw1 + dmu1 * ezg1 + emu1 * nzg1 + 2 * ez1 * nz1
            # write(f, "# order 2\n")
            # write(f, "# z2= $z2 +- $ez2")
            # # write(f, "$nz2  +-  $ez1\n")
            # # println("order 2: ", nz2 + nzg1*mu1 + nzw1 + nz1 * 2 * (-nz1), "  +-  ", ezg1 + ezw1 + ez2 + 2 * abs(z1) * ez1)
        end

    end

end

MC()