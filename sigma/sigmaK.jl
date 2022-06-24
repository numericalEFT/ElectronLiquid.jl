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

const steps = 1e7
const Order = 2

const lgrid = [0, 1]
const Nl = length(lgrid)

const Nk, korder = 4, 4
const minK = 0.2kF

const kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 3kF], [kF,], Nk, minK, korder)
println(kgrid.grid)
println(length(kgrid.grid))
# exit(0)

include("sigma_diagram.jl")

ret = sigmaDiag(Order)
const diagpara = ret[1]
const diag = ret[2]
const root = ret[3]
const extT = ret[4]

function integrand(config)
    order = config.curr
    l = config.var[3][1]
    k = config.var[4][1]
    varK, varT = config.var[1], config.var[2]
    varK.data[1, 1] = kgrid.grid[k]

    ExprTree.evalNaive!(diag[order], varK.data, varT.data, eval)
    w = sum(diag[order].node.current[r] * phase(varT, extT[order][ri], l) for (ri, r) in enumerate(root[order]))
    return w #the current implementation of sigma has an additional minus sign compared to the standard defintion
end

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    l = config.var[3][1]
    k = config.var[4][1]
    # println(config.observable[1][1])
    o = config.curr
    weight = integrand(config)
    config.observable[o, l+1, k] += weight / abs(weight) * factor
end

function MC()
    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kF
    T = MCIntegration.Tau(β, β / 2.0, offset=1)
    T.data[1] = 0.0
    X = MCIntegration.Discrete(lgrid[1], lgrid[end])
    ExtKidx = MCIntegration.Discrete(1, length(kgrid))

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = zeros(ComplexF64, length(dof), Nl, length(kgrid)) # observable for the Fock diagram 

    config = MCIntegration.Configuration(steps, (K, T, X, ExtKidx), dof, obs)
    avg, std = MCIntegration.sample(config, integrand, measure; print=0, Nblock=64, reweight=10000)

    if isnothing(avg) == false
        jldsave("dataK.jld2", order=Order, partition=partition(Order), avg=avg, std=std, kgrid=kgrid)
    end

end

MC()