using Printf, LinearAlgebra
using CompositeGrids
using ElectronGas
using Parameters, Random, DataFrames
using MCIntegration
using Lehmann

using FeynmanDiagram
using StaticArrays
using JLD2

include("./common.jl")

const steps = 4e11

# println(kgrid.grid)
# println(length(kgrid.grid))
# exit(0)

function integrand(config)
    order = config.curr
    l = config.var[3][1]
    k = config.var[4][1]
    varK, varT = config.var[1], config.var[2]
    kgrid = config.para.additional[1]
    varK.data[1, 1] = kgrid.grid[k]

    ExprTree.evalKT!(diag[order], varK.data, varT.data, config.para)
    w = sum(diag[order].node.current[r] * phase(varT, extT[order][ri], l, config.para.β) for (ri, r) in enumerate(root[order]))
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

function MC(para::ParaMC)
    dim, β, kF = para.dim, para.β, para.kF

    ############## external kgrid ###################
    Nk, korder = 4, 4
    minK = 0.2kF
    kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 3kF], [kF,], Nk, minK, korder)
    push!(para.additional, kgrid) # add kgrid to the parameter
    #################################################

    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kF
    T = MCIntegration.Tau(β, β / 2.0, offset=1)
    T.data[1] = 0.0
    X = MCIntegration.Discrete(lgrid[1], lgrid[end])
    ExtKidx = MCIntegration.Discrete(1, length(kgrid))

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = zeros(ComplexF64, length(dof), Nl, length(kgrid)) # observable for the Fock diagram 

    ngb = UEG.neighbor(UEG.partition(Order))
    config = MCIntegration.Configuration((K, T, X, ExtKidx), dof, obs; neighbor=ngb, para=para, reweight_goal=[1.0, 1.0, 1.0, 2.0, 2.0])
    result = MCIntegration.sample(config, integrand, measure; neval=steps, niter=10, print=0, block=16)

    if isnothing(result) == false
        #     jldsave("dataK.jld2", order=Order, partition=partition(Order), avg=avg, std=std, kgrid=kgrid)
        avg, std = result.mean, result.stdev
        println(MCIntegration.summary(result, [o -> real(o[i]) for i in 1:length(dof)]))

        jldopen("dataK.jld2", "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (para, kgrid, avg, std)
        end
    end

end

p = ParaMC(rs=5.0, beta=25.0, Fs=-0.585, order=Order, mass2=0.001)
MC(p)
#p = ParaMC(rs=5.0, beta=100.0, Fs=-0.0, order=Order, mass2=1e-5)
#MC(p)
# p = ParaMC(rs=5.0, beta=100.0, Fs=-1.0, order=Order, mass2=1e-5)
#MC(p)
