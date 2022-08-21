using Revise
using MCIntegration
using FeynmanDiagram
using ElectronGas
using Lehmann
using ElectronLiquid
using LinearAlgebra
using Measurements
using CompositeGrids
using Plots
using LaTeXStrings


const kpool = zeros(3, 2)

##################### propagator and interaction evaluation ##############
function DiagTree.eval(id::BareGreenId, K, extT, varT, p::ParaMC)
    kF, β, me, μ, massratio = p.kF, p.β, p.me, p.μ, p.massratio
    τin, τout = varT[id.extT[1]], varT[id.extT[2]]
    k = norm(K)
    ϵ = k^2 / (2me * massratio) - μ
    # ϵ = kF / me * (k - kF)
    # vF = kF / me
    # ϵ = vF * (K[1] - kF) + (K[2]^2 + K[3]^2) / (2me)

    # if k < 0.4 * kF || k > kF * 1.3
    #     return 0.0
    # end

    τ = τout - τin
    if τ ≈ 0.0
        return Spectral.kernelFermiT(-1e-8, ϵ, β)
    else
        return Spectral.kernelFermiT(τ, ϵ, β)
    end
end

function DiagTree.eval(id::BareInteractionId, K, extT, varT, p::ParaMC)
    dim, e0, ϵ0, mass2 = p.dim, p.e0, p.ϵ0, p.mass2
    qd = sqrt(dot(K, K))
    # qd = sqrt(K[2]^2 + K[3]^2)
    # qd = sqrt(K[1]^2)
    if id.type == Instant
        if interactionTauNum(id.para) == 1
            # return e0^2 / ϵ0 / (dot(K, K) + mass2)
            return UEG.Coulombinstant(qd, p)
        elseif interactionTauNum(id.para) == 2
            # println(id.extT)
            return UEG.interactionStatic(p, qd, varT[id.extT[1]], varT[id.extT[2]])
        else
            error("not implemented!")
        end
    elseif id.type == Dynamic
        return UEG.interactionDynamic(p, qd, varT[id.extT[1]], varT[id.extT[2]])
    else
        error("not implemented!")
    end
end


@inline function phase(varT, extT, l, β)
    tin, tout = varT[extT[1]], varT[extT[2]]
    return exp(1im * π * (2l + 1) / β * (tout - tin))
end

function integrandKW(config)
    para, diag, qgrid, kgrid, ngrid = config.para
    diagram = diag[config.curr]
    weight = diagram.node.current
    object = diagram.node.object
    theta, phi = config.var[1][1], config.var[2][1]
    varT = config.var[3]
    # varK.data[1, 1] = kgrid[k]
    wn = ngrid[config.var[4][1]]
    k = kgrid[config.var[5][1]]
    q = qgrid.grid[config.var[6][1]]
    # q = config.var[6][1]

    kpool[1, 1] = k
    kpool[:, 2] .= [q * sin(theta) * cos(phi) + k, q * sin(theta) * sin(phi), q * cos(theta)]

    ExprTree.evalKT!(diagram, kpool, varT.data, para)
    w = sum(weight[r] * phase(varT, object[r].para.extT, wn, para.β) for r in diagram.root)

    factor = 1.0 / (2π)^3 * q^2 * sin(theta)
    return w * factor #the current implementation of sigma has an additional minus sign compared to the standard defintion
end

function measureKW(config)
    factor = 1.0 / config.reweight[config.curr]
    l = config.var[4][1]
    k = config.var[5][1]
    q = config.var[6][1]
    # println(config.observable[1][1])
    o = config.curr
    weight = integrandKW(config)
    config.observable[o, l, k, q] += weight / abs(weight) * factor
end

function sigmaKW(para::ParaMC, diagram, qgrid;
    kgrid=[para.kF,],
    ngrid=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    kwargs...
)
    UEG.MCinitialize!(para)

    dim, β, kF = para.dim, para.β, para.kF
    partition, diagpara, diag, root, extT = diagram

    # K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF, offset=1)
    # K.data[:, 1] .= 0.0
    # K.data[1, 1] = kF
    Theta = MCIntegration.Continuous(0.0, 1.0 * π, alpha=alpha)
    Phi = MCIntegration.Continuous(0.0, 2.0 * π, alpha=alpha)
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha)
    Q = MCIntegration.Discrete(1, length(qgrid), alpha=alpha)
    # Q = MCIntegration.Continuous(0.0, 6 * para.kF, alpha=alpha)
    T.data[1] = 0.0
    X = MCIntegration.Discrete(1, length(ngrid), alpha=alpha)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid), alpha=alpha)

    dof = [[1, 1, p.totalTauNum - 1, 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = zeros(ComplexF64, length(dof), length(ngrid), length(kgrid), length(qgrid)) # observable for the Fock diagram 

    # if isnothing(neighbor)
    #     neighbor = UEG.neighbor(partition)
    # end
    if isnothing(config)
        config = MCIntegration.Configuration((Theta, Phi, T, X, ExtKidx, Q), dof, obs;
            para=(para, diag, qgrid, kgrid, ngrid),
            kwargs...
            # neighbor=neighbor,
            # reweight_goal=reweight_goal, kwargs...
        )
    end

    result = MCIntegration.sample(config, integrandKW, measureKW; neval=neval, kwargs...)
    # niter=niter, print=print, block=block, kwargs...)

    if isnothing(result) == false

        avg, std = result.mean, result.stdev
        if print >= 0
            MCIntegration.summary(result.config)
            println(MCIntegration.summary(result, [o -> real(o[i]) for i in 1:length(dof)]))
        end

        r = measurement.(real(avg), real(std))
        i = measurement.(imag(avg), imag(std))
        data = Complex.(r, i)
        datadict = Dict{eltype(partition),typeof(data[1, :, :, :])}()
        for i in 1:length(dof)
            datadict[partition[i]] = data[i, :, :, :]
        end
        return datadict, result
    else
        return nothing, nothing
    end
end

const para = UEG.ParaMC(beta=25.0, rs=5.0, mass2=0.00001, isDynamic=true)
const kF = para.kF

const diag = Sigma.diagram(para, [(1, 0, 0),])
# ExprTree.showTree(diag[3][1], diag[3][1].root[1])
# ExprTree.showTree(diag[3][1], diag[3][1].root[2])

const kgrid = [0.0, 0.5kF, kF]

Nk, korder = 4, 4
minK = 0.1 * para.kF
const qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 3kF], [0.0, kF,], Nk, minK, korder)
# const qgrid = para.qgrid

data, result = sigmaKW(para, diag, qgrid; neval=1e6, kgrid=kgrid, ngrid=[-1, 0])

let
    val = data[(1, 0, 0)]
    plot()
    for (ki, k) in enumerate(kgrid)
        println("k/kF=$(k/para.kF)")
        zz = Vector{Measurement{Float64}}()
        for (qi, q) in enumerate(qgrid)
            z = imag(val[2, ki, qi] - val[1, ki, qi]) / (2π / para.β)
            # push!(zz, z.val)
            push!(zz, z)
            println("$(q/para.kF)    $z")
        end
        sw = CompositeGrids.Interp.integrate1D(zz, qgrid)
        println("integration:", sw)
        p = plot!(qgrid.grid / para.kF, zz, label=L"$k/kF = %$(k/para.kF)$, $\int \frac{d Σ_q(k)}{dω}dq  = %$sw$", xlabel=L"$q/kF$", ylabel=L"$\frac{d Σ_q(k)}{dω}$", legend=:topright)
        display(p)
    end
    savefig("sigma_locality.png")
end
