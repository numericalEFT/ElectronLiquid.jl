
function integrandKWV(idx, R, Theta, Phi, varT, N, ExtKidx, config)
    # function integrandKW(idx, varK, varT, config)
    para, diag, extT, kgrid, ngrid, varK = config.userdata
    diagram = diag[idx]
    weight = diagram.node.current
    l = N[1]
    k = ExtKidx[1]
    varK.data[1, 1] = kgrid[k]
    wn = ngrid[l]

    r = R[1]
    θ = Theta[1]
    ϕ = Phi[1]
    # varK[:, i+1] .= [r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ)]
    varK[1, 2] = r * sin(θ) * cos(ϕ)
    varK[2, 2] = r * sin(θ) * sin(ϕ)
    varK[3, 2] = r * cos(θ)
    phifactor *= r^2 * sin(θ)

    ExprTree.evalKT!(diagram, varK.data, varT.data, para)
    w = sum(weight[r] * phase(varT, extT[idx][ri], wn, para.β) for (ri, r) in enumerate(diagram.root))

    loopNum = config.dof[idx][1]
    factor = 1.0 / (2π)^(para.dim * loopNum) * phifactor
    return w * factor #the current implementation of sigma has an additional minus sign compared to the standard defintion
    # return real(w) * factor #the current implementation of sigma has an additional minus sign compared to the standard defintion
end

#for vegas algorithm
function integrandKWV(R, Theta, Phi, varT, N, ExtKidx, config)
    return integrandKWV(1, R, Theta, Phi, varT, N, ExtKidx, config)
end

function KWV(para::ParaMC, diagram;
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
    varK = zeros(para.dim, 128)

    R = MCIntegration.Continuous(0.0, 10.0 * kF; grid=collect(LinRange(0.0, β, 1000)), alpha=alpha)
    Theta = MCIntegration.Continuous(0.0, 1π; grid=collect(LinRange(0.0, β, 1000)), alpha=alpha)
    Phi = MCIntegration.Continuous(0.0, 2π; grid=collect(LinRange(0.0, β, 1000)), alpha=alpha)
    T = MCIntegration.Continuous(0.0, β; grid=collect(LinRange(0.0, β, 1000)), offset=1, alpha=alpha)
    T.data[1] = 0.0
    X = MCIntegration.Discrete(1, length(ngrid), alpha=alpha)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid), alpha=alpha)

    dof = [[p.innerLoopNum, p.innerLoopNum, p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    # observable of sigma diagram of different permutations

    # if isnothing(neighbor)
    #     neighbor = UEG.neighbor(partition)
    # end
    if isnothing(config)
        config = Configuration(;
            var=(R, Theta, Phi, T, X, ExtKidx),
            dof=dof,
            type=ComplexF64, # type of the integrand
            # obs=obs,
            userdata=(para, diag, extT, kgrid, ngrid, varK),
            kwargs...
            # neighbor=neighbor,
            # reweight_goal=reweight_goal, kwargs...
        )
    end

    result = integrate(integrandKWV; config=config, print=print, neval=neval, kwargs...)
    # result = integrate(integrandKW; config=config, print=print, neval=neval, kwargs...)
    # niter=niter, print=print, block=block, kwargs...)

    if isnothing(result) == false

        if print >= 0
            report(result.config)
            println(report(result, o -> first(o)))
            println(result)
        end

        # datadict = Dict{eltype(partition),Complex{Measurement{Float64}}}()
        datadict = Dict{eltype(partition),Any}()

        for o in 1:length(dof)
            avg, std = result.mean[o], result.stdev[o]
            r = measurement.(real(avg), real(std))
            i = measurement.(imag(avg), imag(std))
            data = Complex.(r, i)
            datadict[partition[o]] = data
        end
        return datadict, result
    else
        return nothing, nothing
    end
end