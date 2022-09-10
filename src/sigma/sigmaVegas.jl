
function integrandKWV(idx, vars, config)
    # function integrandKW(idx, varK, varT, config)
    # idx = 1
    K, T, T1, N, ExtKidx = vars
    R, Theta, Phi = K
    Tx, Ty = T
    para, diag, extT, kgrid, ngrid, varK, varT = config.userdata
    diagram = diag[idx]
    weight = diagram.node.current
    loopNum = config.dof[idx][1]
    tauNum = (config.dof[idx][2] + 1) * 2
    l = N[1]
    k = ExtKidx[1]

    varK[1, 1] = kgrid[k]

    # varT[1] = Tx.data[1] + Ty.data[1]
    # varT[2] = Tx.data[1] - Ty.data[1]
    varT[1] = 0.0
    varT[2] = T1.data[1]
    varT[3] = Tx.data[1] + Ty.data[1]
    varT[4] = Tx.data[1] - Ty.data[1]

    for i = 3:tauNum
        if varT[i] < -para.β / 2 || varT[i] > para.β / 2
            return 0.0 + 0.0im
        end
    end
    for i = 3:tauNum
        varT[i] += para.β / 2
    end
    for i = 1:tauNum
        @assert 0 <= varT[i] < para.β "varT[$i] = $(varT[i])"
    end
    # for i = 1:tauNum
    #     varT[i] = T.data[i]
    # end

    wn = ngrid[l]

    phifactor = 1.0
    for i = 1:loopNum
        r = R[i] / (1 - R[i])
        θ = Theta[i]
        ϕ = Phi[i]
        # varK[:, i+1] .= [r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ)]
        varK[1, i+1] = r * sin(θ) * cos(ϕ)
        varK[2, i+1] = r * sin(θ) * sin(ϕ)
        varK[3, i+1] = r * cos(θ)
        phifactor *= r^2 * sin(θ) / (1 - R[i])^2
    end

    ExprTree.evalKT!(diagram, varK, varT, para)
    w = sum(weight[r] * phase(varT, extT[idx][ri], wn, para.β) for (ri, r) in enumerate(diagram.root))

    factor = 1.0 / (2π)^(para.dim * loopNum) * phifactor * (2)^(loopNum - 1)
    return w * factor   #the current implementation of sigma has an additional minus sign compared to the standard defintion
    # return real(w) * factor #the current implementation of sigma has an additional minus sign compared to the standard defintion
end

#for vegas algorithm
function integrandKWV(vars, config)
    return integrandKWV(1, vars, config)::ComplexF64
end

function KWV(para::ParaMC, diagram;
    kgrid=[para.kF,],
    ngrid=[0,],
    neval=1e6, #number of evaluations
    print=-1,
    alpha=3.0, #learning ratio
    config=nothing,
    kwargs...
)
    UEG.MCinitialize!(para)

    dim, β, kF = para.dim, para.β, para.kF
    partition, diagpara, diag, root, extT = diagram
    varK = zeros(para.dim, 128)
    varT = zeros(128)

    R = Continuous(0.0, 1.0; alpha=alpha)
    Theta = Continuous(0.0, 1π; alpha=alpha)
    Phi = Continuous(0.0, 2π; alpha=alpha)
    K = CompositeVar(R, Theta, Phi)
    Tx = Continuous(-β / 2, β / 2; alpha=alpha)
    Ty = Continuous(-β / 2, β / 2; alpha=alpha)
    T = CompositeVar(Tx, Ty)
    T1 = Continuous(0.0, β; alpha=alpha)
    X = Discrete(1, length(ngrid), alpha=alpha)
    ExtKidx = Discrete(1, length(kgrid), alpha=alpha)

    # dof = [[p.innerLoopNum, p.innerLoopNum, p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    dof = [[p.innerLoopNum, Int(p.totalTauNum / 2) - 1, 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    # observable of sigma diagram of different permutations

    # if isnothing(neighbor)
    #     neighbor = UEG.neighbor(partition)
    # end
    if isnothing(config)
        config = Configuration(;
            # var=(R, Theta, Phi, T, X, ExtKidx),
            var=(K, T, T1, X, ExtKidx),
            dof=dof,
            type=ComplexF64, # type of the integrand
            # obs=obs,
            userdata=(para, diag, extT, kgrid, ngrid, varK, varT),
            kwargs...
        )
    end

    result = integrate(integrandKWV; config=config, print=print, neval=neval, kwargs...)
    # result = integrate(integrandKW; config=config, print=print, neval=neval, kwargs...)
    # niter=niter, print=print, block=block, kwargs...)

    if isnothing(result) == false

        if print > -1
            report(result.config)
            println(report(result; pick=o -> first(o)))
        end
        if print > -2
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