"""
    struct GenericPara

The parameters for the sigma.

# Members
- `para`: the parameters for the MC integration
- `kamp`: the amplitude of the external momentum
- `ωn`: vector of the frequency of the external leg
"""
struct GenericPara
    para::ParaMC
    kamp::Float64
    ωn::Vector{Int} #allow measure multiple frequency simultaneously
    function GenericPara(para, kamp, ωn)
        return new(para, kamp, ωn)
    end
end

function integrand_generic(idx, vars, config)
    # function integrandKW(idx, varK, varT, config)
    varK, varT, varX = vars
    diag, extT, paras = config.userdata
    diagram = diag[idx]
    weight = diagram.node.current
    Xidx = varX[1]
    para = paras[Xidx][1]
    varK.data[1, 1] = paras[Xidx][2] # momentum
    wn = paras[Xidx][3] # Matsubara frequency (integer)

    ExprTree.evalKT!(diagram, varK.data, varT.data, para)
    w = sum(weight[r] * phase(varT, extT[idx][ri], wn, para.β) for (ri, r) in enumerate(diagram.root))

    loopNum = config.dof[idx][1]
    factor = 1.0 / (2π)^(para.dim * loopNum)
    return w * factor #the current implementation of sigma has an additional minus sign compared to the standard defintion
    # return real(w) * factor #the current implementation of sigma has an additional minus sign compared to the standard defintion
end

#for vegas algorithm
# function integrandKW(varK, varT, N, ExtKidx, config)
#     return integrandKW(1, varK, varT, N, ExtKidx, config)
# end

function measure_generic(idx, vars, obs, weight, config)
    # l = vars[3][1]  #matsubara frequency
    # k = vars[4][1]  #K
    varK, varT, varX = vars
    Xidx = varX[1]
    obs[idx][Xidx] += weight
end

#for vegas algorithm
# function measureKW(obs, weights, config)
#     l = config.var[3][1]  #matsubara frequency
#     k = config.var[4][1]  #K
#     for o in 1:config.N
#         obs[o][l, k] += weights[o]
#     end
# end

"""
    function Generic(paras::AbstractArray{Tuple{ParaMC,Float64,Int}}, diagram;
        neval=1e6, #number of evaluations
        print=0,
        alpha=3.0, #learning ratio
        config=nothing,
        solver=:mcmc,
        kwargs...
    )

Calculate self-energy with generic set of parameters.

# Arguments
- `paras`: array of tuples of parameters, momentum, and Matsubara frequency
- `diagram`: diagram to be calculated
- `neval`: number of evaluations
- `print`: print level
- `alpha`: learning ratio
- `config`: configuration
- `solver`: solver
- `kwargs`: keyword arguments
"""
function Generic(paras::AbstractArray{Tuple{ParaMC,Float64,Int}}, diagram;
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    kwargs...
)
    # if haskey(kwargs, :solver)
    # @assert kwargs[:solver] == :mcmc "Only :mcmc is supported for Sigma.KW"
    # end
    @assert solver == :mcmc "Only :mcmc is supported for Sigma.Generic"
    dim, β, kF = paras[1].para.dim, paras[1].para.β, paras[1].para.kF
    order = paras[1].para.order
    for p in paras
        para = p.para
        para.isDynamic && UEG.MCinitialize!(para)
        @assert para.dim == dim "All parameters must have the same dimension"
        @assert para.β ≈ β "All parameters must have the same inverse temperature"
        @assert para.kF ≈ kF "All parameters must have the same Fermi momentum"
        @assert para.order == order "All parameters must have the same order"
    end

    partition, diagpara, diag, root, extT = diagram

    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kF
    # T = MCIntegration.Tau(β, β / 2.0, offset=1)
    T = MCIntegration.Continuous(0.0, β; grid=collect(LinRange(0.0, β, 1000)), offset=1, alpha=alpha)
    T.data[1] = 0.0
    X = MCIntegration.Discrete(1, length(paras), alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1] for p in diagpara] # K, T, ExtKidx
    # observable of sigma diagram of different permutations
    obs = [zeros(ComplexF64, size(paras)) for o in 1:length(dof)]

    # if isnothing(neighbor)
    #     neighbor = UEG.neighbor(partition)
    # end
    if isnothing(config)
        config = Configuration(;
            var=(K, T, X),
            dof=dof,
            type=ComplexF64, # type of the integrand
            obs=obs,
            userdata=(diag, extT, paras),
            kwargs...
            # neighbor=neighbor,
            # reweight_goal=reweight_goal, kwargs...
        )
    end

    result = integrate(integrand_generic; config=config, measure=measure_generic, print=print, neval=neval, solver=solver, kwargs...)
    # result = integrate(integrandKW; config=config, print=print, neval=neval, kwargs...)
    # niter=niter, print=print, block=block, kwargs...)

    if isnothing(result) == false

        if print >= 0
            report(result.config)
            println(report(result, pick=o -> first(o)))
            println(result)
        end

        if print >= -2
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


function MC(paras::AbstractArray{Tuple{ParaMC,Float64,Int}}, filename=nothing;
    diagtype=:Parquet,
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    dR=false, # take interaction derivative or not
    reweight_goal::Union{Nothing,AbstractVector{Float64}}=nothing,
    parallel=:nothread,
    kwargs...
)

    dim, β, kF = paras[1][1].dim, paras[1][1].β, paras[1][1].kF
    order = paras[1][1].order

    partition = UEG.partition(_order)
    neighbor = UEG.neighbor(partition)

    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder, vOrder) in partition
            push!(reweight_goal, 4.0^(order + vOrder - 1))
        end
        push!(reweight_goal, 2.0)
    end

    if diagtyep == :Parquet
        diagram = Sigma.diagram(para, partition; dR=dR)
        sigma, result = Sigma.Generic(paras, diagram;
            neighbor=neighbor, reweight_goal=reweight_goal,
            neval=neval, parallel=parallel)
    else
        error("unknown diagrams' generated $diagtype")
    end

    if isnothing(sigma) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (ngrid, kgrid, sigma)
            end
        end
        for (ip, key) in enumerate(partition)
            println("Group ", key)
            @printf("%10s  %10s   %10s   %10s   %10s \n", "q/kF", "real(avg)", "err", "imag(avg)", "err")
            r, i = real(sigma[key]), imag(sigma[key])
            for (in, n) in enumerate(ngrid)
                println("n = $n")
                for (iq, q) in enumerate(kgrid)
                    @printf("%10.6f  %10.6f ± %10.6f   %10.6f ± %10.6f\n", q[1] / kF, r[in, iq].val, r[in, iq].err, i[in, iq].val, i[in, iq].err)
                end
            end
        end
    end

end