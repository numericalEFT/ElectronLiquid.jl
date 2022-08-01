"""
Calculate exchange vertex4 averged on the Fermi surface
"""
# const lgrid = [1, 2]
# const Nl = length(lgrid)

function integrand_l(config)
    kF, β = config.para.kF, config.para.β
    order = config.curr
    x = config.var[3][1]
    varK, varT = config.var[1], config.var[2]

    varK.data[:, 2] = [kF * x, kF * sqrt(1 - x^2), 0.0]

    ExprTree.evalKT!(diag[order], varK.data, varT.data, config.para)
    if !isempty(rootuu[order])
        wuu = sum(diag[order].node.current[root] * phase(varT, extTuu[order][ri], β, isF) for (ri, root) in enumerate(rootuu[order]))
    else
        wuu = 0.0
    end
    if !isempty(rootud[order])
        wud = sum(diag[order].node.current[root] * phase(varT, extTud[order][ri], β, isF) for (ri, root) in enumerate(rootud[order]))
    else
        wud = 0.0
    end
    # println(wuu, ",  ", wud)
    return Weight(wuu / β, wud / β)
end

function measure_l(config)
    factor = 1.0 / config.reweight[config.curr]
    x = config.var[3][1]
    # println(config.observable[1][1])
    if config.curr == 1
        weight = integrand(config)
        config.observable[1, 1] += weight.d / 2 / abs(weight) * factor
        config.observable[1, 2] += weight.e / 2 / abs(weight) * factor
        config.observable[2, 1] += weight.d * x / 2 / abs(weight) * factor
        config.observable[2, 2] += weight.e * x / 2 / abs(weight) * factor
    else
        return
    end
end

function MC(para::ParaMC)
    UEG.MCinitialize!(para)
    dim, β, kF, NF = para.dim, para.β, para.kF, para.NF
    K = MCIntegration.FermiK(para.dim, kF, 0.2 * kF, 10.0 * kF, offset=2)
    K.data[:, 1] .= [kF, 0.0, 0.0]
    T = MCIntegration.Tau(β, β / 2.0)
    X = MCIntegration.Continuous(-1.0, 1.0) #x=cos(θ)

    dof = [[diagpara[o].innerLoopNum, diagpara[o].totalTauNum, 1] for o in 1:para.order] # K, T, ExtKidx
    obs = zeros(Nl, 2) # observable for the Fock diagram 

    config = MCIntegration.Configuration((K, T, X), dof, obs; para=para)
    result = MCIntegration.sample(config, integrand, measure; neval=steps, niter=10, print=0, block=16)

    function info(idx, di)
        return @sprintf("   %8.4f ±%8.4f", avg[idx, di], std[idx, di])
    end

    if isnothing(result) == false
        avg, std = result.mean, result.stdev
        avg *= para.NFstar
        std *= para.NFstar
        N = size(avg)[1]
        grid = lgrid

        println("UpUp ver4: ")
        for li in 1:N
            @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], avg[li, 1], std[li, 1])
        end
        println("UpDown ver4: ")
        for li in 1:N
            @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], avg[li, 2], std[li, 2])
        end

        println("S ver4: ")
        for li in 1:N
            @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], (avg[li, 1] + avg[li, 2]) / 2, (std[li, 1] + std[li, 2]) / 2)
        end
        println("A ver4: ")
        for li in 1:N
            @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], (avg[li, 1] - avg[li, 2]) / 2, (std[li, 1] + std[li, 2]) / 2)
        end

        jldopen("dataF.jld2", "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (para, avg, std)
        end
    end

end

# p = ParaMC(rs=5.0, beta=25.0, Fs=-0.585, order=Order, mass2=0.001)
# MC(p)