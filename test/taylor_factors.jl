using Test
using ElectronLiquid
using FeynmanDiagram
using MCIntegration: CompositeVar, Continuous
import MCIntegration.integrate as mcintegrate

# NOTE: Takes a couple of minutes, suggested to run in parallel.
#       A large runtime is necessary to resolve the bug in (3, 2, 0).

function compare(mean, stdev, expect)
    # println(data, ", ", expect)
    @test isapprox(mean, expect, atol=5 * stdev)
end

function test_diagram(; partitions=[(3, 2, 0)])
    DiagTree.uidreset()
    diag = []
    diagpara = []
    valid_partitions = []
    for p in partitions
        # NOTE: Input partitions are defined using the convention (n_loop, n_lambda, n_mu)
        n_loop, n_lambda, n_mu = p
        n_intn_lines = n_loop - 1
        # Build diagram tree for this partition
        diagparam, d = density_diagram(; n_inner_loop=n_intn_lines)
        # Build tree with counterterms (∂λ(∂μ(DT))) via automatic differentiation
        dp = DiagTree.derivative([d], BareGreenId, n_mu; index=1)
        dpp = DiagTree.derivative(dp, BareInteractionId, n_lambda; index=2)
        # The Taylor expansion should be d^n f(x) / dx^n / n!,
        # so there is a factor of 1/n! for each derivative
        for d in dpp
            d.factor *= 1 / factorial(n_mu) / factorial(n_lambda)
        end
        exprtree = ExprTree.build(dpp, diagparam.dim)
        push!(valid_partitions, p)
        push!(diagpara, diagparam)
        push!(diag, exprtree)
    end
    @assert valid_partitions == partitions
    return diagpara, diag
end

function density_diagram(; n_inner_loop=0)
    # Instantaneous Green's function (occupation number) diagram parameters
    DiagTree.uidreset()
    diagparam = DiagPara{Float64}(;
        type=GreenDiag,
        hasTau=true,
        firstLoopIdx=2,
        innerLoopNum=n_inner_loop,
        firstTauIdx=2,
        totalTauNum=n_inner_loop + 1,
        interaction=[FeynmanDiagram.Interaction(ChargeCharge, Instant)],
        filter=[NoHartree]
    )
    k = DiagTree.getK(diagparam.totalLoopNum, 1)
    extT = (1, 1)
    diagtree = Parquet.green(diagparam, k, extT; name=:nₖ)
    return diagparam, diagtree
end

function integrate_density(
    para::ParaMC,
    diagpara,
    diag;
    alpha=3.0,
    neval=1e7,
    print=-1,
    solver=:vegasmc
)
    @assert all(p.totalTauNum ≤ para.order + 1 for p in diagpara)
    @assert all(length(et.root) == 1 for et in diag)

    totalLoopNums = [p.totalLoopNum for p in diagpara]
    maxloops = maximum(totalLoopNums)
    varK = zeros(3, maxloops)
    (K, T) = density_mc_variables(para, alpha)

    # NOTE: We do not integrate the incoming external time
    # and nₖ is instantaneous, hence n_τ = totalTauNum - 1
    dof = [[p.totalLoopNum, p.totalTauNum - 1] for p in diagpara]
    obs = zeros(length(dof))
    T.data[1] = 0  # τin = 0 (= τout⁺)

    # Phase-space factors
    phase_factors = [para.spin / (2π)^(para.dim * nl) for nl in totalLoopNums]
    prefactors = -phase_factors

    return mcintegrate(
        integrand;
        solver=solver,
        measure=measure,
        neval=neval,
        print=print,
        # MC config kwargs
        userdata=(para, diag, totalLoopNums, prefactors, varK),
        var=(K, T),
        dof=dof,
        obs=obs
    )
end

"""Build variable pools for the density."""
function density_mc_variables(para::UEG.ParaMC, alpha::Float64)
    R = Continuous(0.0, 1.0; alpha=alpha)
    Theta = Continuous(0.0, 1π; alpha=alpha)
    Phi = Continuous(0.0, 2π; alpha=alpha)
    K = CompositeVar(R, Theta, Phi)
    # Offset T pool by 1 for fixed external times (instantaneous Green's function ⟹ τin = τout⁺)
    T = Continuous(0.0, para.β; offset=1, alpha=alpha)
    return (K, T)
end

"""Measurement for multiple diagram trees (counterterm partitions)."""
function measure(vars, obs, weights, config)
    for o in 1:(config.N)
        obs[o] += weights[o]
    end
    return
end

"""Integrand for the total density n."""
function integrand(vars, config)
    # We sample internal momentum/times, and external momentum index
    K, T = vars
    R, Theta, Phi = K
    # Unpack userdata
    para, diag, totalLoopNums, prefactors, varK = config.userdata
    # Evaluate the integrand for each partition
    integrand = Vector(undef, config.N)
    for i in 1:(config.N)
        phifactor = 1.0
        for j in 1:totalLoopNums[i]  # config.dof[i][1]
            r = R[j] / (1 - R[j])
            θ = Theta[j]
            ϕ = Phi[j]
            varK[1, j] = r * sin(θ) * cos(ϕ)
            varK[2, j] = r * sin(θ) * sin(ϕ)
            varK[3, j] = r * cos(θ)
            phifactor *= r^2 * sin(θ) / (1 - R[j])^2
        end
        # Evaluate the density integrand n for this partition
        ExprTree.evalKT!(diag[i], varK, T.data, para)  # (additional = para)
        root = diag[i].root[1]
        weight = diag[i].node.current
        integrand[i] = phifactor * prefactors[i] * weight[root]
    end
    return integrand
end

@testset "Taylor factors" begin
    # Test Taylor factors using total density benchmark
    para = ParaMC(;
        order=5,
        rs=1.0,
        beta=40.0,
        mass2=0.6,
        isDynamic=false
    )
    # Testing (3, 2, 0) partition against benchmark, i.e., the density diagram
    # containing 2 interaction lines, 2 λ derivatives, and 0 μ derivatives
    diagpara, diag = test_diagram()
    res = integrate_density(para, diagpara, diag)
    if !isnothing(res)
        # Compare with (3, 2, 0) density benchmark (given in units of n0)
        expect = 0.02309
        n0 = para.kF^3 / (3 * pi^2)
        println(res.mean / n0, " ± ", res.stdev / n0)
        compare(res.mean / n0, res.stdev / n0, expect)
    end
end
