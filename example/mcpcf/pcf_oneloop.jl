using MCIntegration
using LegendrePolynomials
using LinearAlgebra
using StaticArrays

include("propagators.jl")
using .Propagators
using .Propagators: G0, interaction, response

const iscross = false
const fname = "run/data/PCFdata_30003.jld2"
const steps = 1e7 # 2e8/hr
const ℓ = 0
const θ, rs = 0.1, 3.0
const param = Propagators.Parameter.rydbergUnit(θ, rs, 3)
const α = 0.8
println(param)

function integrand(vars, config)

    ExtT, ExtK, X, T, P, Q = vars
    funcs = config.userdata
    Rt = funcs.Rt
    extT, extK = Rt.mesh[1], Rt.mesh[2]
    param = funcs.param

    t = extT[ExtT[1]]
    k = extK[ExtK[1]]
    x = X[1]
    t3, t4, t1, t2 = T[1], T[2], T[3], T[4]
    p = P[1]
    Qv = SVector{3,Float64}(Q[1], Q[2], Q[3])

    PLX = Pl(x, ℓ)
    q = sqrt(k^2 + p^2 + 2 * k * p * x)
    V = 1.0 / interaction(q, funcs)
    # V = coulomb(q, funcs)[1]
    G1 = G0(t3, p, funcs)
    G021 = G0(t3, -p, funcs)
    G022 = G0(t4, -p, funcs)
    R0 = response(p, funcs) / param.β
    R = response(t3 - t4, p, funcs)
    # R0 = 1.0 / param.β
    # R = 0.0

    result1 = -p^2 / (4π^2) * PLX * V * G1 * (G021 * R0 + G022 * R)

    W = interaction(t, q, funcs) * V
    # W = 0.0
    G21 = G0(t3 - t, -p, funcs)
    G22 = G0(t4 - t, -p, funcs)

    result2 = -p^2 / (4π^2) * PLX * W * G1 * (G21 * R0 + G22 * R)

    if iscross
        Kv = SVector{3,Float64}(k, 0, 0)
        Pv = SVector{3,Float64}(p * x, p * sqrt(1 - x^2), 0)
        kvmq = norm(Kv - Qv)
        pvmq = norm(Pv - Qv)
        V1 = 1.0 / interaction(kvmq, funcs)
        V2 = 1.0 / interaction(pvmq, funcs)
        W1 = interaction(t2, kvmq, funcs) * V1
        W2 = interaction(t - t1, pvmq, funcs) * V2
        V1 = V1 / param.β
        V2 = V2 / param.β

        K1, K2, K3, K4 = norm(Qv), norm(Qv - Pv - Kv), p, -p
        # 1,3 cares if 2 is ins; 2,4 cares if 1 is ins
        # t1--t, t2--0
        Gi1, Gi2 = G0(t, K1, funcs), G0(-t, K2, funcs)
        Gi3, Gi4 = G0(t3 - t, K3, funcs), G0(t4, K4, funcs)
        Gi04 = G0(t3, K4, funcs) # for R0
        Gd1, Gd2 = G0(t1, K1, funcs), G0(t2 - t, K2, funcs)
        Gd3, Gd4 = G0(t3 - t1, K3, funcs), G0(t4 - t2, K4, funcs)
        Gd04 = G0(t3 - t2, K4, funcs) # for R0

        result3 = -p^2 / (2π)^5 * PLX * (
                      V1 * V2 * Gi1 * Gi2 * Gi3 * (Gi4 * R + Gi04 * R0)
                      +
                      W1 * V2 * Gi1 * Gd2 * Gi3 * (Gd4 * R + Gd04 * R0)
                      +
                      V1 * W2 * Gd1 * Gi2 * Gd3 * (Gi4 * R + Gi04 * R0)
                      +
                      W1 * W2 * Gd1 * Gd2 * Gd3 * (Gd4 * R + Gd04 * R0)
                  )
    else
        result3 = 1e-22
    end
    return 1.0, result1, result2, result3
end

function measure(vars, obs, weight, config)
    extt, extk = vars[1], vars[2]
    obs[1][1, extk[1]] += weight[1]
    obs[2][1, extk[1]] += weight[2]
    obs[3][extt[1], extk[1]] += weight[3]
    obs[4][extt[1], extk[1]] += weight[4]
end

function run(steps, param, alg=:vegas)
    println("Prepare propagators")

    mint = 0.001 * param.β
    minK, maxK = 0.001 * sqrt(param.T * param.me), 10param.kF
    order = 6
    rpai, rpat = Propagators.rpa(param; mint=mint, minK=minK, maxK=maxK, order=order)

    mint = 0.001 * param.β
    minK, maxK = 0.1 * sqrt(param.T * param.me), 10param.kF
    order = 3
    Ri, Rt, Rtd = Propagators.loadR(fname, param; mint=mint, minK=minK, maxK=maxK, order=order)
    # Ri, Rt = Propagators.initR(param; mint=mint, minK=minK, maxK=maxK, order=order)
    println(size(Ri))
    println(size(Rt))
    println(size(Rtd))

    # userdata
    funcs = Propagators.Funcs(param, rpai, rpat, Ri, Rt, Rtd)

    println("Prepare variables")
    extT, extK = Rt.mesh[1], Rt.mesh[2]
    T = Continuous(0.0, param.β; alpha=3.0, adapt=true)
    P = Continuous(0.0, maxK; alpha=3.0, adapt=true)
    X = Continuous(-1.0, 1.0; alpha=3.0, adapt=true)
    # Q = FermiK(3, param.kF, 0.2 * param.kF, 10.0 * param.kF)
    Q = Continuous(-maxK, maxK; alpha=3.0, adapt=true)

    ExtT = Discrete(1, length(extT); adapt=false)
    ExtK = Discrete(1, length(extK); adapt=false)

    # ExtT, ExtK, X, T, P, Q
    dof = [
        [0, 1, 0, 0, 0, 0],
        [0, 1, 1, 2, 1, 0],
        [1, 1, 1, 2, 1, 0],
        [1, 1, 1, 4, 1, 3]
    ]
    obs = [zeros(ComplexF64, size(Rt)), zeros(ComplexF64, size(Rt)), zeros(ComplexF64, size(Rt)), zeros(ComplexF64, size(Rt))]

    println("Start")
    result = integrate(integrand; measure=measure, userdata=funcs,
        var=(ExtT, ExtK, X, T, P, Q), dof=dof, obs=obs, solver=alg,
        neval=steps, print=-1, block=8, type=ComplexF64)
    # println(result.mean)
    funcs.Ri.data .= result.mean[1][1, :] .+ result.mean[2][1, :]
    funcs.Rt.data .= result.mean[3] .+ result.mean[4]

    return result, funcs
end

if abspath(PROGRAM_FILE) == @__FILE__
    result, funcs = run(steps, param, :vegasmc)
    Ri, Rt = funcs.Ri, funcs.Rt
    println(sum(abs.(result.mean[4]) / sum(abs.(result.mean[3]))))
    # println(Ri.mesh[1])
    println(real(Ri.data))
    println(real(Rt.data[1, :]))
    # println(result[1][1])
    println("R0=$(real(Propagators.R0(Ri, Rt, param)))")
end