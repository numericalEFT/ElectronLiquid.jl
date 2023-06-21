using MCIntegration
using LegendrePolynomials
using LinearAlgebra
using StaticArrays

include("new_propagators.jl")
using .Propagators
using .Propagators: G0, interaction, response
using .Propagators: diagram_gen, F2F!, add_source
using .Propagators.FeynmanDiagram

using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--cross-diagram", "-x"
        help = "2nd order cross diagram"
        action = :store_true

        "--halfcross-diagram", "-y"
        help = "2nd order half-cross diagram"
        action = :store_true

        "--uid", "-u"
        help = "uid of the job"
        arg_type = Int
        default = 3006

        "--channel", "-l"
        help = "orbital channel"
        arg_type = Int
        default = 0

        "--steps", "-s"
        help = "monte carlo step"
        arg_type = Float64
        default = 1e6

        "--temperature", "-t"
        help = "temperature"
        arg_type = Float64
        default = 0.01

        "--rs", "-r"
        help = "Wigner-Seitz radius"
        arg_type = Float64
        default = 0.3
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
println("Parsed args:")
for (arg, val) in parsed_args
    println("  $arg  =>  $val")
end

const iscross::Bool = parsed_args["cross-diagram"]
const ishalfcross::Bool = parsed_args["halfcross-diagram"]
const uid::Int = parsed_args["uid"]
const fname = "run/data/PCFdata_$uid.jld2"
const savefname = "run/data/mcpcfL1O2X_$uid.jld2"
const steps::Float64 = parsed_args["steps"] # 2e8/hr
const ℓ::Int = parsed_args["channel"]
const θ::Float64, rs::Float64 = parsed_args["temperature"], parsed_args["rs"]
const param = Propagators.Parameter.rydbergUnit(θ, rs, 3)
const α = 0.8
const order = 2
println(param)

function integrand(vars, config)

    ExtT, ExtK, X, T, P, Q = vars
    funcs, diagram, paramc = config.userdata
    dpartition, diagpara, diag, droot, dextT = diagram
    Rt = funcs.Rt
    extT, extK = Rt.mesh[1], Rt.mesh[2]
    param = funcs.param

    t = extT[ExtT[1]]
    k = extK[ExtK[1]]
    x = X[1]
    t1, t2, t3, t4 = T[1], T[2], T[3], T[4]
    p = P[1]
    Qv = SVector{3,Float64}(Q[1], Q[2], Q[3])

    PLX = Pl(x, ℓ)
    # q = sqrt(k^2 + p^2 - 2 * k * p * x)

    K = SMatrix{3,4,Float64}([-p*x -p*sqrt(1 - x^2) 0; -k 0 0; k 0 0; Q[1] Q[2] Q[3]]')
    T = SVector{4,Float64}(t1, t2, t3, t4)

    ExprTree.evalKT!(diag[1], K, T, paramc)
    ExprTree.evalKT!(diag[2], K, T, paramc)

    result = [0.0, 0.0]

    for dorder in 1:1
        weight = diag[dorder].node.current
        for idx in 1:1
            extTu = dextT[idx][dorder]
            factor = -1.0 * p^2 / (2π)^2
            for (ri, r) in enumerate(droot[idx][dorder])
                w = weight[r]
                extT = extTu[ri]
                tInL, tOutL, tInR, tOutR = T[extT[1]], T[extT[2]], T[extT[3]], T[extT[4]]
                F = responsef(tInR - tInL, p, funcs)
                G1 = G0(tOutL - 0, k, funcs)
                G2 = G0(tOutR - t, -k, funcs)
                result[dorder] += factor * G1 * G2 * w * F
            end
        end
    end
    # result ./= 2.0
    return 1.0, result[1], result[2]
end

function measure(vars, obs, weight, config)
    extt, extk = vars[1], vars[2]
    obs[1][1, extk[1]] += weight[1]
    obs[2][extt[1], extk[1]] += weight[2]
    obs[3][extt[1], extk[1]] += weight[3]
end

function run(steps, param, alg=:vegas; order=order)
    println("Prepare propagators")

    paramc, diagram = diagram_gen(param.rs, param.beta; order=order)
    println(paramc)
    # dpartition, diagpara, diag, droot, dextT = diagram

    mint = 0.001 * param.β
    minK, maxK = 0.001 * sqrt(param.T * param.me), 10param.kF
    order = 6
    rpai, rpat = Propagators.rpa(param; mint=mint, minK=minK, maxK=maxK, order=order)

    mint = 0.001 * param.β
    minK, maxK = 0.01 * sqrt(param.T * param.me), 10param.kF
    order = 4
    Ri, Rt, Ft = Propagators.loadR(fname, param; mint=mint, minK=minK, maxK=maxK, order=order)
    # Ri, Rt = Propagators.initR(param; mint=mint, minK=minK, maxK=maxK, order=order)
    # println(size(Ri))
    println(size(Rt))
    println(size(Ft))
    calcF!(Ft, Ri, Rt, param)
    F2F!(Rt, Ft, param) # Rt actually stores Ft histogram
    # userdata
    funcs = Propagators.Funcs(param, rpai, rpat, Ri, Rt, Ft)

    Fw = Propagators.add_source(Rt, param; source=0.0)
    R0 = real(Propagators.R0(Fw, param))
    println("R0=$(R0)")

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
        [1, 1, 0, 0, 0, 0],
        [1, 1, 1, 2, 1, 0],
        [1, 1, 1, 4, 1, 3]
    ]
    obs = [zeros(ComplexF64, size(Rt)), zeros(ComplexF64, size(Rt)), zeros(ComplexF64, size(Rt))]

    println("Start")
    result = integrate(integrand; measure=measure, userdata=(funcs, diagram, paramc),
        var=(ExtT, ExtK, X, T, P, Q), dof=dof, obs=obs, solver=alg,
        neval=steps, print=-1, block=8, type=ComplexF64)
    # println(result.mean)
    # funcs.Ris.data .= result.mean[1][1, :] .+ result.mean[2][1, :]
    funcs.Ris.data .= 0.0
    funcs.Rt.data .= result.mean[2] .+ result.mean[3]

    return result, funcs
end

if abspath(PROGRAM_FILE) == @__FILE__
    result, funcs = run(steps, param, :vegasmc)
    Ri, Rt = funcs.Ris, funcs.Rt
    println(real(Rt.data[1, :]))
    Fw = Propagators.add_source(Rt, param)
    R0 = real(Propagators.R0(Fw, param))
    println("R0=$(R0)")
    # Propagators.save_pcf(savefname, param, funcs.Ris, funcs.Rt; R0=R0, result=result)
end