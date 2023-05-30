using MCIntegration
using ElectronLiquid
using Measurements
using LinearAlgebra
using Lehmann

function _PP(vars, config)
    varK, varT, varX = vars
    R, Theta, Phi = varK
    para, kamp, kamp2 = config.userdata

    t1, t2 = varT[1], varT[2]
    r = R[1] / (1 - R[1])
    θ, ϕ = Theta[1], Phi[1]
    x = varX[1]

    q = [r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ)]

    k1 = [kamp + q[1], q[2], q[3]]
    k2 = [kamp2 * x - q[1], kamp2 * sqrt(1 - x^2) - q[2], -q[3]]

    ek1 = (dot(k1, k1) - para.kF^2) / 2 / para.me
    ek2 = (dot(k2, k2) - para.kF^2) / 2 / para.me

    g1 = Spectral.kernelFermiT(t2, ek1, para.β)
    g2 = Spectral.kernelFermiT(t2, ek2, para.β)
    g3 = Spectral.kernelFermiT(t2 - t1, ek2, para.β)

    qd = sqrt(dot(q, q))
    vq = -UEG.interactionStatic(para, qd, 0.0, t1)
    wq = -UEG.interactionDynamic(para, qd, 0.0, t1)

    factor = r^2 * sin(θ) / (1 - R[1])^2 / (2π)^3 * para.NF

    wud = g1 * (g2 * vq + g3 * wq) * factor * phase((0.0, t2, t1, t2), -1, 0, 0, para.β)
    return wud
end

@inline function phase(extT, ninL, noutL, ninR, β)
    # println(extT)
    tInL, tOutL, tInR, tOutR = extT
    winL, woutL, winR = π * (2ninL + 1) / β, π * (2noutL + 1) / β, π * (2ninR + 1) / β
    woutR = winL + winR - woutL
    return exp(-1im * (tInL * winL - tOutL * woutL + tInR * winR - tOutR * woutR))
end

function vertex3(para::ParaMC;
    neval=1e6, #number of evaluations
    kamp=para.kF,
    kamp2=para.kF,
    config=nothing,
    print=0,
    integrand=_PP,
    kwargs...
)

    UEG.MCinitialize!(para)
    dim, β, kF = para.dim, para.β, para.kF
    R = MCIntegration.Continuous(0.0, 1.0 - 1e-6) # a small cutoff to make sure R is not 1
    Theta = MCIntegration.Continuous(0.0, 1π)
    Phi = MCIntegration.Continuous(0.0, 2π)
    K = CompositeVar(R, Theta, Phi)
    T = MCIntegration.Continuous(0.0, β, offset=1)
    T.data[1] = 0.0
    X = MCIntegration.Continuous(-1.0, 1.0) #x=cos(θ)

    dof = [[1, 2, 1],]

    if isnothing(config)
        config = MCIntegration.Configuration(
            var=(K, T, X),
            dof=dof,
            type=ComplexF64,
            userdata=(para, kamp, kamp2),
            kwargs...
        )
    end

    result = integrate(integrand;
        config=config, neval=neval, print=print, solver=:vegasmc, kwargs...)

    if isnothing(result) == false
        avg, std = result.mean, result.stdev
        r = measurement.(real(avg), real(std))
        i = measurement.(imag(avg), imag(std))
        return Complex(r, i)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    rs = [5.0,]
    mass2 = [0.001,]
    _Fs = [-2.5,]
    beta = [25.0,]
    order = [1,]
    neval = 1e6

    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, _Fs, beta, order)
        para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=1, mass2=_mass2, isDynamic=true, isFock=false)
        result = vertex3(para;
            neval=1e6,
            integrand=_PP,
            print=0
        )
        println(result)
    end

end