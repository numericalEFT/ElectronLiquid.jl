using MCIntegration
using LinearAlgebra
using ElectronGas


function integrand(config)
    para = config.para
    EF, kF, me, e0, NF, β = para.EF, para.kF, para.me, para.e0, para.NF, para.β
    varK = config.var[1]
    varN = config.var[2]
    q = varK[1]
    n = varN[1]
    q = sqrt(dot(q, q))
    # println(q, ", ", n)
    # polar = Polarization.Polarization0_3dZeroTemp_LinearDispersion(q, n, para)
    polar = Polarization.Polarization0_3dZeroTemp(q, n, para)
    # W = 8π / (q^2 - 8π * polar + para.Λs)
    W = 8π / (q^2 + 1.0)
    k = q .+ [kF, 0.0, 0.0]
    wn = 1im * (2 * n) * π / β
    delta = 1im * π / β
    ek = dot(k, k) / (2me) - EF
    Gp = -1.0 / (wn + delta - ek)
    Gm = -1.0 / (wn + delta - ek)
    GG = real(Gp * Gm)
    return GG * W / β / (2π)^3
end

function MC(neval=1e6)
    para = Parameter.rydbergUnit(1.0 / 25, 5.0, Λs=0.01)
    K = FermiK(3, para.kF, para.kF * 0.5, para.kF * 6)
    N = Discrete(-1000, 1000)
    dof = [[1, 1],]
    config = MCIntegration.Configuration((K, N), dof; para=para)
    result = MCIntegration.sample(config, integrand; neval=neval, print=0)
    println(result.mean, "+-", result.stdev)
end

MC()