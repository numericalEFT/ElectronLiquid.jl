using MCIntegration
using ElectronGas
using Measurements
using Plots

const para = Parameter.rydbergUnit(1.0 / 100, 5.0)
const ωp = para.ωp
const kF = para.kF / para.qTF
const EF = para.EF / ωp

const kgrid = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0] .* kF

function integrand(config)
    q = config.var[1][1]
    theta = config.var[2][1]
    ki = config.var[3][1]
    p = kgrid[ki]
    qxy, qz = q * sin(theta), q * cos(theta)
    skq = sqrt(1 + (p - qz)^2 + qxy^2)
    # ek = sqrt(3) / 2 / kF * (q^2 - kF^2)
    ek = para.qTF^2 / (2 * para.me * ωp) * (q^2 - kF^2)
    f = π * q^2 * sin(theta) / skq
    if q < kF
        return f / (skq - ek)^2
    else
        return f / (skq + ek)^2
    end
end

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    ki = config.var[3][1]
    weight = integrand(config)
    config.observable[ki] += weight / abs(weight) * factor
end

function MC(neval=1e6)
    q = MCIntegration.Continuous(0.0, 10.0 * kF)
    theta = MCIntegration.Continuous(0.0, 1π)
    K = MCIntegration.Discrete(1, length(kgrid))
    dof = [[1, 1, 1],]
    obs = zeros(length(kgrid))
    config = MCIntegration.Configuration((q, theta, K), dof, obs)
    result = MCIntegration.sample(config, integrand, measure; neval=neval)

    avg, std = result.mean, result.stdev
    estimate = measurement.(avg, std)
    return estimate
end

estimate = MC()
plot(kgrid ./ kF, estimate)