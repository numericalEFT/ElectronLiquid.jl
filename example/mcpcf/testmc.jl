using MCIntegration

const steps = 1e7

function integrand(vars, config)
    return 1.0, 1.0
end

function measure(vars, obs, weight, config)
    ext1, ext2 = vars[1], vars[2]
    obs[1][ext1[1]] += weight[1]
    obs[2][ext1[1], ext2[1]] += weight[2]
end

function run(steps)
    alg = :vegasmc

    # userdata

    Ext1 = Discrete(1, 10; adapt=false)
    Ext2 = Discrete(1, 10; adapt=false)

    # ExtT, ExtK, X, T, P
    dof = [[1, 0], [1, 1]]
    # obs = [zeros(ComplexF64, (10, 10)), zeros(ComplexF64, (10, 10))]
    obs = [zeros(ComplexF64, 10), zeros(ComplexF64, (10, 10))]

    println("Start")
    result = integrate(integrand; measure=measure,
        var=(Ext1, Ext2), dof=dof, obs=obs, solver=alg,
        neval=steps, print=0, block=8, type=ComplexF64)

    println(result.mean[1])
    println(result.mean[2])
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = run(steps)
end