# include("../common/counterterm.jl")
# include("../common/para_builder.jl")
# using .UEG
# using .CounterTerm
using ElectronLiquid
using Printf
using JLD2

# Zrenorm = false
# Zrenorm = true

rs = [5.0,]
mass2 = [0.01, ]
Fs = [-0.0, ]
beta = [25, ]
order = [2,]
neval = 1e6

# mission = ARGS[1]
# println("mission (Z or K): ", mission)
# exit(0)

filename = "ver3_Z.jld2"

f = jldopen(filename, "r")

# function addbare!(datatuple)
#     para, kgrid, qgrid, nin, nqout, ver3 = datatuple
    # data100 = zeros(Complex{Measurement{Float64}}, 2, length(kgrid), length(qgrid), length(nin), length(nqout))
    # for (qi, q) in enumerate(qgrid)
    #     for (ki, k) in enumerate(kgrid)
    #         Ws, Wa = Ver4.projected_exchange_interaction(0, para, Ver4.exchange_interaction; kamp =k, verbose=0)
    #         Wuu, Wud = Ver4.sa2ud(Ws, Wa)
    #         data100[1, li, ki] = Ws
    #         data100[2, li, ki] = Wa
    #     end
    # end
    # ver4[(1, 0, 0)] = data100
# end

function process(datatuple)
    para, kgrid, qgrid, nin, nqout, ver3 = datatuple
    # addbare!(datatuple)

    mu, sw = CounterTerm.getSigma(para)
    dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)

    _zk = Dict()
    for (p, val) in ver3
        #by definition, z factor expects q=0, w->0 limit
        _zk[p] = real.(val[1, :, 1, 1, 1]+val[2, :, 1, 1, 1])  #only keep the kin dependence
    end

    zk = CounterTerm.chemicalpotential_renormalization(2, _zk, dmu)

    println(dzi)
    println(zk)


    # sumzi = accumulate(+, ver3)
    # # if Zrenorm == false

    # z1 = @. 1.0 / (1.0 + ver3)
    # # println("z without renormalization = $z")
    # # else
    # sumz = accumulate(+, dz)
    # z2 = @. 1.0 + sumz
    # # end
    # # println("z with renormalization = $z")
    # printstyled(UEG.short(para), color=:green)
    # printstyled(@sprintf("%8s   %24s    %24s     %24s\n", "order", "chemical-potential", "no-z-renorm", "z-renorm"), color=:yellow)
    # for o in 1:para.order
    #     @printf("%8d   %24s   %24s     %24s\n", o, "$(dmu[o])", "$(z1[o])", "$(z2[o])")
    # end
end

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=true)
    kF = para.kF
    for key in keys(f)
        if UEG.paraid(f[key][1]) == UEG.paraid(para)
            process(f[key])
        end
    end
end