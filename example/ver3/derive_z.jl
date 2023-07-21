# include("../common/counterterm.jl")
# include("../common/para_builder.jl")
# using .UEG
# using .CounterTerm
using ElectronLiquid
using Printf
using Measurements
using JLD2

# Zrenorm = false
# Zrenorm = true

rs = [5.0,]
mass2 = [0.01,]
Fs = [-0.0,]
beta = [25, 50, 100]
order = [2,]
neval = 1e6

mission = ARGS[1]
println("mission (Z or A): ", mission)

filename = "ver3_$mission.jld2"

f = jldopen(filename, "r")

# function addbare!(datatuple)
#     para, kgrid, qgrid, nin, nqout, ver3 = datatuple
#     data000 = zeros(Complex{Measurement{Float64}}, 2, length(kgrid), length(qgrid), length(nin), length(nqout))
#     for (qi, q) in enumerate(qgrid)
#         for (ki, k) in enumerate(kgrid)
#             for (ni, n) in enumerate(nin)
#                 for (mi, m) in enumerate(nqout)
#                     data100[1, ki, qi, ni, mi] = 1.0
#                     data100[2, ki, qi, ni, mi] = 0.0
#                 end
#             end
#         end
#     end
#     ver3[(0, 0, 0)] = data000
# end

function process(datatuple)
    para, kgrid, qgrid, nin, nqout, ver3 = datatuple
    kF = para.kF
    # addbare!(datatuple)

    mu, sw = CounterTerm.getSigma(para)
    dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)

    _z = Dict()
    for (p, val) in ver3
        #by definition, z factor expects q=0, w->0 limit
        _z[p] = real.(val[1, :, :, 1, 1] + val[2, :, :, 1, 1])  #only keep the kin and qout dependence
    end

    z = CounterTerm.chemicalpotential_renormalization(para.order, _z, dmu)


    z = [-ones(size(z[1])), z...] #add the bare gamma3, which is one
    println(z)

    z = CounterTerm.z_renormalization(para.order + 1, z, dz, 1)
    println(z)

    println(dz)

    for o in 1:para.order
        printstyled("order = $(o)\n", color=:green)
        printstyled(@sprintf("%12s  %12s    %16s\n", "k/kF", "q/kF", "sum"), color=:yellow)
        for (qi, q) in enumerate(qgrid)
            for (ki, k) in enumerate(kgrid)
                s = real(z[o+1][1, ki, qi])
                @printf("%12.6f  %12.6f    %16s\n", k / kF, q / kF, "$s")
            end
        end
    end


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
        loadpara = UEG.ParaMC(key)
        if UEG.paraid(loadpara) == UEG.paraid(para)
            process(f[key])
        end
    end
end