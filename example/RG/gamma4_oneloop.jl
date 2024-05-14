using ElectronLiquid
using ElectronGas
using CompositeGrids
using FiniteDifferences
using GreenFunc
using MCIntegration
using FeynmanDiagram
using JLD2
using Printf

include("gamma4_treelevel.jl")

# rs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
# rs = [1.0, 2.0, 3.0, 4.0, 5.0]
# rs = [6.0, 8.0, 10.0]
rs = [4.0,]
mass2 = [0.001,]
beta = [25.0,]
order = [1,]
neval = 2e8

function print_ver4(ver4, lgrid, Λgrid, para)
    kF = para.kF
    if isnothing(ver4) == false
        for (p, data) in ver4
            printstyled("permutation: $p\n", color=:yellow)
            for (li, l) in enumerate(lgrid)
                printstyled("l = $l\n", color=:green)
                @printf("%12s    %16s    %16s    %16s    %16s    %16s    %16s\n", "k/kF", "uu", "ud", "di", "ex", "symmetric", "asymmetric")
                for (ki, k) in enumerate(Λgrid)
                    factor = 1.0
                    d1, d2 = real(data[1, li, ki]) * factor, real(data[2, li, ki]) * factor
                    s, a = (d1 + d2) / 2.0, (d1 - d2) / 2.0
                    di, ex = (s - a), (a) * 2.0
                    @printf("%12.6f    %16s    %16s    %16s    %16s    %16s    %16s\n", k / kF, "$d1", "$d2", "$di", "$ex", "$s", "$a")
                end
            end
        end

    end
end

for (_rs, _mass2, _beta, _order) in Iterators.product(rs, mass2, beta, order)

    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=-0.0, order=_order, mass2=_mass2, isDynamic=true, isFock=false)
    kF = para.kF
    Λgrid = CompositeGrid.LogDensedGrid(:gauss, [1.0 * para.kF, 100 * para.kF], [para.kF,], 8, 0.01 * para.kF, 8)

    # f = jldopen("data_f.jld2", "r")
    # key = "$(UEG.short(para))"
    # Λgrid, fs, us = f[key]

    _Fs = [KO(para, lambda, lambda; verbose=0)[1] for (li, lambda) in enumerate(Λgrid)]
    println(_Fs)
    # println(fs)

    paras = [UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=_Fs[li], Fa=0.0, order=para.order,
        mass2=para.mass2, isDynamic=true, isFock=false) for li in eachindex(Λgrid)]


    partition = [(2, 0, 0),]
    channel = [
        PHr,
        PHEr,
        PPr
    ]

    neighbor = UEG.neighbor(partition)
    filter = [NoHartree,
        NoBubble,
        Proper
    ]

    n = [-1, 0, 0, -1] # q=0 and w -> 0 to measure F
    lgrid = [0,] # angular momentum

    println("∂Γ/∂f")
    diagram = Ver4.diagram(para, partition; channel=channel, filter=filter, dR=true)
    # ver4_df, result = Ver4.PH_df(paras, diagram;
    #     kamp=Λgrid, n=n, l=lgrid,
    #     neval=neval)
    ver4_df, result = Ver4.PH_df(paras, diagram;
        kamp=[para.kF for l in Λgrid], n=n, l=lgrid,
        neval=neval)
    print_ver4(ver4_df, lgrid, Λgrid, para)

    println("Γ on FS")
    diagram = Ver4.diagram(para, partition; channel=channel, filter=filter, dR=false)
    ver4, result = Ver4.PH(paras[1], diagram;
        kamp=[Λgrid[1],], n=n, l=lgrid,
        neval=neval, print=0)
    print_ver4(ver4, lgrid, [Λgrid[1],], para)

    if isnothing(ver4) == false
        jldopen("ver4_PH.jld2", "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (Λgrid, lgrid, _Fs, ver4, ver4_df)
        end
    end
end