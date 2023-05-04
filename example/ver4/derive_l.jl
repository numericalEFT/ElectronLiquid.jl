using ElectronLiquid
using Measurements
using Printf
using JLD2
using PrettyTables
using ElectronGas

const filename = "ver4_L.jld2"

rs = [1.0,]
mass2 = [0.01,]
Fs = [-0.0, -0.2]
beta = [25.0,]
order = [2,]

function addbare!(para, datatuple)
    kgrid, lgrid, ver4 = datatuple
    data100 = zeros(Complex{Measurement{Float64}}, 2, length(lgrid), length(kgrid))
    for (li, l) in enumerate(lgrid)
        for (ki, k) in enumerate(kgrid)
            Ws, Wa = Ver4.projected_exchange_interaction(0, para, Ver4.exchange_interaction; kamp=k, verbose=0)
            Wuu, Wud = Ver4.sa2ud(Ws, Wa)
            data100[1, li, ki] = Wuu
            data100[2, li, ki] = Wud
        end
    end
    ver4[(1, 0, 0)] = data100
end

# function singularterm(para, kgrid, lgrid)

# end

function process200(para, datatuple)
    # z = CounterTerm.chemicalpotential_renormalization(para.order, sw, dmu)
    key = (2, 0, 0)
    li = 1
    kgrid, lgrid, ver4 = datatuple
    addbare!(para, datatuple)

    # df = CounterTerm.fromFile()
    # mu, sw = CounterTerm.getSigma(df, UEG.paraid(para), para.order)
    # dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)

    sigdyn, sigint = SelfEnergy.G0W0(para.basic, int_type=:ko_const, Fs=-para.Fs)
    z = SelfEnergy.zfactor(para.basic, sigdyn, ngrid=[-1, 0])[1]
    dz = [1 - 1 / z,]
    # println(sw)
    # println(dz)
    # dz =[-0.509, ]
    # println(ver4[(1, 0, 0)][1, 1, 1], ", ", ver4[(1, 0, 0)][2, 1, 1])
    # println((ver4[(1, 0, 0)][1, 1, 1] + ver4[(1, 0, 0)][2, 1, 1]) / 2)
    # println("before: ", (ver4[(2, 0, 0)][1, 1, 1] + ver4[(2, 0, 0)][2, 1, 1]) / 2)
    # println(dz)
    ver4 = CounterTerm.z_renormalization(2, ver4, dz, 2)
    # println("after: ", (ver4[(2, 0)][1, 1, 1] + ver4[(2, 0)][2, 1, 1]) / 2)
    # println(ver4)

    kF = para.kF
    printstyled("l=$(lgrid[li]) for $(UEG.short(para))\n", color=:green)
    for (p, data) in ver4
        # head = ["k/kF", "uu", "ud", "symmetric", "asymmetric"]

        printstyled(@sprintf("%12s    %24s    %24s    %24s    %24s     %24s\n",
                "k/kF", "uu", "ud", "symmetric", "asymmetric", "p: $p"), color=:yellow)
        for (ki, k) in enumerate(kgrid)
            factor = 1.0
            d1, d2 = real(data[1, li, ki]) * factor, real(data[2, li, ki]) * factor
            # s, a = (d1 + d2) / 2.0, (d1 - d2) / 2.0
            s, a = Ver4.ud2sa(d1, d2)
            @printf("%12.6f    %24s    %24s    %24s    %24s\n", k / kF, "$d1", "$d2", "$s", "$a")
        end
    end

    printstyled("summed\n", color=:red)
    printstyled(@sprintf("%12s    %24s    %24s    %24s    %24s\n",
            "k/kF", "uu", "ud", "symmetric", "asymmetric"), color=:yellow)
    for (ki, k) in enumerate(kgrid)
        data = ver4[(1, 0)] .+ ver4[(2, 0)]
        d1, d2 = real(data[1, li, ki]), real(data[2, li, ki])
        # s, a = (d1 + d2) / 2.0, (d1 - d2) / 2.0
        s, a = Ver4.ud2sa(d1, d2)
        @printf("%12.6f    %24s    %24s    %24s    %24s\n", k / kF, "$d1", "$d2", "$s", "$a")
    end

end

if abspath(PROGRAM_FILE) == @__FILE__

    # @assert length(ARGS) >= 1 "One argument for the data file name is required!"
    # filename = ARGS[1]
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    f = jldopen(filename, "r")

    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=true)
        kF = para.kF
        for key in keys(f)
            if key == UEG.short(para)
                process200(para, f[key])
            end
        end
    end

end