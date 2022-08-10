using ElectronLiquid
using Measurements
using Printf
using JLD2

const filename = "ver4_L.jld2"

rs = [5.0,]
mass2 = [0.01, 0.001,]
Fs = [-0.0, ]
beta = [25.0,]
order = [2,]
neval = 1e8

function zfactor(data, β)
    return @. (imag(data[2, 1]) - imag(data[1, 1])) / (2π / β)
end

function mu(data)
    return real(data[1, 1])
end

function process(FileName, isSave)
    f = jldopen(FileName, "r")
    df = CounterTerm.fromFile()
    for key in keys(f)
        println(key)
        value = f[key]
        para, ngrid, kgrid, data = value[1], value[2], value[3], value[4]

        for p in sort([k for k in keys(data)])
            println("$p: μ = $(mu(data[p]))   z = $(zfactor(data[p], para.β))")
        end

        _mu = Dict()
        for (p, val) in data
            _mu[p] = mu(val)
        end
        _z = Dict()
        for (p, val) in data
            _z[p] = zfactor(val, para.β)
        end

        ############# save to csv  #################
        # println(df)
        for o in keys(data)
            # global df
            paraid = UEG.paraid(para)
            df = CounterTerm.appendDict(df, paraid, Dict("order" => o, "μ" => _mu[o].val, "μ.err" => _mu[o].err, "Σw" => _z[o].val, "Σw.err" => _z[o].err); replace=true)
        end
    end
    # println("new dataframe\n$df")
    isSave && CounterTerm.toFile(df)
end

function addbare!(datatuple)
    para, kgrid, lgrid, ver4 = datatuple
    data100 = zeros(Complex{Measurement{Float64}}, 2, length(lgrid), length(kgrid))
    for (li, l) in enumerate(lgrid)
        for (ki, k) in enumerate(kgrid)
            Ws, Wa = Ver4.projected_exchange_interaction(0, para, Ver4.exchange_interaction; kamp =k, verbose=0)
            Wuu, Wud = Ver4.sa2ud(Ws, Wa)
            data100[1, li, ki] = Ws
            data100[2, li, ki] = Wa
        end
    end
    ver4[(1, 0, 0)] = data100
end

function process200(datatuple)
    # z = CounterTerm.chemicalpotential_renormalization(para.order, sw, dmu)
    key = (2, 0, 0)
    li=1
    para, kgrid, lgrid, ver4 = datatuple
    addbare!(datatuple)

    df = CounterTerm.fromFile()
    mu, sw = CounterTerm.getSigma(df, UEG.paraid(para), para.order)
    dmu, dz = CounterTerm.derive_onebody_parameter_from_sigma(para.order, mu, sw, zrenorm=false)

    println(sw)
    println(dz)
    dz =[-0.509, ]
    ver4 = CounterTerm.z_renormalization(2, ver4, dz, 2)

    kF=para.kF
    println(UEG.short(para))
    for (p, data) in ver4
        printstyled("permutation: $p\n", color=:yellow)
        printstyled("l = $(lgrid[li])\n", color=:green)
        @printf("%12s    %16s    %16s    %16s    %16s\n", "k/kF", "uu", "ud", "symmetric", "asymmetric")
        for (ki, k) in enumerate(kgrid)
            factor = 1.0
            d1, d2 = real(data[1, li, ki])*factor, real(data[2, li, ki])*factor
            # s, a = (d1 + d2) / 2.0, (d1 - d2) / 2.0
            s, a = Ver4.ud2sa(d1, d2)
            @printf("%12.6f    %16s    %16s    %16s    %16s\n", k / kF, "$d1", "$d2", "$s", "$a")
        end
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

    for _rs in rs
        for _mass2 in mass2
            for _F in Fs
                for _beta in beta
                    for _order in order

                        para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=true)
                        kF = para.kF
                        for key in keys(f)
                            if UEG.paraid(f[key][1]) == UEG.paraid(para)
                                process200(f[key])
                            end
                        end
                    end
                end
            end
        end
    end

end