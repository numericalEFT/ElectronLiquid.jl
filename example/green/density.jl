using ElectronLiquid
using ElectronGas
using Measurements
using Printf
using JLD2

rs = [1.0,]
mass2 = [0.001,]
Fs = [-0.0,]
beta = [25.0,]
order = [3,]
neval = 1e6

isDynamic = true
isFock = false


const filename = "data_n.jld2"
const parafilename = "para_wn_1minus0.csv"

for (_rs, _F, _beta, _order, _mass2) in Iterators.product(rs, Fs, beta, order, mass2)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, isFock=isFock)

    mu, sw = CounterTerm.getSigma(para, parafile=parafilename)
    dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)

    println("δz = ", dz)
    println("δμ = ", dmu)
    println("δz_inverse = ", dzi)

    # dmu = [-0.61341, -0.08601]

    printstyled("$(UEG.short(para))\n", color=:green)

    f = jldopen(filename, "r")

    kF = para.kF
    for key in keys(f)
        # println(key)
        if UEG.paraid(f[key][1]) == UEG.paraid(para)
            ndict = f[key][2]
            result = CounterTerm.chemicalpotential_renormalization(_order, ndict, dmu)
            println(result)
        end
    end

end