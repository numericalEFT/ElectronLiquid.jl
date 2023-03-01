using JLD2
using Measurements
using ElectronLiquid
using ElectronGas
# using Gaston
using PyCall
pushfirst!(PyVector(pyimport("sys")."path"), ".")

const filename = "data_K.jld2"
const parafilename = "para_wn_1minus0.csv"

function loaddata(para, FileName=filename)
    key = UEG.short(para)
    f = jldopen(FileName, "r")
    # println(key)
    # println(keys(f))
    p, tgrid, kgrid, g = f[key]
    return tgrid, kgrid, g
end

function renormalize(para, sigma, Zrenorm)


    mu, sw = CounterTerm.getSigma(para, parafile=parafilename)
    dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)
    # println(dmu)

    sigma = CounterTerm.chemicalpotential_renormalization(para.order, sigma, dmu)
    return sigma
end

function process(para, Zrenorm)
    dim, β, kF = para.dim, para.β, para.kF

    tgrid, kgrid, data = loaddata(para, filename)

    newdata = Dict()
    kF_label = searchsortedfirst(kgrid, para.kF)
    for key in keys(data)
        # println(key, ": ", data[key][1, kF_label])
        newdata[(key[1] + 1, key[2], key[3])] = data[key]
    end
    # data = mergeInteraction(data)
    # println(size(idata[(1, 0)]))
    return renormalize(para, newdata, Zrenorm), tgrid, kgrid

    # zk = Dict()
    # for (key, val) in idata
    #     zk[key] = zfactor(val, para)
    # end

    # df = fromFile()
    # mu, sw = getSigma(df, UEG.paraid(para), _order - 1)
    # ############ z renormalized  ##########################
    # δμ, δz = derive_onebody_parameter_from_sigma(_order - 1, mu, sw)
    # if Zrenorm
    #     zk = z_renormalization(_order, zk, δz, 1)
    # end
    # zk = chemicalpotential_renormalization(_order, zk, δμ)

    ############ without z renormalized  ##########################
    # δμ, δz = derive_onebody_parameter_from_sigma(_order - 1, mu)
    # zk = chemicalpotential_renormalization(_order, zk, δμ)

    # println(δμ)
    # println(δz)
    # δμ = muCT(df, paraid, _order)
    # println(δμ)
    # println(keys(idata))
    # idata = chemicalpotential_renormalization(_order, idata, δμ)
    # println("zk: ")
    # zk = [zfactor(idata[o]) for o in 1:_order]
    # println(size(zk[1]))
    # δz = zCT(df, paraid, _order)
    # zk[2] = zk[2] - δz[1] * zk[1]
    # return zk, kgrid
end

function plot_k(para, Sw_k, kgrid, Zrenorm)
    println(Sw_k)
    println(size(Sw_k[1]))
    println(kgrid)
    println(length(kgrid))
    dim, β, kF = para.dim, para.β, para.kF
    kF_label = searchsortedfirst(kgrid, kF)
    # zk[1] = zk[1] .- zk[1][kF_label]
    # zk[2] = zk[2] .- zk[2][kF_label]

    plot = pyimport("plot")
    # signal = pyimport("scipy.signal")
    interp = pyimport("scipy.interpolate")

    for o in 1:para.order
        # n0 = @. 1 / (1.0 + exp(para.β * (kgrid^2 / (2 * para.me) - para.EF)))
        # if o == 1
        #     y = n0
        #     e = y .* 0.0
        # else
        #     # println(sum(zk[1:o, 1]))
        #     # zko = 1 ./ (1 .+ sum(zk[1:o]))
        #     # println(length(Sw_k))
        #     zko = sum(Sw_k[1:o-1])
        #     y = [z.val for z in zko[1, :]]
        #     e = [z.err for z in zko[1, :]]
        # end
        zko = sum(Sw_k[1:o])
        y = [z.val for z in zko[1, :]]
        e = [z.err for z in zko[1, :]]
        # println(zk[o])
        plot.plt.errorbar(kgrid / kF, y, yerr=e, color=plot.color[o], capsize=4, fmt="o", markerfacecolor="none", label="Order $o")

        # yfit = signal.savgol_filter(y, 5, 3)
        x = kgrid / kF
        spl = interp.UnivariateSpline(x, y, w=1.0 ./ e)
        yfit = spl(x)
        plot.plt.plot(x, yfit, color=plot.color[o], linestyle="--")
        # println(plt)
        # p = plot(kgrid.grid, zk[o])
        # display(p)
    end
    plot.plt.xlim([kgrid[1] / kF, kgrid[end] / kF])
    plot.plt.xlabel("\$k/k_F\$")
    # plot.plt.ylabel("\$z(k/k_F) = \\left( 1+\\frac{\\partial \\operatorname{Im}\\Sigma(k, i\\omega_0)}{\\partial \\omega}\\right)^{-1}\$")
    if Zrenorm
        plot.plt.ylabel("\$z \\cdot \\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    else
        plot.plt.ylabel("\$\\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    end
    # plot.plt.ylabel("\$z(k/k_F) = \\left( 1+\\frac{\\partial \\operatorname{Im}\\Sigma(k, i\\omega_0)}{\\partial \\omega}\\right)^{-1}\$")
    plot.plt.legend()
    # plot.plt.savefig("sigmaK_rs$(para.rs)_Fs$(para.Fs)_$(para.dim)d.pdf")
    plot.plt.show()
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__

    para = ParaMC(rs=1.0, beta=25.0, Fs=-0.0, order=3, mass2=0.01, isDynamic=true)

    Sw_k, tgrid, kgrid = process(para, false)

    plot_k(para, Sw_k, kgrid, false)
    # readline()
end