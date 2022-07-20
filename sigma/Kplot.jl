using JLD2
using Measurements
# using Gaston
using PyCall
pushfirst!(PyVector(pyimport("sys")."path"), ".")

include("../common/para_builder.jl")
using .UEG
include("../common/counterterm.jl")
using .CounterTerm

const filename = "dataK.jld2"
# const Zrenorm = true    # turn on to renormalize the Z factor
# const Zrenorm = false     # turn off Z renormalization 

function zfactor(idata, para)
    return @. (idata[2, :] - idata[1, :]) / (2π / para.β)
end

function loaddata(para, FileName=filename)
    key = UEG.short(para)
    f = jldopen(FileName, "r")
    # println(key)
    # println(keys(f))
    p, kgrid, avg, std = f[key]
    order = p.order
    _partition = UEG.partition(para.order)
    rdata, idata = Dict(), Dict()
    for (ip, p) in enumerate(_partition)
        rdata[p] = measurement.(real(avg[ip, :, :]), real(std[ip, :, :]))
        idata[p] = measurement.(imag(avg[ip, :, :]), imag(std[ip, :, :]))
    end
    return kgrid, order, _partition, rdata, idata
end

function renormalize(para, sigma, Zrenorm)
    dim, β, kF = para.dim, para.β, para.kF
    df = fromFile()
    mu, sw = getSigma(df, UEG.paraid(para), para.order - 1)
    ############ z renormalized  ##########################
    δμ, δz = derive_onebody_parameter_from_sigma(para.order - 1, mu, sw)
    if Zrenorm
        sigma = z_renormalization(para.order, sigma, δz, 1)
    end
    sigma = chemicalpotential_renormalization(para.order, sigma, δμ)
    return sigma
end

function process(para, Zrenorm)
    dim, β, kF = para.dim, para.β, para.kF

    kgrid, _order, _partition, rdata, idata = loaddata(para, filename)

    rdata = mergeInteraction(rdata)
    idata = mergeInteraction(idata)
    # println(size(idata[(1, 0)]))
    return renormalize(para, rdata, Zrenorm), renormalize(para, idata, Zrenorm), kgrid

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
    return zk, kgrid
end

function plotSw_k(para, Sw_k, kgrid, Zrenorm)
    dim, β, kF = para.dim, para.β, para.kF
    kF_label = searchsortedfirst(kgrid.grid, kF)
    # zk[1] = zk[1] .- zk[1][kF_label]
    # zk[2] = zk[2] .- zk[2][kF_label]

    plot = pyimport("plot")
    # signal = pyimport("scipy.signal")
    interp = pyimport("scipy.interpolate")
    for o in 1:para.order
        # println(sum(zk[1:o, 1]))
        # zko = 1 ./ (1 .+ sum(zk[1:o]))
        zko = Sw_k[o]
        y = [z.val for z in zko]
        e = [z.err for z in zko]
        # println(zk[o])
        plot.plt.errorbar(kgrid.grid / kF, y, yerr=e, color=plot.color[o], capsize=4, fmt="o", markerfacecolor="none", label="Order $o")

        # yfit = signal.savgol_filter(y, 5, 3)
        x = kgrid.grid / kF
        spl = interp.UnivariateSpline(x, y, w=1.0 ./ e)
        yfit = spl(x)
        plot.plt.plot(x, yfit, color=plot.color[o], linestyle="--")
        # println(plt)
        # p = plot(kgrid.grid, zk[o])
        # display(p)
    end
    plot.plt.xlim([kgrid.grid[1] / kF, kgrid.grid[end] / kF])
    plot.plt.xlabel("\$k/k_F\$")
    # plot.plt.ylabel("\$z(k/k_F) = \\left( 1+\\frac{\\partial \\operatorname{Im}\\Sigma(k, i\\omega_0)}{\\partial \\omega}\\right)^{-1}\$")
    if Zrenorm
        plot.plt.ylabel("\$z \\cdot \\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    else
        plot.plt.ylabel("\$\\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    end
    # plot.plt.ylabel("\$z(k/k_F) = \\left( 1+\\frac{\\partial \\operatorname{Im}\\Sigma(k, i\\omega_0)}{\\partial \\omega}\\right)^{-1}\$")
    plot.plt.legend()
    plot.plt.savefig("sigmaK_rs$(para.rs)_Fs$(para.Fs)_$(para.dim)d.pdf")
    plot.plt.show()
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__

    para = ParaMC(rs=5.0, beta=100.0, Fs=-0.0, order=2, mass2=1e-5)

    Sw_k, kgrid, Znorm = process(para, false)

    plotSw_k(para, Sw_k, kgrid, Znorm)
    # readline()
end