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
const Zrenorm = false     # turn off Z renormalization 

const para = ParaMC(rs=5.0, beta=100.0, Fs=-0.0, order=2, mass2=1e-5)

function zfactor(idata, para)
    return @. (idata[2, :] - idata[1, :]) / (2π / para.β)
end

function loaddata(FileName, para)
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

if abspath(PROGRAM_FILE) == @__FILE__

    dim, β, kF = para.dim, para.β, para.kF

    kgrid, _order, _partition, rdata, idata = loaddata(filename, para)

    rdata = mergeInteraction(rdata)
    idata = mergeInteraction(idata)
    println(size(idata[(1, 0)]))

    zk = Dict()
    for (key, val) in idata
        zk[key] = zfactor(val, para)
    end

    df = fromFile()
    mu, sw = getSigma(df, UEG.paraid(para), _order - 1)
    ############ z renormalized  ##########################
    δμ, δz = derive_onebody_parameter_from_sigma(_order - 1, mu, sw)
    if Zrenorm
        zk = z_renormalization(_order, zk, δz, 1)
    end
    zk = chemicalpotential_renormalization(_order, zk, δμ)

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

    kF_label = searchsortedfirst(kgrid.grid, kF)
    # zk[1] = zk[1] .- zk[1][kF_label]
    # zk[2] = zk[2] .- zk[2][kF_label]

    plot = pyimport("plot")
    for o in 1:_order
        # println(sum(zk[1:o, 1]))
        # zko = 1 ./ (1 .+ sum(zk[1:o]))
        zko = zk[o]
        y = [z.val for z in zko]
        e = [z.err for z in zko]
        # println(zk[o])
        plot.plt.errorbar(kgrid.grid / kF, y, yerr=e, color=plot.color[o], label="Order $o")
        # println(plt)
        # p = plot(kgrid.grid, zk[o])
        # display(p)
    end
    plot.plt.xlim([kgrid.grid[1] / kF, kgrid.grid[end] / kF])
    plot.plt.xlabel("\$k/k_F\$")
    plot.plt.ylabel("\$z(k/k_F) = \\left( 1+\\frac{\\partial \\operatorname{Im}\\Sigma(k, i\\omega_0)}{\\partial \\omega}\\right)^{-1}\$")
    plot.plt.legend()
    # plot.plt.savefig("sigmaK_rs5_Fs0_noshift_3d.pdf")
    plot.plt.show()
    # readline()
end