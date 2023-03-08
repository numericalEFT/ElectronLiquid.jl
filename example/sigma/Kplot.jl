using JLD2
using Measurements
# using Gaston
using PyCall
using PyPlot
using ElectronLiquid
using ElectronGas

# pushfirst!(PyVector(pyimport("sys")."path"), ".")

const filename = "data_K.jld2"
const parafilename = "para_wn_1minus0.csv"
# const Zrenorm = true    # turn on to renormalize the Z factor
const Zrenorm = false     # turn off Z renormalization 

# Vibrant qualitative colour scheme from https://personal.sron.nl/~pault/
const cdict = Dict(["blue" => "#0077BB", "cyan" => "#33BBEE", "teal" => "#009988", "orange" => "#EE7733", "red" => "#CC3311", "magenta" => "#EE3377", "grey" => "#BBBBBB"]);

function loaddata(para, FileName=filename)
    key = UEG.short(para)
    f = jldopen(FileName, "r")
    # println(key)
    # println(keys(f))
    p, ngrid, kgrid, sigma = f[key]
    # println(sigma)
    order = p.order
    _partition = UEG.partition(para.order)
    rdata, idata = Dict(), Dict()
    for (ip, p) in enumerate(_partition)
        rdata[p] = real(sigma[p][:, :])
        idata[p] = imag(sigma[p][:, :])
    end
    return ngrid, kgrid, rdata, idata
end

function renormalize(para, sigma, Zrenorm)
    dim, β, kF = para.dim, para.β, para.kF
    mu, sw = CounterTerm.getSigma(para, parafile=parafilename)
    ############ z renormalized  ##########################
    dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)
    println(para.order)
    println(dmu)
    sigma = CounterTerm.chemicalpotential_renormalization(para.order, sigma, dmu)
    return sigma
end

function process(para, Zrenorm)
    dim, β, kF = para.dim, para.β, para.kF

    ngrid, kgrid, rdata, idata = loaddata(para, filename)

    return renormalize(para, rdata, Zrenorm), renormalize(para, idata, Zrenorm), kgrid, ngrid
end

function plotS_k(para, Sw_k, kgrid, Zrenorm)
    dim, β, kF = para.dim, para.β, para.kF
    kF_label = searchsortedfirst(kgrid, kF)
    # zk[1] = zk[1] .- zk[1][kF_label]
    # zk[2] = zk[2] .- zk[2][kF_label]

    # plot = pyimport("plot")
    style = PyPlot.matplotlib."style"
    style.use(["science", "std-colors"])
    color = [cdict["blue"], cdict["red"], "black"]
    #cmap = get_cmap("Paired")
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 16
    rcParams["font.family"] = "Times New Roman"
    figure(figsize=(4, 4))

    function sk(sigma, order, kgrid)
        dk = [(sigma[o][1, 2:end] .- sigma[o][1, 1]) ./ (kgrid[2:end] .^ 2 / (2 * para.me)) for o in 1:order]
        return kgrid[2:end], dk
    end

    # signal = pyimport("scipy.signal")
    interp = pyimport("scipy.interpolate")
    for o in 1:para.order
        # println(sum(zk[1:o, 1]))
        # zko = 1 ./ (1 .+ sum(zk[1:o]))
        # zko = Sw_k[o]
        # y = [z.val for z in zko]
        # e = [z.err for z in zko]
        _kgrid, s_k = sk(Sw_k, o, kgrid)
        z = sum(s_k)
        y = [-z.val for z in z]
        e = [z.err for z in z]
        # println(zk[o])
        # println(length(_kgrid))
        # println(length(y))
        errorbar(_kgrid / kF, y, yerr=e, color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="Order $o")

        # yfit = signal.savgol_filter(y, 5, 3)
        x = _kgrid / kF
        spl = interp.UnivariateSpline(x, y, w=1.0 ./ e)
        yfit = spl(x)
        plot(x, yfit, color=color[o], linestyle="--")
        # println(plt)
        # p = plot(kgrid.grid, zk[o])
        # display(p)
    end
    xlim([kgrid[1] / kF, kgrid[end] / kF])
    xlabel("\$k/k_F\$")
    # plot.plt.ylabel("\$z(k/k_F) = \\left( 1+\\frac{\\partial \\operatorname{Im}\\Sigma(k, i\\omega_0)}{\\partial \\omega}\\right)^{-1}\$")
    if Zrenorm
        ylabel("\$z \\cdot \\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    else
        ylabel("\$\\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    end
    # plot.plt.ylabel("\$z(k/k_F) = \\left( 1+\\frac{\\partial \\operatorname{Im}\\Sigma(k, i\\omega_0)}{\\partial \\omega}\\right)^{-1}\$")
    legend()
    #plot.plt.savefig("sigmaK_rs$(para.rs)_Fs$(para.Fs)_$(para.dim)d.pdf")
    show()
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__

    para = ParaMC(rs=4.0, beta=25.0, Fs=-0.0, order=3, mass2=0.01, isDynamic=true)

    rSw_k, iSw_k, kgrid, ngrid = process(para, false)


    kF_label = searchsortedfirst(kgrid, para.kF)
    for k in keys(rSw_k)
        println(k, ": ", rSw_k[k][1, kF_label])
    end
    plotS_k(para, rSw_k, kgrid, Zrenorm)
    # readline()
end