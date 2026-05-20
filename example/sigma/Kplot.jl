using JLD2
using Measurements
# using Gaston
using PyCall
using PyPlot
using ElectronLiquid
using ElectronGas
using TaylorSeries

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
    p = UEG.ParaMC(key)
    ngrid, kgrid, sigma = f[key]
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

function spline(x, y, e)
    # signal = pyimport("scipy.signal")
    interp = pyimport("scipy.interpolate")
    # yfit = signal.savgol_filter(y, 5, 3)
    w = 1.0 ./ e
    kidx = searchsortedfirst(x, 0.5)
    _x, _y = deepcopy(x[kidx:end]), deepcopy(y[kidx:end])
    _w = 1.0 ./ e[kidx:end]

    #enforce the boundary condition: the derivative at k=0 is zero
    pushfirst!(_x, 0.01)
    pushfirst!(_x, 0.0)
    kidx = searchsortedfirst(_x, 0.7)
    yavr = sum(y[1:kidx] .* w[1:kidx]) / sum(w[1:kidx])
    pushfirst!(_y, yavr)
    pushfirst!(_y, yavr)
    pushfirst!(_w, _w[1] * 10000)
    pushfirst!(_w, _w[1] * 10000)

    # generate knots with spline without constraints
    spl = interp.UnivariateSpline(_x, _y, w=_w, k=3)
    __x = collect(LinRange(0.0, x[end], 100))
    yfit = spl(__x)
    return __x, yfit
end

function plotS_k(para, rSw_k, iSw_k, kgrid, Zrenorm)
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
    figure(figsize=(8, 4))

    function sk(sigma, order, kgrid)
        dk = [(sigma[o][1, 2:end] .- sigma[o][1, 1]) ./ (kgrid[2:end] .^ 2 / (2 * para.me)) for o in 1:order]
        return kgrid[2:end], sum(dk)
    end

    function delta_z(sigma, order, ki)
        # return [(sigma[o][2, ki] .- sigma[o][1, ki]) ./ (2π / para.β) for o in 1:order]
        return [(sigma[o][1, ki]) ./ (π / para.β) for o in 1:order]
    end

    # power series of z-factor
    function sw(sigma, order, kgrid)
        dw = [sum(delta_z(sigma, order, ki)) for ki in 1:length(kgrid)]
        # dw = [(sigma[o][1, :]) ./ (π / para.β) for o in 1:order]
        return kgrid, dw
    end

    # power series of 1/z
    function sw_inv(sigma, order, kgrid)
        dw = []
        for ki in 1:length(kgrid)
            δzi = delta_z(sigma, order, ki)
            zi = Taylor1([1.0, δzi...], order)
            z = 1 / zi
            δz = [getcoeff(z, o) for o in 1:order]
            push!(dw, sum(δz[1:end]))
        end
        return kgrid, dw
    end

    subplot(1, 2, 2)
    for o in 1:para.order
        _kgrid, s_k = sk(rSw_k, o, kgrid)
        z = s_k
        y = [-z.val for z in z]
        e = [z.err for z in z]
        errorbar(_kgrid / kF, y, yerr=e, color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="Order $o")

        _x, _y = spline(_kgrid / kF, y, e)
        plot(_x, _y, color=color[o], linestyle="--")
    end
    xlim([kgrid[1] / kF, 2.0])
    ylim([0.0, 1.0])
    xlabel("\$k/k_F\$")
    # plot.plt.ylabel("\$z(k/k_F) = \\left( 1+\\frac{\\partial \\operatorname{Im}\\Sigma(k, i\\omega_0)}{\\partial \\omega}\\right)^{-1}\$")
    # if Zrenorm
    #     ylabel("\$z \\cdot \\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    # else
    #     ylabel("\$\\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    # end
    # plot.plt.ylabel("\$z(k/k_F) = \\left( 1+\\frac{\\partial \\operatorname{Im}\\Sigma(k, i\\omega_0)}{\\partial \\omega}\\right)^{-1}\$")
    legend()

    subplot(1, 2, 1)
    for o in 1:para.order
        _kgrid, s_w = sw_inv(iSw_k, o, kgrid)
        z = [(1.0 + e) for e in s_w]
        # _kgrid, s_w = sw(iSw_k, o, kgrid)
        # z = [1.0 / (1.0 + e) for e in s_w]
        y = [z.val for z in z]
        e = [z.err for z in z]
        kidx = searchsortedfirst(_kgrid, kF)
        println("order $o at k=$(_kgrid[kidx]/kF): $(z[kidx])")

        errorbar(_kgrid / kF, y, yerr=e, color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="Order $o")

        _x, _y = spline(_kgrid / kF, y, e)
        plot(_x, _y, color=color[o], linestyle="--")
    end
    xlim([kgrid[1] / kF, 2.0])
    ylim([0.0, 1.0])
    xlabel("\$k/k_F\$")
    # if Zrenorm
    #     ylabel("\$z \\cdot \\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    # else
    #     ylabel("\$\\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    # end
    # plot.plt.ylabel("\$z(k/k_F) = \\left( 1+\\frac{\\partial \\operatorname{Im}\\Sigma(k, i\\omega_0)}{\\partial \\omega}\\right)^{-1}\$")
    legend()
    savefig("sigmaK.pdf")
    show()
    # readline()
end

if abspath(PROGRAM_FILE) == @__FILE__

    para = ParaMC(rs=4.0, beta=25.0, Fs=-0.0, order=3, mass2=0.01, isDynamic=true)

    rSw_k, iSw_k, kgrid, ngrid = process(para, false)


    kF_label = searchsortedfirst(kgrid, para.kF)
    for k in keys(rSw_k)
        println(k, ": ", rSw_k[k][1, kF_label])
    end
    plotS_k(para, rSw_k, iSw_k, kgrid, Zrenorm)
    # plotS_k(para, iSw_k, kgrid, Zrenorm)
    # readline()
end