using JLD2
using Measurements
# using Gaston
using PyCall
using PyPlot
using ElectronLiquid
using ElectronGas
using TaylorSeries

dim = 2
isDynamic = false
rs = [1.0]
mass2 = [1.2]
# rs = [0.5]
# mass2 = [2.0]
# mass2 = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6]
Fs = [-0.0,]
beta = [50.0,]
# order = [4,]
order = [3,]
# i_δk = [4,]
i_δk = [1, 2, 3, 4, 5, 6]

const filenameZ = "data_Z.jld2"
const filenameK = "data_K.jld2"
# const filenameK = "data_K_all1.jld2"
const parafilename = "para_wn_1minus0.csv"
const Zrenorm = true     # turn on/off Z renormalization 
# para = ParaMC(rs=1.0, beta=40.0, Fs=-0.0, order=3, mass2=1.0, isDynamic=false)
# Vibrant qualitative colour scheme from https://personal.sron.nl/~pault/
cdict = Dict(["blue" => "#0077BB", "cyan" => "#33BBEE", "teal" => "#009988", "orange" => "#EE7733", "red" => "#CC3311", "magenta" => "#EE3377", "grey" => "#BBBBBB"]);

function zfactor(data, β)
    return @. (imag(data[2, 1]) - imag(data[1, 1])) / (2π / β)
end

function mu(data)
    return real(data[1, 1])
end

function derive_z(datatuple, isSave)
    df = CounterTerm.fromFile(parafilename)
    para, ngrid, kgrid, data = datatuple
    printstyled(UEG.short(para), color=:yellow)
    println()

    # for p in sort([k for k in keys(data)])
    #     println("$p: μ = $(mu(data[p]))   z = $(zfactor(data[p], para.β))")
    # end

    _mu = Dict()
    for (p, val) in data
        _mu[p] = mu(val)
    end
    _z = Dict()
    for (p, val) in data
        _z[p] = zfactor(val, para.β)
    end

    dzi, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _z)
    println("dz: ", dzi)

    z = []
    for i in eachindex(dzi)
        push!(z, 1.0 / (1.0 + sum(dzi[1:i])))
    end
    println("Zfactor:", z)

    ############# save to csv  #################
    # println(df)
    if isSave
        for o in keys(data)
            # println(o)
            # global df
            paraid = UEG.paraid(para)
            df = CounterTerm.appendDict(df, paraid, Dict("order" => o, "μ" => _mu[o].val, "μ.err" => _mu[o].err, "Σw" => _z[o].val, "Σw.err" => _z[o].err); replace=true)
        end
        CounterTerm.toFile(df, parafilename)
    end
    # println("new dataframe\n$df")
    # isSave && CounterTerm.toFile(df, parafilename)

    # return z
end

function loaddata(para, FileName=filename)
    key = UEG.short(para)
    f = jldopen(FileName, "r")
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

function renormalize(para, sigma)
    mu, sw = CounterTerm.getSigma(para, parafile=parafilename)
    ############ z renormalized  ##########################
    dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)
    sigma = CounterTerm.chemicalpotential_renormalization(para.order, sigma, dmu)
    return sigma
end

function process(para, Zrenorm, filename, i_δk)
    ngrid, kgrid, rdata, idata = loaddata(para, filename)
    kF_label = searchsortedfirst(kgrid, para.kF)
    # rSw_k, iSw_k = renormalize(para, rdata), renormalize(para, idata)

    mu, sw = CounterTerm.getSigma(para, parafile=parafilename)
    δsn, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)
    rSw_k = CounterTerm.chemicalpotential_renormalization(para.order, rdata, dmu)
    iSw_k = CounterTerm.chemicalpotential_renormalization(para.order, idata, dmu)

    _dmn = [Dict() for i in eachindex(i_δk)]
    if Zrenorm
        for (p, val) in rdata
            for (i, i_dk) in enumerate(i_δk)
                _dmn[i][p] = (val[1, kF_label+i_dk] - val[1, kF_label-i_dk]) / (kgrid[kF_label+i_dk] - kgrid[kF_label-i_dk]) * para.me / para.kF
            end
        end
        for (i, i_dk) in enumerate(i_δk)
            dmn, dμn, dm = CounterTerm.sigmaCT(para.order, mu, _dmn[i])
            δMn = [dmn[1], (dmn[1])^2 + dmn[2], (dmn[1])^3 + 2dmn[1] * dmn[2] + dmn[3]]
            δrn = [δMn[1] + δsn[1], δMn[2] + δMn[1] * δsn[1] + δsn[2], δMn[3] + δMn[2] * δsn[1] + δMn[1] * δsn[2] + δsn[3]]
            # δMn = [dmn[1], (dmn[1])^2 + dmn[2], (dmn[1])^3 + 2dmn[1] * dmn[2] + dmn[3],
            # (dmn[1])^4 + 3(dmn[1])^2 * dmn[2] + 2 * dmn[1] * dmn[3] + dmn[4]]
            # δrn = [δMn[1] + δsn[1], δMn[2] + δMn[1] * δsn[1] + δsn[2], δMn[3] + δMn[2] * δsn[1] + δMn[1] * δsn[2] + δsn[3],
            # δMn[4] + δMn[3] * δsn[1] + δMn[2] * δsn[2] + δMn[1] * δsn[3] + δsn[4]]
            println("δk/kF = ", (kgrid[kF_label+i_dk] - kgrid[kF_label]) / para.kF)
            for j in eachindex(δrn)
                # push!(m_eff, 1.0 / (1.0 + sum(dmi[1:i])*z_factor[i])) 
                println("order $j: δm_$j = $(dmn[j])  (m*/m)_$j =", 1.0 + sum(δrn[1:j]))
            end
        end
    else
        rS_dk = []
        for val in rSw_k
            push!(rS_dk, @. (val[1, :] - val[1, kF_label]) / (kgrid - kgrid[kF_label]) * para.me / para.kF)
            dm = []
            for o in 1:para.order
                y = sum([v.val for v in rS_dk[j]] for j in 1:o)
                e = [v.err for v in rS_dk[o]]
                push!(dm, measurement((y[kF_label+i_δk] + y[kF_label-i_δk]) / 2, (e[kF_label+i_δk] + e[kF_label-i_δk]) / 2 + abs(y[kF_label+i_δk] - y[kF_label-i_δk]) / 2))
            end
        end
        meff = @. 1.0 / (z * (1.0 - dm))
        println(meff)
    end

    return rSw_k, iSw_k, kgrid
    # return meff
    # return renormalize(para, rdata, Zrenorm), renormalize(para, idata, Zrenorm), kgrid, ngrid
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

function plotS_k(para, rSw_k, iSw_k, kgrid)
    dim, β, kF = para.dim, para.β, para.kF
    kF_label = searchsortedfirst(kgrid, kF)

    # plot = pyimport("plot")
    style = PyPlot.matplotlib."style"
    style.use(["science", "std-colors"])
    color = [cdict["blue"], cdict["red"], cdict["cyan"], "black"]
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
        return [(sigma[o][2, ki] .- sigma[o][1, ki]) ./ (2π / para.β) for o in 1:order]
        # return [(sigma[o][1, ki]) ./ (π / para.β) for o in 1:order]
    end

    # power series of z-factor
    function sw(sigma, order, kgrid)
        println(sigma[order])
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
        y = [z.val for z in z]
        e = [z.err for z in z]
        errorbar(_kgrid / kF, y, yerr=e, color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="Order $o")

        _x, _y = spline(_kgrid / kF, y, e)
        plot(_x, _y, color=color[o], linestyle="--")
    end
    xlim([kgrid[1] / kF, 2.2])
    ylim([-0.26, -0.06])
    xlabel("\$k/k_F\$")
    ylabel(L"$(\Sigma_{k, i\omega_0}-\Sigma_{0, i\omega_0})/(k^2/2m)$")
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
        # _kgrid, s_w = sw_inv(iSw_k, o, kgrid)
        # z = [(1.0 + e) for e in s_w]
        _kgrid, s_w = sw(iSw_k, o, kgrid)
        z = s_w
        # z = [1.0 / (1.0 + e) for e in s_w]
        y = [-z.val for z in z]
        e = [z.err for z in z]
        kidx = searchsortedfirst(_kgrid, kF)
        println("order $o at k=$(_kgrid[kidx]/kF): $(z[kidx])")

        errorbar(_kgrid / kF, y, yerr=e, color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="Order $o")

        _x, _y = spline(_kgrid / kF, y, e)
        plot(_x, _y, color=color[o], linestyle="--")
    end
    xlim([kgrid[1] / kF, 2.0])
    # ylim([0.0, 1.0])
    xlabel("\$k/k_F\$")
    ylabel("\$\\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    # if Zrenorm
    #     ylabel("\$z \\cdot \\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    # else
    #     ylabel("\$\\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    # end
    # plot.plt.ylabel("\$z(k/k_F) = \\left( 1+\\frac{\\partial \\operatorname{Im}\\Sigma(k, i\\omega_0)}{\\partial \\omega}\\right)^{-1}\$")
    legend()
    savefig("sigmaK_rs$(para.rs)_Fs$(para.Fs)_lam$(para.mass2)_$(para.dim)d.pdf")
    # show()
    # readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    f = jldopen(filenameZ, "r")

    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim)
        for key in keys(f)
            if UEG.paraid(f[key][1]) == UEG.paraid(para)
                # derive_z(f[key], isSave)
                rSw_k, iSw_k, kgrid = process(para, Zrenorm, filenameK, i_δk)
                plotS_k(para, rSw_k, iSw_k, kgrid)
            end
        end
    end

    # para = UEG.ParaMC(rs=0.5, beta=50.0, Fs=-0.0, order=4, mass2=1.2, isDynamic=isDynamic, dim=dim)
end