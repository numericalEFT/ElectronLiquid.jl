using JLD2
using Measurements
# using Gaston
using PyCall
using PyPlot
using ElectronLiquid
using ElectronGas
using TaylorSeries
using DelimitedFiles

dim = 2
isDynamic = false
# rs = [3.0]
# rs = [1.0, 2.0, 3.0]
# mass2 = [1.2]
rs = [0.5]
# mass2 = [2.0, 3.0, 4.0]
# mass2 = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4]
# mass2 = [0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.0, 3.2, 3.5, 3.6, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
# mass2 = [4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
mass2 = [4.0]
Fs = [-0.0,]
beta = [50.0,]
order = [4,]
# order = [3,]
# i_δk = [1, 2, 3, 4]
idx_δk = [1, 2, 3]
# i_δk = [1, 2, 3, 4, 5, 6]
δks = [0.015, 0.03, 0.05, 0.1]
# δks = [0.03, 0.05, 0.1]

const Zrenorm = true     # turn on/off Z renormalization 
const filenameZ = "data_Z.jld2"
# const filenameZ = "data_Z_all3.jld2"
# const filenameZ = "data_Z_large.jld2"
# const filenameZ = "data_Z_scan.jld2"
const filenameK = "data_K.jld2"
# const filenameK = "data_K_all3.jld2"
# const filenameK = "data_K_large.jld2"
# const filenameK = "data_K_scan.jld2"
const parafilename = "para_wn_1minus0.csv"
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

function process(para, Zrenorm, filename)
    println("rs=$(para.rs), beta=$(para.beta), mass2=$(para.mass2)")

    ngrid, kgrid, rdata, idata = loaddata(para, filename)
    kF_label = searchsortedfirst(kgrid, para.kF)
    # rSw_k, iSw_k = renormalize(para, rdata), renormalize(para, idata)

    i_δk = collect(1:Int((length(kgrid) - 1) / 2))
    i_δk != idx_δk && (i_δk = idx_δk .+ 1)
    println(kgrid)
    # println(ngrid)
    println("Re part:")
    for (p, val) in rdata
        println(p, val)
    end

    println("Im part:")
    for (p, val) in idata
        println(p, val)
    end

    mu, sw = CounterTerm.getSigma(para, parafile=parafilename)
    δsn, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)
    rSw_k = CounterTerm.chemicalpotential_renormalization(para.order, rdata, dmu)
    iSw_k = CounterTerm.chemicalpotential_renormalization(para.order, idata, dmu)
    # println(rSw_k)

    _dmn = [Dict() for i in eachindex(i_δk)]
    println(δsn)
    meff = Matrix(undef, length(i_δk), para.order)
    δks = []
    if Zrenorm
        for (p, val) in rdata
            for (i, i_dk) in enumerate(i_δk)
                _dmn[i][p] = (val[1, kF_label+i_dk] - val[1, kF_label-i_dk]) / (kgrid[kF_label+i_dk] - kgrid[kF_label-i_dk]) * para.me / para.kF
            end
        end
        for (i, i_dk) in enumerate(i_δk)
            # meff = []
            dmn, dμn, dm = CounterTerm.sigmaCT(para.order, mu, _dmn[i])
            # δMn = [dmn[1], (dmn[1])^2 + dmn[2], (dmn[1])^3 + 2dmn[1] * dmn[2] + dmn[3]]
            # δrn = [δMn[1] + δsn[1], δMn[2] + δMn[1] * δsn[1] + δsn[2], δMn[3] + δMn[2] * δsn[1] + δMn[1] * δsn[2] + δsn[3]]
            δMn = [dmn[1], (dmn[1])^2 + dmn[2], (dmn[1])^3 + 2dmn[1] * dmn[2] + dmn[3],
                (dmn[1])^4 + 3(dmn[1])^2 * dmn[2] + 2 * dmn[1] * dmn[3] + dmn[4]]
            δrn = [δMn[1] + δsn[1], δMn[2] + δMn[1] * δsn[1] + δsn[2], δMn[3] + δMn[2] * δsn[1] + δMn[1] * δsn[2] + δsn[3],
                δMn[4] + δMn[3] * δsn[1] + δMn[2] * δsn[2] + δMn[1] * δsn[3] + δsn[4]]
            δk_kF = (kgrid[kF_label+i_dk] - kgrid[kF_label]) / para.kF
            println("δk/kF = ", δk_kF)
            for j in eachindex(δrn)
                # push!(m_eff, 1.0 / (1.0 + sum(dmi[1:i])*z_factor[i])) 
                _meff = 1.0 + sum(δrn[1:j])
                println("order $j: δm_$j = $(dmn[j])  (m*/m)_$j =", _meff)
                # push!(meff, _meff)
                meff[i, j] = _meff
            end
            push!(δks, δk_kF)
            # push!(data, δk_kF, para.rs, para.beta, para.mass2, para.order)
            # fname = "./meff_2D.txt"
            # open(fname, "a+") do io
            #     writedlm(io, data, ' ')
            # end
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

    return rSw_k, iSw_k, kgrid, meff, i_δk
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

function plotS_k(para, rSw_k, iSw_k, kgrid, Zrenorm=false)
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

    function delta_z(sigma, order, ki)
        return [(sigma[o][2, ki] .- sigma[o][1, ki]) ./ (2π / para.β) for o in 1:order]
        # return [(sigma[o][1, ki]) ./ (π / para.β) for o in 1:order]
    end

    # power series of z-factor
    function z_series(iSigma, order)
        δzi = delta_z(iSigma, order, kF_label)
        z = 1 / Taylor1([1.0, δzi...], order)
        return [getcoeff(z, o) for o in 0:order]
    end

    function sk(sigma, MaxOrder, kgrid)
        dk = [(sigma[o][1, 2:end] .- sigma[o][1, 1]) ./ (kgrid[2:end] .^ 2 / (2 * para.me)) for o in 1:MaxOrder]
        # return kgrid[2:end], sum(dk)
        return kgrid[2:end], dk
    end

    # power series of z-factor
    function sw(sigma, MaxOrder, kgrid)
        # dw = [sum(delta_z(sigma, order, ki)) for ki in 1:length(kgrid)]
        # dw = [(sigma[o][1, :]) ./ (π / para.β) for o in 1:order]
        dw = [(sigma[o][2, :] .- sigma[o][1, :]) ./ (2π / para.β) for o in 1:MaxOrder]
        return kgrid, dw
    end

    function Zrenormalize(Zseries, f_k, MaxOrder, kgrid)
        # println("zseris:", Zseries)
        # z = Taylor1(Zseries, order)
        # println("z:", z)
        fk = hcat(f_k...)
        # println(f_k)
        fk_renorm = []
        for ki in eachindex(kgrid)
            fki = Taylor1(fk[ki, :], MaxOrder)
            fk_renorm_series = [getcoeff(fki * Zseries, o) for o in 0:MaxOrder-1]
            # push!(fk_renorm, sum(fk_renorm_series[1:end]))
            push!(fk_renorm, fk_renorm_series)
            if ki == kF_label
                println("fki:", fki)
                # println("fk_renorm:", fk_renorm_series)
            end
        end
        # println("fk_renorm:", fk_renorm)
        fk_renorm = hcat(fk_renorm...)
        # println(fk_renorm)
        return [fk_renorm[o, :] for o in 1:MaxOrder]
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

    if Zrenorm
        # Zseries = z_series(iSw_k, para.order)
        Zseries = Taylor1([1.0, delta_z(iSw_k, para.order, kF_label)...], para.order)
        println("Zfactor:", Zseries)
    end

    _kgrid, s_k = sk(rSw_k, para.order, kgrid)
    println(s_k)
    Zrenorm && (s_k = Zrenormalize(Zseries, s_k, para.order, _kgrid))
    println(s_k)
    for o in 1:para.order
        z = sum(s_k[1:o])
        # _kgrid, s_k = sk(rSw_k, o, kgrid)
        # if Zrenorm
        #     Zseries = z_series(iSw_k, o)
        #     s_k = Zrenormalize(Zseries, s_k, o, _kgrid)
        # else
        #     s_k = sum(s_k)
        # end
        y = [z.val for z in z]
        e = [z.err for z in z]
        errorbar(_kgrid / kF, y, yerr=e, color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="Order $o")

        _x, _y = spline(_kgrid / kF, y, e)
        plot(_x, _y, color=color[o], linestyle="--")
    end
    xlim([kgrid[1] / kF, 2.2])
    # ylim([-0.26, -0.06])
    xlabel("\$k/k_F\$")
    # plot.plt.ylabel("\$z(k/k_F) = \\left( 1+\\frac{\\partial \\operatorname{Im}\\Sigma(k, i\\omega_0)}{\\partial \\omega}\\right)^{-1}\$")
    if Zrenorm
        ylabel("\$z \\cdot (\\Sigma_{k, i\\omega_0}-\\Sigma_{0, i\\omega_0})/(k^2/2m)\$")
    else
        ylabel("\$(\\Sigma_{k, i\\omega_0}-\\Sigma_{0, i\\omega_0})/(k^2/2m)\$")
    end
    # plot.plt.ylabel("\$z(k/k_F) = \\left( 1+\\frac{\\partial \\operatorname{Im}\\Sigma(k, i\\omega_0)}{\\partial \\omega}\\right)^{-1}\$")
    legend()

    subplot(1, 2, 1)
    _kgrid, s_w = sw(iSw_k, para.order, kgrid)
    Zrenorm && (s_w = Zrenormalize(Zseries, s_w, para.order, _kgrid))
    for o in 1:para.order
        # _kgrid, s_w = sw_inv(iSw_k, o, kgrid)
        # z = [(1.0 + e) for e in s_w]
        # _kgrid, s_w = sw(iSw_k, o, kgrid)
        # if Zrenorm
        #     Zseries = z_series(iSw_k, o)
        #     s_w = Zrenormalize(Zseries, s_w, o, _kgrid)
        # else
        #     s_w = sum(s_w)
        # end
        z = sum(s_w[1:o])

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
    # ylabel("\$\\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    if Zrenorm
        ylabel("\$z \\cdot \\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    else
        ylabel("\$\\frac{\\partial \\Sigma(k, i\\omega_0)}{\\partial i \\omega}\$")
    end
    # plot.plt.ylabel("\$z(k/k_F) = \\left( 1+\\frac{\\partial \\operatorname{Im}\\Sigma(k, i\\omega_0)}{\\partial \\omega}\\right)^{-1}\$")
    legend()
    if Zrenorm
        savefig("sigmaK_rs$(para.rs)_Fs$(para.Fs)_lam$(para.mass2)_$(para.dim)d_Zrenorm.pdf")
    else
        savefig("sigmaK_rs$(para.rs)_Fs$(para.Fs)_lam$(para.mass2)_$(para.dim)d.pdf")
    end
    # show()
    # readline()
end

function plot_optimal(meff; mass2=mass2, idx_dk=length(idx_δk))
    style = PyPlot.matplotlib."style"
    style.use(["science", "std-colors"])
    color = [cdict["blue"], cdict["red"], cdict["cyan"], "black"]
    #cmap = get_cmap("Paired")
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 16
    rcParams["font.family"] = "Times New Roman"
    figure(figsize=(6, 4))

    for o in 1:order[1]
        y = [z.val for z in meff[idx_dk, o, :]]
        yerr = [z.err for z in meff[idx_dk, o, :]]
        # y = [z.val for z in meff[end, o, :]]
        # yerr = [z.err for z in meff[end, o, :]]
        errorbar(mass2, y, yerr=yerr, color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="Order $o")
    end
    xlabel("Mass2")
    ylabel("\$m^*/m\$")
    # ylim(0.9, 1.2)
    legend()
    title("rs=$(rs[1])")
    # savefig("meffvslam_rs$(rs[1])_$(dim)d_comput100_v1.pdf")
    savefig("meffvslam_rs$(rs[1])_$(dim)d_CCQ_scan.pdf")
end

function plot_convergence(meff, dks; idx_mass2=1)
    style = PyPlot.matplotlib."style"
    style.use(["science", "std-colors"])
    color = [cdict["blue"], cdict["red"], cdict["teal"], "black"]
    #cmap = get_cmap("Paired")
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 16
    rcParams["font.family"] = "Times New Roman"
    figure(figsize=(6, 4))

    x = [i for i in 1:order[1]]
    for (idx, dk) in enumerate(dks)
        y = [z.val for z in meff[idx, :, idx_mass2]]
        yerr = [z.err for z in meff[idx, :, idx_mass2]]
        errorbar(x, y, yerr=yerr, color=color[idx], capsize=4, fmt="o", markerfacecolor="none", label="$dk")
    end
    r_s, lam = rs[1], mass2[idx_mass2]
    xlabel("Order")
    ylabel("\$m^*/m\$")
    legend(title="dk/kF")
    title("rs=$(r_s), mass2=$lam")
    # savefig("meffvslam_rs$(rs[1])_$(dim)d_comput100_v1.pdf")
    savefig("meffvsOrder_rs$(r_s)_lam$(lam)_$(dim)d_CCQ.pdf")
end

if abspath(PROGRAM_FILE) == @__FILE__
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    f = jldopen(filenameZ, "r")

    effmass = []
    mass2_existed = []
    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim)
        println("### Mass2:", _mass2)
        for key in keys(f)
            if UEG.paraid(f[key][1]) == UEG.paraid(para)
                derive_z(f[key], isSave)
                rSw_k, iSw_k, kgrid, meff, i_δk = process(para, Zrenorm, filenameK)
                # plotS_k(para, rSw_k, iSw_k, kgrid, true)
                push!(effmass, meff)
                push!(mass2_existed, _mass2)
            end
        end
    end
    # effmass = hcat(effmass...)
    effmass = cat(effmass..., dims=3)
    # plot_optimal(effmass, mass2=mass2_existed)
    # plot_convergence(effmass, δks)

    # para = UEG.ParaMC(rs=0.5, beta=50.0, Fs=-0.0, order=4, mass2=1.2, isDynamic=isDynamic, dim=dim)
end