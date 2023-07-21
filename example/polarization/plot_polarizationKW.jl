using CodecZlib
using ElectronGas
using ElectronLiquid
using JLD2
using Measurements
using PyCall
using PyPlot

# # For style "science"
# @pyimport scienceplots

# For spline interpolations
interp = pyimport("scipy.interpolate")

const filename = "data_polarization_ChargeCharge_KW.jld2"
const parafilename = "para_wn_1minus0.csv"
# const parafilename = "para.csv"

# Vibrant qualitative colour scheme from https://personal.sron.nl/~pault/
const cdict = Dict(["blue" => "#0077BB", "cyan" => "#33BBEE", "teal" => "#009988", "orange" => "#EE7733", "red" => "#CC3311", "magenta" => "#EE3377", "grey" => "#BBBBBB"]);

function loaddata(para, FileName)
    key = UEG.short(para)
    f = jldopen(FileName, "r")
    ngrid, kgrid, polarization = f[key]

    # NOTE: Only a subset of partitions in UEG.partition are measured for Π; the rest are zero 
    partitions = UEG.partition(para.order)
    nonzero_partitions = keys(polarization)
    @assert nonzero_partitions ⊆ partitions

    # Set all vanishing partitions (those not in keys(polarization)) to zero
    data = Dict()
    for P in partitions
        if P in nonzero_partitions
            data[P] = polarization[P][:, :]
        else
            data[P] = zero(polarization[first(nonzero_partitions)][:, :])
        end
    end
    return ngrid, kgrid, data
end

function renormalize(para, polarization)
    dim, β, kF = para.dim, para.β, para.kF
    mu, sw = CounterTerm.getSigma(para, parafile=parafilename)
    dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)
    println(para.order)
    println(dmu)
    polarization = CounterTerm.chemicalpotential_renormalization(para.order, polarization, dmu)
    return polarization
end

function process(para; filename)
    dim, β, kF = para.dim, para.β, para.kF
    ngrid, kgrid, data = loaddata(para, filename)
    return renormalize(para, data), kgrid, ngrid
end

function spline(x, y, e)
    w = 1.0 ./ e

    # generate knots with spline without constraints
    spl = interp.UnivariateSpline(x, y, w=w, k=3)
    __x = collect(LinRange(0.0, x[end], 100))
    yfit = spl(__x)
    return __x, yfit
end

function plot_polarization_kw_vs_k(pi_w_k, kgrid, ngrid; n=0)

    dim, β, kF = para.dim, para.β, para.kF

    if n ∉ ngrid
        error("Requested frequency point n = $n is not included in the measurement ngrid!")
    end

    # Index of n in ngrid
    iw = searchsortedfirst(ngrid, n)

    function polarization(pi_w_k, order)
        return sum(pi_w_k[o][iw, :] for o in 1:order)
    end

    style = PyPlot.matplotlib."style"
    style.use(["science", "std-colors"])
    color = [cdict["blue"], "green", cdict["orange"], cdict["red"], cdict["magenta"]]
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 16
    rcParams["font.family"] = "Times New Roman"
    figure(figsize=(4, 4))
    for o in 1:para.order
        pi_k = polarization(pi_w_k, o)
        y, e = Measurements.value.(pi_k), Measurements.uncertainty.(pi_k)

        errorbar(kgrid / kF, y, yerr=e; color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="Order $o")

        _x, _y = spline(kgrid / kF, y, e)
        plot(_x, _y, color=color[o], linestyle="--")
    end
    xlim([kgrid[1] / kF, 3.0])
    xlabel("\$q / k_F\$")
    # NOTE: We measure -Π using Negele & Orland conventions
    ylabel("\$-\\Pi(q, i\\omega_{$n})\$")
    legend(; loc="best")
    savefig("polarization_q_iw$n.pdf")
    close("all")
    # show()
end

if abspath(PROGRAM_FILE) == @__FILE__
    para = ParaMC(rs=1.0, beta=40.0, Fs=-0.0, order=4, mass2=1.0, isDynamic=false)
    
    pi_w_k, kgrid, ngrid = process(para; filename=filename)

    # Plot Π(q, iωₙ) vs k at a fixed frequency point iωₙ (fixed n)
    plot_polarization_kw_vs_k(pi_w_k, kgrid, ngrid; n=0)
end
