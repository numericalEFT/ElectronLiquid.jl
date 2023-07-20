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

const filename_inst = "data_instantaneous_polarization_ChargeCharge_K.jld2"
const filename_full = "data_polarization_ChargeCharge_KT.jld2"
const parafilename = "para_wn_1minus0.csv"
# const parafilename = "para.csv"

# Vibrant qualitative colour scheme from https://personal.sron.nl/~pault/
const cdict = Dict(["blue" => "#0077BB", "cyan" => "#33BBEE", "teal" => "#009988", "orange" => "#EE7733", "red" => "#CC3311", "magenta" => "#EE3377", "grey" => "#BBBBBB"]);

"""Returns the Hartree-Fock static structure factor of the UEG, S_HF(q) = -Π₀(q, τ=0) / n₀."""
function fock_static_structure_factor(q, para::ParaMC)
    x = q / para.kF
    if x < 2
        return 3x / 4.0 - x^3 / 16.0
    end
    return 1.0
end

function loaddata(para, FileName)
    key = UEG.short(para)
    f = jldopen(FileName, "r")
    tgrid, kgrid, polarization = f[key]

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
    return tgrid, kgrid, data
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
    tgrid, kgrid, data = loaddata(para, filename)
    return renormalize(para, data), kgrid, tgrid
end

function spline(x, y, e)
    w = 1.0 ./ e

    # generate knots with spline without constraints
    spl = interp.UnivariateSpline(x, y, w=w, k=3)
    __x = collect(LinRange(0.0, x[end], 100))
    yfit = spl(__x)
    return __x, yfit
end

function plot_instantaneous_polarization_vs_k(para)
    pi_tau_k, kgrid, tgrid = process(para; filename=filename_inst)

    dim, β, kF = para.dim, para.β, para.kF
    n0 = kF^3 / 3π^2  # non-interacting density

    # Index where τ ≈ 0 in tgrid
    it0 = searchsortedfirst(tgrid, 0.0)

    # Get exact instantaneous bare polarization Π₀(q, τ = 0)
    s_hf_exact = fock_static_structure_factor.(kgrid, [para])

    function instantaneous_polarization_over_n0(pi_tau_k, order)
        return sum(pi_tau_k[o][it0, :] ./ n0 for o in 1:order)
    end

    style = PyPlot.matplotlib."style"
    style.use(["science", "std-colors"])
    color = [cdict["blue"], "green", cdict["orange"], cdict["red"], cdict["magenta"]]
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 16
    rcParams["font.family"] = "Times New Roman"
    figure(figsize=(4, 4))
    # axvline(2.0; linestyle="--", linewidth=1, color="gray")
    axhline(1.0; linestyle="--", linewidth=1, color="gray")
    plot(kgrid / kF, s_hf_exact; color="black", label="\$S_{\\mathrm{HF}}(q)\$")
    for o in 1:para.order
        pi_k_inst_over_n0 = instantaneous_polarization_over_n0(pi_tau_k, o)
        println(pi_k_inst_over_n0)
        y, e = Measurements.value.(pi_k_inst_over_n0), Measurements.uncertainty.(pi_k_inst_over_n0)

        errorbar(kgrid / kF, y, yerr=e; color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="Order $o")

        _x, _y = spline(kgrid / kF, y, e)
        plot(_x, _y, color=color[o], linestyle="--")
    end
    xlim([kgrid[1] / kF, 3.0])
    xlabel("\$q / k_F\$")
    # NOTE: We measure -Π using Negele & Orland conventions
    ylabel("\$-\\Pi(q, \\tau=0) / n_0\$")
    legend(; loc="best")
    savefig("instantaneous_polarization.pdf")
    close("all")
    # show()
end

function plot_polarization_vs_k(para; tau_over_beta_target=0.0)
    pi_tau_k, kgrid, tgrid = process(para; filename=filename_full)

    dim, β, kF = para.dim, para.β, para.kF
    n0 = kF^3 / 3π^2  # non-interacting density

    # Index where τ ≈ 0 in tgrid
    it = searchsortedfirst(tgrid / β, tau_over_beta_target)
    tau_over_beta_approx = round(tgrid[it] / β; sigdigits=3)

    function polarization(pi_tau_k, order)
        return sum(pi_tau_k[o][it, :] for o in 1:order)
    end

    style = PyPlot.matplotlib."style"
    style.use(["science", "std-colors"])
    color = [cdict["blue"], "green", cdict["orange"], cdict["red"], cdict["magenta"]]
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 16
    rcParams["font.family"] = "Times New Roman"
    figure(figsize=(4, 4))
    for o in 1:para.order
        pi_k = polarization(pi_tau_k, o)
        y, e = Measurements.value.(pi_k), Measurements.uncertainty.(pi_k)

        errorbar(kgrid / kF, y, yerr=e; color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="Order $o")

        _x, _y = spline(kgrid / kF, y, e)
        plot(_x, _y, color=color[o], linestyle="--")
    end
    xlim([kgrid[1] / kF, 2.0])
    xlabel("\$q / k_F\$")
    # NOTE: We measure -Π using Negele & Orland conventions
    ylabel("\$-\\Pi(q, \\tau=$(tau_over_beta_approx)\\beta)\$")
    legend(; loc="best")
    savefig("polarization_tau=$(tau_over_beta_approx)beta.pdf")
    close("all")
    # show()
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Plot instantaneous polarization vs k
    para_inst = ParaMC(rs=1.0, beta=40.0, Fs=-0.0, order=4, mass2=1.0, isDynamic=false)
    plot_instantaneous_polarization_vs_k(para_inst)

    # # Plot polarization vs k at a fixed τ / β
    # τ_over_β = 0.5
    # para_full = ParaMC(rs=1.0, beta=40.0, Fs=-0.0, order=4, mass2=1.0, isDynamic=false)
    # plot_polarization_vs_k(para_full; tau_over_beta_target=τ_over_β)
end
