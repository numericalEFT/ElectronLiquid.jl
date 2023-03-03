# include("../common/counterterm.jl")
# include("../common/para_builder.jl")
# using .UEG
# using .CounterTerm
using ElectronLiquid
# using LsqFit
using Printf
using PyCall
using LinearAlgebra
using Plots

curve_fit = pyimport("scipy.optimize").curve_fit

rs = [1.0, 2.0, 4.0, 5.0, 6.0, 8.0]
#rs = [2.0, 4.0, ]
mass2 = [0.01, 0.003, 0.001, 0.0003, 0.0001]
Fs = [-0.0,]
beta = [25.0,]
order = [3,]

rs_old = [1.0, 2.0, 3.0, 4.0]
z_old = [0.8725, 0.7984, 0.7219, 0.6571]
z_error_old = [0.0002, 0.0002, 0.0002, 0.0002]

nz = zeros(order[1], length(rs))
nz_error = similar(nz)
z = similar(nz)
z_error = similar(nz)


# const filename = "para.csv"
const filename = "para_wn_1minus0.csv"

cache = []

# @. model(x, p) = p[1] * x + p[2]

# function f(x, a, b)
#     return a * x + b
# end

py"""
def f(x, a, b):
    return a * x + b
"""

for (_rs, _F, _beta, _order) in Iterators.product(rs, Fs, beta, order)
    norenorm_z = zeros(_order, length(mass2), 2)
    renorm_z = zeros(_order, length(mass2), 2)
    for (mi, _mass2) in enumerate(mass2)
        para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=true)

        # Zrenorm = false
        # Zrenorm = true

        mu, sw = CounterTerm.getSigma(para, parafile=filename)
        dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)

        println("δz = ", dz)
        println("δμ = ", dmu)
        println("δz_inverse = ", dzi)
        sumzi = accumulate(+, dzi)
        # if Zrenorm == false

        z1 = @. 1.0 / (1.0 + sumzi)
        # z1 = sumzi
        # println("z without renormalization = $z")
        push!(cache, sumzi[3])

        # else
        sumz = accumulate(+, dz)
        z2 = @. 1.0 + sumz
        # end
        # println("z with renormalization = $z")
        printstyled("$(UEG.short(para))\n", color=:green)
        printstyled(@sprintf("%8s   %24s    %24s     %24s\n", "order", "chemical-potential", "no-z-renorm", "z-renorm"), color=:yellow)
        for o in 1:para.order
            norenorm_z[o, mi, 1] = z1[o].val
            norenorm_z[o, mi, 2] = z1[o].err
            renorm_z[o, mi, 1] = z2[o].val
            renorm_z[o, mi, 2] = z2[o].err
            @printf("%8d   %24s   %24s     %24s\n", o, "$(dmu[o])", "$(z1[o])", "$(z2[o])")
        end
    end

    rsi = searchsortedfirst(rs, _rs)
    println("without renormalization fit")
    for o in 1:_order
        x = mass2 .^ 0.5
        popt, pcov = curve_fit(py"f", x, norenorm_z[o, :, 1])
        println(o, " ", popt[1], " +- ", sqrt.(diag(pcov))[1], " ", popt[2], " +- ", sqrt.(diag(pcov))[2])
        nz[o, rsi] = popt[2]
        nz_error[o, rsi] = sqrt.(diag(pcov))[2]
    end
    println("with renormalization fit")
    for o in 1:_order
        x = mass2 .^ 0.5
        popt, pcov = curve_fit(py"f", x, renorm_z[o, :, 1])
        println(o, " ", popt[1], " +- ", sqrt.(diag(pcov))[1], " ", popt[2], " +- ", sqrt.(diag(pcov))[2])
        z[o, rsi] = popt[2]
        z_error[o, rsi] = sqrt.(diag(pcov))[2]
    end
end

p = plot(rs, nz[3, :], yerror=abs.(nz[3, :] .- nz[2, :]) .+ nz_error[3, :], label="without renormalization")
plot!(p, rs, z[3, :], yerror=abs.(z[3, :] .- z[2, :]) .+ z_error[3, :], label="with renormalization")
plot!(p, rs_old, z_old, yerror=z_error_old, label="Previous")
display(p)
readline()

