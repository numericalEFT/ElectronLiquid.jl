# include("../common/counterterm.jl")
# include("../common/para_builder.jl")
# using .UEG
# using .CounterTerm
using ElectronLiquid
using Printf

rs = [5.0,]
mass2 = [0.01, ]
Fs = [-0.0,]
beta = [25.0, 50.0, 100.0, ]
order = [3,]

const filename = "para.csv"

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=true)

    # Zrenorm = false
    # Zrenorm = true

    mu, sw = CounterTerm.getSigma(para, parafile=filename)
    dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)

    println("δz = ", dz)
    println("δz_inverse = ", dzi)
    sumzi = accumulate(+, dzi)
    # if Zrenorm == false

    z1 = @. 1.0 / (1.0 + sumzi)
    # println("z without renormalization = $z")
    # else
    sumz = accumulate(+, dz)
    z2 = @. 1.0 + sumz
    # end
    # println("z with renormalization = $z")
    printstyled("$(UEG.short(para))\n", color=:green)
    printstyled(@sprintf("%8s   %24s    %24s     %24s\n", "order", "chemical-potential", "no-z-renorm", "z-renorm"), color=:yellow)
    for o in 1:para.order
        @printf("%8d   %24s   %24s     %24s\n", o, "$(dmu[o])", "$(z1[o])", "$(z2[o])")
    end
end