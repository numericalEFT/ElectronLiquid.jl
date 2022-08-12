# include("../common/counterterm.jl")
# include("../common/para_builder.jl")
# using .UEG
# using .CounterTerm
using ElectronLiquid
using Printf

para = UEG.ParaMC(
    rs=3.0,
    # Fs=-0.585,
    Fs=-0.0,
    beta=25.0,
    mass2=0.001,
    order=3,
    isDynamic=true,
)
# Zrenorm = false
# Zrenorm = true

mu, sw = CounterTerm.getSigma(para)
dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)

# println("δz = ", dz)
# println("δz_inverse = ", dzi)
sumzi = accumulate(+, dzi)
# if Zrenorm == false

z1 = @. 1.0 / (1.0 + sumzi)
# println("z without renormalization = $z")
# else
sumz = accumulate(+, dz)
z2 = @. 1.0 + sumz
# end
# println("z with renormalization = $z")
printstyled(UEG.short(para), color=:green)
printstyled(@sprintf("%8s   %24s    %24s     %24s\n", "order", "chemical-potential", "no-z-renorm", "z-renorm"), color=:yellow)
for o in 1:para.order
    @printf("%8d   %24s   %24s     %24s\n", o, "$(dmu[o])", "$(z1[o])", "$(z2[o])")
end